#include "patidx.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdint>
#include <chrono>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

/**
 * @brief High-speed builder/converter for CDRPAT01 binary patent index files.
 *
 * Reads:
 *   1. patents.txt (line-separated patent document IDs)
 *   2. occurrences.bin (raw RawOccLog binary stream of (oid, pat_idx, seq_id_no))
 *
 * Emits:
 *   Single memory-mappable .patidx file with prefix-sum OID offsets for O(1) lookups.
 */

static void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << " <db_dir> [patents.txt] [occurrences.bin] [output.patidx]\n"
              << "       " << prog << " --help\n\n"
              << "Arguments:\n"
              << "  db_dir           Base directory containing index_state/ or database files\n"
              << "  patents.txt      Path to patents list (default: <db_dir>/index_state/patents.txt)\n"
              << "  occurrences.bin  Path to raw occurrences (default: <db_dir>/index_state/occurrences.bin)\n"
              << "  output.patidx    Path to output file (default: <db_dir>/wipo.patidx)\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2 || std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h") {
        print_usage(argv[0]);
        return 0;
    }

    std::string db_dir = argv[1];
    std::string patents_path = (argc > 2) ? argv[2] : (db_dir + "/index_state/patents.txt");
    std::string occ_path     = (argc > 3) ? argv[3] : (db_dir + "/index_state/occurrences.bin");
    std::string out_path     = (argc > 4) ? argv[4] : (db_dir + "/wipo.patidx");

    auto t0 = std::chrono::steady_clock::now();
    std::cout << "[1/4] Reading patent list from: " << patents_path << "..." << std::endl;

    std::ifstream pat_file(patents_path);
    if (!pat_file) {
        std::cerr << "Error: Failed to open patents file: " << patents_path << std::endl;
        return 1;
    }

    std::vector<std::string> patents;
    std::string line;
    while (std::getline(pat_file, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        patents.push_back(line);
    }
    pat_file.close();

    uint32_t num_patents = static_cast<uint32_t>(patents.size());
    std::cout << "  Loaded " << num_patents << " patents." << std::endl;

    // Build StringPool and PatentMeta
    std::vector<char> str_pool;
    str_pool.push_back('\0'); // offset 0 is empty string

    std::vector<SPatIdxPatentMeta> pat_meta(num_patents);
    for (uint32_t i = 0; i < num_patents; ++i) {
        uint32_t offset = static_cast<uint32_t>(str_pool.size());
        for (char c : patents[i]) {
            str_pool.push_back(c);
        }
        str_pool.push_back('\0');
        pat_meta[i].doc_id_offset = offset;
        pat_meta[i].is_patent = 1;
    }
    uint64_t pat_strings_size = str_pool.size();

    std::cout << "[2/4] Counting occurrences by OID from: " << occ_path << "..." << std::endl;
    int occ_fd = open(occ_path.c_str(), O_RDONLY);
    if (occ_fd < 0) {
        std::cerr << "Error: Failed to open occurrences file: " << occ_path << std::endl;
        return 1;
    }
    struct stat occ_st;
    fstat(occ_fd, &occ_st);
    uint64_t total_occ = occ_st.st_size / sizeof(RawOccLog);
    std::cout << "  Total occurrences: " << total_occ << " (" << (occ_st.st_size / 1024 / 1024) << " MB)" << std::endl;

    const RawOccLog* raw_occs = reinterpret_cast<const RawOccLog*>(
        mmap(nullptr, occ_st.st_size, PROT_READ, MAP_SHARED, occ_fd, 0)
    );
    if (raw_occs == MAP_FAILED) {
        std::cerr << "Error: Failed to mmap occurrences file: " << occ_path << std::endl;
        close(occ_fd);
        return 1;
    }

    // Determine maximum OID
    uint32_t max_oid = 0;
    for (uint64_t i = 0; i < total_occ; ++i) {
        if (raw_occs[i].oid > max_oid) {
            max_oid = raw_occs[i].oid;
        }
    }
    uint32_t num_oids = max_oid + 1;
    std::cout << "  Max OID: " << max_oid << " (num_oids = " << num_oids << ")" << std::endl;

    // Count occurrences per OID
    std::vector<uint32_t> counts(num_oids + 1, 0);
    for (uint64_t i = 0; i < total_occ; ++i) {
        counts[raw_occs[i].oid]++;
    }

    std::cout << "[3/4] Computing prefix sums for OID offsets..." << std::endl;
    std::vector<uint64_t> oid_offsets(num_oids + 1, 0);
    uint64_t running = 0;
    for (uint32_t i = 0; i < num_oids; ++i) {
        oid_offsets[i] = running;
        running += counts[i];
    }
    oid_offsets[num_oids] = running;
    if (running != total_occ) {
        std::cerr << "Error: Mismatch: running = " << running << ", total_occ = " << total_occ << std::endl;
        return 1;
    }

    std::cout << "[4/4] Placing sorted occurrences..." << std::endl;
    std::vector<uint64_t> write_pos = oid_offsets; // copy initial offsets
    std::vector<SPatIdxOccurrence> sorted_occs(total_occ);

    for (uint64_t i = 0; i < total_occ; ++i) {
        uint32_t oid = raw_occs[i].oid;
        uint64_t pos = write_pos[oid]++;
        sorted_occs[pos].pat_idx = raw_occs[i].pat_idx;
        sorted_occs[pos].seq_id_no = raw_occs[i].seq_id_no;
    }

    munmap(const_cast<RawOccLog*>(raw_occs), occ_st.st_size);
    close(occ_fd);

    // Write binary CDRPAT01 file
    std::cout << "Writing output index to: " << out_path << "..." << std::endl;
    std::ofstream out(out_path, std::ios::binary);
    if (!out) {
        std::cerr << "Error: Failed to create output file: " << out_path << std::endl;
        return 1;
    }

    SPatIdxHeader hdr;
    std::memset(&hdr, 0, sizeof(hdr));
    std::memcpy(hdr.magic, "CDRPAT01", 8);
    hdr.version = 1;
    hdr.num_patents = num_patents;
    hdr.num_oids = num_oids;
    hdr.reserved = 0;
    hdr.total_occurrences = total_occ;
    hdr.pat_strings_size = pat_strings_size;
    hdr.reserved2 = 0;

    out.write(reinterpret_cast<const char*>(&hdr), sizeof(hdr));
    out.write(reinterpret_cast<const char*>(pat_meta.data()), pat_meta.size() * sizeof(SPatIdxPatentMeta));
    out.write(reinterpret_cast<const char*>(oid_offsets.data()), oid_offsets.size() * sizeof(uint64_t));
    out.write(reinterpret_cast<const char*>(sorted_occs.data()), sorted_occs.size() * sizeof(SPatIdxOccurrence));
    out.write(str_pool.data(), str_pool.size());
    out.close();

    auto t1 = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(t1 - t0).count();
    struct stat out_st;
    stat(out_path.c_str(), &out_st);
    std::cout << "[DONE] Successfully generated " << out_path << " (" 
              << (out_st.st_size / 1024 / 1024) << " MB) in " << elapsed << " seconds." << std::endl;

    return 0;
}
