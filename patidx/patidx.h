#ifndef CDR_PATIDX_H
#define CDR_PATIDX_H

#include <cstdint>

#pragma pack(push, 1)

/**
 * @brief Fixed 48-byte header for CDRPAT01 binary patent index files.
 *
 * File layout:
 *   [0 .. 47]                            : SPatIdxHeader (48 bytes)
 *   [48 .. 48 + P*8 - 1]                 : SPatIdxPatentMeta[P] (P * 8 bytes)
 *   [.. + (N+1)*8 - 1]                   : uint64_t oid_offsets[N + 1] ((N + 1) * 8 bytes)
 *   [.. + M*8 - 1]                       : SPatIdxOccurrence[M] (M * 8 bytes)
 *   [.. + S - 1]                         : char string_pool[S] (null-terminated strings)
 */
struct SPatIdxHeader {
    char     magic[8];          // Magic bytes: "CDRPAT01"
    uint32_t version;           // Format version: 1
    uint32_t num_patents;       // Total unique patent documents (P)
    uint32_t num_oids;          // Total sequence OIDs in database (N)
    uint32_t reserved;          // Reserved padding (0)
    uint64_t total_occurrences; // Total sequence-to-patent occurrences (M)
    uint64_t pat_strings_size;  // Size of null-terminated string pool in bytes (S)
    uint64_t reserved2;         // Reserved padding (0)
}; // 48 bytes

/**
 * @brief Metadata record for a patent document.
 */
struct SPatIdxPatentMeta {
    uint32_t doc_id_offset;     // Byte offset into string pool for document ID
    uint32_t is_patent;         // 1 if official patent document, 0 otherwise
}; // 8 bytes

/**
 * @brief Occurrence entry mapping an OID to a patent document and sequence number.
 */
struct SPatIdxOccurrence {
    uint32_t pat_idx;           // Index into patent metadata array (0 .. P-1)
    uint32_t seq_id_no;         // Sequence number within patent (<212> SEQ ID NO)
}; // 8 bytes

/**
 * @brief Intermediate raw occurrence record logged during extraction.
 */
struct RawOccLog {
    uint32_t oid;               // BLAST sequence OID
    uint32_t pat_idx;           // Patent index
    uint32_t seq_id_no;         // Sequence number in patent
}; // 12 bytes

#pragma pack(pop)

#endif // CDR_PATIDX_H
