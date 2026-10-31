# NCBI BLAST+ Patched: Patent Sequence, HSP Filtering & Antibody CDR Search

High-performance C++ modifications, post-traceback identity/coverage filters, unified antibody/TCR CDR search engine, and ClickHouse storage architecture built on top of **NCBI BLAST+ 2.17.0**.

---

## Architecture & Patch Overview

This repository provides three modular, stackable patch layers that can be compiled into NCBI BLAST+ (`blastp`, `blastn`, `tblastn`, `blastx`, `tblastx`, `makeblastdb`):

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                       NCBI BLAST+ 2.17.0 Source Tree                        │
└──────────────────────────────────────┬──────────────────────────────────────┘
                                       │
                                       ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ Patch 1: ncbi_blast_2.17.0_patch (Base Infrastructure)                     │
│  • IPC progress tracking pipe (-progress_pipe, -progress_key)               │
│  • Subject length boundary filters (-subject_length_min/max)                │
│  • Middle line formatting field (midline) in tabular -outfmt 7              │
└──────────────────────────────────────┬──────────────────────────────────────┘
                                       │
                                       ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ Patch 2: ncbi_blast_2.17.0_filter_patch (Post-Traceback HSP Filters)        │
│  • 10 CLI parameters for per-HSP identity and coverage cutoffs              │
│  • min/max align identity, query identity, subject identity                 │
│  • min/max query coverage, subject coverage                                 │
│  • Injected after traceback convergence; works across all 5 BLAST programs  │
└──────────────────────────────────────┬──────────────────────────────────────┘
                                       │
                                       ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ Patch 3: ncbi_blast_2.17.0_cdr_patch (Antibody / TCR CDR Search Engine)     │
│  • Multi-CDR simultaneous search (up to 16 query slots) in blastp & blastn  │
│  • Exact mismatch verification with uniform IUPAC score matrix              │
│  • Dynamic CDR chaining with distance rules (-dist12, -dist23, -any_order)  │
│  • Native patent resolution via memory-mapped CDRPAT01 binary index         │
│  • 32-column compact streaming output (-cdr_outfmt compact)                 │
└──────────────────────────────────────┬──────────────────────────────────────┘
                                       │
                                       ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ ClickHouse Columnar Analytics & Storage Engine                              │
│  • High-throughput batch ingestion (10k rows/block)                         │
│  • Real-time 7-sector Venn Diagram bitmask breakdown (O(1))                 │
│  • Patent ranking prioritized by 100% identity and multi-CDR combinations   │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Repository Structure

```text
.
├── README.md                            # Main project overview and documentation
├── README_PATCHED_BLAST.md              # Detailed testing & filter verification report
├── build_patched_blast.sh               # Automated one-step download, patch, and build script
│
├── ncbi_blast_2.17.0_patch              # Patch 1: Progress pipe, subject limits, midline
├── ncbi_blast_2.17.0_filter_patch       # Patch 2: 10 post-traceback HSP filters
├── ncbi_blast_2.17.0_cdr_patch          # Patch 3: Unified Antibody CDR Search & Patent Engine
│
├── patidx/                              # Patent Index (CDRPAT01) tools & specs
│   ├── patidx.h                         # C++ header for binary index layout
│   ├── convert_nucl_patidx.cpp          # Fast C++ converter for prefix-sum .patidx index
│   └── build_wipo_blastdb.py            # Deduplication & BLAST DB builder script
│
├── clickhouse/                          # ClickHouse analytical schema & documentation
│   ├── init.sql                         # Production DDL for searches, hits, and metadata
│   └── clickhouse_storage_spec.md       # Full 32-column spec, bitmasks, and query recipes
│
└── test_data/                           # Verified test queries
    ├── test_hcdrs.fa                    # Heavy chain CDR1 + CDR2 + CDR3 (Trastuzumab)
    ├── test_hcdrs_plus_chain.fa         # 3 CDRs + full Heavy Chain sequence
    ├── test_lcdrs.fa                    # Light chain CDR1 + CDR2 + CDR3
    ├── test_hcdr3_only.fa               # Isolated HCDR3 loop
    └── test_nucl_primers.fa             # Nucleotide CDR primer queries
```

---

## 1. Quickstart & Build Instructions

### Prerequisites
- Linux or macOS
- GCC & G++ (support for C++17 / C99)
- Make, Patch, `curl` or `wget`

### Automated Build (All 3 Patches)
Run the automated build script:

```bash
./build_patched_blast.sh --jobs 8
```

The script will automatically:
1. Download `ncbi-blast-2.17.0+-src.tar.gz` from the NCBI official repository.
2. Unpack and sequentially apply `ncbi_blast_2.17.0_patch`, `ncbi_blast_2.17.0_filter_patch`, and `ncbi_blast_2.17.0_cdr_patch`.
3. Configure the minimal required projects (`app/blast`, `app/blastdb`) with `--with-mt --without-debug --without-gui --without-vdb`.
4. Compile release binaries (`ReleaseMT/bin/blastp`, `blastn`, `makeblastdb`, etc.).

### Build Flags
- `--dest <dir>`: Custom build destination directory (default: `./blast-build`).
- `--jobs <N>`: Parallel compilation threads (default: `4`).
- `--skip-download`: Use existing local tarball if present.
- `--no-patch2`: Omit HSP identity/coverage filters.
- `--no-cdr`: Omit antibody CDR engine.

---

## 2. Command-Line Parameters

### A. Post-Traceback HSP Filtering Parameters (Patch 2)
Filtering cutoffs are applied after traceback is complete, eliminating low-identity or low-coverage hits before formatting:

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-min_align_identity` | Real 0..100 | 0 | Minimum % identity across aligned HSP segment (`num_ident / align_len * 100`) |
| `-max_align_identity` | Real 0..100 | 0 | Maximum % identity across aligned HSP segment (useful for excluding self-matches) |
| `-min_query_identity` | Real 0..100 | 0 | Minimum % identity relative to full query length |
| `-max_query_identity` | Real 0..100 | 0 | Maximum % identity relative to full query length |
| `-min_subject_identity`| Real 0..100 | 0 | Minimum % identity relative to subject length |
| `-max_subject_identity`| Real 0..100 | 0 | Maximum % identity relative to subject length |
| `-min_query_coverage` | Real 0..100 | 0 | Minimum % query length covered by alignment |
| `-max_query_coverage` | Real 0..100 | 0 | Maximum % query length covered by alignment |
| `-min_subject_coverage`| Real 0..100 | 0 | Minimum % subject length covered by alignment |
| `-max_subject_coverage`| Real 0..100 | 0 | Maximum % subject length covered by alignment |
| `-subject_length_min` | Integer | 0 | Reject subjects with sequence length below cutoff |
| `-subject_length_max` | Integer | 0 | Reject subjects with sequence length above cutoff |

### B. Antibody / CDR Search Parameters (Patch 3)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-cdr_mode` | Flag | Off | Activates multi-query CDR and antibody alignment engine |
| `-cdr_outfmt` | String / Int | `full` | Output format: `blast` (0), `full` (1), `tsv` (2), `json` (3), `compact` / `clickhouse` (4) |
| `-cdr_chain_min_ident` | Real 0..100 | 70.0 | Minimum % identity required for the optional 4th full Chain / V-Region query |
| `-cdr_max_mismatches` | String | `2,2,2` | Comma-separated max mismatches per CDR slot (e.g. `1,2,1` or single integer `2`) |
| `-cdr_dist` | String | `10:25,30:50` | Allowed inter-CDR framework distances in aa/nt (`10:25,30:50` for protein; `30:75,90:150` for nucleotide) |
| `-cdr_any_order` | Flag | Off | Allow CDR1, CDR2, CDR3 to occur in arbitrary order along the subject sequence |
| `-append` (`makeblastdb`) | Flag | Off | Append new sequences to an existing BLAST database in-place |

---

## 3. Usage Examples

### Standard Patent BLAST with Strict Quality Thresholds
Search a protein query against patent database, retaining only hits with $\ge 80\%$ identity and $\ge 70\%$ query coverage:
```bash
blastp \
  -query test_data/test_hcdrs_plus_chain.fa \
  -db /data/wipo_prot \
  -evalue 1e6 \
  -max_target_seqs 100000 \
  -min_align_identity 80 \
  -min_query_coverage 70 \
  -outfmt "7 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
  -out results_standard.tsv
```

### Multi-CDR Antibody Search with Compact Streaming
Search 3 CDRs (HCDR1, HCDR2, HCDR3) simultaneously, enforcing chaining distance constraints and outputting the 32-column compact format:
```bash
blastp \
  -query test_data/test_hcdrs.fa \
  -db /data/wipo_prot \
  -cdr_mode \
  -cdr_outfmt compact \
  -out results_cdr_compact.tsv
```

---

## 4. Patent Index Format (`CDRPAT01`)

To eliminate the latency of external database lookups, the patched engine reads a binary memory-mapped index file (`<db>.patidx`) directly in C++:

- **Header (48 bytes)**: Magic bytes `CDRPAT01`, format version, patent count ($P$), OID count ($N$), total occurrences ($M$), and string pool size ($S$).
- **Patent Metadata Array**: `SPatIdxPatentMeta[P]` (8 bytes per patent: string pool offset and patent flag).
- **OID Offsets Array**: `uint64_t[N + 1]` prefix sums enabling $O(1)$ lookup from any BLAST OID to its occurrence slice.
- **Occurrences Array**: `SPatIdxOccurrence[M]` (8 bytes per entry: `pat_idx` + `seq_id_no`).
- **String Pool**: Null-terminated ASCII patent publication numbers (e.g. `WO2024123456\0`).

See [`patidx/patidx.h`](patidx/patidx.h) and [`patidx/convert_nucl_patidx.cpp`](patidx/convert_nucl_patidx.cpp) for implementation details.

---

## 5. ClickHouse Result Storage

The output of `-cdr_outfmt compact` streams directly into ClickHouse.

```bash
# Direct batch ingestion into ClickHouse (prepending search_id UUID and skipping '#' header lines)
SEARCH_ID="6b4766a3-fc88-4200-9833-f5e51ca2b259"
grep -v '^#' results_cdr_compact.tsv | awk -v sid="$SEARCH_ID" 'BEGIN{OFS="\t"} {print sid, $0}' | \
  clickhouse-client --query "INSERT INTO wipo.search_hits (search_id, chain_type, oid, subj_len, is_patent, pat_idx, doc_id, seq_id_no, seq_venn_mask, seq_exact_mask, pat_venn_mask, pat_exact_mask, cdr_tier, pat_tier, pat_match_type, seq_combo_score, pat_combo_score, total_cdr_mm, cdr1_start, cdr1_end, cdr1_mm, cdr2_start, cdr2_end, cdr2_mm, cdr3_start, cdr3_end, cdr3_mm, dist_12, dist_23, chain_pct_ident, chain_mm, chain_start, chain_end) FORMAT TabSeparated"
```

For table definitions, bitmask specifications, and analytical query recipes, see:
- [`clickhouse/init.sql`](clickhouse/init.sql)
- [`clickhouse/clickhouse_storage_spec.md`](clickhouse/clickhouse_storage_spec.md)

---

## 6. License & Attribution

NCBI BLAST+ is in the public domain and distributed by NCBI/NLM/NIH. The custom patches, CDR algorithms, binary indexing formats, and ClickHouse schemas in this repository are developed for patent sequence intelligence.
