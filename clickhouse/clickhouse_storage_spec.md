# ClickHouse Result Storage & Analytical Specification

This document defines the result storage format, table schema, compression strategy, and analytical query recipes used for storing and querying NCBI BLAST+ results (both Standard Patent BLAST and Antibody CDR Search) in ClickHouse.

---

## 1. Architecture Overview

```
                      ┌──────────────────────────────────────┐
                      │  Patched blastp / blastn Executable  │
                      │  -cdr_mode -cdr_outfmt compact       │
                      └──────────────────┬───────────────────┘
                                         │
                                         │ High-throughput TSV stream
                                         │ (32 columns per hit)
                                         ▼
                      ┌──────────────────────────────────────┐
                      │    Batch Ingestion Buffer (10k rows) │
                      └──────────────────┬───────────────────┘
                                         │
                                         │ Native columnar block insert
                                         ▼
                      ┌──────────────────────────────────────┐
                      │      ClickHouse: wipo.search_hits    │
                      │      MergeTree Engine                │
                      │      ORDER BY (search_id, venn_mask) │
                      └──────────────────┬───────────────────┘
                                         │
                 ┌───────────────────────┴───────────────────────┐
                 ▼                                               ▼
   Instant Venn Aggregations (7 Sectors)         Ranked Patent Hits & Schematics
   O(1) - O(log N) through indexed masks         Best hit per patent with 100% hits first
```

Rather than storing heavyweight pairwise alignments or verbose text files, the engine outputs compact, fixed-width columnar rows. Patent identifiers are pre-resolved directly inside the C++ engine via memory-mapped `.patidx` files, completely eliminating external lookup overhead.

---

## 2. Compact Output Format Specification (32 Columns)

When `-cdr_outfmt compact` is active, the patched BLAST executable outputs tab-delimited records with 32 fields (Matching the header `#chain_type	oid	subj_len	...	chain_end`):

| Col # | Field Name | Type | Description |
|:---:|---|---|---|
| 1 | `chain_type` | String | Chain classification: `'H'` (Heavy), `'L'` (Light), `'NT'`, `'C'` |
| 2 | `oid` | UInt32 | Internal BLAST sequence ordinal ID in database |
| 3 | `subj_len` | UInt32 | Subject sequence length (amino acids or nucleotides) |
| 4 | `is_patent` | UInt8 | `1` if sequence originates from patent literature, `0` otherwise |
| 5 | `pat_idx` | UInt32 | Patent index within `.patidx` string pool |
| 6 | `doc_id` | String | Patent publication number (e.g. `WO2024123456`) |
| 7 | `seq_id_no` | UInt32 | Sequence identifier in patent (`<212> SEQ ID NO`) |
| 8 | `seq_venn_mask` | UInt8 | Hit-level Venn bitmask (1=CDR1, 2=CDR2, 4=CDR3, 8=Chain) |
| 9 | `seq_exact_mask`| UInt8 | Exact 0-mismatch match bitmask |
| 10 | `pat_venn_mask` | UInt8 | Patent-level cumulative Venn bitmask |
| 11 | `pat_exact_mask`| UInt8 | Patent-level cumulative exact bitmask |
| 12 | `cdr_tier` | UInt8 | CDR match tier (1: CDR1+2+3, 2: CDR3+1, 3: CDR3, 4: Other) |
| 13 | `pat_tier` | UInt8 | Patent-level best tier |
| 14 | `pat_match_type`| UInt8 | Canonical match type category |
| 15 | `seq_combo_score`| UInt32| Sequence combination score |
| 16 | `pat_combo_score`| UInt32| Cumulative patent combination score |
| 17 | `total_cdr_mm` | Int16 | Total mismatches across all aligned CDRs |
| 18 | `cdr1_start` | UInt32 | 1-based start coordinate of CDR1 on subject (0 if unaligned) |
| 19 | `cdr1_end` | UInt32 | 1-based end coordinate of CDR1 on subject (0 if unaligned) |
| 20 | `cdr1_mm` | Int8 | Number of mismatches in CDR1 (-1 if unaligned) |
| 21 | `cdr2_start` | UInt32 | 1-based start coordinate of CDR2 on subject (0 if unaligned) |
| 22 | `cdr2_end` | UInt32 | 1-based end coordinate of CDR2 on subject (0 if unaligned) |
| 23 | `cdr2_mm` | Int8 | Number of mismatches in CDR2 (-1 if unaligned) |
| 24 | `cdr3_start` | UInt32 | 1-based start coordinate of CDR3 on subject (0 if unaligned) |
| 25 | `cdr3_end` | UInt32 | 1-based end coordinate of CDR3 on subject (0 if unaligned) |
| 26 | `cdr3_mm` | Int8 | Number of mismatches in CDR3 (-1 if unaligned) |
| 27 | `dist_12` | Int16 | Framework amino acid distance between CDR1 end and CDR2 start |
| 28 | `dist_23` | Int16 | Framework amino acid distance between CDR2 end and CDR3 start |
| 29 | `chain_pct_ident`| Float32| Alignment % identity for chain / standard BLAST hit |
| 30 | `chain_mm` | Int16 | Total mismatch count in chain alignment |
| 31 | `chain_start` | UInt32 | 1-based start coordinate on subject sequence |
| 32 | `chain_end` | UInt32 | 1-based end coordinate on subject sequence |

When ingesting standard BLAST searches, 6 additional canonical alignment metrics are populated:
`qstart`, `qend`, `evalue`, `bitscore`, `qcov` (% query coverage), `scov` (% subject coverage).

---

## 3. Venn Diagram Bitmask Encoding

To support real-time 3-way Venn diagram visualization (CDR1, CDR2, CDR3) with zero recomputation, hits are encoded with a 4-bit bitmask:

```
Bit 0 (value 1): CDR1 matched within threshold
Bit 1 (value 2): CDR2 matched within threshold
Bit 2 (value 4): CDR3 matched within threshold
Bit 3 (value 8): Chain alignment matched within threshold
```

### The 7 Venn Sectors

| Bitmask Value | Binary | Matched Regions | Description |
|:---:|:---:|---|---|
| **7** | `0111` | CDR1 + CDR2 + CDR3 | **Triple match**: Complete antigen-binding site |
| **6** | `0110` | CDR2 + CDR3 | Dual match (CDR2 & CDR3) |
| **5** | `0101` | CDR1 + CDR3 | Dual match (CDR1 & CDR3) |
| **3** | `0011` | CDR1 + CDR2 | Dual match (CDR1 & CDR2) |
| **4** | `0100` | CDR3 only | Isolated hypervariable loop |
| **2** | `0010` | CDR2 only | Isolated CDR2 loop |
| **1** | `0001` | CDR1 only | Isolated CDR1 loop |

---

## 4. ClickHouse Schema

```sql
CREATE DATABASE IF NOT EXISTS wipo;

-- Table for search hits
CREATE TABLE IF NOT EXISTS wipo.search_hits (
    search_id UUID,
    chain_type LowCardinality(String),
    oid UInt32,
    subj_len UInt32,
    is_patent UInt8,
    pat_idx UInt32,
    doc_id LowCardinality(String),
    seq_id_no UInt32,

    seq_venn_mask UInt8,
    seq_exact_mask UInt8,
    pat_venn_mask UInt8,
    pat_exact_mask UInt8,
    cdr_tier UInt8,
    pat_tier UInt8,
    pat_match_type UInt8,
    seq_combo_score UInt32,
    pat_combo_score UInt32,
    total_cdr_mm Int16,

    cdr1_start UInt32 DEFAULT 0,
    cdr1_end UInt32 DEFAULT 0,
    cdr1_mm Int8 DEFAULT -1,
    cdr2_start UInt32 DEFAULT 0,
    cdr2_end UInt32 DEFAULT 0,
    cdr2_mm Int8 DEFAULT -1,
    cdr3_start UInt32 DEFAULT 0,
    cdr3_end UInt32 DEFAULT 0,
    cdr3_mm Int8 DEFAULT -1,
    dist_12 Int16 DEFAULT -1,
    dist_23 Int16 DEFAULT -1,

    chain_pct_ident Float32 DEFAULT 0.0,
    chain_mm Int16 DEFAULT 0,
    chain_start UInt32 DEFAULT 0,
    chain_end UInt32 DEFAULT 0,
    qstart UInt32 DEFAULT 0,
    qend UInt32 DEFAULT 0,
    evalue Float64 DEFAULT 0.0,
    bitscore Float32 DEFAULT 0.0,
    qcov Float32 DEFAULT 0.0,
    scov Float32 DEFAULT 0.0
) ENGINE = MergeTree()
ORDER BY (search_id, seq_venn_mask, chain_pct_ident, subj_len, doc_id);
```

### Primary Key & Ordering Rationales
- `search_id`: The partition/lead filter for every query.
- `seq_venn_mask`: Enables instant Venn breakdown queries without reading unneeded hit data.
- `chain_pct_ident`: Allows push-down filtering for minimum identity sliders (e.g. `chain_pct_ident >= 85.0`).
- `subj_len`: Accelerates sequence length range filters (`subj_len BETWEEN 100 AND 150`).
- `doc_id`: Facilitates patent-level `GROUP BY doc_id` aggregation.

---

## 5. Key Analytical Query Recipes

### A. 7-Sector Venn Breakdown (Milliseconds on Millions of Rows)
```sql
SELECT
    bitAnd(seq_venn_mask, 7) AS sector,
    count() AS hit_count,
    uniqExact(doc_id) AS patent_count
FROM wipo.search_hits
WHERE search_id = {sid:UUID}
GROUP BY sector
ORDER BY sector DESC;
```

### B. Top Patents Ranked by Best Hit (100% Matches First)
```sql
SELECT
    doc_id,
    max(chain_pct_ident) AS max_ident,
    max(qcov) AS max_qcov,
    min(evalue) AS min_evalue,
    max(bitscore) AS max_score,
    max(bitAnd(seq_venn_mask, 7)) AS max_venn,
    count() AS total_hits
FROM wipo.search_hits
WHERE search_id = {sid:UUID}
GROUP BY doc_id
ORDER BY
    max_ident DESC,
    max_venn DESC,
    max_score DESC
LIMIT 50;
```

### C. Identity Distribution Histogram (10-Bucket Steps)
```sql
SELECT
    intDiv(toUInt32(chain_pct_ident), 10) * 10 AS bucket,
    count() AS hit_count
FROM wipo.search_hits
WHERE search_id = {sid:UUID}
GROUP BY bucket
ORDER BY bucket ASC;
```

### D. Coordinate Retrieval for Schematic Alignment Visualizer
```sql
SELECT
    doc_id,
    seq_id_no,
    subj_len,
    chain_pct_ident,
    qstart,
    qend,
    chain_start,
    chain_end,
    cdr1_start, cdr1_end, cdr1_mm,
    cdr2_start, cdr2_end, cdr2_mm,
    cdr3_start, cdr3_end, cdr3_mm
FROM wipo.search_hits
WHERE search_id = {sid:UUID} AND doc_id = {target_doc:String}
ORDER BY chain_pct_ident DESC, bitscore DESC;
```
