-- =============================================================================
-- WIPO Patent Sequence & CDR Search Schema for ClickHouse
-- =============================================================================

CREATE DATABASE IF NOT EXISTS wipo;

-- 1. Table for registered searches (Standard BLAST and Antibody CDR)
CREATE TABLE IF NOT EXISTS wipo.searches (
    search_id UUID,
    name String,
    search_type LowCardinality(String) DEFAULT 'antibody',  -- 'standard' or 'antibody'
    chain_type LowCardinality(String),                      -- 'H', 'L', 'H+L', or 'NT'
    query_chain String,
    cdr1_seq String DEFAULT '',
    cdr2_seq String DEFAULT '',
    cdr3_seq String DEFAULT '',
    max_mm_cdr1 UInt8 DEFAULT 0,
    max_mm_cdr2 UInt8 DEFAULT 0,
    max_mm_cdr3 UInt8 DEFAULT 0,
    min_chain_ident Float32 DEFAULT 75.0,
    created_at DateTime DEFAULT now(),
    status Enum8('PENDING'=1, 'RUNNING'=2, 'COMPLETED'=3, 'FAILED'=4),
    execution_time_sec Float32 DEFAULT 0.0,
    total_hits UInt32 DEFAULT 0,
    unique_seqs UInt32 DEFAULT 0,
    unique_patents UInt32 DEFAULT 0,
    error_message String DEFAULT '',
    db_type LowCardinality(String) DEFAULT 'prot',          -- 'prot' or 'nucl'
    search_mode LowCardinality(String) DEFAULT 'cdr_prot',  -- 'blastp', 'blastn', 'cdr_prot', etc.
    evalue Float64 DEFAULT 1000000.0,
    max_target_seqs UInt32 DEFAULT 100000,
    min_query_cov Float32 DEFAULT 75.0,
    min_subj_len Int32 DEFAULT -1,
    max_subj_len Int32 DEFAULT -1
) ENGINE = ReplacingMergeTree()
PRIMARY KEY (search_id);

-- 2. Table for search hits (denormalized, optimized for ultra-fast aggregation)
CREATE TABLE IF NOT EXISTS wipo.search_hits (
    search_id UUID,
    chain_type LowCardinality(String),
    oid UInt32,
    subj_len UInt32,
    is_patent UInt8,
    pat_idx UInt32,
    doc_id LowCardinality(String),
    seq_id_no UInt32,

    -- Venn classification & scoring
    seq_venn_mask UInt8,        -- Bitmask: 1=CDR1, 2=CDR2, 4=CDR3, 8=Chain
    seq_exact_mask UInt8,       -- Bitmask of exact 0-mismatch CDR matches
    pat_venn_mask UInt8,        -- Patent-level aggregated Venn bitmask
    pat_exact_mask UInt8,       -- Patent-level aggregated exact bitmask
    cdr_tier UInt8,             -- 1: all 3 CDRs, 2: CDR3 + 1, 3: CDR3 only, 4: other
    pat_tier UInt8,
    pat_match_type UInt8,
    seq_combo_score UInt32,
    pat_combo_score UInt32,
    total_cdr_mm Int16,

    -- Individual CDR alignment coordinates & mismatches
    cdr1_start UInt32 DEFAULT 0,
    cdr1_end UInt32 DEFAULT 0,
    cdr1_mm Int8 DEFAULT -1,
    cdr2_start UInt32 DEFAULT 0,
    cdr2_end UInt32 DEFAULT 0,
    cdr2_mm Int8 DEFAULT -1,
    cdr3_start UInt32 DEFAULT 0,
    cdr3_end UInt32 DEFAULT 0,
    cdr3_mm Int8 DEFAULT -1,
    dist_12 Int16 DEFAULT -1,   -- Framework distance between CDR1 and CDR2
    dist_23 Int16 DEFAULT -1,   -- Framework distance between CDR2 and CDR3

    -- Full chain / standard BLAST alignment coordinates & metrics
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

-- 3. Table for Patent Metadata (Applicants, Title, Publication Year & Date)
CREATE TABLE IF NOT EXISTS wipo.patents_meta (
    pat_idx UInt32,
    doc_id LowCardinality(String),
    pub_year UInt16,
    applicant LowCardinality(String),
    title String,
    pub_date Nullable(Date)
) ENGINE = ReplacingMergeTree()
PRIMARY KEY (doc_id);
