#!/usr/bin/env python3
"""
Incremental WIPO ST.26 Protein Extractor, 256-Shard Sequence Deduplicator,
and Binary Patent Index (.patidx) Builder for BLAST+ CDR Patent Search.

Pipeline:
  Stage 1: Parse new WIPO ST.26 .zip archives, deduplicate sequences within the batch in RAM
           using BLAKE2b-128, and partition unique batch hashes into 256 buckets (00..ff).
  Stage 2: Iterate b = 0..255, load only shard_{b:02x}.bin from disk, resolve existing oids
           or allocate new sequential oids, append new unique sequences to FASTA, and append
           (oid, pat_idx, seq_id_no) to occurrences.bin.
  Stage 3: Build mmap-ready binary Patent Index (<db>.patidx) and run makeblastdb.
"""
import argparse
import hashlib
import json
import re
import shutil
import struct
import subprocess
import sys
import time
import unicodedata
import xml.etree.ElementTree as ET
import zipfile
from array import array
from pathlib import Path

WIPO_ROOT = Path("/Users/vladimir/Projects/AirFlow/data/wipo")
OUT_DIR = Path("/Users/vladimir/Projects/BLAST+CDR/wipo_2026_db")
STATE_DIR = OUT_DIR / "index_state"
SHARDS_DIR = STATE_DIR / "shards"
FASTA_PATH = OUT_DIR / "wipo_2026_prot.fa"
DB_PREFIX = OUT_DIR / "wipo_2026_prot"
PATIDX_PATH = OUT_DIR / "wipo_2026_prot.patidx"
MAKEBLASTDB = Path("/Users/vladimir/Projects/BLAST+CDR/ncbi-blast-2.17.0+-src/c++/build_min/bin/makeblastdb")

CLEAN_WS = re.compile(r"\s+")
VALID_AA = re.compile(r"[^ACDEFGHIKLMNPQRSTVWYBZXUOJ*]")

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "ASX": "B", "GLX": "Z", "XLE": "J", "SEC": "U", "PYL": "O",
    "XAA": "X", "UNK": "X", "TER": "*", "STP": "*",
}
AA3_RE = re.compile(r"\b([A-Za-z]{3})\b")
ST25_TAG_RE = re.compile(r"^\s*<(\d{3})>\s*(.*?)\s*$")
ST25_INT_RE = re.compile(r"\d+")


def parse_st25_proteins(text: str) -> list[tuple[int, str]]:
    """
    Parse WIPO ST.25 sequence listing (.txt / .app / .seq) and extract (seq_id_no, aa_seq)
    for all <212> PRT entries.
    """
    seqs: list[tuple[int, str]] = []
    cur_id = 0
    cur_type = ""
    in_400 = False
    buf_400: list[str] = []

    def flush():
        nonlocal in_400, buf_400
        if in_400 and cur_type == "PRT" and buf_400:
            raw_block = " ".join(buf_400)
            if raw_block.strip() != "000":
                tokens = AA3_RE.findall(raw_block)
                aa_chars = [AA3_TO_1[t.upper()] for t in tokens if t.upper() in AA3_TO_1]
                if aa_chars:
                    seqs.append((cur_id, "".join(aa_chars)))
        in_400 = False
        buf_400 = []

    for line in text.splitlines():
        if "<" in line:
            m = ST25_TAG_RE.match(line)
            if m:
                code, val = m.group(1), m.group(2)
                if code == "210":
                    flush()
                    m_id = ST25_INT_RE.search(val)
                    cur_id = int(m_id.group(0)) if m_id else 0
                    cur_type = ""
                elif code == "212":
                    cur_type = val.upper().strip()
                elif code == "400":
                    in_400 = True
                    buf_400 = []
                else:
                    if in_400:
                        flush()
                continue
        if in_400 and cur_type == "PRT":
            buf_400.append(line)
    flush()
    return seqs

SHARD_REC_STRUCT = struct.Struct("<16sII")  # hash128 (16B), oid (4B), seq_len (4B) = 24 bytes
OCC_LOG_STRUCT = struct.Struct("<III")      # oid (4B), pat_idx (4B), seq_id_no (4B) = 12 bytes
PATIDX_HDR_STRUCT = struct.Struct("<8sIIIIQQQ")  # 48 bytes
PATIDX_PAT_STRUCT = struct.Struct("<II")         # doc_id_offset (4B), is_patent (4B) = 8 bytes
PATIDX_OCC_STRUCT = struct.Struct("<II")         # pat_idx (4B), seq_id_no (4B) = 8 bytes


def load_state():
    SHARDS_DIR.mkdir(parents=True, exist_ok=True)
    ingested_file = STATE_DIR / "ingested_zips.txt"
    patents_file = STATE_DIR / "patents.jsonl"
    meta_file = STATE_DIR / "meta.json"

    ingested_zips: set[str] = set()
    if ingested_file.exists():
        ingested_zips = set(ingested_file.read_text(encoding="utf-8").splitlines())

    patents: list[dict] = []
    doc_to_pat_idx: dict[str, int] = {}
    if patents_file.exists():
        with open(patents_file, "r", encoding="utf-8") as f:
            for line in f:
                if line.strip():
                    raw_rec = json.loads(line)
                    rec = {
                        "pat_idx": int(raw_rec["pat_idx"]),
                        "doc_id": str(raw_rec.get("doc_id", "")),
                        "is_patent": int(raw_rec.get("is_patent", 1)),
                    }
                    doc_to_pat_idx[rec["doc_id"]] = rec["pat_idx"]
                    patents.append(rec)

    next_oid = 0
    total_occ = 0
    if meta_file.exists():
        meta = json.loads(meta_file.read_text(encoding="utf-8"))
        next_oid = meta.get("next_oid", 0)
        total_occ = meta.get("total_occ", 0)

    return ingested_zips, patents, doc_to_pat_idx, next_oid, total_occ


def save_meta(next_oid: int, total_occ: int):
    meta_file = STATE_DIR / "meta.json"
    meta_file.write_text(
        json.dumps({"next_oid": next_oid, "total_occ": total_occ}, indent=2),
        encoding="utf-8",
    )


def ingest_batch(zip_files: list[Path],
                 ingested_zips: set[str],
                 patents: list[dict],
                 doc_to_pat_idx: dict[str, int],
                 next_oid: int,
                 total_occ: int) -> tuple[int, int, int, int]:
    """
    Stage 1 (In-RAM Intra-Batch Deduplication & 256-Bucket Partitioning) +
    Stage 2 (Sequential 256-Shard Lookup & Incremental Update).
    """
    t0 = time.time()
    new_zips = [z for z in zip_files if str(z) not in ingested_zips]
    if not new_zips:
        print("No new .zip archives to ingest (all already in ingested_zips.txt).")
        return 0, 0, next_oid, total_occ

    print(f"[Stage 1] Parsing {len(new_zips)} new WIPO ST.26 .zip archives & deduplicating batch in RAM...")
    # Stage 1 structures in RAM for the incoming batch:
    # batch_seqs[h128] = clean_seq_str
    # batch_occs[h128] = list of (pat_idx, seq_id_no)
    # buckets[0..255] = list of unique h128 in this batch with h128[0] == b
    batch_seqs: dict[bytes, str] = {}
    batch_occs: dict[bytes, list[tuple[int, int]]] = {}
    buckets: list[list[bytes]] = [[] for _ in range(256)]

    new_patents: list[dict] = []
    newly_ingested_zips: list[str] = []
    raw_aa_count = 0
    skipped_short = 0

    def add_protein_seq(pat_idx: int, seq_id_no: int, raw_seq: str):
        nonlocal raw_aa_count, skipped_short
        clean_seq = VALID_AA.sub("", raw_seq.strip().upper())
        if len(clean_seq) < 4:
            skipped_short += 1
            return
        raw_aa_count += 1
        h128 = hashlib.blake2b(clean_seq.encode("ascii"), digest_size=16).digest()
        if h128 not in batch_seqs:
            batch_seqs[h128] = clean_seq
            batch_occs[h128] = [(pat_idx, seq_id_no)]
            buckets[h128[0]].append(h128)
        else:
            batch_occs[h128].append((pat_idx, seq_id_no))

    for zpath in new_zips:
        wo_dir = zpath.parent.name

        try:
            with zipfile.ZipFile(zpath, "r") as zf:
                for member in zf.namelist():
                    m_low = member.lower()
                    if m_low.endswith((".tif", ".tiff", ".pdf", ".png", ".jpg", ".jpeg", ".lst")):
                        continue
                    if not m_low.endswith((".xml", ".txt", ".app", ".seq")):
                        continue

                    doc_id = wo_dir
                    if doc_id not in doc_to_pat_idx:
                        pat_idx = len(patents)
                        pat_rec = {
                            "pat_idx": pat_idx,
                            "doc_id": doc_id,
                            "is_patent": 1,
                        }
                        doc_to_pat_idx[doc_id] = pat_idx
                        patents.append(pat_rec)
                        new_patents.append(pat_rec)
                    else:
                        pat_idx = doc_to_pat_idx[doc_id]

                    info = zf.getinfo(member)
                    if m_low.endswith(".xml") and info.file_size > 100_000_000:
                        with zf.open(member, "r") as xf:
                            for _event, elem in ET.iterparse(xf, events=("end",)):
                                if elem.tag == "SequenceData":
                                    seq_id_str = elem.attrib.get("sequenceIDNumber", "0")
                                    try:
                                        seq_id_no = int(seq_id_str)
                                    except ValueError:
                                        seq_id_no = 0
                                    moltype = (elem.findtext(".//INSDSeq_moltype") or "").strip().upper()
                                    if moltype == "AA":
                                        raw_seq = elem.findtext(".//INSDSeq_sequence") or ""
                                        add_protein_seq(pat_idx, seq_id_no, raw_seq)
                                    elem.clear()
                        continue

                    raw_bytes = zf.read(member)
                    if m_low.endswith(".xml") or raw_bytes.lstrip()[:5].startswith((b"<?xml", b"<!DOC", b"<ST26")):
                        root = ET.fromstring(raw_bytes)
                        for seq_elem in root.iter("SequenceData"):
                            seq_id_str = seq_elem.attrib.get("sequenceIDNumber", "0")
                            try:
                                seq_id_no = int(seq_id_str)
                            except ValueError:
                                seq_id_no = 0
                            moltype = (seq_elem.findtext(".//INSDSeq_moltype") or "").strip().upper()
                            if moltype != "AA":
                                continue
                            raw_seq = seq_elem.findtext(".//INSDSeq_sequence") or ""
                            add_protein_seq(pat_idx, seq_id_no, raw_seq)
                    else:
                        text = raw_bytes.decode("latin-1", errors="ignore")
                        if "<ST26SequenceListing" in text[:1000]:
                            root = ET.fromstring(raw_bytes)
                            for seq_elem in root.iter("SequenceData"):
                                seq_id_str = seq_elem.attrib.get("sequenceIDNumber", "0")
                                try:
                                    seq_id_no = int(seq_id_str)
                                except ValueError:
                                    seq_id_no = 0
                                moltype = (seq_elem.findtext(".//INSDSeq_moltype") or "").strip().upper()
                                if moltype != "AA":
                                    continue
                                raw_seq = seq_elem.findtext(".//INSDSeq_sequence") or ""
                                add_protein_seq(pat_idx, seq_id_no, raw_seq)
                        else:
                            for seq_id_no, prt_seq in parse_st25_proteins(text):
                                add_protein_seq(pat_idx, seq_id_no, prt_seq)
            newly_ingested_zips.append(str(zpath))
        except Exception as exc:
            print(f"[WARN] Failed to parse {zpath}: {exc}", file=sys.stderr)
            newly_ingested_zips.append(str(zpath))

    t1 = time.time()
    print(
        f"[Stage 1 Done in {t1 - t0:.2f}s] Raw AA sequences in batch: {raw_aa_count:,} -> "
        f"Unique in batch: {len(batch_seqs):,} "
        f"(intra-batch dupes removed: {raw_aa_count - len(batch_seqs):,})"
    )

    # Stage 2: Compare bucket by bucket (00..ff) against on-disk shards
    print("[Stage 2] Comparing 256 buckets against on-disk hash shards (shard_00.bin .. shard_ff.bin)...")
    new_unique_seqs: list[tuple[int, str]] = []
    occ_log_path = STATE_DIR / "occurrences.bin"
    reused_existing = 0

    with open(occ_log_path, "ab") as occ_f:
        for b in range(256):
            b_hashes = buckets[b]
            if not b_hashes:
                continue

            shard_path = SHARDS_DIR / f"shard_{b:02x}.bin"
            shard_map: dict[bytes, int] = {}
            if shard_path.exists():
                raw_shard = shard_path.read_bytes()
                for h_key, oid_val, _slen in SHARD_REC_STRUCT.iter_unpack(raw_shard):
                    shard_map[h_key] = oid_val

            new_shard_bytes = bytearray()
            occ_chunk = bytearray()

            for h128 in b_hashes:
                existing_oid = shard_map.get(h128)
                if existing_oid is not None:
                    oid = existing_oid
                    reused_existing += 1
                else:
                    oid = next_oid
                    next_oid += 1
                    seq_str = batch_seqs[h128]
                    shard_map[h128] = oid
                    new_shard_bytes.extend(SHARD_REC_STRUCT.pack(h128, oid, len(seq_str)))
                    new_unique_seqs.append((oid, seq_str))

                for pat_idx, seq_id_no in batch_occs[h128]:
                    occ_chunk.extend(OCC_LOG_STRUCT.pack(oid, pat_idx, seq_id_no))
                    total_occ += 1

            if new_shard_bytes:
                with open(shard_path, "ab") as sf:
                    sf.write(new_shard_bytes)
            if occ_chunk:
                occ_f.write(occ_chunk)

    # Sort newly assigned unique sequences by oid and append to FASTA
    new_unique_seqs.sort(key=lambda x: x[0])
    FASTA_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(FASTA_PATH, "a", encoding="utf-8") as fa_f:
        for oid, seq_str in new_unique_seqs:
            fa_f.write(f">USEQ_{oid} len={len(seq_str)}\n{seq_str}\n")

    batch_new_fa = STATE_DIR / "new_batch.fa"
    if new_unique_seqs:
        with open(batch_new_fa, "w", encoding="utf-8") as bfa:
            for oid, seq_str in new_unique_seqs:
                bfa.write(f">USEQ_{oid} len={len(seq_str)}\n{seq_str}\n")
    elif batch_new_fa.exists():
        batch_new_fa.unlink()

    # Persist updated patent metadata, ingested_zips, and meta counters
    with open(STATE_DIR / "patents.jsonl", "w", encoding="utf-8") as pf:
        for p in patents:
            pf.write(json.dumps(p) + "\n")

    with open(STATE_DIR / "ingested_zips.txt", "a", encoding="utf-8") as zf_out:
        for zp in newly_ingested_zips:
            zf_out.write(zp + "\n")

    save_meta(next_oid, total_occ)
    t2 = time.time()
    print(
        f"[Stage 2 Done in {t2 - t1:.2f}s] Matched existing DB oids: {reused_existing:,}, "
        f"Added new unique sequences: {len(new_unique_seqs):,} "
        f"(Total DB unique oids: {next_oid:,}, Total patent occurrences: {total_occ:,})"
    )
    return len(newly_ingested_zips), raw_aa_count, next_oid, total_occ, len(new_unique_seqs)


def build_patidx_file(patents: list[dict], num_oids: int, total_occ: int):
    """
    Stage 3A: Compile occurrences.bin and patents list into the flat 8-byte-aligned
    mmap-ready binary file wipo_2026_prot.patidx using O(N) counting sort by oid.
    """
    t0 = time.time()
    print(f"[Stage 3A] Building mmap-ready binary Patent Index {PATIDX_PATH} ...")

    # 1. Build StringPool and PatentMetaEntry table (24 bytes per patent)
    str_pool = bytearray(b"\x00")  # offset 0 is empty string ""
    str_cache: dict[str, int] = {"": 0}

    def intern_str(s: str) -> int:
        if not s:
            return 0
        off = str_cache.get(s)
        if off is not None:
            return off
        off = len(str_pool)
        str_pool.extend(s.encode("ascii", errors="ignore") + b"\x00")
        str_cache[s] = off
        return off

    pat_entries_bytes = bytearray()
    for p in patents:
        doc_off = intern_str(p.get("doc_id", ""))
        pat_entries_bytes.extend(
            PATIDX_PAT_STRUCT.pack(
                doc_off,
                int(p.get("is_patent", 1)),
            )
        )

    # 2. Group occurrences.bin by oid using O(N) counting sort
    occ_log_path = STATE_DIR / "occurrences.bin"
    raw_occ = occ_log_path.read_bytes()
    actual_occ = len(raw_occ) // OCC_LOG_STRUCT.size
    assert actual_occ == total_occ, f"Expected {total_occ} occurrences, found {actual_occ}"

    counts = array("I", [0]) * (num_oids + 1)
    for oid, _pat_idx, _seq_id_no in OCC_LOG_STRUCT.iter_unpack(raw_occ):
        counts[oid] += 1

    # Compute prefix sums for OidOffsets
    oid_offsets = array("Q", [0]) * (num_oids + 1)
    running = 0
    for i in range(num_oids):
        oid_offsets[i] = running
        running += counts[i]
    oid_offsets[num_oids] = running

    # Place 8-byte (pat_idx, seq_id_no) records into sorted_occs buffer
    write_pos = array("I", oid_offsets[:-1])
    sorted_occs = bytearray(actual_occ * PATIDX_OCC_STRUCT.size)
    for oid, pat_idx, seq_id_no in OCC_LOG_STRUCT.iter_unpack(raw_occ):
        pos = write_pos[oid]
        write_pos[oid] = pos + 1
        PATIDX_OCC_STRUCT.pack_into(sorted_occs, pos * 8, pat_idx, seq_id_no)

    # 3. Write header + PatentTable + OidOffsets + Occurrences + StringPool
    hdr_bytes = PATIDX_HDR_STRUCT.pack(
        b"CDRPAT01",
        1,                  # version
        len(patents),       # num_patents
        num_oids,           # num_oids
        0,                  # reserved
        actual_occ,         # total_occurrences
        len(str_pool),      # pat_strings_size
        0,                  # reserved2
    )

    with open(PATIDX_PATH, "wb") as out_f:
        out_f.write(hdr_bytes)
        out_f.write(pat_entries_bytes)
        if sys.byteorder != "little":
            oid_offsets.byteswap()
        out_f.write(oid_offsets.tobytes())
        out_f.write(sorted_occs)
        out_f.write(str_pool)

    t1 = time.time()
    size_mb = PATIDX_PATH.stat().st_size / (1024 * 1024)
    print(
        f"[Stage 3A Done in {t1 - t0:.2f}s] Wrote {PATIDX_PATH} ({size_mb:.2f} MB): "
        f"{len(patents):,} patents, {num_oids:,} unique oids, {actual_occ:,} occurrences."
    )


def main():
    parser = argparse.ArgumentParser(description="Incremental WIPO BLAST DB + .patidx builder")
    parser.add_argument(
        "--years",
        type=str,
        default="",
        help="Comma-separated list of years (e.g. 2020,2019,2018) or empty for all years",
    )
    parser.add_argument(
        "--weeks",
        type=str,
        default="",
        help="Comma-separated list of week folders (e.g. 0813,0820) or empty for all weeks",
    )
    parser.add_argument(
        "--reset",
        action="store_true",
        help="Reset existing deduplication state and rebuild from scratch",
    )
    parser.add_argument(
        "--rebuild-patidx",
        action="store_true",
        help="Rebuild wipo_2026_prot.patidx and clean patents.jsonl without re-parsing zips",
    )
    args = parser.parse_args()

    if args.reset:
        print("Resetting index state and FASTA...")
        if STATE_DIR.exists():
            shutil.rmtree(STATE_DIR)
        if FASTA_PATH.exists():
            FASTA_PATH.unlink()
        if PATIDX_PATH.exists():
            PATIDX_PATH.unlink()

    ingested_zips, patents, doc_to_pat_idx, next_oid, total_occ = load_state()

    if args.rebuild_patidx:
        with open(STATE_DIR / "patents.jsonl", "w", encoding="utf-8") as pf:
            for p in patents:
                pf.write(json.dumps(p) + "\n")
        build_patidx_file(patents, next_oid, total_occ)
        return

    if args.years:
        year_set = [y.strip() for y in args.years.split(",") if y.strip()]
        zip_files = []
        for y in sorted(year_set):
            ydir = WIPO_ROOT / y
            if ydir.exists():
                zip_files.extend(sorted(ydir.rglob("*.zip")))
    elif args.weeks:
        week_set = [w.strip() for w in args.weeks.split(",") if w.strip()]
        zip_files = []
        for w in sorted(week_set):
            wdir = WIPO_ROOT / w
            if wdir.exists():
                zip_files.extend(sorted(wdir.rglob("*.zip")))
    else:
        zip_files = sorted(WIPO_ROOT.rglob("*.zip"))

    new_zips_cnt, _raw_cnt, next_oid, total_occ, new_unique_cnt = ingest_batch(
        zip_files, ingested_zips, patents, doc_to_pat_idx, next_oid, total_occ
    )

    if new_zips_cnt > 0 or not PATIDX_PATH.exists():
        build_patidx_file(patents, next_oid, total_occ)
        pin_file = Path(str(DB_PREFIX) + ".pin")
        batch_new_fa = STATE_DIR / "new_batch.fa"
        if pin_file.exists() and batch_new_fa.exists() and new_unique_cnt > 0:
            cmd = [
                str(MAKEBLASTDB),
                "-in", str(batch_new_fa),
                "-dbtype", "prot",
                "-title", f"WIPO_PCT_Proteins_Unified_2010_2026 ({next_oid} unique seqs, {len(patents)} patents)",
                "-out", str(DB_PREFIX),
                "-append",
            ]
            print("[Stage 3B] Running makeblastdb (in-place append mode):", " ".join(cmd))
        else:
            cmd = [
                str(MAKEBLASTDB),
                "-in", str(FASTA_PATH),
                "-dbtype", "prot",
                "-title", f"WIPO_PCT_Proteins_Unified_2010_2026 ({next_oid} unique seqs, {len(patents)} patents)",
                "-out", str(DB_PREFIX),
            ]
            print("[Stage 3B] Running makeblastdb (full build mode):", " ".join(cmd))
        t_db = time.time()
        subprocess.run(cmd, check=True)
        print(f"[Stage 3B Done in {time.time() - t_db:.2f}s] Updated BLAST DB at {DB_PREFIX}")
        kmer_idx = Path(str(DB_PREFIX) + ".kmeridx")
        if kmer_idx.exists():
            kmer_idx.unlink()
            print(f"Removed outdated {kmer_idx} (will auto-rebuild on next search).")


if __name__ == "__main__":
    main()
