#!/usr/bin/env bash
# =============================================================================
# build_patched_blast.sh
#
# Downloads NCBI BLAST+ 2.17.0 source, applies two custom patches, and
# compiles the BLAST suite (blastn, blastp, blastx, tblastn, tblastx,
# makeblastdb, blastdbcmd).
#
# Usage:
#   ./build_patched_blast.sh [OPTIONS]
#
# Options:
#   --dest DIR       Directory to unpack and build in  (default: ./blast-build)
#   --jobs N         Parallel make jobs                (default: 4)
#   --skip-download  Reuse existing tarball if present
#   --no-patch2      Apply only the base patch (skip filter params patch)
#   --help           Show this help
#
# Required files (expected alongside this script):
#   ncbi_blast_2.17.0_patch         Base patch (progress, subject length, midline)
#   ncbi_blast_2.17.0_filter_patch  Filter patch (min/max identity/coverage params)
#
# Requirements:
#   gcc, g++, make, wget or curl, patch
# =============================================================================

set -euo pipefail

# ──────────────────────────────────────────────────────────────────────────────
# Defaults
# ──────────────────────────────────────────────────────────────────────────────
DEST_DIR="./blast-build"
JOBS=4
SKIP_DOWNLOAD=0
NO_PATCH2=0

SOURCE_VERSION="2.17.0+"
SOURCE_TARBALL="ncbi-blast-${SOURCE_VERSION}-src.tar.gz"
SOURCE_URL="https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/${SOURCE_TARBALL}"
SOURCE_DIR="ncbi-blast-${SOURCE_VERSION}-src"

PATCH1="ncbi_blast_2.17.0_patch"
PATCH2="ncbi_blast_2.17.0_filter_patch"

# ──────────────────────────────────────────────────────────────────────────────
# Colours
# ──────────────────────────────────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; BOLD='\033[1m'; NC='\033[0m'

log_info()  { echo -e "${CYAN}[INFO]${NC}  $*"; }
log_ok()    { echo -e "${GREEN}[OK]${NC}    $*"; }
log_warn()  { echo -e "${YELLOW}[WARN]${NC}  $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $*" >&2; }
log_step()  { echo -e "\n${BOLD}══════════════════════════════════════════${NC}"; \
              echo -e "${BOLD} $*${NC}"; \
              echo -e "${BOLD}══════════════════════════════════════════${NC}"; }

# ──────────────────────────────────────────────────────────────────────────────
# Argument parsing
# ──────────────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  sed -n '2,30p' "$0" | grep -E "^#" | sed 's/^# \?//'
  exit 0
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dest)          DEST_DIR="$2"; shift 2 ;;
    --jobs)          JOBS="$2";     shift 2 ;;
    --skip-download) SKIP_DOWNLOAD=1; shift ;;
    --no-patch2)     NO_PATCH2=1;   shift ;;
    --help|-h)       usage ;;
    *) log_error "Unknown option: $1"; usage ;;
  esac
done

# ──────────────────────────────────────────────────────────────────────────────
# Preflight checks
# ──────────────────────────────────────────────────────────────────────────────
log_step "Preflight checks"

check_cmd() {
  if command -v "$1" &>/dev/null; then
    log_ok "$1 found: $(command -v "$1")"
  else
    log_error "Required command not found: $1"
    exit 1
  fi
}

check_cmd gcc
check_cmd g++
check_cmd make
check_cmd patch

# wget or curl
if command -v wget &>/dev/null; then
  DOWNLOADER="wget"
elif command -v curl &>/dev/null; then
  DOWNLOADER="curl"
else
  log_error "Neither wget nor curl found. Install one to proceed."
  exit 1
fi
log_ok "Downloader: $DOWNLOADER"

# Patch files
for pf in "$SCRIPT_DIR/$PATCH1"; do
  if [[ ! -f "$pf" ]]; then
    log_error "Patch file not found: $pf"
    log_error "Ensure patch files are in the same directory as this script."
    exit 1
  fi
  log_ok "Patch file found: $(basename "$pf")"
done

if [[ $NO_PATCH2 -eq 0 ]]; then
  if [[ ! -f "$SCRIPT_DIR/$PATCH2" ]]; then
    log_error "Patch file not found: $SCRIPT_DIR/$PATCH2"
    exit 1
  fi
  log_ok "Patch file found: $PATCH2"
fi

# Create destination directory
mkdir -p "$DEST_DIR"
DEST_DIR="$(cd "$DEST_DIR" && pwd)"
log_info "Build destination: $DEST_DIR"

# ──────────────────────────────────────────────────────────────────────────────
# Download
# ──────────────────────────────────────────────────────────────────────────────
log_step "Download source tarball"

TARBALL_PATH="$DEST_DIR/$SOURCE_TARBALL"

if [[ $SKIP_DOWNLOAD -eq 1 && -f "$TARBALL_PATH" ]]; then
  log_warn "--skip-download set and tarball exists, skipping download."
else
  log_info "Downloading $SOURCE_URL ..."
  if [[ "$DOWNLOADER" == "wget" ]]; then
    wget -q --show-progress -O "$TARBALL_PATH" "$SOURCE_URL"
  else
    curl -L --progress-bar -o "$TARBALL_PATH" "$SOURCE_URL"
  fi
  log_ok "Downloaded: $TARBALL_PATH ($(du -h "$TARBALL_PATH" | cut -f1))"
fi

# ──────────────────────────────────────────────────────────────────────────────
# Unpack
# ──────────────────────────────────────────────────────────────────────────────
log_step "Unpack source"

SRC_PATH="$DEST_DIR/$SOURCE_DIR"

if [[ -d "$SRC_PATH" ]]; then
  log_warn "Source directory already exists: $SRC_PATH"
  log_warn "Removing and re-extracting..."
  rm -rf "$SRC_PATH"
fi

log_info "Extracting $SOURCE_TARBALL ..."
tar -xzf "$TARBALL_PATH" -C "$DEST_DIR"
log_ok "Extracted to: $SRC_PATH"

# ──────────────────────────────────────────────────────────────────────────────
# Apply patches
# ──────────────────────────────────────────────────────────────────────────────
log_step "Apply patches"

cd "$DEST_DIR"

log_info "Applying patch 1: $PATCH1 (progress tracking, subject length, midline field)..."
if patch -p0 < "$SCRIPT_DIR/$PATCH1"; then
  log_ok "Patch 1 applied successfully."
else
  log_error "Patch 1 failed. Check for conflicts."
  exit 1
fi

if [[ $NO_PATCH2 -eq 0 ]]; then
  log_info "Applying patch 2: $PATCH2 (min/max identity and coverage filters)..."
  if patch -p0 < "$SCRIPT_DIR/$PATCH2"; then
    log_ok "Patch 2 applied successfully."
  else
    log_error "Patch 2 failed. Check for conflicts."
    exit 1
  fi
else
  log_warn "--no-patch2 set, skipping filter parameter patch."
fi

# ──────────────────────────────────────────────────────────────────────────────
# Configure
# ──────────────────────────────────────────────────────────────────────────────
log_step "Configure (ReleaseMT)"

BUILD_DIR="$SRC_PATH/c++/ReleaseMT/build"
cd "$SRC_PATH/c++"

log_info "Running ./configure ..."
./configure \
  --with-mt \
  --without-debug \
  --without-gui \
  --without-vdb \
  --with-projects="$SRC_PATH/c++/projects.txt" \
  CC=/usr/bin/gcc \
  CXX=/usr/bin/g++ \
  2>&1 | tee "$DEST_DIR/configure.log" | grep -E "^configure:|error:|warning: " || true

if [[ ! -d "$BUILD_DIR" ]]; then
  log_error "Configure failed — build directory not created. See $DEST_DIR/configure.log"
  exit 1
fi
log_ok "Configure complete. Build directory: $BUILD_DIR"

# ──────────────────────────────────────────────────────────────────────────────
# Compile
# ──────────────────────────────────────────────────────────────────────────────
log_step "Compile (make -j${JOBS})"

cd "$BUILD_DIR"

log_info "Running make -j${JOBS} all_r ... (this takes 10-20 minutes)"
START_TS=$(date +%s)

make -j"$JOBS" all_r 2>&1 | tee "$DEST_DIR/make.log"
MAKE_STATUS=${PIPESTATUS[0]}

if [[ $MAKE_STATUS -ne 0 ]]; then
  log_error "Compilation failed (exit $MAKE_STATUS). See $DEST_DIR/make.log"
  exit 1
fi

# Verify key binaries exist
BIN_DIR="$SRC_PATH/c++/ReleaseMT/bin"
log_step "Verify output binaries"

EXPECTED_BINS=(blastn blastp blastx tblastn tblastx makeblastdb blastdbcmd)
ALL_OK=1
for bin in "${EXPECTED_BINS[@]}"; do
  if [[ -x "$BIN_DIR/$bin" ]]; then
    SIZE=$(du -h "$BIN_DIR/$bin" | cut -f1)
    log_ok "$bin  ($SIZE)"
  else
    log_error "MISSING: $bin"
    ALL_OK=0
  fi
done

END_TS=$(date +%s)
ELAPSED=$(( END_TS - START_TS ))
ELAPSED_MIN=$(( ELAPSED / 60 ))
ELAPSED_SEC=$(( ELAPSED % 60 ))

echo ""
if [[ $ALL_OK -eq 1 ]]; then
  log_ok "Build complete in ${ELAPSED_MIN}m ${ELAPSED_SEC}s."
  echo ""
  echo -e "${GREEN}${BOLD}Binaries are in:${NC} $BIN_DIR"
  echo ""
  echo "  Quick test:"
  echo "    $BIN_DIR/blastn -help | grep -A3 'Filtering options'"
  echo "    $BIN_DIR/blastp -help | grep -A3 'Filtering options'"
  echo ""
else
  log_error "Some binaries are missing. Check $DEST_DIR/make.log for errors."
  exit 1
fi
