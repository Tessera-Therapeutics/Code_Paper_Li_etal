#!/bin/bash
# Mosdepth coverage pipeline for CNV analysis
# Runs mosdepth (500kb bins) on control + all samples.
# Li et al. 2026 - Gene Writing with engineered retrotransposons

set -euo pipefail

# Load configuration
if [ -f config.local.sh ]; then
    source config.local.sh
elif [ -f ../config.local.sh ]; then
    source ../config.local.sh
else
    echo "Warning: config.local.sh not found, using defaults from config.example.sh"
    source config.example.sh
fi

# ── Configuration ────────────────────────────────────────────────────────────
OUT_DIR="${DATA_DIR}/mosdepth"
BIN_SIZE=500000
S3_BAMS="${S3_INPUT_BUCKET}"

CONTROL_BAM="${CONTROL_BAM}"

# ── Helpers ──────────────────────────────────────────────────────────────────
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

check_disk() {
    local avail_gb
    avail_gb=$(df -BG "${DATA_DIR}" | awk 'NR==2 {gsub("G",""); print $4}')
    if [[ ${avail_gb} -lt 150 ]]; then
        log "ERROR: Less than 150GB free (${avail_gb}GB). Free space before continuing."
        exit 1
    fi
    log "Disk check OK: ${avail_gb}GB available"
}

run_mosdepth() {
    local sample="$1"
    local bam="$2"
    local prefix="${OUT_DIR}/${sample}"

    if [[ -f "${prefix}.regions.bed.gz" ]]; then
        log "mosdepth output already exists — skipping: ${prefix}.regions.bed.gz"
        return
    fi

    log "Running mosdepth for ${sample}..."
    mosdepth \
        --no-per-base \
        --by "${BIN_SIZE}" \
        --threads "${THREADS}" \
        "${prefix}" \
        "${bam}"
    log "mosdepth complete: ${sample}"
}

# ── Setup ────────────────────────────────────────────────────────────────────
mkdir -p "${OUT_DIR}"
log "=== RW-DA-004 mosdepth pipeline ==="
log "Output dir: ${OUT_DIR}"
log "Bin size: ${BIN_SIZE} bp"

# ── C8 control (local) ───────────────────────────────────────────────────────
log "========================================"
log "Processing: GM25256-C8 (control)"
log "========================================"

if [[ ! -f "${CONTROL_BAM}" ]]; then
    log "ERROR: C8 BAM not found at ${CONTROL_BAM}"
    exit 1
fi
if [[ ! -f "${CONTROL_BAM}.bai" ]]; then
    log "Indexing C8 BAM..."
    samtools index "${CONTROL_BAM}"
fi

run_mosdepth "GM25256-C8" "${CONTROL_BAM}"

# ── B samples from S3 ────────────────────────────────────────────────────────
declare -A S3_BAM_PATHS=(
    [GM25256-B1]="${S3_BAMS}/GM25256-B1/B1.sorted.bam"
    [GM25256-B2]="${S3_BAMS}/GM25256-B2.sorted.bam"
    [GM25256-B3]="${S3_BAMS}/GM25256-B3/B3.sorted.bam"
    [GM25256-B4]="${S3_BAMS}/GM25256-B4/B4.sorted.bam"
    [GM25256-B5]="${S3_BAMS}/GM25256-B5/B5.sorted.bam"
    [GM25256-B6]="${S3_BAMS}/GM25256-B6/B6.sorted.bam"
)

SAMPLES=(GM25256-B1 GM25256-B2 GM25256-B3 GM25256-B4 GM25256-B5 GM25256-B6)

for SAMPLE in "${SAMPLES[@]}"; do
    log "========================================"
    log "Processing: ${SAMPLE}"
    log "========================================"

    BAM="${DATA_DIR}/${SAMPLE}.sorted.bam"

    if [[ -f "${OUT_DIR}/${SAMPLE}.regions.bed.gz" ]]; then
        log "mosdepth output already exists — skipping download: ${SAMPLE}"
    else
        check_disk

        log "Downloading BAM from S3..."
        aws s3 cp "${S3_BAM_PATHS[${SAMPLE}]}"     "${BAM}"
        aws s3 cp "${S3_BAM_PATHS[${SAMPLE}]}.bai" "${BAM}.bai"
        log "Download complete."

        run_mosdepth "${SAMPLE}" "${BAM}"

        log "Removing local BAM..."
        rm -f "${BAM}" "${BAM}.bai"
        log "Cleanup done."
    fi

    log "Uploading mosdepth results to S3..."
    aws s3 cp "${OUT_DIR}/${SAMPLE}.regions.bed.gz"       "${S3_RESULTS}/${SAMPLE}.regions.bed.gz"
    aws s3 cp "${OUT_DIR}/${SAMPLE}.mosdepth.summary.txt" "${S3_RESULTS}/${SAMPLE}.mosdepth.summary.txt"
    log "Upload complete: ${S3_RESULTS}/${SAMPLE}"

    log "Done: ${SAMPLE}"
done

# Upload C8 results too
log "Uploading C8 mosdepth results to S3..."
aws s3 cp "${OUT_DIR}/GM25256-C8.regions.bed.gz"       "${S3_RESULTS}/GM25256-C8.regions.bed.gz"
aws s3 cp "${OUT_DIR}/GM25256-C8.mosdepth.summary.txt" "${S3_RESULTS}/GM25256-C8.mosdepth.summary.txt"

log "========================================"
log "All samples complete. Results in: ${S3_RESULTS}"
log "========================================"
