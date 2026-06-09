#!/bin/bash
# RW-DA-004 — DRAGEN somatic SV + CNV pipeline
# Processes all 6 clones (B1–B6) vs C8 sequentially.
# Run on the DRAGEN EC2 instance.

set -euo pipefail
trap '' HUP  # ignore hangup — keeps running if SSH session drops

# ── Configuration ────────────────────────────────────────────────────────────
DATA_DIR="/home/ec2-user/wgs/data"
REF_FASTA="${DATA_DIR}/ref/Homo_sapiens_assembly38.fasta"
REF_DIR="/staging/ref/hg38"
THREADS=8
S3_RESULTS="s3://compbio-discovery-shared/nonLTR/wgs/wgs_sv"

# C8 control BAM — adjust path if C8 is in a subdirectory
CONTROL_BAM="${DATA_DIR}/GM25256-C8.sorted.bam"

# ── Sample definitions ───────────────────────────────────────────────────────
declare -A S3_PATHS=(
    [GM25256-B1]="s3://fvl58-compbio/archive/baltshuler/baltshuler-wgs-01/scratch/from-s3/rt-platform/rt-wgs/2021-08_rt-wgs_EXP21001214/GM25256-B1/GM25256-B1.name-sorted.cram"
    [GM25256-B2]="s3://fvl58-compbio/analyses/baltshuler/wgs-ec2-storage-to-sort/wgs-02/from-s3/2021-08_rt-wgs_EXP21001214/GM25256-B2/GM25256-B2.name-sorted.cram"
    [GM25256-B3]="s3://fvl58-compbio/analyses/baltshuler/wgs-ec2-storage-to-sort/wgs-03/from-s3/rt-platform/rt-wgs/2021-08_rt-wgs_EXP21001214/GM25256-B3/GM25256-B3.name-sorted.cram"
    [GM25256-B4]="s3://fvl58-compbio/analyses/baltshuler/wgs-ec2-storage-to-sort/wgs-04/from-s3/rt-platform/rt-wgs/2021-08_rt-wgs_EXP21001214/GM25256-B4/GM25256-B4.name-sorted.cram"
    [GM25256-B5]="s3://fvl58-compbio/analyses/baltshuler/wgs-ec2-storage-to-sort/wgs-05/from-s3/rt-platform/rt-wgs/2021-08_rt-wgs_EXP21001214/GM25256-B5/GM25256-B5.name-sorted.cram"
    [GM25256-B6]="s3://fvl58-compbio/analyses/baltshuler/wgs-ec2-storage-to-sort/wgs-06/from-s3/rt-platform/rt-wgs/2021-08_rt-wgs_EXP21001214/GM25256-B6/GM25256-B6.name-sorted.cram"
)

SAMPLES=(GM25256-B3 GM25256-B4 GM25256-B5 GM25256-B6)
# B1 and B2 BAMs already in S3 — skipped
# SAMPLES=(GM25256-B1 GM25256-B2 GM25256-B3 GM25256-B4 GM25256-B5 GM25256-B6)

# ── Helpers ──────────────────────────────────────────────────────────────────
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

check_disk() {
    local avail_gb
    avail_gb=$(df -BG "${DATA_DIR}" | awk 'NR==2 {gsub("G",""); print $4}')
    if [[ ${avail_gb} -lt 80 ]]; then
        log "ERROR: Less than 80GB free (${avail_gb}GB). Free space before continuing."
        exit 1
    fi
    log "Disk check OK: ${avail_gb}GB available"
}

# ── Pre-flight checks ────────────────────────────────────────────────────────
log "=== RW-DA-004 DRAGEN pipeline ==="
log "Samples: ${SAMPLES[*]}"

if [[ ! -f "${CONTROL_BAM}" ]]; then
    log "ERROR: Control BAM not found at ${CONTROL_BAM}"
    exit 1
fi
if [[ ! -f "${CONTROL_BAM}.bai" ]]; then
    log "Control BAM index missing — indexing now..."
    samtools index "${CONTROL_BAM}"
fi
log "Control BAM OK: ${CONTROL_BAM}"

# ── Per-sample pipeline ──────────────────────────────────────────────────────
for SAMPLE in "${SAMPLES[@]}"; do
    log "========================================"
    log "Processing: ${SAMPLE}"
    log "========================================"

    SAMPLE_DIR="${DATA_DIR}/${SAMPLE}"
    CRAM="${SAMPLE_DIR}/${SAMPLE}.name-sorted.cram"
    BAM="${SAMPLE_DIR}/${SAMPLE}.sorted.bam"
    DRAGEN_OUT="${SAMPLE_DIR}/dragen_sv_out"
    PREFIX="${SAMPLE##GM25256-}_vs_C8"   # e.g. B1_vs_C8

    mkdir -p "${SAMPLE_DIR}" "${DRAGEN_OUT}"

    # Step 1 — Download CRAM
    if [[ -f "${BAM}" && -f "${BAM}.bai" ]]; then
        log "BAM already exists and indexed — skipping download+sort: ${BAM}"
    else
        check_disk

        log "Downloading CRAM from S3..."
        aws s3 cp "${S3_PATHS[${SAMPLE}]}" "${CRAM}"
        log "Download complete: ${CRAM}"

        # Step 2 — Coordinate sort + convert to BAM
        log "Sorting and converting to BAM..."
        samtools sort -@ "${THREADS}" \
            --reference "${REF_FASTA}" \
            -O BAM \
            -o "${BAM}" \
            "${CRAM}" \
            2>"${SAMPLE_DIR}/sort.log"
        log "Sort complete. Removing CRAM..."
        rm "${CRAM}"

        # Step 3 — Index BAM
        log "Indexing BAM..."
        samtools index "${BAM}"
        log "Index complete."
    fi

    # Step 4 — DRAGEN somatic SV
    if [[ -f "${DRAGEN_OUT}/${PREFIX}.sv.vcf.gz" ]]; then
        log "DRAGEN output already exists — skipping: ${DRAGEN_OUT}/${PREFIX}.sv.vcf.gz"
    else
        log "Running DRAGEN somatic SV..."
        dragen \
            --enable-sv true \
            --enable-map-align false \
            --tumor-bam-input "${BAM}" \
            --bam-input "${CONTROL_BAM}" \
            --ref-dir "${REF_DIR}" \
            --output-directory "${DRAGEN_OUT}" \
            --output-file-prefix "${PREFIX}" \
            > "${SAMPLE_DIR}/dragen_sv.log" 2>&1
        log "DRAGEN complete: ${SAMPLE}"
    fi

    # Step 5 — Upload DRAGEN results to S3
    log "Uploading DRAGEN results to S3..."
    aws s3 cp "${DRAGEN_OUT}/" "${S3_RESULTS}/${SAMPLE}/" --recursive
    log "Upload complete: ${S3_RESULTS}/${SAMPLE}/"

    # Step 6 — Upload BAM + index to S3 alongside DRAGEN output
    log "Uploading BAM and index to S3..."
    aws s3 cp "${BAM}"      "${S3_RESULTS}/${SAMPLE}/${SAMPLE##GM25256-}.sorted.bam"
    aws s3 cp "${BAM}.bai"  "${S3_RESULTS}/${SAMPLE}/${SAMPLE##GM25256-}.sorted.bam.bai"
    log "BAM upload complete."

    # Step 7 — Remove large local files to free disk space
    log "Removing local BAM and DRAGEN output..."
    rm -f "${BAM}" "${BAM}.bai"
    rm -rf "${DRAGEN_OUT}"
    log "Cleanup complete for ${SAMPLE}."

    log "Done: ${SAMPLE}"
done

log "========================================"
log "All samples complete."
log "========================================"
