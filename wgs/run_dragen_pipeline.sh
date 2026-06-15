#!/bin/bash
# DRAGEN somatic SV + CNV pipeline
# Processes edited clones vs control for structural variant calling
# Run on a DRAGEN-enabled system (e.g., AWS F1 instance)

set -euo pipefail
trap '' HUP  # ignore hangup — keeps running if SSH session drops

# ==============================================================================
# LOAD CONFIGURATION
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

if [ -f "${REPO_ROOT}/config.local.sh" ]; then
    echo "[INFO] Loading config from config.local.sh"
    source "${REPO_ROOT}/config.local.sh"
elif [ -f "${REPO_ROOT}/config.example.sh" ]; then
    echo "[WARN] config.local.sh not found, using config.example.sh"
    echo "[WARN] Copy config.example.sh to config.local.sh and customize for your system"
    source "${REPO_ROOT}/config.example.sh"
else
    echo "[ERROR] No configuration file found"
    exit 1
fi

# ==============================================================================
# VALIDATE CONFIGURATION
# ==============================================================================

if [ -z "${DATA_DIR}" ] || [ -z "${REF_FASTA}" ] || [ -z "${REF_DIR}" ]; then
    echo "[ERROR] Required configuration variables not set"
    echo "       Please set DATA_DIR, REF_FASTA, and REF_DIR in config.local.sh"
    exit 1
fi

if [ -z "${S3_RESULTS}" ]; then
    echo "[ERROR] S3_RESULTS not set in config.local.sh"
    echo "       Set this to your S3 bucket for results, e.g.:"
    echo "       S3_RESULTS=\"s3://your-bucket/wgs/results\""
    exit 1
fi

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

check_disk() {
    local avail_gb
    avail_gb=$(df -BG "${DATA_DIR}" | awk 'NR==2 {gsub("G",""); print $4}')
    if [[ ${avail_gb} -lt 80 ]]; then
        log "ERROR: Less than 80GB free (${avail_gb}GB). Free space before continuing."
        exit 1
    fi
    log "Disk check OK: ${avail_gb}GB available"
}

# ==============================================================================
# PRE-FLIGHT CHECKS
# ==============================================================================

log "=== DRAGEN SV Calling Pipeline ==="
log "Samples: ${SAMPLES[*]}"

# Check control BAM exists
if [[ ! -f "${CONTROL_BAM}" ]]; then
    log "ERROR: Control BAM not found at ${CONTROL_BAM}"
    log "       Update CONTROL_BAM in config.local.sh"
    exit 1
fi

# Check or create index
if [[ ! -f "${CONTROL_BAM}.bai" ]]; then
    log "Control BAM index missing — indexing now..."
    samtools index "${CONTROL_BAM}"
fi
log "Control BAM OK: ${CONTROL_BAM}"

# Check reference files
if [[ ! -f "${REF_FASTA}" ]]; then
    log "ERROR: Reference FASTA not found: ${REF_FASTA}"
    log "       See metadata/README.md for download instructions"
    exit 1
fi

# ==============================================================================
# PER-SAMPLE PIPELINE
# ==============================================================================

for SAMPLE in "${SAMPLES[@]}"; do
    log "========================================"
    log "Processing: ${SAMPLE}"
    log "========================================"

    SAMPLE_DIR="${DATA_DIR}/${SAMPLE}"
    DRAGEN_OUT="${SAMPLE_DIR}/dragen_sv_out"

    # Clean sample name for output prefix (remove any prefix like "GM25256-")
    PREFIX="${SAMPLE##*-}_vs_control"

    mkdir -p "${SAMPLE_DIR}" "${DRAGEN_OUT}"

    # ── Step 1: Get CRAM/BAM ─────────────────────────────────────────────────

    # Try to find local CRAM first
    CRAM="${SAMPLE_DIR}/${SAMPLE}.name-sorted.cram"
    BAM="${SAMPLE_DIR}/${SAMPLE}.sorted.bam"

    if [[ -f "${BAM}" && -f "${BAM}.bai" ]]; then
        log "BAM already exists and indexed — skipping download+sort: ${BAM}"
    else
        check_disk

        # Check if we have a local CRAM path defined
        if [[ -v LOCAL_CRAM_PATHS[${SAMPLE}] ]]; then
            log "Using local CRAM: ${LOCAL_CRAM_PATHS[${SAMPLE}]}"
            CRAM="${LOCAL_CRAM_PATHS[${SAMPLE}]}"
        # Check if we need to download from S3
        elif [[ -v S3_PATHS[${SAMPLE}] ]]; then
            log "Downloading CRAM from S3..."
            aws s3 cp "${S3_PATHS[${SAMPLE}]}" "${CRAM}"
            log "Download complete: ${CRAM}"
        else
            log "ERROR: No path defined for sample ${SAMPLE}"
            log "       Update S3_PATHS or LOCAL_CRAM_PATHS in config.local.sh"
            exit 1
        fi

        # ── Step 2: Sort and convert to BAM ──────────────────────────────────
        log "Sorting and converting to BAM..."
        samtools sort -@ "${THREADS}" \
            --reference "${REF_FASTA}" \
            -O BAM \
            -o "${BAM}" \
            "${CRAM}" \
            2>"${SAMPLE_DIR}/sort.log"

        log "Sort complete. Removing CRAM..."
        if [[ "${CRAM}" == "${SAMPLE_DIR}"* ]]; then
            # Only remove if it's in our sample dir (not a local file elsewhere)
            rm "${CRAM}"
        fi

        # ── Step 3: Index BAM ────────────────────────────────────────────────
        log "Indexing BAM..."
        samtools index "${BAM}"
        log "Index complete."
    fi

    # ── Step 4: Run DRAGEN somatic SV ────────────────────────────────────────
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

    # ── Step 5: Upload results to S3 (if S3_RESULTS is set) ─────────────────
    if [[ -n "${S3_RESULTS}" ]]; then
        log "Uploading DRAGEN results to S3..."
        aws s3 cp "${DRAGEN_OUT}/" "${S3_RESULTS}/${SAMPLE}/" --recursive
        log "Upload complete: ${S3_RESULTS}/${SAMPLE}/"

        # ── Step 6: Upload BAM + index ───────────────────────────────────────
        log "Uploading BAM and index to S3..."
        aws s3 cp "${BAM}"      "${S3_RESULTS}/${SAMPLE}/${SAMPLE}.sorted.bam"
        aws s3 cp "${BAM}.bai"  "${S3_RESULTS}/${SAMPLE}/${SAMPLE}.sorted.bam.bai"
        log "BAM upload complete."
    else
        log "S3_RESULTS not set — skipping upload"
    fi

    # ── Step 7: Cleanup (optional) ───────────────────────────────────────────
    # Uncomment to remove large local files after upload
    # log "Removing local BAM and DRAGEN output..."
    # rm -f "${BAM}" "${BAM}.bai"
    # rm -rf "${DRAGEN_OUT}"
    # log "Cleanup complete for ${SAMPLE}."

    log "Done: ${SAMPLE}"
done

log "========================================"
log "All samples complete."
log "========================================"
