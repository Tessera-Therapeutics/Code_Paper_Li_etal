#!/bin/bash
# Download only the first 5 MB of each BAM (enough for the header)
# and print @RG lines to detect library/flowcell batch differences.

set -euo pipefail

TMP=$(mktemp -d)
trap "rm -rf ${TMP}" EXIT

declare -A BAMS=(
    [B1]="s3://compbio-discovery-shared/nonLTR/wgs/wgs_sv/GM25256-B1/B1.sorted.bam"
    [B2]="s3://compbio-discovery-shared/nonLTR/wgs/wgs_bam_files/GM25256-B2.sorted.bam"
    [B3]="s3://compbio-discovery-shared/nonLTR/wgs/wgs_sv/GM25256-B3/B3.sorted.bam"
    [B4]="s3://compbio-discovery-shared/nonLTR/wgs/wgs_sv/GM25256-B4/B4.sorted.bam"
    [B5]="s3://compbio-discovery-shared/nonLTR/wgs/wgs_sv/GM25256-B5/B5.sorted.bam"
    [B6]="s3://compbio-discovery-shared/nonLTR/wgs/wgs_sv/GM25256-B6/B6.sorted.bam"
)

for S in B1 B2 B3 B4 B5 B6; do
    echo "=== GM25256-${S} ==="
    head_bam="${TMP}/${S}.head.bam"
    aws s3api get-object \
        --bucket "$(echo ${BAMS[$S]} | awk -F/ '{print $3}')" \
        --key    "$(echo ${BAMS[$S]} | cut -d/ -f4-)" \
        --range  "bytes=0-5242880" \
        "${head_bam}" > /dev/null 2>&1
    samtools view -H "${head_bam}" 2>/dev/null | grep "^@RG" || echo "  (no @RG found — header may need larger range)"
    echo
done
