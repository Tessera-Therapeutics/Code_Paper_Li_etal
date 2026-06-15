#!/bin/bash
# Compute GC content per 500kb bin from hg38 reference FASTA.
# Uses only samtools + python3 (no bedtools needed).
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

REF="${REF_FASTA}"
BIN_SIZE=500000
GC_OUT="/tmp/hg38_500kb_gc.tsv"
S3_OUT="${S3_RESULTS}/hg38_500kb_gc.tsv"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

log "Checking tools..."
for tool in samtools python3 aws; do
    command -v ${tool} >/dev/null || { log "ERROR: ${tool} not found"; exit 1; }
done

log "Indexing reference (if needed)..."
[[ -f "${REF}.fai" ]] || samtools faidx "${REF}"

log "Computing GC content per ${BIN_SIZE}bp bin..."
python3 - "${REF}" "${BIN_SIZE}" "${GC_OUT}" <<'EOF'
import sys, gzip

ref_path  = sys.argv[1]
bin_size  = int(sys.argv[2])
out_path  = sys.argv[3]

# read FAI to get chrom sizes
chroms = []
with open(ref_path + ".fai") as f:
    for line in f:
        parts = line.split("\t")
        chroms.append((parts[0], int(parts[1])))

open_fn = gzip.open if ref_path.endswith(".gz") else open

count = 0
with open_fn(ref_path, "rt") as fa, open(out_path, "w") as out:
    out.write("chrom\tstart\tend\tgc_frac\tn_frac\n")
    chrom_name = None
    seq_buf = []

    def flush_chrom(name, seq):
        global count
        chrom_len = len(seq)
        for start in range(0, chrom_len, bin_size):
            end = min(start + bin_size, chrom_len)
            window = seq[start:end].upper()
            gc = window.count("G") + window.count("C")
            at = window.count("A") + window.count("T")
            acgt = gc + at
            win_len = end - start
            n_frac = (win_len - acgt) / win_len if win_len > 0 else 1.0
            frac = gc / acgt if acgt > 0 else 0.0
            out.write(f"{name}\t{start}\t{end}\t{frac:.6f}\t{n_frac:.6f}\n")
            count += 1

    for line in fa:
        line = line.rstrip()
        if line.startswith(">"):
            if chrom_name is not None:
                flush_chrom(chrom_name, "".join(seq_buf))
            chrom_name = line[1:].split()[0]
            seq_buf = []
        else:
            seq_buf.append(line)
    if chrom_name is not None:
        flush_chrom(chrom_name, "".join(seq_buf))

print(f"Written {count} bins to {out_path}", flush=True)
EOF

log "GC content computed: ${GC_OUT} ($(tail -n +2 ${GC_OUT} | wc -l) bins)"

log "Uploading to S3..."
aws s3 cp "${GC_OUT}" "${S3_OUT}"
log "Done: ${S3_OUT}"
