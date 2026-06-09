# DRAGEN Somatic SV Analysis — GM25256 B1–B6 vs C8

## Goal
Detect genome→genome translocations (secondary genomic instability) caused by vector insertion in experiment samples B1–B6, using untransduced C8 as control. Vector locus (chr6:~73520000–73522000) to be masked from results.

## DRAGEN Instance
- Separate EC2 instance (no Claude installed — commands run manually)
- **DRAGEN version:** 4.3.13
- **Reference build:** hg38 (Homo_sapiens_assembly38)
  - FASTA: `/home/ec2-user/wgs/data/ref/Homo_sapiens_assembly38.fasta`
  - Hash table: `/staging/ref/hg38/`

## Sample Directories
All sample data lives under `/home/ec2-user/wgs/data/<sample>/`

| Role | Sample | BAM path |
|------|--------|---------|
| Control (normal) | GM25256-C8 | `/home/ec2-user/wgs/data/GM25256-C8/GM25256-C8.sorted.bam` |
| Experiment | GM25256-B1 | `/home/ec2-user/wgs/data/GM25256-B1/GM25256-B1.sorted.bam` |
| Experiment | GM25256-B2 | `/home/ec2-user/wgs/data/GM25256-B2/GM25256-B2.sorted.bam` ✅ done |
| Experiment | GM25256-B3 | `/home/ec2-user/wgs/data/GM25256-B3/GM25256-B3.sorted.bam` |
| Experiment | GM25256-B4 | `/home/ec2-user/wgs/data/GM25256-B4/GM25256-B4.sorted.bam` |
| Experiment | GM25256-B5 | `/home/ec2-user/wgs/data/GM25256-B5/GM25256-B5.sorted.bam` |
| Experiment | GM25256-B6 | `/home/ec2-user/wgs/data/GM25256-B6/GM25256-B6.sorted.bam` |

See [cram_locations.md](cram_locations.md) for S3 paths to all original CRAMs.

## Per-Sample Pipeline

Process one sample at a time due to disk space constraints. Delete CRAM after BAM is confirmed good before starting the next sample.

### Step 1 — Download CRAM from S3
```bash
aws s3 cp <S3_PATH> /home/ec2-user/wgs/data/<SAMPLE>/<SAMPLE>.name-sorted.cram
```

### Step 2 — Coordinate sort + convert to BAM (deletes CRAM on success)
```bash
nohup samtools sort -@ 8 \
  --reference /home/ec2-user/wgs/data/ref/Homo_sapiens_assembly38.fasta \
  -O BAM \
  -o /home/ec2-user/wgs/data/<SAMPLE>/<SAMPLE>.sorted.bam \
  /home/ec2-user/wgs/data/<SAMPLE>/<SAMPLE>.name-sorted.cram \
  > /home/ec2-user/wgs/data/<SAMPLE>/sort.log 2>&1 && \
rm /home/ec2-user/wgs/data/<SAMPLE>/<SAMPLE>.name-sorted.cram &
```

### Step 3 — Index BAM
```bash
samtools index /home/ec2-user/wgs/data/<SAMPLE>/<SAMPLE>.sorted.bam
```

### Step 4 — Run DRAGEN somatic SV
```bash
mkdir -p /home/ec2-user/wgs/data/<SAMPLE>/dragen_sv_out

nohup dragen \
  --enable-sv true \
  --enable-map-align false \
  --tumor-bam-input /home/ec2-user/wgs/data/<SAMPLE>/<SAMPLE>.sorted.bam \
  --bam-input /home/ec2-user/wgs/data/GM25256-C8/GM25256-C8.sorted.bam \
  --ref-dir /staging/ref/hg38 \
  --output-directory /home/ec2-user/wgs/data/<SAMPLE>/dragen_sv_out \
  --output-file-prefix <SAMPLE>_vs_C8 \
  > /home/ec2-user/wgs/data/<SAMPLE>/dragen_sv.log 2>&1 &
```

Monitor: `tail -f /home/ec2-user/wgs/data/<SAMPLE>/dragen_sv.log`

## B2 Results (completed 2026-04-15)

DRAGEN output: `/home/ec2-user/wgs/data/dragen_sv_out2/B2_vs_C8.sv.vcf.gz`

**SV Summary (PASS):** 18 total — 3 deletions, 7 insertions, 0 duplications, 8 BND pairs

**PASS BND calls:**
- chr6:73520491–73521211 hub → chr1, chr2, chr18, chr19 — **vector artefact** (chr6 locus is human promoter in vector)
- chr3:127792419 ↔ chr3:127793865 — **somatic** (0 alt reads in C8, 3–4 in B2)
- chr15:64080412 ↔ chr15:64081828 — **somatic** (0 alt reads in C8, 4–6 in B2)
- chr15:64080696 ↔ chr15:64081842 — **likely artefact** (1 alt read in both B2 and C8)

## After All Samples Complete
1. Filter each VCF: `bcftools view -i 'INFO/SVTYPE="BND" && FILTER="PASS"'`
2. Exclude chr6:73520000–73522000 (vector promoter artefact)
3. Compare across B1–B6 — recurring BNDs at same locus are higher confidence
4. Annotate against gene bodies, fragile sites, repeat elements
5. Filter by read support (PR/SR fields)
