# WGS Data Processing — Methods & Materials

Draft text for paper M&M section. Covers CRAM→BAM conversion, somatic SV calling (DRAGEN), and genome-wide coverage analysis (mosdepth).

---

## Whole-genome sequencing data processing

Whole-genome sequencing (WGS) was performed on six single-cell-derived edited clones (GM25256-B1 through B6) and one unedited parental control (GM25256-C8) as part of experiment EXP21001214. Raw sequencing data were stored as name-sorted CRAM files aligned to the human reference genome (GRCh38/hg38, Homo_sapiens_assembly38).

CRAM files were downloaded from cloud storage and converted to coordinate-sorted BAM format using `samtools sort` (version 1.x, 8 threads) with the hg38 reference FASTA for CRAM decoding. BAM files were indexed with `samtools index`.

---

## Somatic structural variant calling

Somatic structural variants (SVs) were called for each edited clone against the unedited control (C8) using the DRAGEN Bio-IT Platform (Illumina, version 4.3.13) running on a dedicated AWS F1 FPGA instance. DRAGEN was run in SV-only mode (alignment disabled, pre-aligned BAMs provided) using the hg38 hash table. Each clone BAM was supplied as the tumor input and the C8 BAM as the matched normal input:

```
dragen \
  --enable-sv true \
  --enable-map-align false \
  --tumor-bam-input <clone>.sorted.bam \
  --bam-input GM25256-C8.sorted.bam \
  --ref-dir <hg38_hash_table> \
  --output-file-prefix <clone>_vs_C8
```

SV calls were filtered to PASS-tier breakend (BND) events. Calls overlapping the vector insertion locus on chromosome 6 (chr6:73,520,000–73,522,000) were excluded as vector-derived artefacts arising from the human promoter sequence present in the vector construct.

---

## Genome-wide copy-number coverage analysis

To assess genome-wide copy-number stability in edited clones, read depth was computed in 500 kb non-overlapping windows across the genome using mosdepth (version 0.3.10) with the `--no-per-base` flag to limit output to windowed depth:

```
mosdepth \
  --no-per-base \
  --by 500000 \
  <sample_prefix> \
  <sample>.sorted.bam
```

GC content and N-base fraction per 500 kb bin were computed directly from the GRCh38 reference FASTA. Analysis was restricted to autosomes (chr1–22) and chrX. The following bins were excluded to remove systematic artefacts:

1. **Assembly gaps**: bins with >10% ambiguous (N) bases.
2. **Pericentromeric regions**: bins within ±3 Mb of UCSC-defined centromere coordinates, which contain alpha-satellite and other repetitive sequences where reads map ambiguously.
3. **Depth outliers**: bins with raw depth above the per-sample 99th percentile, which correspond to segmental duplications and other high-copy loci where multi-mapping reads produce artefactually elevated coverage.

After filtering, sample-specific GC bias was corrected using a piecewise-median approach: bins were stratified into 50 equal-width GC bins, the median depth per GC bin was computed, and each bin's depth was divided by the median depth of its GC stratum and rescaled to the genome-wide median. For each clone, log2 ratios of GC-corrected clone depth relative to GC-corrected C8 control depth were computed per 500 kb bin. Each clone's log2-ratio distribution was then median-centred to remove residual library-level offsets. Reference lines at log2 ratio 0 (diploid) and ±1 (2-fold change) are indicated on the genome-wide plots. Analysis was performed in Python using pandas, numpy, scipy, and matplotlib.

---

## Notes (not for paper)

- B1 and B2 CRAMs were processed prior to the automated pipeline; their BAMs were archived to S3 separately (`wgs_bam_files/`).
- B3–B6 were processed by the automated pipeline (`run_dragen_pipeline.sh`) which downloads CRAM, sorts, indexes, runs DRAGEN SV, uploads BAM + DRAGEN output to S3, then removes local files.
- mosdepth was run on a separate EC2 instance with sufficient disk space (`run_mosdepth.sh`); results uploaded to `s3://compbio-discovery-shared/nonLTR/wgs/mosdepth_500kb/`.
- samtools version: confirm on the DRAGEN instance with `samtools --version`.
