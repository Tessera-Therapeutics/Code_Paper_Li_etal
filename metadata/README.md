# Reference Data

This directory is for reference genome files used in the analyses.

## Required Files

Download these hg38/GRCh38 reference files and place them in `metadata/hg38/`:

### 1. Reference FASTA

```bash
mkdir -p hg38
cd hg38

# Download reference genome
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai
```

### 2. Gene Annotations (GTF)

```bash
# Download UCSC RefSeq annotations
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz

# Sort and compress
gunzip hg38.ncbiRefSeq.gtf.gz
grep -v "^#" hg38.ncbiRefSeq.gtf | sort -k1,1 -k4,4n | bgzip > hg38.ncbiRefSeq.sorted.transcript.gtf.gz
tabix -p gff hg38.ncbiRefSeq.sorted.transcript.gtf.gz
```

### 3. TSS (Transcription Start Sites)

```bash
# Extract TSS from GTF
zcat hg38.ncbiRefSeq.sorted.transcript.gtf.gz | \
  awk '$3=="transcript"' | \
  awk 'BEGIN{OFS="\t"} {if($7=="+") print $1,$4-1,$4,$10,$6,$7; else print $1,$5-1,$5,$10,$6,$7}' | \
  sed 's/"//g; s/;//g' | \
  sort -k1,1 -k2,2n | uniq | \
  bgzip > hg38.ncbiRefSeq.sorted.TSS2.bed.gz
```

### 4. Exons

```bash
# Extract exons from GTF
zcat hg38.ncbiRefSeq.sorted.transcript.gtf.gz | \
  awk '$3=="exon"' | \
  bgzip > hg38.ncbiRefSeq.sorted_slice_exons.gtf.gz
```

## File Sizes

- Reference FASTA: ~3 GB
- GTF annotations: ~50 MB
- TSS BED: ~10 MB
- Exons GTF: ~30 MB

**Total**: ~3-4 GB

## Note

These files are too large for git and are excluded via `.gitignore`. Download them separately for each system where you run the analyses.

## Sample Metadata

User-provided sample metadata files should go in `metadata/sample_info/`:
- Master sample tracking (Excel/CSV)
- Clone annotations
- Integration site references

## Public Databases

Cancer gene databases can be placed in `metadata/public_databases/`:
- MSK Cancer Gene List (download from oncokb.org)
- CancerMine annotations (download from bionlp.bcgsc.ca/cancermine/)
- UniProt oncogene/tumor suppressor lists

See analysis scripts for expected file formats.
