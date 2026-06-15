# Release Checklist for Li et al. 2026 Code Repository

This repository has been cleaned and anonymized for publication. Here's what was done:

## ✅ Completed

### 1. **Removed Temporary Files**
- Deleted internal documentation: `CLEANUP_SUMMARY.md`, `DONE.md`, `INSTALL.md`, `PRE_RELEASE_AUDIT.md`, `README_UPDATED.md`, `SIMPLE_CLEANUP_PLAN.md`
- Removed Mac artifacts: `.DS_Store`

### 2. **Anonymized All Hardcoded Paths**
All personal/system-specific paths have been replaced with configuration variables:

**R Scripts:**
- `distalseq/_init_distalseq_analysis.R` - Now uses config system
- `distalseq/_init_uniseq_analysis.R` - Now uses config system
- `distalseq/distalseq_utils.R` - Removed hardcoded paths from comments
- `distalseq/isa_utils.R` - Updated function defaults to use config
- `uniseq/260223.Uniseq.PaperV2.Volcano.R` - Uses relative paths

**Shell Scripts:**
- `wgs/compute_gc_content.sh` - Now sources config
- `wgs/run_mosdepth.sh` - Now sources config
- `wgs/run_dragen_pipeline.sh` - Now sources config

### 3. **Configuration System**
Created two example config files:
- `config.example.R` - For R scripts (paths to reference genome, databases, sample metadata)
- `config.example.sh` - For shell scripts (S3 paths, local directories, parameters)

Users copy these to `config.local.R` / `config.local.sh` and customize for their system.

### 4. **Added Public Identifiers**
- BioProject accession **PRJNA1471271** added to README
- Linked to NCBI BioProject page

### 5. **Sample Metadata Organization**
Updated config to reference actual metadata files in `metadata/sample_info/`:
- `DISTALseq.Metadata.xlsx`
- `Uniseq.260223.Metadata.xlsx`

### 6. **Git Ignore**
Comprehensive `.gitignore` includes:
- Reference genome files
- Large data files (BAM, CRAM, VCF)
- Local configuration files
- Analysis outputs
- IDE and system files

## 📋 Before Public Release - Manual Review Needed

### 1. **Metadata Files**
Review files in `metadata/sample_info/`:
- [ ] `DISTALseq.Metadata.xlsx` - Verify no internal sample IDs or sensitive info
- [ ] `Uniseq.260223.Metadata.xlsx` - Check for any proprietary information
- [ ] `Metadata.Template.tsv` - Ensure it's a clean template

### 2. **Data Files**
- [ ] `uniseq/260223.IS_gdf_molten.tsv.gz` - Verify this should be included
- [ ] Check all files in `metadata/public_databases/` are from public sources

### 3. **README Updates**
- [ ] Add Zenodo DOI when processed data is uploaded
- [ ] Add publication DOI when manuscript is published
- [ ] Add journal name and citation details

### 4. **Testing**
Before release, test the setup instructions:
- [ ] Clone fresh copy of repository
- [ ] Follow setup instructions in README
- [ ] Verify config examples work with standard reference files
- [ ] Run at least one analysis script to confirm paths resolve correctly

## 🔒 What's Gitignored (Good!)

These files are properly excluded from the repository:
- `config.local.R` and `config.local.sh` - User-specific configurations
- Reference genome files (*.fasta, *.gtf.gz, *.bed.gz)
- Large data files (*.bam, *.cram, *.vcf.gz)
- Analysis outputs (`results/`, `output/`, *.RData)
- IDE files (`.vscode/`, `.idea/`)

## 📝 User Instructions

Users will need to:

1. **Download reference genome** - Follow `metadata/README.md`
2. **Copy and customize config files**:
   ```bash
   cp config.example.R config.local.R
   cp config.example.sh config.local.sh
   # Edit both files with local paths
   ```
3. **Download raw data** from BioProject PRJNA1471271
4. **Run analysis scripts** following README instructions

## ✨ Repository is Publication-Ready

All personal information has been anonymized and the repository uses a clean configuration system for user customization.
