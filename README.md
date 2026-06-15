# Code for Li et al., 2026

Analysis code supporting the manuscript:

**"In vivo Gene Writing with engineered retrotransposons"**  
Li et al., 2026  
*Link to publication to be added upon release.*

## Repository Contents

This repository contains computational analysis code for three sequencing methods:

```
Code_Paper_Li_etal/
├── distalseq/      # DISTAL-seq integration site analysis (R)
├── uniseq/         # Uni-seq analysis (R)
├── wgs/            # Whole-genome sequencing analysis (bash/Python)
└── metadata/       # Reference genome files go here (see metadata/README.md)
```

## Methods Overview

- **DISTAL-seq** - Long-read sequencing for precise integration site mapping and vector characterization
- **Uni-seq** - Short-read amplicon sequencing for quantitative integration site profiling  
- **WGS** - Whole-genome sequencing for off-target analysis and copy number variation detection

Each directory contains analysis scripts used to generate results in the manuscript.

## Setup

### Requirements

- **R** (≥ 4.0) with packages: dplyr, ggplot2, ggrepel, pheatmap, openxlsx, reshape2, vegan
- **Python** (≥ 3.8) with: pysam, pandas, boto3
- **Tools**: samtools, bedtools, mosdepth
- **Reference**: hg38/GRCh38 genome

### Quick Start

1. **Download reference genome files**
   ```bash
   cd metadata
   # See metadata/README.md for download instructions
   ```

2. **Update file paths in analysis scripts**
   - Look for `# UPDATE PATH:` or `# USER CONFIGURATION` comments
   - Point to your reference genome files and sample metadata

3. **Run analyses**
   ```r
   # DISTAL-seq
   setwd("distalseq")
   source("_init_distalseq_analysis.R")
   # Then run analysis scripts in date order

   # Uni-seq
   setwd("uniseq")
   source("260223.Uniseq.PaperV2.Volcano.R")
   ```

   ```bash
   # WGS
   cd wgs
   # Update configuration in run_dragen_pipeline.sh
   bash run_dragen_pipeline.sh
   ```

## Data Availability

- **Raw sequencing data**: NCBI BioProject [PRJNA1471271](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1471271)
- **Processed integration sites**: [Zenodo DOI to be added]
- **Reference genome**: hg38/GRCh38 - download from UCSC or Broad Institute (see [metadata/README.md](metadata/README.md))

## Citation

If you use this code, please cite:

```
Li et al. (2026). "In vivo Gene Writing with engineered retrotransposons."
[Journal]. DOI: [to be added]
```

## License

MIT License - Copyright (c) 2026 Tessera Therapeutics

See [LICENSE](LICENSE) file for details.

## Notes

- Analysis scripts are provided as-is to document methods used in the manuscript
- File paths need to be updated for your system (see comments in scripts)
- For detailed WGS methods, see `wgs/methods_wgs.md`

## Contact

For questions, see the manuscript for corresponding author contact information.
