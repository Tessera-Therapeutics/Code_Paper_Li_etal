#!/usr/bin/env python3
"""
Configuration file for Li et al. 2026 Python analysis scripts
Copy this file to config.local.py and customize for your system
config.local.py is gitignored and won't be tracked
"""

import os
from pathlib import Path

# ==============================================================================
# PROJECT PATHS
# ==============================================================================

# Base directory (automatically detected as the repo root)
PROJECT_ROOT = Path(__file__).parent.absolute()

# Metadata directory
METADATA_DIR = PROJECT_ROOT / "metadata"

# ==============================================================================
# REFERENCE GENOME FILES (hg38/GRCh38)
# ==============================================================================
# Download instructions in metadata/README.md

REF_GENOME_DIR = METADATA_DIR / "hg38"

# Reference FASTA (for CRAM decoding and sequence extraction)
REF_FASTA = REF_GENOME_DIR / "Homo_sapiens_assembly38.fasta"

# Gene annotations (GTF format)
REF_GTF = REF_GENOME_DIR / "genes.gtf.gz"

# ==============================================================================
# SAMPLE METADATA
# ==============================================================================

SAMPLE_INFO_DIR = METADATA_DIR / "sample_info"

# scRNA-seq sample metadata
SCRNASEQ_SAMPLES = SAMPLE_INFO_DIR / "scRNAseq.Samples.csv"

# ==============================================================================
# scRNA-seq ANALYSIS PARAMETERS
# ==============================================================================

SCRNASEQ_PARAMS = {
    # Data input paths (set to your cellranger output locations)
    "cellranger_base_path": "/path/to/your/cellranger/outputs",

    # Quality control thresholds
    "min_genes_per_cell": 200,
    "min_counts_per_cell": 1000,
    "max_pct_mt": 10,  # Maximum mitochondrial percentage
    "min_cells_per_gene_pct": 0.03,  # Gene must be in at least 3% of cells

    # InferCNV parameters
    "infercnv_window_size": 100,
    "infercnv_step": 1,
    "infercnv_chunksize": 5000,
    "infercnv_n_jobs": 64,

    # Aneuploidy calling
    "abs_cnv_threshold": 0,
    "pct_aberrant_chr_threshold": 0.8,

    # Reference and case treatments
    "reference_treatment": "Activated-Tcells",
    "case_treatments": [
        "Etop 500 nM",
        "Etop 250nM",
        "Nucleofection-only",
        "Vingi EN- RT- 100ng",
        "Vingi WT 100ng",
        "Cas9-B2M-TRAC-DKO",
    ],
}

# ==============================================================================
# OUTPUT DIRECTORIES
# ==============================================================================

# Output directories (will be created if they don't exist)
OUTPUT_DIR = PROJECT_ROOT / "output"
FIGURES_DIR = OUTPUT_DIR / "figures"
TABLES_DIR = OUTPUT_DIR / "tables"

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

def check_config():
    """Check if required files exist"""
    required_files = [
        SCRNASEQ_SAMPLES,
    ]

    missing = [f for f in required_files if not f.exists()]

    if missing:
        print("WARNING: Missing required files:")
        for f in missing:
            print(f"  - {f}")
        print("\nSee metadata/README.md for setup instructions.")
        return False

    print("✓ All required metadata files found")
    return True


def setup_output_dirs():
    """Create output directories if they don't exist"""
    dirs = [OUTPUT_DIR, FIGURES_DIR, TABLES_DIR]
    for d in dirs:
        if not d.exists():
            d.mkdir(parents=True)
            print(f"Created directory: {d}")


# ==============================================================================
# USAGE
# ==============================================================================
#
# In your analysis scripts:
#
# import sys
# from pathlib import Path
#
# # Load configuration
# if Path("config.local.py").exists():
#     import config.local as config
# else:
#     import config.example as config
#     print("Warning: Using example config. Copy config.example.py to config.local.py and customize.")
#
# # Check that files exist
# config.check_config()
#
# # Create output directories
# config.setup_output_dirs()
#