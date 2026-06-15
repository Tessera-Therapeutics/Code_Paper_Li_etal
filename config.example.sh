#!/bin/bash
# Configuration for WGS analysis scripts
# Copy this file to config.local.sh and customize for your system
# config.local.sh is gitignored and won't be tracked

# ==============================================================================
# AWS S3 CONFIGURATION
# ==============================================================================
# IMPORTANT: Set these to YOUR S3 bucket locations

# S3 bucket for input CRAMs (if downloading from S3)
# Leave empty if using local files
S3_INPUT_BUCKET=""

# S3 bucket for analysis results
# Example: S3_RESULTS="s3://your-bucket/wgs/results"
S3_RESULTS=""

# ==============================================================================
# REFERENCE GENOME
# ==============================================================================

# Local data directory
DATA_DIR="/path/to/your/wgs/data"

# Reference FASTA for hg38
REF_FASTA="${DATA_DIR}/ref/Homo_sapiens_assembly38.fasta"

# DRAGEN reference hash table directory (if using DRAGEN)
REF_DIR="/path/to/dragen/ref/hg38"

# ==============================================================================
# SAMPLE PATHS
# ==============================================================================
# Define paths to your CRAM files

# Control sample BAM/CRAM
CONTROL_BAM="${DATA_DIR}/control-sample.sorted.bam"

# Sample CRAMs - customize this associative array for your samples
declare -A S3_PATHS=(
    [sample-B1]="s3://your-bucket/path/to/sample-B1.cram"
    [sample-B2]="s3://your-bucket/path/to/sample-B2.cram"
    # Add more samples as needed
)

# Or if using local files:
# declare -A LOCAL_CRAM_PATHS=(
#     [sample-B1]="/path/to/sample-B1.cram"
#     [sample-B2]="/path/to/sample-B2.cram"
# )

# ==============================================================================
# ANALYSIS PARAMETERS
# ==============================================================================

# Number of threads for samtools/processing
THREADS=8

# Samples to process (space-separated list)
SAMPLES=(sample-B1 sample-B2)

# ==============================================================================
# CHIMERIC READ DETECTION PARAMETERS
# ==============================================================================

# Minimum mapping quality
MIN_MAPQ=30

# Minimum alignment length
MIN_ALN_LEN=20

# Minimum junction supporting reads
MIN_JUNCTION_READS=4

# Context bases around breakpoint
CONTEXT_BP=30

# ==============================================================================
# USAGE
# ==============================================================================
#
# In your shell scripts, source this config:
#
# if [ -f config.local.sh ]; then
#     source config.local.sh
# else
#     echo "Warning: config.local.sh not found, using defaults from config.example.sh"
#     source config.example.sh
# fi
#
