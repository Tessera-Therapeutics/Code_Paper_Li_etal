# Configuration file for Li et al. 2026 analysis scripts
# Copy this file to config.local.R and customize for your system
# config.local.R is gitignored and won't be tracked

# ==============================================================================
# PROJECT PATHS
# ==============================================================================

# Base directory (automatically detected as the repo root)
PROJECT_ROOT <- here::here()

# Metadata directory
METADATA_DIR <- file.path(PROJECT_ROOT, "metadata")

# ==============================================================================
# REFERENCE GENOME FILES (hg38/GRCh38)
# ==============================================================================
# Download instructions in metadata/README.md

REF_GENOME_DIR <- file.path(METADATA_DIR, "hg38")

# Reference FASTA (for CRAM decoding and sequence extraction)
REF_FASTA <- file.path(REF_GENOME_DIR, "Homo_sapiens_assembly38.fasta")

# Gene annotations (GTF format)
REF_GTF_TRANSCRIPT <- file.path(REF_GENOME_DIR, "hg38.ncbiRefSeq.sorted.transcript.gtf.gz")

# TSS (Transcription Start Sites)
REF_TSS_BED <- file.path(REF_GENOME_DIR, "hg38.ncbiRefSeq.sorted.TSS2.bed.gz")

# Exon annotations
REF_GTF_EXONS <- file.path(REF_GENOME_DIR, "hg38.ncbiRefSeq.sorted_slice_exons.gtf.gz")

# Low mappability regions (optional - for filtering)
REF_BLACKLIST <- file.path(REF_GENOME_DIR, "hg38.ucsc.encode_grc.lowmappability.2024.sort.gtf")

# ==============================================================================
# PUBLIC DATABASE FILES
# ==============================================================================

DB_DIR <- file.path(METADATA_DIR, "public_databases")

# MSK Cancer Gene List
MSK_CANCER_GENES <- file.path(DB_DIR, "msk.cancerGeneList.tsv")

# CancerMine database
CANCERMINE_DB <- file.path(DB_DIR, "cancermine_collated.tsv")

# UniProt oncogene/tumor suppressor lists
UNIPROT_ONCOGENE <- file.path(DB_DIR, "uniprot-Proto-oncogene.tsv")
UNIPROT_TUMOR_SUPPRESSOR <- file.path(DB_DIR, "uniprot-Tumor-suppressor.tsv")

# ==============================================================================
# SAMPLE METADATA
# ==============================================================================

SAMPLE_INFO_DIR <- file.path(METADATA_DIR, "sample_info")

# DISTAL-seq master file
DISTALSEQ_MASTER_FILE <- file.path(SAMPLE_INFO_DIR, "DISTALseq.Metadata.xlsx")

# Uni-seq master file
UNISEQ_MASTER_FILE <- file.path(SAMPLE_INFO_DIR, "Uniseq.260223.Metadata.xlsx")

# ==============================================================================
# ANALYSIS PARAMETERS
# ==============================================================================

# DISTAL-seq parameters
DISTALSEQ_PARAMS <- list(
  threshold_IS_span = 30,
  min_read_len = 1000,
  min_target_aln_len = 100,
  chromosomes_to_use = c(paste0("chr", 1:22), "chrX", "chrY")
)

# Uni-seq parameters
UNISEQ_PARAMS <- list(
  threshold_IS_span = 30,
  chromosomes_to_use = c(paste0("chr", 1:22), "chrX", "chrY")
)

# WGS parameters
WGS_PARAMS <- list(
  threads = 8,
  min_mapq = 30,
  bin_size = 500000  # 500kb for CNV analysis
)

# ==============================================================================
# OUTPUT DIRECTORIES
# ==============================================================================

# Output directories (will be created if they don't exist)
OUTPUT_DIR <- file.path(PROJECT_ROOT, "output")
FIGURES_DIR <- file.path(OUTPUT_DIR, "figures")
TABLES_DIR <- file.path(OUTPUT_DIR, "tables")

# ==============================================================================
# VISUAL THEME
# ==============================================================================

# Tessera color schema for plots
TESSERA_COLORS <- data.frame(
  txblue = "#9999FF",
  txred = "#FF0066",
  txpurple = "#333399",
  txgray = "#999999",
  orange = "#FF9933",
  red = "#FF3333",
  azure = "#3399FF"
)

# Sample order for plots
SAMPLE_ORDER <- c("RTE25", "RTE3", "CR1", "Vingi")

SAMPLE_LABELS <- c(
  "RTE25" = "RTE-25_Lmi",
  "RTE3"  = "RTE-3_BF",
  "CR1"   = "CR1-1_PH",
  "Vingi" = "Vingi-1_Acar"
)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# Function to check if required files exist
check_config <- function() {
  required_files <- c(
    REF_FASTA,
    REF_GTF_TRANSCRIPT,
    REF_TSS_BED,
    REF_GTF_EXONS
  )

  missing <- required_files[!file.exists(required_files)]

  if (length(missing) > 0) {
    warning("Missing required reference files:\n",
            paste0("  - ", missing, collapse = "\n"),
            "\n\nSee metadata/README.md for download instructions.")
    return(FALSE)
  }

  message("✓ All required reference files found")
  return(TRUE)
}

# Create output directories if they don't exist
setup_output_dirs <- function() {
  dirs <- c(OUTPUT_DIR, FIGURES_DIR, TABLES_DIR)
  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE)
      message("Created directory: ", d)
    }
  }
}

# ==============================================================================
# USAGE
# ==============================================================================
#
# In your analysis scripts:
#
# # Load configuration
# if (file.exists("config.local.R")) {
#   source("config.local.R")
# } else {
#   source("config.example.R")
#   warning("Using example config. Copy config.example.R to config.local.R and customize.")
# }
#
# # Check that files exist
# check_config()
#
# # Create output directories
# setup_output_dirs()
#
