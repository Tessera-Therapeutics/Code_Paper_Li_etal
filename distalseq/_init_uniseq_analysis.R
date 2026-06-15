###### Init UNI-seq analysis ########
# Configuration and setup for Uni-seq analysis pipeline
# Li et al. 2026 - Gene Writing with engineered retrotransposons

# ==============================================================================
# LOAD CONFIGURATION
# ==============================================================================

# Load project configuration
if (file.exists("config.local.R")) {
  source("config.local.R")
} else if (file.exists("../config.local.R")) {
  source("../config.local.R")
} else {
  source("config.example.R")
  warning("Using example config. Copy config.example.R to config.local.R and customize.")
}

# Check that required files exist
check_config()

# ==============================================================================
# LOAD REQUIRED LIBRARIES
# ==============================================================================

required_packages <- c("openxlsx", "dplyr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    stop("Required package '", pkg, "' is not installed. Install with: install.packages('", pkg, "')")
  }
}

# ==============================================================================
# PROJECT SETTINGS
# ==============================================================================

uni_id_cols <- c("chr", "integration_locus", "integration_strand", "GeneName", "GeneStrand", "GeneDistance")
old_id_cols <- c("refIS_chr", "refIS_locus", "integration_strand", "feature_geneName", "feature_strand", "feature_distance")

# Color schemas (from config)
tessera_color_schema <- TESSERA_COLORS

paper_barplot_colors <- data.frame("ppgray" = "#a4a4a4",
                                   "ppblue" = "#3c3c93",
                                   "ppwater" = "#44acc3",
                                   "ppolivegreen" = "#9bbb56"
                                   )

vector_rename <- data.frame("Vector" = c("CR1", "RTE1", "RTE3", "RTE25", "Vingi", "Lenti", "Random"),
                            "ExtendedVector" = c("CR1-1_PH", "RTE-1_MD", "RTE-3_BF", "RTE-25_Lmi", "Vingi-1_Acar", "Lenti", "Random") )

# Load utility functions
source("uniseq_utils.R")

# Load Master file (from config)
uniseq_master_df <- read.xlsx(xlsxFile = UNISEQ_MASTER_FILE, sheet = "uniseq_sample_list")
vector_df <- read.xlsx(xlsxFile = UNISEQ_MASTER_FILE, sheet = "vector")

# Reference genome files (from config)
target_genome_fa_file <- REF_FASTA
target_genome_GTF_file <- REF_GTF_TRANSCRIPT
target_genome_TSS_file <- REF_TSS_BED
target_genome_exons_file <- REF_GTF_EXONS
target_genome_blacklist_mappability_GTF <- REF_BLACKLIST

# Integration site span threshold (from config)
threshold_IS_span <- UNISEQ_PARAMS$threshold_IS_span

# ==============================================================================
# LOAD PUBLIC DATABASES
# ==============================================================================

# MSK Cancer Gene List
msk_db <- read.csv(MSK_CANCER_GENES, header=TRUE, fill=T, sep='\t', check.names = FALSE)
msk_db$GeneName <- msk_db$`Hugo Symbol` # just add a name that I am used to merge
msk_db <- msk_db[which((msk_db$`Gene Type` %in% c("ONCOGENE", "ONCOGENE_AND_TSG", "TSG")) & 
                         msk_db$`# of occurrence within resources (Column J-P)` > 1),] # avoid unknown genes and false positive / weak genes

