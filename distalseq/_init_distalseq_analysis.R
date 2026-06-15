###### Init DISTAL-seq analysis ########
# Configuration and setup for DISTAL-seq analysis pipeline
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

ProjectID <- "DS"

# Column name mappings
id_cols <- c("chr", "integration_locus", "integration_strand", "GeneName", "GeneStrand", "GeneDistance")
old_id_cols <- c("refIS_chr", "refIS_locus", "integration_strand", "feature_geneName", "feature_strand", "feature_distance")

# Color schema for plots (from config)
tessera_color_schema <- TESSERA_COLORS

# Chromosomes to analyze (from config)
chromosomes_to_use <- DISTALSEQ_PARAMS$chromosomes_to_use

# Integration site span threshold (from config)
threshold_IS_span <- DISTALSEQ_PARAMS$threshold_IS_span

# ddPCR faceting order
faceting_order_ddpcr <- c("ddPCR_Set1", "ddPCR_Set2", "ddPCR_NoPrimer")

# ==============================================================================
# LOAD METADATA FILES
# ==============================================================================

message("[INFO] Loading metadata files...")

# Clone annotations (if applicable)
if (file.exists(CLONE_ANNOTATIONS)) {
  dict_clones <- read.xlsx(xlsxFile = CLONE_ANNOTATIONS)
  message("[INFO] Loaded clone annotations: ", CLONE_ANNOTATIONS)
} else {
  warning("[WARN] Clone annotations file not found: ", CLONE_ANNOTATIONS)
  dict_clones <- NULL
}

# Master sample metadata file
if (!file.exists(DISTALSEQ_MASTER_FILE)) {
  stop("[ERROR] Master metadata file not found: ", DISTALSEQ_MASTER_FILE,
       "\n       See metadata/sample_info/README.md for format specifications.")
}

master_df <- read.xlsx(xlsxFile = DISTALSEQ_MASTER_FILE, sheet = "sample_list")
message("[INFO] Loaded ", nrow(master_df), " samples from master file")

# Load vector annotations (if sheet exists)
if ("vector" %in% getSheetNames(DISTALSEQ_MASTER_FILE)) {
  vector_df <- read.xlsx(xlsxFile = DISTALSEQ_MASTER_FILE, sheet = "vector", startRow = 2)
  master_df <- merge(x = master_df, y = vector_df, by = c("VectorName_BAM"), all.x = TRUE)
  message("[INFO] Merged vector annotations")
}

# Load flanking sequence info (if sheet exists)
if ("IS_Sequence" %in% getSheetNames(DISTALSEQ_MASTER_FILE)) {
  flankseq_df <- read.xlsx(xlsxFile = DISTALSEQ_MASTER_FILE, sheet = "IS_Sequence", startRow = 1)
  master_df <- merge(x = master_df, y = flankseq_df, by = c("CommonPrefixAllFiles"), all.x = TRUE)
  message("[INFO] Merged flanking sequence info")
}

# ==============================================================================
# REFERENCE GENOME FILES
# ==============================================================================

# These are loaded from config.R
target_genome_GTF_file <- REF_GTF_TRANSCRIPT
target_genome_TSS_file <- REF_TSS_BED
target_genome_exons_file <- REF_GTF_EXONS
target_genome_blacklist_mappability_GTF <- REF_BLACKLIST

message("[INFO] Reference genome files:")
message("       GTF:   ", basename(target_genome_GTF_file))
message("       TSS:   ", basename(target_genome_TSS_file))
message("       Exons: ", basename(target_genome_exons_file))

# ==============================================================================
# PUBLIC DATABASES
# ==============================================================================

# Load cancer gene annotations
if (file.exists(MSK_CANCER_GENES)) {
  msk_db <- read.csv(MSK_CANCER_GENES, header = TRUE, fill = TRUE, sep = '\t', check.names = FALSE)
  msk_db$GeneName <- msk_db$`Hugo Symbol`

  # Filter for high-confidence cancer genes
  msk_db <- msk_db[which(
    (msk_db$`Gene Type` %in% c("ONCOGENE", "ONCOGENE_AND_TSG", "TSG")) &
    msk_db$`# of occurrence within resources (Column J-P)` > 1
  ), ]

  message("[INFO] Loaded ", nrow(msk_db), " cancer genes from MSK database")
} else {
  warning("[WARN] MSK cancer gene database not found: ", MSK_CANCER_GENES)
  msk_db <- NULL
}

# ==============================================================================
# TEMPLATE FOR BAM RE-TAGGING (if using clone-specific BAMs)
# ==============================================================================

template_string_retag <- "#!/bin/bash
cd 'TOKEN_FILEFOLDER'

# Generate sub-BAM with specific reads
picard FilterSamReads \\
  I=TOKEN_INPUTRAWBAM \\
  O=TOKEN_OUTPUT_SINGLE_IS_BAM \\
  READ_LIST_FILE=TOKEN_READLIST \\
  FILTER=includeReadList

samtools index TOKEN_OUTPUT_SINGLE_IS_BAM

# Re-tag BAM with clone-specific read groups
picard AddOrReplaceReadGroups \\
  I=TOKEN_OUTPUT_SINGLE_IS_BAM \\
  O=TOKEN_FINAL_SINGLE_IS \\
  RGID='TOKEN_KNOWNIS_ID' \\
  RGLB='TOKEN_CLONEID' \\
  RGPL='TOKEN_SAMPLEID' \\
  RGPU='TOKEN_CLONEID' \\
  RGSM='TOKEN_VECTORID'

samtools index TOKEN_FINAL_SINGLE_IS

# Move outputs to final location
mv TOKEN_OUTPUT_SINGLE_IS_BAM \\
   TOKEN_OUTPUT_SINGLE_IS_BAM.bai \\
   TOKEN_READLIST \\
   TOKEN_FINAL_SINGLE_IS \\
   TOKEN_FINAL_SINGLE_IS.bai \\
   TOKEN_OUTPUTFOLDER
"

# ==============================================================================
# PAF ALIGNMENT COLUMN MAPPING (if using PAF format)
# ==============================================================================

map_pafcols_to_mycols <- data.frame(
  mycols = c("chr", "start", "end", "ratio_bp_aligned_on_raw",
             "query_start", "query_end", "strand", "name", "read_len"),
  pafcols = c("tname", "tstart", "tend", "mapq",
              "qstart", "qend", "strand", "qname", "qlen")
)

# ==============================================================================
# SETUP COMPLETE
# ==============================================================================

message("[INFO] DISTAL-seq initialization complete")
message("[INFO] Ready to process ", sum(master_df$ToProcess, na.rm = TRUE), " samples")
