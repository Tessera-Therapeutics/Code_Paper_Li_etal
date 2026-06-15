#' @author Andrea 
#' @note analysis of DISTAL-seq AND Uni-seq results for the Paper
#' @date March 30, 2026

library(dplyr)
# library(ISAnalytics)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(gplots)
# library(xlsx)
library(ggplot2)
library(scales) 
library(splines)
library(gridExtra)
library(stringr)
# library(sqldf)
# library(plyr)
library(psych)
library(reshape2)
# library(rGREAT)
# library(trackViewer)
# library(GenomicAlignments)
library(openxlsx)
library(Hmisc)
# library(circlize)
library(ggridges)
library(ggpubr)
library(rstatix)
library(ggbreak)
# library(webr)
library(vegan)
library(parallel)
# library(wordcloud)
# library(pafr)
library(ggExtra)

##### =============================================================== #####
##### ------------- Input and global functions ---------------------- #####
##### =============================================================== #####
# Load utility functions
source("distalseq_utils.R")
source("_init_distalseq_analysis.R")
source("isa_utils.R")
source("uniseq_utils.R")
source("_init_uniseq_analysis.R")

source_folder <- "analyses/260331_Paper_DistalUni/"
dir.create(file.path(source_folder), showWarnings = FALSE)
analysis_folder_date <- "260331"
# analysis_prefix <- ".RTE3"
chromosomes_to_use <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")

threshold_IS_span <- 30
pipeline_steps_order <- c("WithoutVector", "WithVector", "VectorOnly", "VectorOnlyPassingP5filter", "VectorOnlyPassingMinLenFilter", "Chimeric", "ChimeraPassingP5filter", "ChimericPassingMinLenFilter", "ChimericPassingProperFilter")
pipeline_steps_order_chimeric <- c("WithVector", "Chimeric", "ChimeraPassingP5filter", "ChimericPassingMinLenFilter", "ChimericPassingProperFilter")

sample_order  <- c("RTE25", "RTE3", "CR1", "Vingi")
sample_labels <- c(
  "RTE25" = "RTE-25_Lmi",
  "RTE3"  = "RTE-3_BF",
  "CR1"   = "CR1-1_PH",
  "Vingi" = "Vingi-1_Acar"
)


##### =============================================================== #####
##### ------------- DISTALSeq                        ---------------- #####
##### =============================================================== #####

##### =============================================================== #####
##### ------------- import data from output Pipeline ---------------- #####
##### =============================================================== #####

# get all results you need to parse
full_df <- NULL
import_only_proper <- F

# data
for (filename in master_df[which(master_df$ToProcess == T), "FileName"]) {
  mat_df <- parseDistalseqPipeTSVOutput(filename = filename, 
                                        master_df = master_df,
                                        import_only_proper = import_only_proper, 
                                        rownames_as_id = F, 
                                        threshold_IS_span = threshold_IS_span)
  if (!is.null(dim(mat_df))) {
    # create a full data df
    if (!is.null(dim(full_df))) {
      #take only valid cols... unfortunately this is needed since several runs showed different col names
      shared_cols <- intersect(colnames(full_df), colnames(mat_df))
      # message(paste0("[AP]\tNumber of shared cols: ", length(shared_cols), ", initially full_df has ", length(colnames(full_df)), " cols, mat_df instead ", length(colnames(mat_df)), "\n\nlist: ", paste(shared_cols, collapse = ' - '), "\n\n"))
      # bind rows
      full_df <- rbind(full_df[shared_cols], mat_df[shared_cols])
    } else {
      full_df <- mat_df
    } # if (length(full_df) > 0)
  } # if (!(dim(mat_df) == NULL))
  # mat_df <- NULL
} # for (filename in master_df[which(master_df$ToProcess == T), "FileName"])

# get filtering steps stats from these reads
filtering_steps_stats <- getPoolStats_FilteringSteps(df = full_df, perc_scale100 = T, incremental_delta_perc = T)



# get all results you need to parse in the vector
vect_df <- NULL
# data
for (filename in master_df[which(master_df$ToProcess == T), "FileName"]) {
  vmat_df <- parseDistalseqPipeTSVOutput_VectorOnly(filename = filename, 
                                                    master_df = master_df,
                                                    rownames_as_id = F)
  if (!is.null(dim(mat_df))) {
    # create a full data df
    if (!is.null(dim(full_df))) {
      # bind rows
      vect_df <- rbind(vect_df, vmat_df)
    } else {
      vect_df <- vmat_df
    } # if (length(full_df) > 0)
  } # if (!(dim(mat_df) == NULL))
} # for (filename in master_df[which(master_df$ToProcess == T), "FileName"])

# get filtering steps stats from these reads
filtering_steps_stats_vector <- getPoolStats_FilteringSteps_VectorOnly(df = vect_df, perc_scale100 = T)


# pipe steps stats
stats_df <- NULL
for (filename in master_df[which(master_df$ToProcess == T), "FileName"][1:8]) {
  message(paste0("[AP]\tRunning the file ", filename))
  # get stats
  mat_stats <- parseDistalseqPipeStatsFiles(filename = filename, master_df = master_df, molten = T, perc_scale100 = T)
  if (!is.null(dim(stats_df))) {
    stats_df <- rbind(stats_df, mat_stats)
  } else {
    stats_df <- mat_stats
  }
  # mat_stats <- NULL
} # for (filename in master_df[which(master_df$ToProcess == T), "FileName"])





###############################################################
# Filter only valid reads
###############################################################
#### query DB ### do summary stats
# full_df_proper_readlenabovemin <- full_df[which(full_df$Proper == "True" & full_df$AboveMinReadLen == T),]
full_df_proper_readlenabovemin <- full_df[which(full_df$Proper == "True" & full_df$AboveMinReadLen == T & full_df$TargetGenomeAlmSize_AboveMinLen == T),]

# write.table(x = full_df_proper_readlenabovemin, 
#             file = gzfile(paste(source_folder, analysis_folder_date, ".Source.full_df_proper_readlenabovemin.tsv.gz", sep = "")), 
#             sep = "\t", quote = FALSE, row.names = FALSE, col.names = T, na = '')


###############################################################
# Find unique target genome IS
###############################################################
full_df_proper_readlenabovemin$integration_locus <- as.numeric(full_df_proper_readlenabovemin$integration_locus)
full_df_proper_readlenabovemin_sorted <- full_df_proper_readlenabovemin %>% 
  arrange(targetRegion_chr, integration_locus, integration_strand)

# get unique IS list with counts
ispileup <- full_df_proper_readlenabovemin_sorted %>% 
  # group_by(targetRegion_chr, integration_locus, targetRegion_strand, .drop = FALSE) %>%
  group_by(targetRegion_chr, integration_locus, integration_strand, .drop = FALSE) %>%
  dplyr::count() 

# find genomic intervals
# # OLD:: is_id_col_sorted <- findGenomicIntervals(df = ispileup[1000:2000,], show_status_bar = T)
# the parallelized version
is_list <- splitDfByChr(df = ispileup, chr_colname = "targetRegion_chr")
cl <- makeCluster(8)  # Crea un cluster con 4 thread
is_id_col_sorted_pres <- parLapply(cl, is_list, findGenomicIntervals)
stopCluster(cl)  # Ferma il cluster quando hai finito
is_id_col_sorted <- do.call(rbind, is_id_col_sorted_pres)


# add all info
ispileup_ext <- cbind(ispileup, is_id_col_sorted)
# get peak (reference IS) by interval
ispileup_ext_stats <- ispileup_ext %>% 
  group_by(RefID) %>%
  arrange(desc(n), .by_group = TRUE) %>%
  mutate(RefID_maxspan = max(RefID_span), RefID_endspan = max(RefID_end), RefID_peakReads = max(n), RefID_nStrand = n_distinct(RefID_strand)) %>% 
  filter(row_number()==1)
ispileup_ext_stats <- ispileup_ext_stats %>% 
  arrange(targetRegion_chr, integration_locus, integration_strand)

# TODO: IS stats such as: max range, len, etc.

# write the bed files to run bedtools
# ispileup_ext_stats$FooStrand <- "+"
write.table(x = ispileup_ext_stats[c("Ref_chr", "integration_locus", "integration_locus", "RefID", "RefID_span", "RefID_strand")], 
            file = paste(source_folder, analysis_folder_date, ".is_refmax.bed", sep = ""), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = F, na = '')

write.table(x = full_df_proper_readlenabovemin_sorted[c("targetRegion_chr", "integration_locus", "integration_locus", "Read_Name", "Num_TotalAln", "targetRegion_strand")], 
            file = paste(source_folder, analysis_folder_date, ".all_reads.bed", sep = ""), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = F, na = '')

# run bedtools such this:
# bedtools closest -b 230907.is_refmax.bed -a 230907.all_reads.bed -D a -t first > 230907.closestRefIS.bed
# bedtools closest -b ../../metadata/hg38/hg38.ncbiRefSeq.sorted.slice_exon.gtf -a 230907.all_reads.bed -D a -t first > 230907.all_reads.annohg38.ucsc.bed ## NOT bedtools closest -b ../../metadata/hg38/ucsc.hg38.refGene.2309.exons.gtf -a 230907.all_reads.bed -D a -t first > 230907.all_reads.annohg38.ucsc.bed 
system(command = paste0("bedtools closest -b '", 
                        paste(source_folder, analysis_folder_date, ".is_refmax.bed", sep = ""), 
                        "' -a ", 
                        paste(source_folder, analysis_folder_date, ".all_reads.bed", sep = ""), 
                        " -D a -t first > ", paste(source_folder, analysis_folder_date, ".closestRefIS.bed", sep = "")))

system(command = paste0("bedtools closest -b '", 
                        target_genome_GTF_file, 
                        "' -a ", 
                        paste(source_folder, analysis_folder_date, ".is_refmax.bed", sep = ""), 
                        " -D a -t first > ", paste(source_folder, analysis_folder_date, ".is_refmax.annohg38.ucsc.bed", sep = "")))

# read the annotation of all the reads
full_df_proper_readlenabovemin_annotatedIS <- read.csv(file = paste(source_folder, analysis_folder_date, ".closestRefIS.bed", sep = ""), 
                                                       header=F, fill=T, sep='\t', check.names = FALSE, 
                                                       na.strings = c("NONE", "NA", "NULL", "NaN", "ND", ""))
names(full_df_proper_readlenabovemin_annotatedIS) <- c("targetRegion_chr", "integration_locus", "human_IS", "Read_Name", "Num_TotalAln", "targetRegion_strand", "refIS_chr", "refIS_locus", "refIS_end", "RefID", "RefID_v2", "refIS_strand", "DistanceToRef")
# rownames(full_df_proper_readlenabovemin_annotatedIS) <- full_df_proper_readlenabovemin_annotatedIS$Read_Name

# combine annotatinos with previous data
full_df_proper_readlenabovemin_annotatedIS_ext <- cbind(full_df_proper_readlenabovemin_sorted, 
                                                        full_df_proper_readlenabovemin_annotatedIS[setdiff(colnames(full_df_proper_readlenabovemin_annotatedIS), c("targetRegion_chr", "integration_locus", "human_IS", "Read_Name", "Num_TotalAln", "targetRegion_strand", "RefID_v2", "RefID_v3", "feature_dot", "feature_dot2"))])


# read annotation of unique IS
ispileup_ext_stats_exon <- read.csv(file = paste(source_folder, analysis_folder_date, ".is_refmax.annohg38.ucsc.bed", sep = ""), 
                                    header=F, fill=T, sep='\t', check.names = FALSE, 
                                    na.strings = c("NONE", "NA", "NULL", "NaN", "ND", ""))
names(ispileup_ext_stats_exon) <- c("Ref_chr", "IS_base", "integration_locus_copy", "RefID", "RefID_span", "IS_strand", 
                                    "feature_chr", "feature_source", "feature_type", "feature_start", "feature_end", "feature_dot", "feature_strand", "feature_dot2", "feature_string", "feature_distance")
ispileup_ext_stats_exon$feature_geneName <- apply(ispileup_ext_stats_exon[c("feature_string", "feature_distance")], 1, function(x) {
  strsplit( strsplit(x[1], ';', fixed = T)[[1]][1] , ' ', fixed = T)[[1]][2]
} )

# combine all (not merge since now IDs are not unique for a pipe problem)
full_df_proper_readlenabovemin_annotatedIS_full <- merge(x = full_df_proper_readlenabovemin_annotatedIS_ext,
                                                         y = ispileup_ext_stats_exon,
                                                         by = c("RefID"), all.x = T)

# add a single read count and write this data
full_df_proper_readlenabovemin_annotatedIS_full$SC <- 1
write.table(x = full_df_proper_readlenabovemin_annotatedIS_full, 
            file = gzfile(paste(source_folder, analysis_folder_date, ".full_df_proper_readlenabovemin_annotatedIS_full.tsv.gz", sep = "")), 
            sep = "\t", quote = FALSE, row.names = T, col.names = T, na = '')

full_df_proper_readlenabovemin_annotatedIS_full <- classifyCompleteIntegrationByIS(df = full_df_proper_readlenabovemin_annotatedIS_full, 
                                                                                   group_by_cols = c("RefID", "Complete"),
                                                                                   complete_colname = "Complete", 
                                                                                   complete_positive_val = "True", 
                                                                                   complete_negative_val = "False", 
                                                                                   output_plot_file_prefix = paste(source_folder, analysis_folder_date, ".QC.Complete_vs_Truncated.MixedISCases", sep = ""),
                                                                                   height=7, width=9)

IS_matrix_D <- dcast(data = full_df_proper_readlenabovemin_annotatedIS_full, 
                   # refIS_chr + refIS_locus + integration_strand + feature_geneName + feature_strand + feature_distance ~ SampleID, 
                   refIS_chr + refIS_locus + integration_strand + feature_geneName + feature_strand + feature_distance + CompleteIntegration ~ SampleID, 
                   value.var = "SC", fun.aggregate = sum)
# adjust colnames
# id_cols <- c("chr", "integration_locus", "integration_strand", "GeneName", "GeneStrand", "GeneDistance")
id_cols <- c("chr", "integration_locus", "integration_strand", "GeneName", "GeneStrand", "GeneDistance", "CompleteIntegration")
old_id_cols <- c("refIS_chr", "refIS_locus", "integration_strand", "feature_geneName", "feature_strand", "feature_distance", "CompleteIntegration")
old_id_cols %in% colnames(IS_matrix_D)
names(IS_matrix_D) <- c(id_cols, setdiff(colnames(IS_matrix_D), old_id_cols))
IS_matrix_D <- rownamesAsIS(df = IS_matrix_D)

# keep only some chrs
IS_matrix_D <- IS_matrix_D[which(IS_matrix_D$chr %in% chromosomes_to_use),]
# write the matrix
write.table(x = IS_matrix_D, 
            file = paste(source_folder, analysis_folder_date, ".Dis", ".IS_matrix_D.tsv", sep = ""), 
            sep = "\t", quote = FALSE, row.names = T, col.names = T, na = '')
# do stats
IS_matrix_D_dataonly <- IS_matrix_D[setdiff(colnames(IS_matrix_D), id_cols)]
IS_matrix_D_stats <- data.frame(
  "nIS" = apply(IS_matrix_D_dataonly, 2, function(x) {length(x[x>0])}),
  "nReads" = apply(IS_matrix_D_dataonly, 2, function(x) {sum(x, na.rm = T)}),
  "Hindex" = diversity(t(IS_matrix_D_dataonly))
)
IS_matrix_D_stats$SampleID <- rownames(IS_matrix_D_stats)
IS_matrix_D_stats_meta <- merge(x = IS_matrix_D_stats, y = master_df[c("GroupName", "SampleID", "SampleName", "Replica", "PoolID", "Transgene", "SampleType")], by = c("SampleID"), all.x = T)
rownames(IS_matrix_D_stats_meta) <- IS_matrix_D_stats_meta$SampleID
# IS_matrix_D_stats_meta <- merge(x = IS_matrix_D_stats_meta, y = stats_df[c("Total_Reads", "Vector_Reads", "TargetGenome_Reads", "PassingP5filter_Reads", "SampleID")], by = c("SampleID"), all.x = T)
# rownames(IS_matrix_D_stats_meta) <- IS_matrix_D_stats_meta$SampleID
write.xlsx(x = IS_matrix_D_stats_meta, 
           file = paste(source_folder, analysis_folder_date, ".Dis", ".DescriptiveStats.IS.stats.xlsx", sep = ""), 
           rowNames = T)

# master_df$gdf <- apply(master_df[c("SampleName", "SampleType")], 1, function(x) {paste0(x[1], "_", x[2], collapse = "_")})
master_df$gdf <- master_df$Vector
slice_master_metadata <- master_df[which(master_df$ToProcess == TRUE),]
rownames(slice_master_metadata) <- slice_master_metadata$SampleID
IS_gdf_dataonly <- aggregateDfColumnsByName(df = IS_matrix_D, 
                                            metadata_df = slice_master_metadata, 
                                            key_field = c("gdf"),
                                            # key_field = c("Vector"),
                                            # key_field = c("SampleName"), 
                                            starting_data_col_index = length(id_cols)+1, 
                                            number_of_last_cols_to_remove = 0)
IS_gdf <- cbind(IS_matrix_D[id_cols], IS_gdf_dataonly)
# check. IS sharing
IS_gdf_sharing <- as.data.frame(getSharedISnumber(df = IS_gdf_dataonly, compact = T, left_to_rigth_reading_output = T))
rownames(IS_gdf_sharing) <- colnames(IS_gdf_dataonly)
write.xlsx(x = IS_gdf_sharing, 
           file = paste(source_folder, analysis_folder_date, ".Dis", ".DescriptiveStats.IS_gdf.sharing_N.xlsx", sep = ""), 
           rowNames = T)
# do stats
IS_gdf_stats <- data.frame(
  "nIS" = apply(IS_gdf_dataonly, 2, function(x) {length(x[x>0])}),
  "nReads" = apply(IS_gdf_dataonly, 2, function(x) {sum(x, na.rm = T)}),
  "Hindex" = diversity(t(IS_gdf_dataonly))
)
IS_gdf_stats$SampleID <- rownames(IS_gdf_stats)
IS_gdf_no0 <- IS_gdf
IS_gdf_no0[IS_gdf_no0==0] <- NA
IS_gdf_molten <- melt(data = IS_gdf_no0, 
                      id.vars = id_cols,
                      variable.name = "gdf", 
                      na.rm = T, 
                      value.name = "SC")
write.xlsx(x = IS_gdf_stats, 
           file = paste(source_folder, analysis_folder_date, ".Dis", ".DescriptiveStats.IS_gdf.stats.xlsx", sep = ""), 
           rowNames = T)

# compute top genes
topGenes <- getTopHitGenes(df = IS_gdf, 
                           data_column_to_use = setdiff(colnames(IS_gdf), id_cols), 
                           group_by_colums = "GeneName", 
                           maxGenesToReturn = 50)
write.xlsx(x = topGenes, 
           file = paste(source_folder, analysis_folder_date, ".Dis", ".DescriptiveStats.IS_gdf.topGenes.xlsx", sep = ""), 
           rowNames = F)


##### =============================================================== #####
##### ------------- Uniseq                           ---------------- #####
##### =============================================================== #####

##### =============================================================== #####
##### ------------- import data from output Pipeline ---------------- #####
##### =============================================================== #####

# get all results you need to parse
uni_full_df <- NULL
# data
for (datafolder in uniseq_master_df[which(uniseq_master_df$ToProcess == T), "RootFolder"]) {
  mat_df <- readUniseqBed(rootfolder_id = datafolder, uniseq_master_df = uniseq_master_df)
  if (!is.null(dim(mat_df))) {
    # create a full data df
    if (!is.null(dim(uni_full_df))) {
      uni_full_df <- rbind(uni_full_df, mat_df)
    } else {
      uni_full_df <- mat_df
    } # if (length(full_df) > 0)
  } # if (!(dim(mat_df) == NULL))
} # for (filename in master_df[which(master_df$ToProcess == T), "FileName"])


###############################################################
# Find unique target genome IS
###############################################################
# keep only some chrs
uni_full_df <- uni_full_df[which(uni_full_df$chr %in% chromosomes_to_use),]
uni_IS_matrix <- dcast(data = uni_full_df, 
                   chr + integration_locus + integration_strand ~ SampleID,
                   value.var = "SC", fun.aggregate = sum)
rownames(uni_IS_matrix) <- apply(uni_IS_matrix[c("chr", "integration_locus", "integration_strand")], 1, function(x) {
  paste0(x[1], "_", as.character(as.numeric(x[2])), "_", x[3])})
# create a new metadata file
sample_metadata <- unique(uni_full_df[c("SampleID", "SampleName", "Vector", "PoolID")])
rownames(sample_metadata) <- sample_metadata$SampleID
# annotate IS
write.table(x = uni_IS_matrix[c("chr", "integration_locus", "integration_locus", "integration_locus", "integration_locus", "integration_strand")], 
            file = paste(source_folder, analysis_folder_date, ".Uni", ".uni_IS_matrix_perSample.bed", sep = ""), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = F, na = '')
system(command = paste0("bedtools closest -b '", 
                        target_genome_GTF_file, 
                        "' -a ", 
                        paste(source_folder, analysis_folder_date, ".Uni", ".uni_IS_matrix_perSample.bed", sep = ""), 
                        " -D a -t first > ", paste(source_folder, analysis_folder_date, ".Uni", ".uni_IS_matrix_perSample.annohg38.ucsc.bed", sep = "")))
# read annotation of unique IS
uni_IS_matrix_bedannotation <- read.csv(file = paste(source_folder, analysis_folder_date, ".Uni", ".uni_IS_matrix_perSample.annohg38.ucsc.bed", sep = ""), 
                                    header=F, fill=T, sep='\t', check.names = FALSE, 
                                    na.strings = c("NONE", "NA", "NULL", "NaN", "ND", ""))
names(uni_IS_matrix_bedannotation) <- c("chr", "integration_locus", "integration_locus_copy", "integration_locus_copy2", "integration_locus_copy3",
                                    "integration_strand", "feature_chr", "feature_source", "feature_type", "feature_start", "feature_end", "feature_dot", "GeneStrand", "feature_dot2", "feature_string", "GeneDistance")
rownames(uni_IS_matrix_bedannotation) <- apply(uni_IS_matrix_bedannotation[c("chr", "integration_locus", "integration_strand")], 1, function(x) {
  paste0(x[1], "_", as.character(as.numeric(x[2])), "_", x[3])})
uni_IS_matrix_bedannotation$GeneName <- apply(uni_IS_matrix_bedannotation[c("feature_string", "GeneDistance")], 1, function(x) {
  strsplit( strsplit(x[1], ';', fixed = T)[[1]][1] , ' ', fixed = T)[[1]][2]
} )

# combine all (not merge since now IDs are not unique for a pipe problem)
uni_IS_matrix_annotated <- merge(x = uni_IS_matrix, 
                             y = uni_IS_matrix_bedannotation[c("GeneName", "GeneStrand", "GeneDistance")], by = 0)
rownames(uni_IS_matrix_annotated) <- uni_IS_matrix_annotated$Row.names
uni_IS_matrix_annotated <- uni_IS_matrix_annotated[-1]
uni_IS_matrix_annotated <- uni_IS_matrix_annotated[c(uni_id_cols, setdiff(colnames(uni_IS_matrix_annotated), id_cols))]

# do stats by sample
uni_IS_matrix_dataonly <- uni_IS_matrix_annotated[setdiff(colnames(uni_IS_matrix_annotated), id_cols)]
uni_IS_matrix_stats <- data.frame(
  "nIS" = apply(uni_IS_matrix_dataonly, 2, function(x) {length(x[x>0])}),
  "nReads" = apply(uni_IS_matrix_dataonly, 2, function(x) {sum(x, na.rm = T)}),
  "Hindex" = diversity(t(uni_IS_matrix_dataonly))
)
uni_IS_matrix_stats <- merge(x = uni_IS_matrix_stats, 
                         y = sample_metadata, 
                         by = 0)
rownames(uni_IS_matrix_stats) <- uni_IS_matrix_stats$Row.names
uni_IS_matrix_stats <- uni_IS_matrix_stats[-1]
write.xlsx(x = uni_IS_matrix_stats, 
           file = paste(source_folder, analysis_folder_date, ".Uni", ".DescriptiveStats.uni_IS_matrix_bySample.stats.xlsx", sep = ""), 
           rowNames = T)

# aggregate samples by vector
uni_IS_gdf_dataonly <- aggregateDfColumnsByName(df = uni_IS_matrix_annotated, 
                                            metadata_df = sample_metadata, 
                                            key_field = c("Vector"),
                                            # key_field = c("SampleName"), 
                                            starting_data_col_index = length(id_cols)+1, 
                                            number_of_last_cols_to_remove = 0)
uni_IS_gdf <- cbind(uni_IS_matrix_annotated[uni_id_cols], uni_IS_gdf_dataonly)
uni_IS_gdf <- rownamesAsIS(df = uni_IS_gdf, id_cols = uni_id_cols)
write.table(x = uni_IS_gdf, 
            file = paste(source_folder, analysis_folder_date, ".Uni", ".IS_gdf.tsv", sep = ""), 
            sep = "\t", quote = FALSE, row.names = F, col.names = T, na = '')


# check. IS sharing
uni_IS_gdf_sharing <- as.data.frame(getSharedISnumber(df = uni_IS_gdf_dataonly, compact = T, left_to_rigth_reading_output = T))
rownames(uni_IS_gdf_sharing) <- colnames(uni_IS_gdf_dataonly)
write.xlsx(x = uni_IS_gdf_sharing, 
           file = paste(source_folder, analysis_folder_date, ".Uni", ".DescriptiveStats.IS_gdf.sharing_N.xlsx", sep = ""), 
           rowNames = T)
# and plot sharing
pheatmap(uni_IS_gdf_sharing, scale = "none", 
         height = 4.5, width = 6, 
         display_numbers = T, fontsize_number = 7, 
         number_format = "%d",
         filename = paste(source_folder, analysis_folder_date, ".DescriptiveStats.IS.shared.numbers.png", sep = ""),
         main = paste0("Shared ISs among samples"),
         cluster_rows=FALSE, cluster_cols=FALSE,
         color=colorRampPalette(c("white", "orange", "red"))(50)
         # annotation_row = uni_IS_matrix_stats_meta[c("SampleType", "GroupName", "SampleName")]
         # cutree_rows = 2
         # cutree_cols = 2
)


# do stats
uni_IS_gdf_stats <- data.frame(
  "nIS" = apply(uni_IS_gdf_dataonly, 2, function(x) {length(x[x>0])}),
  "nReads" = apply(uni_IS_gdf_dataonly, 2, function(x) {sum(x, na.rm = T)}),
  "Hindex" = diversity(t(uni_IS_gdf_dataonly))
)
uni_IS_gdf_stats$SampleID <- rownames(uni_IS_gdf_stats)
write.xlsx(x = uni_IS_gdf_stats, 
           file = paste(source_folder, analysis_folder_date, ".Uni", ".DescriptiveStats.IS_gdf.stats.xlsx", sep = ""), 
           rowNames = T)
# melt data
uni_IS_gdf_no0 <- uni_IS_gdf
uni_IS_gdf_no0[uni_IS_gdf_no0==0] <- NA
uni_IS_gdf_molten <- melt(data = uni_IS_gdf_no0, 
                      id.vars = uni_id_cols,
                      variable.name = "SampleID", 
                      na.rm = T, 
                      value.name = "SC")
uni_IS_gdf_molten$GeneDistance <- ifelse(is.na(uni_IS_gdf_molten$GeneDistance), 0, uni_IS_gdf_molten$GeneDistance)
write.table(x = uni_IS_gdf_molten, 
            file = gzfile(paste(source_folder, analysis_folder_date, ".Uni", ".uni_IS_gdf_molten.tsv.gz", sep = "")), 
            sep = "\t", quote = FALSE, row.names = F, col.names = T, na = '')


## allIS df
uni_IS_gdf_molten$Assay <- "Uniseq"
uni_IS_gdf_molten$CompleteIntegration <- NA
allIS <- uni_IS_gdf_molten


# check IS distance ---- !
uni_IS_gdf_distance <- uni_IS_gdf
uni_IS_gdf_distance <- uni_IS_gdf_distance %>% 
  arrange(chr, integration_locus)
delta <- c(0)
for (i in seq(2, nrow(uni_IS_gdf_distance))) {
  ifelse(uni_IS_gdf_distance[i,"chr"] == uni_IS_gdf_distance[i-1, "chr"], 
         delta <- c(delta, uni_IS_gdf_distance[i, "integration_locus"] - uni_IS_gdf_distance[i-1, "integration_locus"] ),
         delta <- c(delta, 0))
}
uni_IS_gdf_distance$delta <- delta
uni_IS_refIS <- findGenomicIntervals(df = uni_IS_gdf_distance, threshold_IS_span = 30, show_status_bar = T, in_chr_colname = "chr", in_locus_colname = "integration_locus", in_strand_colname = "integration_strand")
# evaluate the bias
# dim(uni_IS_gdf_distance[which(uni_IS_gdf_distance$delta < 5),])
dim(uni_IS_gdf_distance[which(uni_IS_gdf_distance$delta < threshold_IS_span),])
dim(uni_IS_gdf_distance[which(!(uni_IS_gdf_distance$integration_strand %in% c("+", "-"))),])






IS_gdf_molten$Assay <- "Distalseq"
IS_gdf_molten$SampleID <- IS_gdf_molten$gdf
IS_gdf_molten$SampleID <- gsub("_Pos", "", IS_gdf_molten$SampleID)



uniseq_allIS <- allIS
# add these ISs in the allIS matrix
allIS <- rbind(allIS, IS_gdf_molten[colnames(allIS)])
allIS <- allIS[which(allIS$SampleID %in% sample_order),]

allIS$GeneDistance <- ifelse(is.na(allIS$GeneDistance), 0, allIS$GeneDistance) # double check this
# filter weird cased like chr3:93470483 (chr3_93470483_-_RNU6-488P)
allIS <- allIS[which(allIS$integration_strand %in% c("+", "-")),]

write.table(x = allIS, 
            file = gzfile(paste(source_folder, analysis_folder_date, ".AllIS_combined", ".molten.tsv.gz", sep = "")), 
            sep = "\t", quote = FALSE, row.names = F, col.names = T, na = '')


# Add all the annotations
# allIS <- read.csv(paste(source_folder, analysis_folder_date, ".AllIS_combined", ".molten.tsv.gz", sep = ""), 
#                   header=TRUE, fill=T, sep='\t', check.names = FALSE)

allIS_mat <- dcast(data = allIS, 
                   # chr + integration_locus + integration_strand + GeneName + GeneStrand + GeneDistance ~ SampleID + Assay, ## with GeneDistance DOES NOT WORK! to be fized
                   chr + integration_locus + integration_strand ~ SampleID + Assay,
                   value.var = "SC", fun.aggregate = mean)


allIS_fullanno <- 
  annotateISMatrix(df = allIS_mat, 
                   cols_chr_start_end_strand = c("chr", "integration_locus", "integration_locus", "integration_locus", "integration_locus", "integration_strand"), 
                   id_cols = c("chr", "integration_locus", "integration_strand"), 
                   features_to_annotate_gtf_bed = target_genome_TSS_file, 
                   df_bedfile_towrite = paste(source_folder, analysis_folder_date, ".IS_matrix.TSS.bed", sep = ""),
                   df_bedfile_annotated = paste(source_folder, analysis_folder_date, ".IS_matrix.annotated.TSS.bed", sep = ""),
                   feature_output_names = c("TSSName", "TSSStrand", "TSSDistance"), ref_bed_format = T
  )
allIS_fullanno <- 
  annotateISMatrix(df = allIS_fullanno, 
                   cols_chr_start_end_strand = c("chr", "integration_locus", "integration_locus", "integration_locus", "integration_locus", "integration_strand"), 
                   id_cols = c("chr", "integration_locus", "integration_strand"), 
                   features_to_annotate_gtf_bed = target_genome_exons_file, 
                   df_bedfile_towrite = paste(source_folder, analysis_folder_date, ".IS_matrix.exon.bed", sep = ""),
                   df_bedfile_annotated = paste(source_folder, analysis_folder_date, ".IS_matrix.annotated.exon.bed", sep = ""),
                   feature_output_names = c("ExonName", "ExonStrand", "ExonDistance"), ref_bed_format = F
  )
allIS_fullanno <- 
  annotateISMatrix(df = allIS_fullanno, 
                   cols_chr_start_end_strand = c("chr", "integration_locus", "integration_locus", "integration_locus", "integration_locus", "integration_strand"), 
                   id_cols = c("chr", "integration_locus", "integration_strand"), 
                   features_to_annotate_gtf_bed = target_genome_GTF_file, 
                   df_bedfile_towrite = paste(source_folder, analysis_folder_date, ".IS_matrix.transc.bed", sep = ""),
                   df_bedfile_annotated = paste(source_folder, analysis_folder_date, ".IS_matrix.annotated.transc.bed", sep = ""), ref_bed_format = F
  )

# and heir classes
allIS_fullanno <- uniqueAnnotation(df = allIS_fullanno)

allIS_fullanno_molten <- melt(data = allIS_fullanno, 
                          id.vars = c("chr", "integration_locus", "integration_strand", "GeneName", "GeneStrand", "GeneDistance", "ExonName", "ExonStrand", "ExonDistance", "TSSName", "TSSStrand", "TSSDistance", "AnnotationClass"),
                          variable.name = "SampleID", 
                          na.rm = T, 
                          value.name = "SC")

allIS_fullanno_molten$Assay <- sub(".*_(Distalseq|Uniseq)$", "\\1", allIS_fullanno_molten$SampleID)
allIS_fullanno_molten$SampleAssayID <- allIS_fullanno_molten$SampleID
allIS_fullanno_molten$SampleID <- sub("_(Distalseq|Uniseq)$", "", allIS_fullanno_molten$SampleAssayID)

write.table(x = allIS_fullanno_molten, 
            file = gzfile(paste(source_folder, analysis_folder_date, ".AllIS_combined.annotated.molten.tsv.gz", sep = "")), 
            sep = "\t", quote = FALSE, row.names = F, col.names = T, na = '')

# make annotation label stats
allIS_fullanno_molten_summarylabel <- allIS_fullanno_molten %>% 
  dplyr::group_by(Assay, SampleID, AnnotationClass) %>% 
  dplyr::mutate(nIS = n()) %>%
  get_summary_stats(nIS, type = "full") 
allIS_fullanno_molten_summarylabel <- allIS_fullanno_molten_summarylabel %>% 
  dplyr::group_by(Assay, SampleID) %>% 
  dplyr::mutate(SumIS = sum(n))
allIS_fullanno_molten_summarylabel$GroupPerc <- allIS_fullanno_molten_summarylabel$n / allIS_fullanno_molten_summarylabel$SumIS
write.xlsx(x = allIS_fullanno_molten_summarylabel, 
           file = paste(source_folder, analysis_folder_date, ".AllIS_combined.annotated.molten.labeled.summary.xlsx", sep = ""), 
           rowNames = T)


# visual check
plot_stackedbar_annotationlabel_byassay <- 
  ggplot(allIS_fullanno_molten_summarylabel, 
         aes(x = SampleID, y = GroupPerc, color = AnnotationClass, fill = AnnotationClass)
  ) +
  # scale_fill_manual(values=wes_palette(n=4, name="Royal1")) +
  # scale_color_manual(values=wes_palette(n=4, name="Royal1")) +
  scale_color_manual(values = c(t(paper_barplot_colors[1,]))) +
  scale_fill_manual(values = c(t(paper_barplot_colors[1,]))) +
  geom_bar(position = position_stack(reverse = TRUE), stat = "identity") +
  # stat_summary(fun.data = n_fun, geom = "text", hjust = 0.9, vjust = 0.5, size = 5, angle = 90) +
  # geom_jitter(size = 1, alpha = 0.3) +
  facet_wrap(. ~ Assay, labeller = label_wrap_gen(width=6), scales = "free_x") +
  scale_y_continuous(labels = scales::percent) +
  geom_text(
    aes(label = n),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    # size = 4,
    color = "white", 
    angle = 0
  ) +
  theme_bw() +
  theme(strip.text.y = element_text(size = 14, colour = "blue", angle = 0)) +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(strip.text.x = element_text(size = 14, colour = "darkblue", angle = 0), 
        strip.text.y = element_text(size = 14, colour = "darkred", angle = 270)) +
  # theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=14, angle = 45, hjust=0.95,vjust=0.95), 
        axis.text.y = element_text(size=16), 
        axis.title = element_text(size=16), plot.title = element_text(size=22)) +
  labs(title = paste0("IS annotation"), 
       x = "Sample", y = "Annotated features [%]", color = "Feature", fill = "Feature",
       subtitle = paste0("UCSC RefGene.") ) 


plot_stackedbar_annotationlabel_byassay
pdf(file = paste(source_folder, analysis_folder_date, ".AllIS", ".FeatureLabel.stackedbar.pdf", sep = ""), height=6, width=7)
plot(plot_stackedbar_annotationlabel_byassay)
dev.off()
png(file = paste(source_folder, analysis_folder_date, ".AllIS", ".FeatureLabel.stackedbar.png", sep = ""), height=6, width=7, units = "in", res = 300)
plot(plot_stackedbar_annotationlabel_byassay)
dev.off()


# -------------------------------------------------------
# Colors
# -------------------------------------------------------
annot_colors <- c(
  "TSS_10000"  = "#E63946",
  "Exon"       = "#F4A261",
  "Intron"     = "#457B9D",
  "Intergenic" = "#A8DADC"
)

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(openxlsx)
library(irr)
library(patchwork)


# -------------------------------------------------------
# Fisher's exact test with equalized N per SampleID
# sample min(n_Distalseq, n_Uniseq) from each assay
# -------------------------------------------------------
set.seed(42)
n_boot   <- 100
n_sample <- 2000

# Compute min N per SampleID
n_per_sample <- allIS_fullanno_molten %>%
  group_by(SampleID, Assay) %>%
  summarise(n_total = n(), .groups = "drop") %>%
  pivot_wider(names_from = Assay, values_from = n_total) %>%
  mutate(n_equalized = pmin(Distalseq, Uniseq))

message("[AP]\tEqualized N per SampleID:")
print(n_per_sample)

# Subsample IS to equalized N and compute counts per AnnotationClass
counts_eq <- do.call(rbind, lapply(unique(allIS_fullanno_molten$SampleID), function(sid) {
  n_eq <- n_per_sample[n_per_sample$SampleID == sid, "n_equalized"]$n_equalized
  do.call(rbind, lapply(c("Distalseq", "Uniseq"), function(assay) {
    sub_df  <- allIS_fullanno_molten[allIS_fullanno_molten$SampleID == sid &
                                       allIS_fullanno_molten$Assay    == assay, ]
    sampled <- sub_df[sample(nrow(sub_df), size = n_eq, replace = FALSE), ]
    sampled %>%
      group_by(AnnotationClass) %>%
      summarise(n_IS = n(), .groups = "drop") %>%
      mutate(SampleID = sid, Assay = assay, n_total = n_eq)
  }))
}))

# Run Fisher per SampleID x AnnotationClass on equalized counts
fisher_eq <- do.call(rbind, lapply(unique(counts_eq$SampleID), function(sid) {
  do.call(rbind, lapply(unique(counts_eq$AnnotationClass), function(cat) {
    dis <- counts_eq[counts_eq$SampleID == sid & counts_eq$Assay == "Distalseq", ]
    uni <- counts_eq[counts_eq$SampleID == sid & counts_eq$Assay == "Uniseq",    ]
    
    n_cat_dis  <- dis[dis$AnnotationClass == cat, "n_IS"]$n_IS
    n_cat_uni  <- uni[uni$AnnotationClass == cat, "n_IS"]$n_IS
    n_rest_dis <- unique(dis$n_total) - n_cat_dis
    n_rest_uni <- unique(uni$n_total) - n_cat_uni
    
    mat <- matrix(c(n_cat_dis, n_rest_dis,
                    n_cat_uni, n_rest_uni),
                  nrow = 2,
                  dimnames = list(c(cat, "Other"),
                                  c("Distalseq", "Uniseq")))
    ft <- fisher.test(mat)
    
    data.frame(
      SampleID        = sid,
      AnnotationClass = cat,
      n_eq            = unique(dis$n_total),
      n_cat_Distalseq = n_cat_dis,
      n_cat_Uniseq    = n_cat_uni,
      prop_Distalseq  = n_cat_dis / unique(dis$n_total),
      prop_Uniseq     = n_cat_uni / unique(uni$n_total),
      OR              = ft$estimate,
      log2OR          = log2(ft$estimate),
      CI_lower        = ft$conf.int[1],
      CI_upper        = ft$conf.int[2],
      p_value         = ft$p.value
    )
  }))
}))

rownames(fisher_eq) <- NULL
fisher_eq$p_adj     <- p.adjust(fisher_eq$p_value, method = "bonferroni")

# OR threshold: flag as biologically relevant only if |log2(OR)| > or_threshold
# OR > 1.5 or OR < 0.67 -> log2(OR) > 0.58
# OR > 1.25 or OR < 0.80 -> log2(OR) > 0.32
or_threshold <- 0.32  # corresponds to OR > 1.25 or OR < 0.80

fisher_eq$sig_label <- ifelse(fisher_eq$p_adj < 0.05 & abs(fisher_eq$log2OR) > or_threshold,
                              ifelse(fisher_eq$p_adj < 0.001, "***",
                                     ifelse(fisher_eq$p_adj < 0.01,  "**", "*")),
                              "ns")

fisher_eq$AnnotationClass <- factor(fisher_eq$AnnotationClass, levels = names(annot_colors))
fisher_eq$SampleID        <- factor(fisher_eq$SampleID, levels = c("CR1", "RTE3", "RTE25", "Vingi"))

message("[AP]\tEqualized Fisher results:")
print(fisher_eq)

# -------------------------------------------------------
# Forest plot (with sample order and labels)
# -------------------------------------------------------
fisher_eq$SampleID <- factor(fisher_eq$SampleID,
                             levels  = rev(sample_order),
                             labels  = rev(sample_labels))

plot_forest_eq <- ggplot(fisher_eq,
                         aes(x = OR, y = SampleID)) +
                         # aes(x = OR, y = SampleID, color = AnnotationClass)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.8) +
  # geom_vline(xintercept = c(1/1.25, 1.25), linetype = "dotted",
  #            color = "orange", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper),
                 height = 0.25, linewidth = 0.8) +
  geom_point(size = 3.5) +
  geom_text(aes(label = sig_label, x = CI_upper),
            hjust = -0.3, size = 4, fontface = "bold", color = "black") +
  scale_color_manual(values = annot_colors, guide = "none") +
  scale_x_log10(breaks = c(0.5, 0.67, 0.8, 1, 1.25, 1.5, 2), 
                limits = c(0.6, 1.6)) +
  # facet_wrap(~ AnnotationClass, nrow = 1, scales = "free_x") +
  facet_wrap(AnnotationClass ~ ., ncol = 1) +
  labs(
       # title    = "Fisher's exact test: Distalseq vs Uniseq (equalized N)",
       # subtitle = paste0("Significance requires p_adj < 0.05 AND |log2(OR)| > ",
       #                   or_threshold, " (OR > 1.25 or OR < 0.80).\nOrange lines = OR thresholds"),
       x        = "Odds Ratio (log scale)",
       y        = NULL) +
  theme_bw(base_size = 13) +
  theme(plot.title    = element_text(size = 15, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey40"),
        strip.text    = element_text(size = 12, face = "bold"),
        axis.text.y   = element_text(size = 12))

for (ext in c("pdf", "png")) {
  ggsave(filename = paste0(source_folder, analysis_folder_date,
                           ".Comparison.DistalUni.GenomicAnnotation.FisherForest.Equalized.", ext),
         plot = plot_forest_eq, width = 13, height = 4,
         dpi = 300, device = ext)
}
for (ext in c("pdf", "png", "svg")) {
  
  device_fun <- switch(
    ext,
    pdf = cairo_pdf,
    png = "png",
    svg = svglite::svglite
  )
  
  ggsave(
    filename = paste0(
      source_folder,
      analysis_folder_date,
      ".Comparison.DistalUni.GenomicAnnotation.FisherForest.Equalized.vertical.",
      ext
    ),
    plot = plot_forest_eq,
    width = 3.5,
    height = 7,
    dpi = 300,
    device = device_fun
  )
}

# -------------------------------------------------------
# Save xlsx
# -------------------------------------------------------
wb <- createWorkbook()

addWorksheet(wb, "Fisher_Equalized")
writeDataTable(wb, "Fisher_Equalized", fisher_eq, rowNames = FALSE)

addWorksheet(wb, "N_Equalized")
writeDataTable(wb, "N_Equalized", n_per_sample, rowNames = FALSE)

saveWorkbook(wb, file = paste0(source_folder, analysis_folder_date,
                               ".Comparison.DistalUni.GenomicAnnotation.Fisher.Equalized.xlsx"), overwrite = TRUE)

message("[AP]\tEqualized Fisher analysis complete. All outputs saved.")










