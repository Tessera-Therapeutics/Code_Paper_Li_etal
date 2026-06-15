#' @author Andrea 
#' @note analysis of UNI-seq results for the paper Revision in Cell
#' @date 24 Feb 2026
#' 

library(patchwork)
library(simplermarkdown)
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
library(circlize)
library(ggridges)
library(ggpubr)
library(rstatix)
library(ggbreak)
library(webr)
library(vegan)
# library(parallel)
# library(wordcloud)
library(ggplot2)
library(ggseqlogo)
library(wesanderson)

##### =============================================================== #####
##### ------------- Input and global functions ---------------------- #####
##### =============================================================== #####
# Load utility functions and initialization
source("../distalseq/uniseq_utils.R")
source("../distalseq/_init_uniseq_analysis.R")

source_folder <- "analyses/260223/"
dir.create(file.path(source_folder), showWarnings = FALSE)
analysis_folder_date <- "260223"
# analysis_prefix <- ".RTE3"
chromosomes_to_use <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
N_topgenes <- 50
significance_threshold_minus_log_p <- -log(0.05, base = 10)

subsampling_n_is <- 25000 # to be used for random and subsampling

# faceting_order_vector <- c("RTE1", "Lenti", "Random")
# RTE25, RTE3, CR1, Vingi
faceting_order_vector <- c("Random", c("RTE25", "RTE3", "CR1", "Vingi"), "Lenti")
# faceting_order_vector_extended <- c("Random", "RTE-1_MD", "Vingi-1_Acar", "RTE-3_BF", "RTE-25_Lmi", "CR1-1_PH", "Lenti")
faceting_order_vector_extended <- c("Random", "Vingi-1_Acar", "RTE-3_BF", "RTE-25_Lmi", "CR1-1_PH", "Lenti")
samples_withlowreads <- c("Vingi_gDNA2065_3site_S12",
                          "D21_3_2_S8",
                          "RTE3_rep1_12_S12",
                          "RTE3_rep2_12_S12") # post-hoc knowledge, here overimposed



##### =============================================================== #####
##### ------------- import data from output Pipeline ---------------- #####
##### =============================================================== #####
# write metadata file 
write.xlsx(x = uniseq_master_df[which(uniseq_master_df$ToProcess == T), ], 
           file = paste(source_folder, analysis_folder_date, ".Metadata.TRUE.xlsx", sep = ""), 
           rowNames = T)

# get all results you need to parse
full_df <- NULL
# data
for (datafolder in uniseq_master_df[which(uniseq_master_df$ToProcess == T), "RootFolder"]) {
  message(paste("[AP]\tProcessing:", datafolder))
  mat_df <- readUniseqBed(rootfolder_id = datafolder, uniseq_master_df = uniseq_master_df)
  if (!is.null(dim(mat_df))) {
    # create a full data df
    if (!is.null(dim(full_df))) {
      full_df <- rbind(full_df, mat_df)
    } else {
      full_df <- mat_df
    } # if (length(full_df) > 0)
  } # if (!(dim(mat_df) == NULL))
} # for (filename in master_df[which(master_df$ToProcess == T), "FileName"])

full_df <- full_df[which(!(full_df$SampleID %in% samples_withlowreads)),]

###############################################################
# Find unique target genome IS
###############################################################
# keep only some chrs
full_df <- full_df[which(full_df$chr %in% chromosomes_to_use),]
IS_matrix <- dcast(data = full_df, 
                   chr + integration_locus + integration_strand ~ SampleID,
                   value.var = "SC", fun.aggregate = sum)
rownames(IS_matrix) <- apply(IS_matrix[c("chr", "integration_locus", "integration_strand")], 1, function(x) {
  paste0(x[1], "_", as.character(as.numeric(x[2])), "_", x[3])})
# create a new metadata file
sample_metadata <- unique(full_df[c("SampleID", "SampleName", "Vector", "PoolID")])
rownames(sample_metadata) <- sample_metadata$SampleID

# annotate IS
IS_matrix_annotated_all <- 
  annotateISMatrix(df = IS_matrix, 
                   cols_chr_start_end_strand = c("chr", "integration_locus", "integration_locus", "integration_locus", "integration_locus", "integration_strand"), 
                   id_cols = c("chr", "integration_locus", "integration_strand"), 
                   features_to_annotate_gtf_bed = target_genome_blacklist_mappability_GTF, 
                   df_bedfile_towrite = paste(source_folder, analysis_folder_date, ".IS_matrix.bed", sep = ""),
                   df_bedfile_annotated = paste(source_folder, analysis_folder_date, ".IS_matrix.annotated.bed", sep = ""), 
                   feature_output_names = c("MappabilityName", "MappabilityStrand", "MappabilityDistance")
  )
IS_matrix_annotated_all <- 
  annotateISMatrix(df = IS_matrix_annotated_all, 
                   cols_chr_start_end_strand = c("chr", "integration_locus", "integration_locus", "integration_locus", "integration_locus", "integration_strand"), 
                   id_cols = c("chr", "integration_locus", "integration_strand"), 
                   features_to_annotate_gtf_bed = target_genome_TSS_file, 
                   df_bedfile_towrite = paste(source_folder, analysis_folder_date, ".IS_matrix.TSS.bed", sep = ""),
                   df_bedfile_annotated = paste(source_folder, analysis_folder_date, ".IS_matrix.annotated.TSS.bed", sep = ""),
                   feature_output_names = c("TSSName", "TSSStrand", "TSSDistance"), ref_bed_format = T
  )
IS_matrix_annotated_all <- 
  annotateISMatrix(df = IS_matrix_annotated_all, 
                   cols_chr_start_end_strand = c("chr", "integration_locus", "integration_locus", "integration_locus", "integration_locus", "integration_strand"), 
                   id_cols = c("chr", "integration_locus", "integration_strand"), 
                   features_to_annotate_gtf_bed = target_genome_exons_file, 
                   df_bedfile_towrite = paste(source_folder, analysis_folder_date, ".IS_matrix.exon.bed", sep = ""),
                   df_bedfile_annotated = paste(source_folder, analysis_folder_date, ".IS_matrix.annotated.exon.bed", sep = ""),
                   feature_output_names = c("ExonName", "ExonStrand", "ExonDistance"), ref_bed_format = F
  )
IS_matrix_annotated_all <- 
  annotateISMatrix(df = IS_matrix_annotated_all, 
                   cols_chr_start_end_strand = c("chr", "integration_locus", "integration_locus", "integration_locus", "integration_locus", "integration_strand"), 
                   id_cols = c("chr", "integration_locus", "integration_strand"), 
                   features_to_annotate_gtf_bed = target_genome_GTF_file, 
                   df_bedfile_towrite = paste(source_folder, analysis_folder_date, ".IS_matrix.transc.bed", sep = ""),
                   df_bedfile_annotated = paste(source_folder, analysis_folder_date, ".IS_matrix.annotated.transc.bed", sep = ""), ref_bed_format = F
  )

# filter out IS in black list regions
IS_matrix_annotated <- IS_matrix_annotated_all[which(IS_matrix_annotated_all$MappabilityDistance > 0 | IS_matrix_annotated_all$MappabilityDistance < 0),]
# IS_matrix_annotated <- IS_matrix_annotated[setdiff(colnames(IS_matrix_annotated), c("MappabilityName", "MappabilityStrand", "MappabilityDistance"))]
write.xlsx(x = IS_matrix_annotated, 
           file = paste(source_folder, analysis_folder_date, ".IS_matrix.filtered.xlsx", sep = ""), 
           rowNames = T)

# label univocally each annotated IS
IS_matrix_annotated_labeled <- uniqueAnnotation(df = IS_matrix_annotated)
write.xlsx(x = IS_matrix_annotated_labeled, 
           file = paste(source_folder, analysis_folder_date, ".IS_matrix.filtered.labeled.xlsx", sep = ""), 
           rowNames = T)
# annotate cols with vars
uni_id_cols_fulllabels <- c(uni_id_cols, "TSSName", "TSSStrand", "TSSDistance", "ExonName", "ExonStrand", "ExonDistance", "MappabilityName", "MappabilityStrand", "MappabilityDistance", "AnnotationClass")
uni_data_cols <- setdiff(colnames(IS_matrix_annotated_labeled), uni_id_cols_fulllabels)

# test plots
plot_hist_density <- 
  ggplot(IS_matrix_annotated, aes(x = TSSDistance) ) +
  geom_histogram(bins = 500) + 
  scale_x_continuous(limits = c(-50000, 50000))

# reshape data
IS_matrix_annotated_labeled_cast <- melt(data = IS_matrix_annotated_labeled[c(uni_id_cols, "AnnotationClass", uni_data_cols)], 
                                         id.vars = c(uni_id_cols, "AnnotationClass"),
                                         variable.name = "SampleID", 
                                         na.rm = T, 
                                         value.name = "SC")
IS_matrix_annotated_labeled_cast <- IS_matrix_annotated_labeled_cast[which(IS_matrix_annotated_labeled_cast$SC>0),]
IS_matrix_annotated_labeled_cast_summarylabel <- IS_matrix_annotated_labeled_cast %>% 
  dplyr::group_by(SampleID, AnnotationClass) %>% 
  dplyr::mutate(nIS = n()) %>%
  get_summary_stats(nIS, type = "full") 
IS_matrix_annotated_labeled_cast_summarylabel <- IS_matrix_annotated_labeled_cast_summarylabel %>% 
  dplyr::group_by(SampleID) %>% 
  dplyr::mutate(SumIS = sum(n))
IS_matrix_annotated_labeled_cast_summarylabel$GroupPerc <- IS_matrix_annotated_labeled_cast_summarylabel$n / IS_matrix_annotated_labeled_cast_summarylabel$SumIS
# IS_matrix_annotated_labeled_cast_summarylabel <- merge(x = IS_matrix_annotated_labeled_cast_summarylabel, 
#                                                        y = uniseq_master_df, 
#                                                        by = c("SampleID"))
write.xlsx(x = IS_matrix_annotated_labeled_cast_summarylabel, 
           file = paste(source_folder, analysis_folder_date, ".IS_matrix.filtered.labeled.summary.xlsx", sep = ""), 
           rowNames = T)

plot_stackedbar_annotationlabel <- 
  ggplot(IS_matrix_annotated_labeled_cast_summarylabel, 
         aes(x = SampleID, y = GroupPerc, color = AnnotationClass, fill = AnnotationClass)
  ) +
  # scale_fill_manual(values=wes_palette(n=4, name="Royal1")) +
  # scale_color_manual(values=wes_palette(n=4, name="Royal1")) +
  scale_color_manual(values = c(t(paper_barplot_colors[1,]))) +
  scale_fill_manual(values = c(t(paper_barplot_colors[1,]))) +
  geom_bar(position = position_stack(reverse = TRUE), stat = "identity") +
  # stat_summary(fun.data = n_fun, geom = "text", hjust = 0.9, vjust = 0.5, size = 5, angle = 90) +
  # geom_jitter(size = 1, alpha = 0.3) +
  # facet_wrap(. ~ SampleID, labeller = label_wrap_gen(width=6), scales = "free_x") +
  scale_y_continuous(labels = scales::percent) +
  geom_text(
    aes(label = n),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    # size = 4,
    color = "white", 
    angle = 90
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
  labs(title = paste0("Uni-seq - Per sample IS annotation"), 
       x = "Sample", y = "Annotated features [%]", color = "Feature", fill = "Feature",
       subtitle = paste0("UCSC RefGene.") ) 

pdf(file = paste(source_folder, analysis_folder_date, ".Uni", ".FeatureLabel.stackedbar.pdf", sep = ""), height=6, width=19)
plot(plot_stackedbar_annotationlabel)
dev.off()
png(file = paste(source_folder, analysis_folder_date, ".Uni", ".FeatureLabel.stackedbar.png", sep = ""), height=6, width=19, units = "in", res = 300)
plot(plot_stackedbar_annotationlabel)
dev.off()


# do stats by sample
IS_matrix_dataonly <- IS_matrix_annotated[setdiff(uni_data_cols, "AnnotationClass")]
# IS_matrix_dataonly <- IS_matrix_annotated[setdiff(colnames(IS_matrix_annotated), uni_id_cols)]
IS_matrix_stats <- data.frame(
  "nIS" = apply(IS_matrix_dataonly, 2, function(x) {length(x[x>0])}),
  "nReads" = apply(IS_matrix_dataonly, 2, function(x) {sum(x, na.rm = T)}),
  "Hindex" = diversity(t(IS_matrix_dataonly))
)
IS_matrix_stats <- merge(x = IS_matrix_stats, 
                         y = sample_metadata, 
                         by = 0)
rownames(IS_matrix_stats) <- IS_matrix_stats$Row.names
IS_matrix_stats <- IS_matrix_stats[-1]
write.xlsx(x = IS_matrix_stats, 
           file = paste(source_folder, analysis_folder_date, ".DescriptiveStats.IS_matrix_bySample.stats.xlsx", sep = ""), 
           rowNames = T)

# we identify the samples to filter out
samples_unreliable <- rownames(IS_matrix_stats[which(IS_matrix_stats$nReads < 100), ])

# Vingi_gDNA2065_3site_S12
# D21_3_2_S8
# RTE3_rep1_12_S12
# RTE3_rep2_12_S12

# aggregate samples by vector
IS_gdf_dataonly <- aggregateDfColumnsByName(df = IS_matrix_annotated,
# IS_gdf_dataonly <- aggregateDfColumnsByName(df = IS_matrix_annotated[grep("LV|RTE1", colnames(IS_matrix_annotated), value = T)], 
                                            metadata_df = sample_metadata, 
                                            key_field = c("Vector"),
                                            starting_data_col_index = length(uni_id_cols)+1,
                                            # starting_data_col_index = 1, 
                                            number_of_last_cols_to_remove = 0)
IS_gdf <- cbind(IS_matrix_annotated[uni_id_cols], IS_gdf_dataonly)
IS_gdf <- rownamesAsIS(df = IS_gdf, id_cols = uni_id_cols)
write.table(x = IS_gdf, 
            file = paste(source_folder, analysis_folder_date, ".IS_gdf.tsv", sep = ""), 
            sep = "\t", quote = FALSE, row.names = F, col.names = T, na = '')


# check. IS sharing
IS_gdf_sharing <- as.data.frame(getSharedISnumber(df = IS_gdf_dataonly, compact = T, left_to_rigth_reading_output = T))
rownames(IS_gdf_sharing) <- colnames(IS_gdf_dataonly)
write.xlsx(x = IS_gdf_sharing, 
           file = paste(source_folder, analysis_folder_date, ".DescriptiveStats.IS_gdf.sharing_N.xlsx", sep = ""), 
           rowNames = T)
# and plot sharing
pheatmap(IS_gdf_sharing, scale = "none", 
         height = 4.5, width = 6, 
         display_numbers = T, fontsize_number = 7, 
         number_format = "%d",
         filename = paste(source_folder, analysis_folder_date, ".DescriptiveStats.IS.shared.numbers.png", sep = ""),
         main = paste0("Shared ISs among samples"),
         cluster_rows=FALSE, cluster_cols=FALSE,
         color=colorRampPalette(c("white", "orange", "red"))(50)
         # annotation_row = IS_matrix_stats_meta[c("SampleType", "GroupName", "SampleName")]
         # cutree_rows = 2
         # cutree_cols = 2
)

# do stats
IS_gdf_stats <- data.frame(
  "nIS" = apply(IS_gdf_dataonly, 2, function(x) {length(x[x>0])}),
  "nReads" = apply(IS_gdf_dataonly, 2, function(x) {sum(x, na.rm = T)}),
  "Hindex" = diversity(t(IS_gdf_dataonly))
)
IS_gdf_stats$SampleID <- rownames(IS_gdf_stats)
write.xlsx(x = IS_gdf_stats, 
           file = paste(source_folder, analysis_folder_date, ".DescriptiveStats.IS_gdf.stats.xlsx", sep = ""), 
           rowNames = T)
# melt data
IS_gdf_no0 <- IS_gdf
IS_gdf_no0[IS_gdf_no0==0] <- NA
IS_gdf_molten <- melt(data = IS_gdf_no0, 
                      id.vars = uni_id_cols,
                      variable.name = "SampleID", 
                      na.rm = T, 
                      value.name = "SC")
IS_gdf_molten$GeneDistance <- ifelse(is.na(IS_gdf_molten$GeneDistance), 0, IS_gdf_molten$GeneDistance)
write.table(x = IS_gdf_molten, 
            file = gzfile(paste(source_folder, analysis_folder_date, ".IS_gdf_molten.tsv.gz", sep = "")), 
            sep = "\t", quote = FALSE, row.names = F, col.names = T, na = '')

# compute top genes
N_topgenes <- 20
# vector_topGenes_bySample <- getTopHitGenes(df = IS_matrix_annotated, 
#                            data_column_to_use = setdiff(colnames(IS_matrix_annotated), uni_id_cols), 
#                            group_by_colums = "GeneName", 
#                            maxGenesToReturn = N_topgenes)

vector_topGenes_byVector <- getTopHitGenes(df = IS_gdf, 
                           data_column_to_use = setdiff(colnames(IS_gdf), uni_id_cols), 
                           group_by_colums = "GeneName", 
                           maxGenesToReturn = N_topgenes)

# vector_topGenes_bySample_list <- vector_topGenes_bySample[grep("Gene", colnames(vector_topGenes_bySample), value = T)]
vector_topGenes_byVector_list <- vector_topGenes_byVector[grep("Gene", colnames(vector_topGenes_byVector), value = T)]
# write.xlsx(x = vector_topGenes_bySample_list, 
#            file = paste(source_folder, analysis_folder_date, ".Top", as.character(N_topgenes), "Genes.UpSet.bySample.xlsx", sep = ""), 
#            rowNames = T)
write.xlsx(x = vector_topGenes_byVector_list, 
           file = paste(source_folder, analysis_folder_date, ".Top", as.character(N_topgenes), "Genes.UpSet.byVector.xlsx", sep = ""), 
           rowNames = F)

###############################################################
## BED FILES
###############################################################
# subsampling_n_is <- min(IS_gdf_stats$nIS) - 50
# subsampling_n_is <- 3500
human_cytoband <- read.cytoband(species = "hg38")$df
cytoband_rbind <- rbind(human_cytoband)
cytoband <- read.cytoband(cytoband_rbind)
cytoband_df = cytoband$df
chromosome = cytoband$chromosome
xrange = cytoband$chr.len
isa_mat_becols <- c("chr", "integration_locus", "integration_locus", "GeneName", "SC", "integration_strand")
bed_colnames <- c("chr","start","end","value1", "score", "strand")

###############################################################
## FEATURE ANNOTATION AND RANDOM 
###############################################################
# bed_random <- generateRandomBed(nr = subsampling_n_is, species = "hg38", fun = function(k) sample(letters, k, replace = TRUE))
bed_random <- generateRandomBed(nr = (subsampling_n_is*1.5), species = "hg38", fun = function(k) sample(letters, k, replace = TRUE))
bed_random$end <- bed_random$start
bed_random$score <- sample(seq(1, 30), nrow(bed_random), replace = T)
bed_random$strand <- sample(c("+", "-"), nrow(bed_random), replace = T)
# names(bed_random) <- bed_colnames
write.table(x = bed_random, 
            file = paste(source_folder, analysis_folder_date, ".IS_matrix.subsampling_", subsampling_n_is, ".Random.bed", sep = ""), 
            sep = "\t", col.names = F, row.names = F, na = '', quote = F)
# annotate IS
system(command = paste0("bedtools sort -i '", 
                        paste(source_folder, analysis_folder_date, ".IS_matrix.subsampling_", subsampling_n_is, ".Random.bed", sep = ""), "' > '",
                        paste(source_folder, analysis_folder_date, ".IS_matrix.subsampling_", subsampling_n_is, ".Random.sorted.bed", sep = ""), "'") )
system(command = paste0("bedtools closest -b '", 
                        target_genome_GTF_file, 
                        "' -a ", 
                        paste(source_folder, analysis_folder_date, ".IS_matrix.subsampling_", subsampling_n_is, ".Random.sorted.bed", sep = ""), 
                        " -D ref -t first > ", paste(source_folder, analysis_folder_date, ".IS_matrix.subsampling_", subsampling_n_is, ".Random.sorted.annohg38.ucsc.bed", sep = "") ))
# read annotation of unique IS
IS_random <- read.csv(file = paste(source_folder, analysis_folder_date, ".IS_matrix.subsampling_", subsampling_n_is, ".Random.sorted.annohg38.ucsc.bed", sep = ""),
                      header=F, fill=T, sep='\t', check.names = FALSE, 
                      na.strings = c("NONE", "NA", "NULL", "NaN", "ND", ""))
names(IS_random) <- c("chr", "integration_locus", "integration_locus_copy", "integration_locus_copy2", "integration_locus_copy3",
                                    "integration_strand", "feature_chr", "feature_source", "feature_type", "feature_start", "feature_end", "feature_dot", "GeneStrand", "feature_dot2", "feature_string", "GeneDistance")
rownames(IS_random) <- apply(IS_random[c("chr", "integration_locus", "integration_strand")], 1, function(x) {
  paste0(x[1], "_", as.character(as.numeric(x[2])), "_", x[3])})
IS_random$GeneName <- apply(IS_random[c("feature_string", "GeneDistance")], 1, function(x) {
  strsplit( strsplit(x[1], ';', fixed = T)[[1]][1] , ' ', fixed = T)[[1]][2]
} )

# allIS_random <- IS_random[uni_id_cols]
allIS_random <- IS_random[c("chr", "integration_locus", "integration_strand")]
# annotate IS
allIS_random_withmappab <- 
  annotateISMatrix(df = allIS_random, 
                   cols_chr_start_end_strand = c("chr", "integration_locus", "integration_locus", "integration_locus", "integration_locus", "integration_strand"), 
                   id_cols = c("chr", "integration_locus", "integration_strand"), 
                   features_to_annotate_gtf_bed = target_genome_blacklist_mappability_GTF, 
                   df_bedfile_towrite = paste(source_folder, analysis_folder_date, ".IS_randomM.bed", sep = ""),
                   df_bedfile_annotated = paste(source_folder, analysis_folder_date, ".IS_randomM.annotated.bed", sep = ""), 
                   feature_output_names = c("MappabilityName", "MappabilityStrand", "MappabilityDistance")
  )
# allIS_random <- allIS_random_withmappab[which((allIS_random_withmappab$MappabilityDistance > 0 | allIS_random_withmappab$MappabilityDistance < 0)), uni_id_cols]
allIS_random <- allIS_random_withmappab[which((allIS_random_withmappab$MappabilityDistance > 0 | allIS_random_withmappab$MappabilityDistance < 0)), ]

# add the other annotations
allIS_random_all <- 
  annotateISMatrix(df = allIS_random, 
                   cols_chr_start_end_strand = c("chr", "integration_locus", "integration_locus", "integration_locus", "integration_locus", "integration_strand"), 
                   id_cols = c("chr", "integration_locus", "integration_strand"), 
                   features_to_annotate_gtf_bed = target_genome_TSS_file, 
                   df_bedfile_towrite = paste(source_folder, analysis_folder_date, ".IS_random.TSS.bed", sep = ""),
                   df_bedfile_annotated = paste(source_folder, analysis_folder_date, ".IS_random.annotated.TSS.bed", sep = ""),
                   feature_output_names = c("TSSName", "TSSStrand", "TSSDistance"), ref_bed_format = T
  )
allIS_random_all <- 
  annotateISMatrix(df = allIS_random_all, 
                   cols_chr_start_end_strand = c("chr", "integration_locus", "integration_locus", "integration_locus", "integration_locus", "integration_strand"), 
                   id_cols = c("chr", "integration_locus", "integration_strand"), 
                   features_to_annotate_gtf_bed = target_genome_exons_file, 
                   df_bedfile_towrite = paste(source_folder, analysis_folder_date, ".IS_random.exon.bed", sep = ""),
                   df_bedfile_annotated = paste(source_folder, analysis_folder_date, ".IS_random.annotated.exon.bed", sep = ""),
                   feature_output_names = c("ExonName", "ExonStrand", "ExonDistance"), ref_bed_format = F
  )
allIS_random_all <- 
  annotateISMatrix(df = allIS_random_all, 
                   cols_chr_start_end_strand = c("chr", "integration_locus", "integration_locus", "integration_locus", "integration_locus", "integration_strand"), 
                   id_cols = c("chr", "integration_locus", "integration_strand"), 
                   features_to_annotate_gtf_bed = target_genome_GTF_file, 
                   df_bedfile_towrite = paste(source_folder, analysis_folder_date, ".IS_random.transc.bed", sep = ""),
                   df_bedfile_annotated = paste(source_folder, analysis_folder_date, ".IS_random.annotated.transc.bed", sep = ""), ref_bed_format = F
  )
# need to reorder the cols!!! precise ordering!
allIS_random_all <- allIS_random_all[uni_id_cols_fulllabels[1:(length(uni_id_cols_fulllabels)-1)]]
allIS_random_all$SampleID <- "Random"
allIS_random_all$SC <- 1
allIS_random_all <- allIS_random_all[which(!(allIS_random_all$GeneName %in% c("LOC124903223"))),]
allIS_random_all_labeled <- uniqueAnnotation(df = allIS_random_all)
write.xlsx(x = allIS_random_all_labeled, 
           file = paste(source_folder, analysis_folder_date, ".IS_random.filtered.labeled.xlsx", sep = ""), 
           rowNames = T)

allIS_random_all_labeled_summarylabel <- allIS_random_all_labeled %>% 
  dplyr::group_by(SampleID, AnnotationClass) %>% 
  dplyr::mutate(nIS = n()) %>%
  get_summary_stats(nIS, type = "full") 
allIS_random_all_labeled_summarylabel <- allIS_random_all_labeled_summarylabel %>% 
  dplyr::group_by(SampleID) %>% 
  dplyr::mutate(SumIS = sum(n))
allIS_random_all_labeled_summarylabel$GroupPerc <- allIS_random_all_labeled_summarylabel$n / allIS_random_all_labeled_summarylabel$SumIS
allIS_random_all_labeled_summarylabel$Vector <- "Random"
write.xlsx(x = allIS_random_all_labeled_summarylabel, 
           file = paste(source_folder, analysis_folder_date, ".IS_random.filtered.labeled.summary.xlsx", sep = ""), 
           rowNames = T)

# merge info
IS_matrix_annotated_labeled_cast_summarylabel <- merge(x = IS_matrix_annotated_labeled_cast_summarylabel, y = IS_matrix_stats[c("SampleID", "Vector")], by = c("SampleID"))
# combine all data
shared_cols <- intersect(colnames(IS_matrix_annotated_labeled_cast_summarylabel), colnames(allIS_random_all_labeled_summarylabel))
all_is_labeled_summary <- rbind(IS_matrix_annotated_labeled_cast_summarylabel[shared_cols], allIS_random_all_labeled_summarylabel[shared_cols])
# plot it by sample
plot_stackedbar_annotationlabel_all <- 
  ggplot(all_is_labeled_summary, 
         aes(x = SampleID, y = GroupPerc, color = AnnotationClass, fill = AnnotationClass)
  ) +
  geom_bar(position = position_stack(reverse = TRUE), stat = "identity") +
  # stat_summary(fun.data = n_fun, geom = "text", hjust = 0.9, vjust = 0.5, size = 5, angle = 90) +
  # geom_jitter(size = 1, alpha = 0.3) +
  facet_grid(. ~ Vector, labeller = label_wrap_gen(width=6), space = "free_x", scales = "free") +
  scale_y_continuous(labels = scales::percent) +
  geom_text(
    aes(label = n),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    # size = 4,
    color = "white", 
    angle = 90
  ) + 
  theme_bw() +
  theme(strip.text.y = element_text(size = 14, colour = "blue", angle = 0)) +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(strip.text.x = element_text(size = 14, colour = "darkblue", angle = 0), 
        strip.text.y = element_text(size = 14, colour = "darkred", angle = 270)) +
  # theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=10, angle = 45, hjust=0.95,vjust=0.95), 
        axis.text.y = element_text(size=16), 
        axis.title = element_text(size=16), plot.title = element_text(size=22)) +
  labs(title = paste0("Uni-seq - Per sample IS annotation"), 
       x = "Sample", y = "Annotated features [%]", color = "Feature", fill = "Feature",
       subtitle = paste0("UCSC RefGene.") ) 

pdf(file = paste(source_folder, analysis_folder_date, ".Uni", ".FeatureLabel.stackedbar.all.pdf", sep = ""), height=6, width=20)
plot(plot_stackedbar_annotationlabel_all)
dev.off()
png(file = paste(source_folder, analysis_folder_date, ".Uni", ".FeatureLabel.stackedbar.all.png", sep = ""), height=6, width=20, units = "in", res = 300)
plot(plot_stackedbar_annotationlabel_all)
dev.off()

# do stats avg per vector
all_is_labeled_summary_avg <- all_is_labeled_summary %>% 
  dplyr::group_by(Vector, AnnotationClass) %>% 
  get_summary_stats(GroupPerc, type = "full") 
all_is_labeled_summary_avg$AnnotationClass <- ifelse(all_is_labeled_summary_avg$AnnotationClass == "TSS_10000", "Upstr.TSS", all_is_labeled_summary_avg$AnnotationClass)
# all_is_labeled_summary_avg$Vector <- ifelse(all_is_labeled_summary_avg$Vector == "RTE1", "GW", all_is_labeled_summary_avg$Vector)
# all_is_labeled_summary_avg$Vector <- ifelse(all_is_labeled_summary_avg$Vector == "Lenti", "LV", all_is_labeled_summary_avg$Vector)
all_is_labeled_summary_avg$Vector_Order <- factor(all_is_labeled_summary_avg$Vector, levels = faceting_order_vector)

# faceting_order_ucscannotation <- c("Upstr.TSS", "Exon", "Intron", "Intergenic")
faceting_order_ucscannotation <- c("Intergenic", "Upstr.TSS", "Exon", "Intron")
all_is_labeled_summary_avg$AnnotationClass_Order <- factor(all_is_labeled_summary_avg$AnnotationClass, levels = faceting_order_ucscannotation)
# # faceting_order_vector <- c("RTE1", "Lenti", "Random")
# faceting_order_vector <- c("GW", "Random", "LV")
# all_is_labeled_summary_avg$Vector_Order <- factor(all_is_labeled_summary_avg$Vector, levels = faceting_order_vector)
all_is_labeled_summary_avg <- all_is_labeled_summary_avg[which(!is.na(all_is_labeled_summary_avg$Vector_Order)),]
all_is_labeled_summary_avg <- merge(x = all_is_labeled_summary_avg, y = vector_rename, by = "Vector", all.x=T)
faceting_order_vector_extended <- c("Random", "CR1-1_PH", "RTE-3_BF", "RTE-25_Lmi", "Vingi-1_Acar", "Lenti")
all_is_labeled_summary_avg$Vector_Order_Extended <- factor(all_is_labeled_summary_avg$ExtendedVector, levels = faceting_order_vector_extended)

write.xlsx(x = all_is_labeled_summary_avg, 
           file = paste(source_folder, analysis_folder_date, ".Uni", ".FeatureLabel.stackedbar.avgVector.xlsx", sep = ""), 
           rowNames = T)

# plot it by vector
plot_stackedbar_annotationlabel_avg <- 
  ggplot(all_is_labeled_summary_avg, 
         aes(x = Vector_Order_Extended, y = median, color = AnnotationClass_Order, fill = AnnotationClass_Order)
  ) +
  geom_bar(position = position_stack(reverse = TRUE), stat = "identity") +
  scale_color_manual(values = c(t(tessera_color_schema[1,]))) +
  scale_fill_manual(values = c(t(tessera_color_schema[1,]))) +
  # scale_color_manual(values = c(t(paper_barplot_colors[1,]))) +
  # scale_fill_manual(values = c(t(paper_barplot_colors[1,]))) +
  # stat_summary(fun.data = n_fun, geom = "text", hjust = 0.9, vjust = 0.5, size = 5, angle = 90) +
  # geom_jitter(size = 1, alpha = 0.3) +
  # facet_grid(. ~ Vector, labeller = label_wrap_gen(width=6), space = "free_x", scales = "free") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(strip.text.y = element_text(size = 14, colour = "blue", angle = 0)) +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(strip.text.x = element_text(size = 14, colour = "darkblue", angle = 0), 
        strip.text.y = element_text(size = 14, colour = "darkred", angle = 270)) +
  # theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(legend.direction = "vertical", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=12, 
                                   hjust=0.95,vjust=0.95,
                                   angle = 30), 
        axis.text.y = element_text(size=16), 
        axis.title = element_text(size=16), plot.title = element_text(size=22)) +
  labs(title = paste0("IS annotation"), 
       x = "Sample", y = "Annotated features [%]", color = "Feature", fill = "Feature",
       subtitle = paste0("UCSC RefGene.") ) 

pdf(file = paste(source_folder, analysis_folder_date, ".Uni", ".FeatureLabel.stackedbar.avgVector.pdf", sep = ""), height=5.8, width=4)
plot(plot_stackedbar_annotationlabel_avg)
dev.off()
png(file = paste(source_folder, analysis_folder_date, ".Uni", ".FeatureLabel.stackedbar.avgVector.png", sep = ""), height=5.8, width=4, units = "in", res = 300)
plot(plot_stackedbar_annotationlabel_avg)
dev.off()



###############################################################
## ANNOTATIONS
###############################################################
#### feature annotations
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggupset)
library(ggimage)

#### RTE1
peak <- readPeakFile(peakfile = paste(source_folder, analysis_folder_date, ".IS_matrix.subsampling_", subsampling_n_is, ".RTE1.bed", sep = ""))
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)
# upset plot
peakAnno <- annotatePeak(peak = paste(source_folder, analysis_folder_date, ".IS_matrix.subsampling_", subsampling_n_is, ".RTE1.bed", sep = ""), 
                         tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
# vennpie(peakAnno)
# upsetplot(peakAnno)
# upsetplot(peakAnno, vennpie=TRUE)
# plot
pdf(file = paste(source_folder, analysis_folder_date, ".IS_matrix.subsampling_", subsampling_n_is, ".RTE1.features.pdf", sep = ""), height=5, width=8)
upsetplot(peakAnno)
upsetplot(peakAnno, vennpie=TRUE) # + grid.text("RTE1",x = 0.1, y=0.1, gp=gpar(fontsize=20))
vennpie(peakAnno)
dev.off()

#### LV
peak <- readPeakFile(peakfile = paste(source_folder, analysis_folder_date, ".IS_matrix.subsampling_", subsampling_n_is, ".LV.bed", sep = ""))
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)
# upset plot
peakAnno <- annotatePeak(peak = paste(source_folder, analysis_folder_date, ".IS_matrix.subsampling_", subsampling_n_is, ".LV.bed", sep = ""), 
                         tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
# plot
pdf(file = paste(source_folder, analysis_folder_date, ".IS_matrix.subsampling_", subsampling_n_is, ".LV.features.pdf", sep = ""), height=5, width=8)
upsetplot(peakAnno)
upsetplot(peakAnno, vennpie=TRUE) # + grid.text("RTE1",x = 0.1, y=0.1, gp=gpar(fontsize=20))
vennpie(peakAnno)
dev.off()

#### Random
peak <- readPeakFile(peakfile = paste(source_folder, analysis_folder_date, ".IS_matrix.subsampling_", subsampling_n_is, ".Random.bed", sep = ""))
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)
# upset plot
peakAnno <- annotatePeak(peak = paste(source_folder, analysis_folder_date, ".IS_matrix.subsampling_", subsampling_n_is, ".Random.bed", sep = ""), 
                         tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
# plot
pdf(file = paste(source_folder, analysis_folder_date, ".IS_matrix.subsampling_", subsampling_n_is, ".Random.features.pdf", sep = ""), height=5, width=8)
upsetplot(peakAnno)
upsetplot(peakAnno, vennpie=TRUE) # + grid.text("RTE1",x = 0.1, y=0.1, gp=gpar(fontsize=20))
vennpie(peakAnno)
dev.off()


###############################################################
## Exon counts
###############################################################
IS_gdf_molten_withrandom <- rbind(IS_gdf_molten, allIS_random_all[colnames(IS_gdf_molten)])
IS_gdf_molten_withrandom_oncotsg <- merge(x = IS_gdf_molten_withrandom, y = msk_db, by = c("GeneName"), all.x = T)
IS_gdf_molten_withrandom_oncotsg$KnownGeneClass <- ifelse(!is.na(IS_gdf_molten_withrandom_oncotsg$`Gene Type`), 
                                        (ifelse((IS_gdf_molten_withrandom_oncotsg$`Gene Type`) %in% c("ONCOGENE", "ONCOGENE_AND_TSG"), 
                                                "OncoGene", "TumSuppressor")), 
                                        "Other")


IS_gdf_molten_withrandom_oncotsg_desc <- IS_gdf_molten_withrandom_oncotsg %>%
  dplyr::group_by(SampleID, KnownGeneClass) %>%
  dplyr::summarise(nIS = n(), .groups = "drop") %>%
  dplyr::group_by(SampleID) %>%
  dplyr::mutate(freq = nIS / sum(nIS) * 100) 

IS_gdf_molten_withrandom_oncotsg$KnownGeneClassBinary <- ifelse(!is.na(IS_gdf_molten_withrandom_oncotsg$`Gene Type`), 
                                                                "OG-TSG", 
                                                                "Other")
IS_gdf_molten_withrandom_oncotsg_desc_simplified <- IS_gdf_molten_withrandom_oncotsg %>%
  dplyr::group_by(SampleID, KnownGeneClassBinary) %>%
  dplyr::summarise(nIS = n(), .groups = "drop") %>%
  dplyr::group_by(SampleID) %>%
  dplyr::mutate(freq = nIS / sum(nIS) * 100) 

# get Rnadom as the reference
random_data <- IS_gdf_molten_withrandom_oncotsg_desc_simplified %>%
  dplyr::filter(SampleID == "Random") %>%
  dplyr::select(KnownGeneClassBinary, nIS) %>%  # <-- solo queste due colonne
  tidyr::pivot_wider(names_from = KnownGeneClassBinary, values_from = nIS, values_fill = 0)

# Fisher test for each SampleID vs Random
fisher_results <- IS_gdf_molten_withrandom_oncotsg_desc_simplified %>%
  dplyr::filter(SampleID != "Random") %>%
  dplyr::group_by(SampleID) %>%
  dplyr::group_map(~ {
    
    sample_data <- .x %>%
      dplyr::select(KnownGeneClassBinary, nIS) %>%  # <-- solo queste due colonne
      tidyr::pivot_wider(names_from = KnownGeneClassBinary, values_from = nIS, values_fill = 0)
    
    if (!"OG-TSG" %in% names(sample_data)) sample_data$`OG-TSG` <- 0
    if (!"Other"  %in% names(sample_data)) sample_data$Other     <- 0
    if (!"OG-TSG" %in% names(random_data)) random_data$`OG-TSG` <- 0
    if (!"Other"  %in% names(random_data)) random_data$Other     <- 0
    
    contingency_table <- matrix(
      c(sample_data$`OG-TSG`, sample_data$Other,
        random_data$`OG-TSG`, random_data$Other),
      nrow = 2, byrow = TRUE,
      dimnames = list(
        c(.y$SampleID, "Random"),
        c("OG-TSG", "Other")
      )
    )
    
    ft <- fisher.test(contingency_table)
    
    tibble::tibble(
      SampleID   = .y$SampleID,
      p.value    = ft$p.value,
      odds_ratio = ft$estimate,
      conf.low   = ft$conf.int[1],
      conf.high  = ft$conf.int[2]
    )
  }, .keep = TRUE) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(p.adj = p.adjust(p.value, method = "BH"))

fisher_results

# but I need consistent numbers, so:
# 1. subsample all samples
# 2. run test
# 3. repeat N times
# 4. get summary

# subsampling and N repliates
set.seed(42)
n_subsample <- 3500
n_reps <- 100

replicate_fisher <- function(rep_id) {
  
  # Sottocampiona ogni SampleID a n_subsample
  subsampled <- IS_gdf_molten_withrandom_oncotsg %>%
    dplyr::group_by(SampleID) %>%
    dplyr::slice_sample(n = n_subsample, replace = FALSE) %>%
    dplyr::ungroup()
  
  # Ricalcola conteggi binari
  counts <- subsampled %>%
    dplyr::mutate(KnownGeneClassBinary = ifelse(!is.na(`Gene Type`), "OG-TSG", "Other")) %>%
    dplyr::group_by(SampleID, KnownGeneClassBinary) %>%
    dplyr::summarise(nIS = n(), .groups = "drop")
  
  # Estrai Random
  random_data <- counts %>%
    dplyr::filter(SampleID == "Random") %>%
    dplyr::select(KnownGeneClassBinary, nIS) %>%
    tidyr::pivot_wider(names_from = KnownGeneClassBinary, values_from = nIS, values_fill = 0)
  
  # Fisher per ogni sample vs Random
  results <- counts %>%
    dplyr::filter(SampleID != "Random") %>%
    dplyr::group_by(SampleID) %>%
    dplyr::group_map(~ {
      
      sample_data <- .x %>%
        dplyr::select(KnownGeneClassBinary, nIS) %>%
        tidyr::pivot_wider(names_from = KnownGeneClassBinary, values_from = nIS, values_fill = 0)
      
      if (!"OG-TSG" %in% names(sample_data)) sample_data$`OG-TSG` <- 0
      if (!"Other"  %in% names(sample_data)) sample_data$Other     <- 0
      if (!"OG-TSG" %in% names(random_data)) random_data$`OG-TSG` <- 0
      if (!"Other"  %in% names(random_data)) random_data$Other     <- 0
      
      contingency_table <- matrix(
        c(sample_data$`OG-TSG`, sample_data$Other,
          random_data$`OG-TSG`, random_data$Other),
        nrow = 2, byrow = TRUE,
        dimnames = list(c(.y$SampleID, "Random"), c("OG-TSG", "Other"))
      )
      
      ft <- fisher.test(contingency_table)
      
      tibble::tibble(
        SampleID   = .y$SampleID,
        p.value    = ft$p.value,
        odds_ratio = ft$estimate,
        conf.low   = ft$conf.int[1],
        conf.high  = ft$conf.int[2]
      )
    }, .keep = TRUE) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(
      p.adj  = p.adjust(p.value, method = "BH"),
      rep_id = rep_id
    )
  
  return(results)
}

# run N replicates
all_reps <- purrr::map_dfr(1:n_reps, replicate_fisher)

# Riepilogo: per ogni SampleID, in quante repliche p.adj < 0.05?
summary_reps <- all_reps %>%
  dplyr::group_by(SampleID) %>%
  dplyr::summarise(
    n_significant    = sum(p.adj < 0.05),
    mean_odds_ratio  = mean(odds_ratio),
    sd_odds_ratio    = sd(odds_ratio),
    mean_p.adj       = mean(p.adj)
  )

summary_reps

# now, do it against BUT no longer vs Random rather vs Lenti
set.seed(45)
n_subsample <- 3500
n_reps <- 100

replicate_fisher_lenti <- function(rep_id) {
  
  # Subsample each SampleID to n_subsample
  subsampled <- IS_gdf_molten_withrandom_oncotsg %>%
    dplyr::group_by(SampleID) %>%
    dplyr::slice_sample(n = n_subsample, replace = FALSE) %>%
    dplyr::ungroup()
  
  # Recompute binary counts
  counts <- subsampled %>%
    dplyr::mutate(KnownGeneClassBinary = ifelse(!is.na(`Gene Type`), "OG-TSG", "Other")) %>%
    dplyr::group_by(SampleID, KnownGeneClassBinary) %>%
    dplyr::summarise(nIS = n(), .groups = "drop")
  
  # Extract Lenti as reference
  lenti_data <- counts %>%
    dplyr::filter(SampleID == "Lenti") %>%
    dplyr::select(KnownGeneClassBinary, nIS) %>%
    tidyr::pivot_wider(names_from = KnownGeneClassBinary, values_from = nIS, values_fill = 0)
  
  # Fisher test for each sample vs Lenti (excluding Random and Lenti)
  results <- counts %>%
    dplyr::filter(!SampleID %in% c("Lenti", "Random")) %>%
    dplyr::group_by(SampleID) %>%
    dplyr::group_map(~ {
      
      sample_data <- .x %>%
        dplyr::select(KnownGeneClassBinary, nIS) %>%
        tidyr::pivot_wider(names_from = KnownGeneClassBinary, values_from = nIS, values_fill = 0)
      
      if (!"OG-TSG" %in% names(sample_data)) sample_data$`OG-TSG` <- 0
      if (!"Other"  %in% names(sample_data)) sample_data$Other     <- 0
      if (!"OG-TSG" %in% names(lenti_data))  lenti_data$`OG-TSG`  <- 0
      if (!"Other"  %in% names(lenti_data))  lenti_data$Other      <- 0
      
      contingency_table <- matrix(
        c(sample_data$`OG-TSG`, sample_data$Other,
          lenti_data$`OG-TSG`,  lenti_data$Other),
        nrow = 2, byrow = TRUE,
        dimnames = list(c(.y$SampleID, "Lenti"), c("OG-TSG", "Other"))
      )
      
      ft <- fisher.test(contingency_table)
      
      tibble::tibble(
        SampleID   = .y$SampleID,
        p.value    = ft$p.value,
        odds_ratio = ft$estimate,
        conf.low   = ft$conf.int[1],
        conf.high  = ft$conf.int[2]
      )
    }, .keep = TRUE) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(
      p.adj  = p.adjust(p.value, method = "BH"),
      rep_id = rep_id
    )
  
  return(results)
}

# Run N replicates
all_reps_lenti <- purrr::map_dfr(1:n_reps, replicate_fisher_lenti)

# Summary: how many replicates are significant per SampleID?
summary_reps_lenti <- all_reps_lenti %>%
  dplyr::group_by(SampleID) %>%
  dplyr::summarise(
    n_significant   = sum(p.adj < 0.05),
    mean_odds_ratio = mean(odds_ratio),
    sd_odds_ratio   = sd(odds_ratio),
    mean_p.adj      = mean(p.adj)
  )

summary_reps_lenti

# write results in an excel file
write.xlsx(x = summary_reps, 
           file = paste(source_folder, analysis_folder_date, ".OG-TSG.vsRandom.summary_reps.xlsx", sep = ""))
write.xlsx(x = summary_reps_lenti, 
           file = paste(source_folder, analysis_folder_date, ".OG-TSG.vsRandom.summary_reps_lenti.xlsx", sep = ""))


# ── Plot Data ──────────────────────────────────────────────────────────────────────

# Desired display order (top to bottom on y-axis)
sample_order  <- c("RTE25", "RTE3", "CR1", "Vingi", "Lenti")
sample_labels <- c(
  "RTE25" = "RTE-25_Lmi",
  "RTE3"  = "RTE-3_BF",
  "CR1"   = "CR1-1_PH",
  "Vingi" = "Vingi-1_Acar",
  "Lenti" = "Lenti"
)

# Results vs Random (reference = Random) — 100 replicates
vs_random <- summary_reps %>%
  mutate(
    ci_lo = mean_odds_ratio - 1.96 * sd_odds_ratio,
    ci_hi = mean_odds_ratio + 1.96 * sd_odds_ratio,
    n_sig = n_significant,
  )

# Results vs Lenti (reference = Lenti) — 100 replicates
vs_lenti <- summary_reps_lenti %>%
  mutate(
    ci_lo = mean_odds_ratio - 1.96 * sd_odds_ratio,
    ci_hi = mean_odds_ratio + 1.96 * sd_odds_ratio,
    n_sig = n_significant
  )

# ── Factor levels: rev() so that first in order = top of y-axis ──────────────

vs_random_ordered <- vs_random %>%
  mutate(
    SampleID  = factor(SampleID, levels = rev(sample_order)),
    dot_color = case_when(
      SampleID == "Lenti" ~ "#c0392b",
      TRUE                ~ "#9ca3af"
    ),
    # Significance thresholds calibrated on 100 replicates
    sig_label = case_when(
      n_sig == 100  ~ "***",
      n_sig >= 95   ~ "**",
      n_sig >= 80   ~ "*",
      TRUE          ~ "ns"
    )
  )

vs_lenti_ordered <- vs_lenti %>%
  mutate(
    SampleID  = factor(SampleID, levels = rev(sample_order[sample_order != "Lenti"])),
    dot_color = "#1a1a2e",
    sig_label = "***"
  )

# ── Shared publication theme ──────────────────────────────────────────────────

pub_theme <- theme_classic(base_size = 10, base_family = "Helvetica") +
  theme(
    axis.text          = element_text(size = 9, color = "black"),
    axis.title         = element_text(size = 9, color = "black"),
    axis.line          = element_line(color = "black", linewidth = 0.4),
    axis.ticks         = element_line(color = "black", linewidth = 0.4),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title         = element_text(size = 10, face = "bold", color = "black", margin = margin(b = 4)),
    plot.subtitle      = element_text(size = 8, color = "grey40", margin = margin(b = 10)),
    legend.position    = "none",
    plot.margin        = margin(10, 16, 10, 10)
  )

# ── Panel A: vs Random ───────────────────────────────────────────────────────

pA <- ggplot(vs_random_ordered, aes(y = SampleID)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.4) +
  # geom_segment(aes(x = ci_lo, xend = ci_hi, yend = SampleID, color = dot_color),
  #              linewidth = 0.6) +
  geom_segment(aes(x = ci_lo, xend = ci_hi, yend = SampleID),
               color = "black", linewidth = 0.6) +
  # geom_point(aes(x = mean_odds_ratio, color = dot_color,
  #                shape = ifelse(SampleID == "Lenti", 16, 16),
  #                size  = ifelse(SampleID == "Lenti", 3, 3))) +
  geom_point(aes(x = mean_odds_ratio,
                 shape = ifelse(SampleID == "Lenti", 16, 16),
                 size  = ifelse(SampleID == "Lenti", 3, 3)),
             color = "black") +
  # geom_text(aes(x = ci_hi + 0.1, label = sig_label, color = dot_color),
  #           hjust = 0, size = 3, fontface = "bold") +
  geom_text(aes(x = ci_hi + 0.1, label = sig_label),
            color = "black", hjust = 0, size = 3, fontface = "bold") +
  scale_color_identity() +
  scale_size_identity() +
  scale_shape_identity() +
  scale_y_discrete(labels = sample_labels) +
  # scale_x_continuous(
  #   limits = c(0.5, 4.3),
  #   breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5),
  #   expand = expansion(mult = c(0.02, 0.1))
  # ) +  
  scale_x_log10(
    # limits = c(0.5, 4.3),
    # breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5),
    # expand = expansion(mult = c(0.02, 0.1))
  ) +
  labs(
    title    = "OG-TSG enrichment vs Random",
    # subtitle = paste0("Mean OR ± 95% CI across N=", n_reps, " replicates, subsampling at ", n_subsample, "; p.adj < 0.05."),
    x        = "Odds Ratio",
    y        = NULL
  ) +
  # theme_bw() +
  pub_theme

pAv2 <- ggplot(vs_random_ordered, aes(x = SampleID)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_segment(aes(y = ci_lo, yend = ci_hi, xend = SampleID),
               color = "black", linewidth = 0.6) +
  geom_point(aes(y = mean_odds_ratio,
                 shape = ifelse(SampleID == "Lenti", 16, 16),
                 size  = ifelse(SampleID == "Lenti", 3, 3)),
             color = "black") +
  geom_text(aes(y = ci_hi, label = sig_label),
            color = "black", vjust = -0.6, hjust = 0.5, size = 3, fontface = "bold") +
  scale_x_discrete(labels = sample_labels,
                   limits = rev(levels(vs_random_ordered$SampleID))) + 
  scale_color_identity() +
  scale_size_identity() +
  scale_shape_identity() +
  # scale_x_discrete(labels = sample_labels) +
  scale_y_log10(limits = c(0.5, 4), breaks = c(0.5, 0.75, 1, 1.5, 2, 3, 4)) +
  labs(
    title = "OG-TSG enrichment vs Random",
    y     = "Odds Ratio",
    x     = NULL
  ) +
  pub_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ── Panel B: vs Lenti ────────────────────────────────────────────────────────

pB <- ggplot(vs_lenti_ordered, aes(y = SampleID)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  annotate("rect", xmin = -Inf, xmax = 1, ymin = -Inf, ymax = Inf,
           fill = "#fee2e2", alpha = 0.35) +
  geom_segment(aes(x = ci_lo, xend = ci_hi, yend = SampleID),
               color = "#1a1a2e", linewidth = 0.6) +
  geom_point(aes(x = mean_odds_ratio), color = "#1a1a2e", size = 2.8, shape = 16) +
  geom_text(aes(x = ci_hi + 0.015, label = sig_label),
            hjust = 0, size = 3, color = "#1a1a2e", fontface = "bold") +
  scale_y_discrete(labels = sample_labels) +
  scale_x_continuous(
    limits = c(0.25, 1.3),
    breaks = c(0.25, 0.5, 0.75, 1.0, 1.25),
    expand = expansion(mult = c(0.02, 0.1))
  ) +
  labs(
    title    = "OG-TSG enrichment vs Lenti",
    subtitle = paste0("All vectors significantly depleted (N=", n_reps," replicates, p.adj < 0.05)"),
    x        = "Odds Ratio",
    y        = NULL
  ) +
  pub_theme

# ── Combine and save ─────────────────────────────────────────────────────────

pdf(file = paste(source_folder, analysis_folder_date, ".BIO.OG-TSG.vsRandom.pdf", sep = ""), height=5, width=6)
plot(pA)
dev.off()
png(file = paste(source_folder, analysis_folder_date, ".BIO.OG-TSG.vsRandom.png", sep = ""), height=5, width=6, units = "in", res = 300)
plot(pA)
dev.off()

pdf(file = paste(source_folder, analysis_folder_date, ".BIO.OG-TSG.vsRandom.v2.pdf", sep = ""), height=4, width=5)
plot(pA)
dev.off()
png(file = paste(source_folder, analysis_folder_date, ".BIO.OG-TSG.vsRandom.v2.png", sep = ""), height=4, width=5, units = "in", res = 300)
plot(pA)
dev.off()

pdf(file = paste(source_folder, analysis_folder_date, ".BIO.OG-TSG.vsRandom.v2.new.pdf", sep = ""), height=4, width=3.5)
plot(pAv2)
dev.off()
png(file = paste(source_folder, analysis_folder_date, ".BIO.OG-TSG.vsRandom.v2.new.png", sep = ""), height=4, width=3.5, units = "in", res = 300)
plot(pAv2)
dev.off()

pdf(file = paste(source_folder, analysis_folder_date, ".BIO.OG-TSG.vsLenti.pdf", sep = ""), height=5, width=6)
plot(pB)
dev.off()
png(file = paste(source_folder, analysis_folder_date, ".BIO.OG-TSG.vsLenti.png", sep = ""), height=5, width=6, units = "in", res = 300)
plot(pB)
dev.off()

 
# fig <- pA | pB
# 
# ggsave(
#   filename = "/mnt/user-data/outputs/forest_plot_figure.pdf",
#   plot     = fig,
#   width    = 7.5,
#   height   = 3.2,
#   units    = "in",
#   device   = cairo_pdf
# )
# 
# ggsave(
#   filename = "/mnt/user-data/outputs/forest_plot_figure.png",
#   plot     = fig,
#   width    = 7.5,
#   height   = 3.2,
#   units    = "in",
#   dpi      = 300
# )
# 
# message("Saved: forest_plot_figure.pdf and .png")

library(AnnotationHub)
library(GenomicRanges)
library(GenomeInfoDb)
library(ensembldb)
library(regioneR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(purrr)

# ── 1. Load Ensembl gene annotation for hg38 ─────────────────────────────────

ah <- AnnotationHub()
ensdb    <- ah[["AH113665"]]  # replace with correct ID from query()
genes_gr <- genes(ensdb)

seqlevelsStyle(genes_gr) <- "UCSC"
genes_gr <- keepStandardChromosomes(genes_gr, pruning.mode = "coarse")

# ── 2. Define OG/TSG gene list from OncoKB ───────────────────────────────────

oncotsg_genes <- IS_gdf_molten_withrandom_oncotsg %>%
  dplyr::filter(!is.na(`Gene Type`)) %>%
  dplyr::pull(GeneName) %>%
  unique()

message("N OG/TSG genes from OncoKB: ", length(oncotsg_genes))

oncotsg_gr <- genes_gr[genes_gr$gene_name %in% oncotsg_genes]

message("N OG/TSG genes mapped to Ensembl: ", length(oncotsg_gr))
message("N OG/TSG genes NOT found in Ensembl: ",
        sum(!oncotsg_genes %in% genes_gr$gene_name))

# ── 3. Convert IS to GRanges ─────────────────────────────────────────────────

is_gr <- GRanges(
  seqnames = IS_gdf_molten_withrandom_oncotsg$chr,
  ranges   = IRanges(
    start = IS_gdf_molten_withrandom_oncotsg$integration_locus,
    width = 1
  ),
  SampleID = IS_gdf_molten_withrandom_oncotsg$SampleID
)

# ── 4. Harmonize chromosomes ──────────────────────────────────────────────────

seqlevelsStyle(is_gr)    <- "UCSC"
seqlevelsStyle(oncotsg_gr) <- "UCSC"

is_gr      <- keepStandardChromosomes(is_gr,      pruning.mode = "coarse")
oncotsg_gr <- keepStandardChromosomes(oncotsg_gr, pruning.mode = "coarse")

shared_chroms <- intersect(seqlevels(is_gr), seqlevels(oncotsg_gr))
is_gr         <- keepSeqlevels(is_gr,      shared_chroms, pruning.mode = "coarse")
oncotsg_gr    <- keepSeqlevels(oncotsg_gr, shared_chroms, pruning.mode = "coarse")

# ── 5. Associate each IS to nearest OG/TSG gene within 1kb ───────────────────

max_dist <- 1000L  # 1kb threshold

# distanceToNearest returns distance to closest gene
nearest_hits <- distanceToNearest(is_gr, oncotsg_gr, ignore.strand = TRUE)

# IS is associated to OG/TSG only if nearest gene is within 1kb
in_oncotsg_vec <- rep(FALSE, length(is_gr))
within_1kb     <- mcols(nearest_hits)$distance <= max_dist
in_oncotsg_vec[queryHits(nearest_hits)[within_1kb]] <- TRUE

message("IS associated to OG/TSG within 1kb: ", sum(in_oncotsg_vec))
message("IS NOT associated to any OG/TSG within 1kb: ", sum(!in_oncotsg_vec))

# Attach annotation to is_gr
is_gr$in_oncotsg <- in_oncotsg_vec

# Sanity check per SampleID
is_df_full <- data.frame(
  SampleID   = as.character(is_gr$SampleID),
  in_oncotsg = is_gr$in_oncotsg
)

message("Annotation check (full dataset):")
print(table(is_df_full$in_oncotsg, is_df_full$SampleID))

# ── 6. Subsample each SampleID to 3,500 IS ───────────────────────────────────

n_subsample <- 3500
set.seed(42)

is_gr_subsampled <- lapply(unique(as.character(is_gr$SampleID)), function(sid) {
  
  gr_sid <- is_gr[is_gr$SampleID == sid]
  
  if (length(gr_sid) > n_subsample) {
    gr_sid <- gr_sid[sample(length(gr_sid), n_subsample)]
    message("Subsampled ", sid, " to ", n_subsample, " IS")
  } else {
    message(sid, " has ", length(gr_sid), " IS — no subsampling needed")
  }
  
  gr_sid
  
}) %>% do.call(c, .)

message("\nIS counts after subsampling:")
print(table(is_gr_subsampled$SampleID))

# Sanity check after subsampling
is_df <- data.frame(
  SampleID   = as.character(is_gr_subsampled$SampleID),
  in_oncotsg = is_gr_subsampled$in_oncotsg
)

message("Annotation check (after subsampling):")
print(table(is_df$in_oncotsg, is_df$SampleID))

# ── 7. Fisher test vs Random — 100 subsampling replicates ────────────────────

n_reps <- 100

replicate_fisher_oncotsg <- function(rep_id) {
  
  # Subsample from full annotated is_gr (annotation already done)
  sub_df <- as.data.frame(is_gr) %>%
    dplyr::mutate(row_idx = dplyr::row_number()) %>%
    dplyr::group_by(SampleID) %>%
    dplyr::slice_sample(n = n_subsample, replace = FALSE) %>%
    dplyr::ungroup()
  
  # Counts per SampleID
  counts <- sub_df %>%
    dplyr::group_by(SampleID) %>%
    dplyr::summarise(
      n_in  = sum(in_oncotsg),
      n_out = sum(!in_oncotsg),
      .groups = "drop"
    )
  
  random_row <- counts %>% dplyr::filter(SampleID == "Random")
  
  # Fisher test for each sample vs Random
  counts %>%
    dplyr::filter(SampleID != "Random") %>%
    dplyr::group_by(SampleID) %>%
    dplyr::group_map(~ {
      
      ct <- matrix(
        c(.x$n_in,         .x$n_out,
          random_row$n_in,  random_row$n_out),
        nrow = 2, byrow = TRUE,
        dimnames = list(c(.y$SampleID, "Random"), c("n_in", "n_out"))
      )
      
      ft <- fisher.test(ct)
      
      tibble::tibble(
        SampleID   = .y$SampleID,
        p.value    = ft$p.value,
        odds_ratio = ft$estimate,
        conf.low   = ft$conf.int[1],
        conf.high  = ft$conf.int[2],
        rep_id     = rep_id
      )
    }, .keep = TRUE) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(p.adj = p.adjust(p.value, method = "BH"))
}

# Test single replicate first
message("Testing single replicate...")
test_rep <- replicate_fisher_oncotsg(1)
print(test_rep)

# Run all 100 replicates
message("Running 100 replicates...")
all_reps_oncotsg <- purrr::map_dfr(1:n_reps, replicate_fisher_oncotsg)

summary_oncotsg_vs_random <- all_reps_oncotsg %>%
  dplyr::group_by(SampleID) %>%
  dplyr::summarise(
    n_significant   = sum(p.adj < 0.05),
    mean_odds_ratio = mean(odds_ratio),
    sd_odds_ratio   = sd(odds_ratio),
    mean_p.adj      = mean(p.adj),
    .groups         = "drop"
  )

print(summary_oncotsg_vs_random)

# ── 8. Permutation test vs genome background ──────────────────────────────────

genome_hg38  <- BSgenome.Hsapiens.UCSC.hg38
sample_ids   <- unique(as.character(is_gr_subsampled$SampleID))
sample_ids   <- sample_ids[sample_ids != "Random"]
is_random_gr <- is_gr_subsampled[is_gr_subsampled$SampleID == "Random"]

# Extend OG/TSG genes by 1kb on each side to match the 1kb association rule
oncotsg_gr_extended <- resize(oncotsg_gr,
                              width = width(oncotsg_gr) + 2000L,
                              fix   = "center")

perm_results_oncotsg <- lapply(sample_ids, function(sid) {
  
  message("Running permutation test for: ", sid)
  
  is_sample <- is_gr_subsampled[is_gr_subsampled$SampleID == sid]
  
  tryCatch({
    
    pt <- permTest(
      A                  = is_sample,
      B                  = oncotsg_gr_extended,
      randomize.function = randomizeRegions,
      genome             = genome_hg38,
      evaluate.function  = numOverlaps,
      ntimes             = 1000,
      mc.cores           = 1,
      verbose            = FALSE
    )
    
    tibble::tibble(
      SampleID   = sid,
      observed   = pt$numOverlaps$observed,
      expected   = mean(pt$numOverlaps$permuted, na.rm = TRUE),
      sd_perm    = sd(pt$numOverlaps$permuted,   na.rm = TRUE),
      zscore     = pt$numOverlaps$zscore,
      p.value    = pt$numOverlaps$pval,
      enrichment = pt$numOverlaps$observed / mean(pt$numOverlaps$permuted, na.rm = TRUE)
    )
    
  }, error = function(e) {
    message("Error for ", sid, ": ", e$message)
    tibble::tibble(
      SampleID   = sid,
      observed   = NA_real_, expected   = NA_real_,
      sd_perm    = NA_real_, zscore     = NA_real_,
      p.value    = NA_real_, enrichment = NA_real_
    )
  })
})

# Self-check: Random vs genome
message("Running self-check for Random...")

pt_random_oncotsg <- tryCatch({
  permTest(
    A                  = is_random_gr,
    B                  = oncotsg_gr_extended,
    randomize.function = randomizeRegions,
    genome             = genome_hg38,
    evaluate.function  = numOverlaps,
    ntimes             = 1000,
    mc.cores           = 1,
    verbose            = FALSE
  )
}, error = function(e) {
  message("Self-check error: ", e$message)
  NULL
})

# Combine permutation results
results_oncotsg <- dplyr::bind_rows(perm_results_oncotsg) %>%
  dplyr::mutate(p.adj = p.adjust(p.value, method = "BH"))

if (!is.null(pt_random_oncotsg)) {
  results_oncotsg <- dplyr::bind_rows(
    results_oncotsg,
    tibble::tibble(
      SampleID   = "Random",
      observed   = pt_random_oncotsg$numOverlaps$observed,
      expected   = mean(pt_random_oncotsg$numOverlaps$permuted, na.rm = TRUE),
      sd_perm    = sd(pt_random_oncotsg$numOverlaps$permuted,   na.rm = TRUE),
      zscore     = pt_random_oncotsg$numOverlaps$zscore,
      p.value    = pt_random_oncotsg$numOverlaps$pval,
      enrichment = pt_random_oncotsg$numOverlaps$observed / mean(pt_random_oncotsg$numOverlaps$permuted, na.rm = TRUE),
      p.adj      = NA_real_
    )
  )
}

print(results_oncotsg)

### plot results

# ── Sample order and labels ───────────────────────────────────────────────────

sample_order  <- c("RTE25", "RTE3", "CR1", "Vingi", "Lenti")
sample_labels <- c(
  "RTE25" = "RTE-25_Lmi",
  "RTE3"  = "RTE-3_BF",
  "CR1"   = "CR1-1_PH",
  "Vingi" = "Vingi-1_Acar",
  "Lenti" = "Lenti"
)

# ── Prepare Panel B data from summary_oncotsg_vs_random ──────────────────────

pB_data <- summary_oncotsg_vs_random %>%
  mutate(
    ci_lo     = mean_odds_ratio - 1.96 * sd_odds_ratio,
    ci_hi     = mean_odds_ratio + 1.96 * sd_odds_ratio,
    SampleID  = factor(SampleID, levels = rev(sample_order)),
    dot_color = case_when(
      SampleID == "Lenti" ~ "#c0392b",
      TRUE                ~ "#1a1a2e"
    ),
    sig_label = case_when(
      n_significant == 100 ~ "***",
      n_significant >= 95  ~ "**",
      n_significant >= 80  ~ "*",
      TRUE                 ~ "ns"
    )
  )

# ── Prepare Panel C data from summary_reps_lenti ─────────────────────────────
# summary_reps_lenti contains OR of each vector vs Lenti (reference)
# OR < 1 = depleted vs Lenti

pC_data <- summary_reps_lenti %>%
  mutate(
    ci_lo     = mean_odds_ratio - 1.96 * sd_odds_ratio,
    ci_hi     = mean_odds_ratio + 1.96 * sd_odds_ratio,
    # Exclude Lenti itself (it's the reference)
    SampleID  = factor(SampleID,
                       levels = rev(sample_order[sample_order != "Lenti"])),
    dot_color = "#1a1a2e",
    sig_label = case_when(
      n_significant == 100 ~ "***",
      n_significant >= 95  ~ "**",
      n_significant >= 80  ~ "*",
      TRUE                 ~ "ns"
    )
  ) %>%
  dplyr::filter(!is.na(SampleID))

# ── Shared publication theme ──────────────────────────────────────────────────

pub_theme <- theme_classic(base_size = 10, base_family = "Helvetica") +
  theme(
    axis.text          = element_text(size = 9, color = "black"),
    axis.title         = element_text(size = 9, color = "black"),
    axis.line          = element_line(color = "black", linewidth = 0.4),
    axis.ticks         = element_line(color = "black", linewidth = 0.4),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title         = element_text(size = 10, face = "bold", color = "black",
                                      margin = margin(b = 4)),
    plot.subtitle      = element_text(size = 8, color = "grey40",
                                      margin = margin(b = 10)),
    legend.position    = "none",
    plot.margin        = margin(10, 16, 10, 10)
  )


# ── Panel B: OR vs Random IS ──────────────────────────────────────────────────

x_max_B <- max(pB_data$ci_hi, na.rm = TRUE) * 1.12

pB <- ggplot(pB_data, aes(y = SampleID)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50",
             linewidth = 0.4) +
  annotate("rect", xmin = 1, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#dbeafe", alpha = 0.35) +
  geom_segment(aes(x = ci_lo, xend = ci_hi, yend = SampleID,
                   color = dot_color), linewidth = 0.6) +
  geom_point(aes(x = mean_odds_ratio, color = dot_color,
                 shape = ifelse(SampleID == "Lenti", 18, 16),
                 size  = ifelse(SampleID == "Lenti", 4, 2.8))) +
  geom_text(aes(x = ci_hi + x_max_B * 0.03,
                label = sig_label, color = dot_color),
            hjust = 0, size = 3, fontface = "bold") +
  scale_color_identity() +
  scale_size_identity() +
  scale_shape_identity() +
  scale_y_discrete(labels = sample_labels) +
  scale_x_continuous(
    limits = c(NA, x_max_B),
    expand = expansion(mult = c(0.05, 0.1))
  ) +
  labs(
    title    = "OG-TSG enrichment vs Random IS",
    subtitle = "Mean OR \u00b1 95% CI across 100 subsampling replicates (n = 3,500, BH correction)",
    x        = "Odds Ratio",
    y        = NULL
  ) +
  pub_theme

# ── Panel C: OR vs Lenti ──────────────────────────────────────────────────────

x_max_C <- max(pC_data$ci_hi, na.rm = TRUE) * 1.12

pC <- ggplot(pC_data, aes(y = SampleID)) +
  # Green background first
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#d1fae5", alpha = 0.4) +
  # Vertical line on top
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey30",
             linewidth = 0.6) +
  geom_segment(aes(x = ci_lo, xend = ci_hi, yend = SampleID),
               color = "#1a1a2e", linewidth = 0.6) +
  geom_point(aes(x = mean_odds_ratio), color = "#1a1a2e",
             size = 2.8, shape = 16) +
  geom_text(aes(x = ci_hi + x_max_C * 0.03, label = sig_label),
            hjust = 0, size = 3, color = "#1a1a2e", fontface = "bold") +
  scale_color_identity() +
  scale_y_discrete(labels = sample_labels) +
  scale_x_continuous(
    limits = c(NA, x_max_C),
    expand = expansion(mult = c(0.05, 0.1))
  ) +
  labs(
    title    = "OG-TSG enrichment vs Lenti",
    subtitle = "Mean OR \u00b1 95% CI across 100 subsampling replicates (n = 3,500, BH correction)",
    x        = "Odds Ratio (ref = Lenti)",
    y        = NULL
  ) +
  pub_theme +
  # Override panel background to allow annotate rect to show through
  theme(
    panel.background = element_rect(fill = NA, color = NA),
    plot.background  = element_rect(fill = NA, color = NA)
  )

# ── Combine and save ──────────────────────────────────────────────────────────


pdf(file = paste(source_folder, analysis_folder_date, ".BIO.OG-TSG.Ensmbl.vsRandom.pdf", sep = ""), height=5, width=6)
plot(pB)
dev.off()
png(file = paste(source_folder, analysis_folder_date, ".BIO.OG-TSG.Ensmbl.vsRandom.png", sep = ""), height=5, width=6, units = "in", res = 300)
plot(pB)
dev.off()

pdf(file = paste(source_folder, analysis_folder_date, ".BIO.OG-TSG.Ensmbl.vsLenti.pdf", sep = ""), height=5, width=6)
plot(pC)
dev.off()
png(file = paste(source_folder, analysis_folder_date, ".BIO.OG-TSG.Ensmbl.vsLenti.png", sep = ""), height=5, width=6, units = "in", res = 300)
plot(pC)
dev.off()


fig <- pA | pB | pC
# paste(source_folder, analysis_folder_date, ".BIO.OpenChromatin.GW.pdf", sep = "")
ggsave(
  filename = paste(source_folder, analysis_folder_date, ".BIO.OG-TSG.forest.panel.pdf", sep = ""),
  plot     = fig,
  width    = 12.0,
  height   = 3.6,
  units    = "in",
  device   = cairo_pdf
)

ggsave(
  filename = paste(source_folder, analysis_folder_date, ".BIO.OG-TSG.forest.panel.png", sep = ""),
  plot     = fig,
  width    = 14.0,
  height   = 3.6,
  units    = "in",
  dpi      = 300
)

message("Saved: forest_plot_oncotsg_3panels.pdf and .png")

###############################################################
## Annotation enrichment, open chromatin
###############################################################
library(AnnotationHub)
library(GenomicRanges)
library(GenomeInfoDb)
library(regioneR)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)

# ── 1. Load HSC chromatin states from Roadmap Epigenomics (E035) ─────────────

ah <- AnnotationHub()

chrom_states <- ah[["AH46890"]]

open_states <- c(
  "Active TSS",
  "Flanking Active TSS",
  "Genic enhancers",
  "Enhancers",
  "Bivalent/Poised TSS",
  "Flanking Bivalent TSS/Enh"
)

hsc_open <- chrom_states[chrom_states$name %in% open_states]

# ── 2. Liftover hsc_open from hg19 to hg38 ───────────────────────────────────

chain <- ah[["AH14150"]]  # hg19 -> hg38

hsc_open_hg38 <- unlist(liftOver(hsc_open, chain))
hsc_open_hg38 <- keepStandardChromosomes(hsc_open_hg38, pruning.mode = "coarse")

message("N open chromatin regions after liftover: ", length(hsc_open_hg38))
message("Total coverage bp: ", sum(width(reduce(hsc_open_hg38))))

# ── 3. Convert IS to GRanges ─────────────────────────────────────────────────

is_gr <- GRanges(
  seqnames = IS_gdf_molten_withrandom_oncotsg$chr,
  ranges   = IRanges(
    start = IS_gdf_molten_withrandom_oncotsg$integration_locus,
    width = 1
  ),
  SampleID = IS_gdf_molten_withrandom_oncotsg$SampleID
)

# ── 4. Harmonize chromosome names and keep only standard chromosomes ──────────

seqlevelsStyle(is_gr)        <- "UCSC"
seqlevelsStyle(hsc_open_hg38) <- "UCSC"

is_gr        <- keepStandardChromosomes(is_gr,        pruning.mode = "coarse")
hsc_open_hg38 <- keepStandardChromosomes(hsc_open_hg38, pruning.mode = "coarse")

shared_chroms <- intersect(seqlevels(is_gr), seqlevels(hsc_open_hg38))
is_gr         <- keepSeqlevels(is_gr,         shared_chroms, pruning.mode = "coarse")
hsc_open_hg38 <- keepSeqlevels(hsc_open_hg38, shared_chroms, pruning.mode = "coarse")

# ── 5. Subsample each SampleID to 3,500 IS and apply ±100bp window ───────────

n_subsample <- 3500
set.seed(42)

is_gr_subsampled <- lapply(unique(as.character(is_gr$SampleID)), function(sid) {
  
  gr_sid <- is_gr[is_gr$SampleID == sid]
  
  if (length(gr_sid) > n_subsample) {
    gr_sid <- gr_sid[sample(length(gr_sid), n_subsample)]
    message("Subsampled ", sid, " to ", n_subsample, " IS")
  } else {
    message(sid, " has ", length(gr_sid), " IS — no subsampling needed")
  }
  
  gr_sid
  
}) %>% do.call(c, .)

is_gr_window <- resize(is_gr_subsampled, width = 200, fix = "center")  # ±100bp

message("\nIS counts after subsampling:")
print(table(is_gr_window$SampleID))

# ── 6. Define Random IS window as resampling universe ────────────────────────

is_random_window <- is_gr_window[is_gr_window$SampleID == "Random"]

# ── 7. Permutation test per SampleID using randomizeRegions ──────────────────

genome_hg38  <- BSgenome.Hsapiens.UCSC.hg38
sample_ids   <- unique(as.character(is_gr_window$SampleID))
sample_ids   <- sample_ids[sample_ids != "Random"]

perm_results <- lapply(sample_ids, function(sid) {
  
  message("Running permutation test for: ", sid)
  
  is_sample <- is_gr_window[is_gr_window$SampleID == sid]
  
  tryCatch({
    
    pt <- permTest(
      A                  = is_sample,
      B                  = hsc_open_hg38,
      randomize.function = randomizeRegions,
      genome             = genome_hg38,
      evaluate.function  = numOverlaps,
      ntimes             = 1000,
      mc.cores           = 1,
      verbose            = FALSE
    )
    
    tibble::tibble(
      SampleID   = sid,
      observed   = pt$numOverlaps$observed,
      expected   = mean(pt$numOverlaps$permuted, na.rm = TRUE),
      sd_perm    = sd(pt$numOverlaps$permuted,   na.rm = TRUE),
      zscore     = pt$numOverlaps$zscore,
      p.value    = pt$numOverlaps$pval,
      enrichment = pt$numOverlaps$observed / mean(pt$numOverlaps$permuted, na.rm = TRUE)
    )
    
  }, error = function(e) {
    message("Error for ", sid, ": ", e$message)
    tibble::tibble(
      SampleID   = sid,
      observed   = NA_real_, expected   = NA_real_,
      sd_perm    = NA_real_, zscore     = NA_real_,
      p.value    = NA_real_, enrichment = NA_real_
    )
  })
})

# ── 8. Self-check: Random vs Random (enrichment should be ~1) ────────────────

message("Running self-check for Random...")

pt_random <- tryCatch({
  
  permTest(
    A                  = is_random_window,
    B                  = hsc_open_hg38,
    randomize.function = randomizeRegions,
    genome             = genome_hg38,
    evaluate.function  = numOverlaps,
    ntimes             = 1000,
    mc.cores           = 1,
    verbose            = FALSE
  )
  
}, error = function(e) {
  message("Self-check error: ", e$message)
  NULL
})

# ── 9. Combine and print results ─────────────────────────────────────────────

results_open <- dplyr::bind_rows(perm_results) %>%
  dplyr::mutate(p.adj = p.adjust(p.value, method = "BH"))

if (!is.null(pt_random)) {
  results_open <- dplyr::bind_rows(
    results_open,
    tibble::tibble(
      SampleID   = "Random",
      observed   = pt_random$numOverlaps$observed,
      expected   = mean(pt_random$numOverlaps$permuted, na.rm = TRUE),
      sd_perm    = sd(pt_random$numOverlaps$permuted,   na.rm = TRUE),
      zscore     = pt_random$numOverlaps$zscore,
      p.value    = pt_random$numOverlaps$pval,
      enrichment = pt_random$numOverlaps$observed / mean(pt_random$numOverlaps$permuted, na.rm = TRUE),
      p.adj      = NA_real_
    )
  )
}

print(results_open)



# ── Annotate each IS as in/out open chromatin ─────────────────────────────────

hits <- findOverlaps(is_gr_window, hsc_open_hg38)

is_gr_window$in_open <- seq_len(length(is_gr_window)) %in% queryHits(hits)

# ── Correct annotation using a dataframe directly ────────────────────────────

is_df <- as.data.frame(is_gr_window) %>%
  dplyr::mutate(is_index = dplyr::row_number())

# Get overlapping indices
hits <- findOverlaps(is_gr_window, hsc_open_hg38)
open_indices <- unique(queryHits(hits))

# Annotate in dataframe
is_df$in_open <- is_df$is_index %in% open_indices

# Verify
table(is_df$in_open, is_df$SampleID)

# Count per SampleID
open_counts <- is_df %>%
  dplyr::group_by(SampleID) %>%
  dplyr::summarise(
    in_open  = sum(in_open),
    out_open = sum(!in_open),
    total    = n(),
    pct_open = round(in_open / total * 100, 1),
    .groups  = "drop"
  )

print(open_counts)

# ── Fisher test vs Random with 100 subsampling replicates ────────────────────

set.seed(42)
n_reps      <- 100
n_subsample <- 3500

replicate_fisher_open <- function(rep_id) {
  
  sub_df <- as.data.frame(is_gr) %>%
    dplyr::group_by(SampleID) %>%
    dplyr::slice_sample(n = n_subsample, replace = FALSE) %>%
    dplyr::ungroup()
  
  sub_gr <- GRanges(
    seqnames = sub_df$seqnames,
    ranges   = IRanges(start = sub_df$start, width = 1),
    SampleID = sub_df$SampleID
  )
  sub_gr <- resize(sub_gr, width = 200, fix = "center")
  sub_gr <- keepSeqlevels(sub_gr, shared_chroms, pruning.mode = "coarse")
  
  hits_sub    <- findOverlaps(sub_gr, hsc_open_hg38)
  open_idx    <- unique(queryHits(hits_sub))
  in_open_vec <- seq_len(length(sub_gr)) %in% open_idx
  
  sub_ann <- data.frame(
    SampleID = as.character(sub_gr$SampleID),
    in_open  = in_open_vec
  )
  
  # Use unambiguous column names: n_in and n_out
  counts <- sub_ann %>%
    dplyr::group_by(SampleID) %>%
    dplyr::summarise(
      n_in  = sum(in_open),
      n_out = sum(!in_open),
      .groups = "drop"
    )
  
  random_row <- counts %>% dplyr::filter(SampleID == "Random")
  
  counts %>%
    dplyr::filter(SampleID != "Random") %>%
    dplyr::group_by(SampleID) %>%
    dplyr::group_map(~ {
      
      ct <- matrix(
        c(.x$n_in,          .x$n_out,
          random_row$n_in,   random_row$n_out),
        nrow = 2, byrow = TRUE,
        dimnames = list(c(.y$SampleID, "Random"), c("n_in", "n_out"))
      )
      
      ft <- fisher.test(ct)
      
      tibble::tibble(
        SampleID   = .y$SampleID,
        p.value    = ft$p.value,
        odds_ratio = ft$estimate,
        conf.low   = ft$conf.int[1],
        conf.high  = ft$conf.int[2],
        rep_id     = rep_id
      )
    }, .keep = TRUE) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(p.adj = p.adjust(p.value, method = "BH"))
}

# Test single replicate
test_rep <- replicate_fisher_open(1)
print(test_rep)

# Run 100 replicates
all_reps_open <- purrr::map_dfr(1:n_reps, replicate_fisher_open)

# Summary across replicates
summary_open_vs_random <- all_reps_open %>%
  dplyr::group_by(SampleID) %>%
  dplyr::summarise(
    n_significant   = sum(p.adj < 0.05),
    mean_odds_ratio = mean(odds_ratio),
    sd_odds_ratio   = sd(odds_ratio),
    mean_p.adj      = mean(p.adj),
    .groups         = "drop"
  )

print(summary_open_vs_random)


# ── Sample order and labels ───────────────────────────────────────────────────

sample_order  <- c("RTE25", "RTE3", "CR1", "Vingi", "Lenti")
sample_labels <- c(
  "RTE25" = "RTE-25_Lmi",
  "RTE3"  = "RTE-3_BF",
  "CR1"   = "CR1-1_PH",
  "Vingi" = "Vingi-1_Acar",
  "Lenti" = "Lenti"
)

# ── Prepare Panel A data from results_open ────────────────────────────────────

pA_data <- results_open %>%
  mutate(
    SampleID  = factor(SampleID, levels = rev(c(sample_order, "Random"))),
    dot_color = case_when(
      SampleID == "Random" ~ "#9ca3af",
      SampleID == "Lenti"  ~ "#c0392b",
      TRUE                 ~ "#1a1a2e"
    ),
    sig_label = case_when(
      SampleID == "Random"               ~ "",
      p.value < 0.001                    ~ "***",
      p.value < 0.01                     ~ "**",
      p.value < 0.05                     ~ "*",
      TRUE                               ~ "ns"
    )
  )

# ── Prepare Panel B data from summary_open_vs_random ─────────────────────────

pB_data <- summary_open_vs_random %>%
  mutate(
    ci_lo     = mean_odds_ratio - 1.96 * sd_odds_ratio,
    ci_hi     = mean_odds_ratio + 1.96 * sd_odds_ratio,
    SampleID  = factor(SampleID, levels = rev(sample_order)),
    dot_color = case_when(
      SampleID == "Lenti" ~ "#c0392b",
      TRUE                ~ "#1a1a2e"
    ),
    sig_label = case_when(
      n_significant == 100 ~ "***",
      n_significant >= 95  ~ "**",
      n_significant >= 80  ~ "*",
      TRUE                 ~ "ns"
    )
  )

# ── Shared publication theme ──────────────────────────────────────────────────

pub_theme <- theme_classic(base_size = 10, base_family = "Helvetica") +
  theme(
    axis.text          = element_text(size = 9, color = "black"),
    axis.title         = element_text(size = 9, color = "black"),
    axis.line          = element_line(color = "black", linewidth = 0.4),
    axis.ticks         = element_line(color = "black", linewidth = 0.4),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title         = element_text(size = 10, face = "bold", color = "black", margin = margin(b = 4)),
    plot.subtitle      = element_text(size = 8, color = "grey40", margin = margin(b = 10)),
    legend.position    = "none",
    plot.margin        = margin(10, 16, 10, 10)
  )

# ── Panel A: enrichment vs genome background ──────────────────────────────────

x_max_A <- max(pA_data$enrichment, na.rm = TRUE) * 1.15

pA <- ggplot(pA_data, aes(y = SampleID)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_segment(aes(x = 1, xend = enrichment, yend = SampleID, color = dot_color),
               linewidth = 0.6) +
  geom_point(aes(x = enrichment, color = dot_color,
                 shape = ifelse(SampleID == "Lenti", 18,
                                ifelse(SampleID == "Random", 15, 16)),
                 size  = ifelse(SampleID == "Lenti", 4, 2.8))) +
  geom_text(aes(x = enrichment + (x_max_A * 0.02), label = sig_label, color = dot_color),
            hjust = 0, size = 3, fontface = "bold") +
  geom_hline(yintercept = 1.5, linetype = "dotted", color = "grey70", linewidth = 0.3) +
  scale_color_identity() +
  scale_size_identity() +
  scale_shape_identity() +
  scale_y_discrete(labels = c(sample_labels, "Random" = "Random")) +
  scale_x_continuous(
    limits = c(0.8, x_max_A),
    expand = expansion(mult = c(0.02, 0.1))
  ) +
  labs(
    title    = "Open chromatin enrichment vs genome background",
    subtitle = "Enrichment = observed / expected overlaps (permutation test, n = 1,000, \u00b1100bp window)",
    x        = "Enrichment (obs / exp)",
    y        = NULL
  ) +
  pub_theme

# ── Panel B: OR vs Random IS ──────────────────────────────────────────────────

x_max_B <- max(pB_data$ci_hi, na.rm = TRUE) * 1.12

pB <- ggplot(pB_data, aes(y = SampleID)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  annotate("rect", xmin = 1, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#dbeafe", alpha = 0.35) +
  geom_segment(aes(x = ci_lo, xend = ci_hi, yend = SampleID, color = dot_color),
               linewidth = 0.6) +
  geom_point(aes(x = mean_odds_ratio, color = dot_color,
                 shape = ifelse(SampleID == "Lenti", 18, 16),
                 size  = ifelse(SampleID == "Lenti", 4, 2.8))) +
  geom_text(aes(x = ci_hi + (x_max_B * 0.02), label = sig_label, color = dot_color),
            hjust = 0, size = 3, fontface = "bold") +
  scale_color_identity() +
  scale_size_identity() +
  scale_shape_identity() +
  scale_y_discrete(labels = sample_labels) +
  scale_x_continuous(
    limits = c(NA, x_max_B),
    expand = expansion(mult = c(0.05, 0.1))
  ) +
  labs(
    title    = "Open chromatin enrichment vs Random IS",
    subtitle = "Mean OR \u00b1 95% CI across 100 subsampling replicates (n = 3,500, BH correction)",
    x        = "Odds Ratio",
    y        = NULL
  ) +
  pub_theme

# ── Combine and save ──────────────────────────────────────────────────────────

pdf(file = paste(source_folder, analysis_folder_date, ".BIO.OpenChromatin.GW.pdf", sep = ""), height=5, width=6)
plot(pA)
dev.off()
png(file = paste(source_folder, analysis_folder_date, ".BIO.OpenChromatin.GW.png", sep = ""), height=5, width=6, units = "in", res = 300)
plot(pA)
dev.off()

pdf(file = paste(source_folder, analysis_folder_date, ".BIO.OpenChromatin.vsRandom.pdf", sep = ""), height=5, width=6)
plot(pB)
dev.off()
png(file = paste(source_folder, analysis_folder_date, ".BIO.OpenChromatin.vsRandom.png", sep = ""), height=5, width=6, units = "in", res = 300)
plot(pB)
dev.off()


# fig <- pA | pB
# 
# ggsave(
#   filename = "forest_plot_openchromatin.pdf",
#   plot     = fig,
#   width    = 8.0,
#   height   = 3.6,
#   units    = "in",
#   device   = cairo_pdf
# )
# 
# ggsave(
#   filename = "forest_plot_openchromatin.png",
#   plot     = fig,
#   width    = 8.0,
#   height   = 3.6,
#   units    = "in",
#   dpi      = 300
# )
# 
# message("Saved: forest_plot_openchromatin.pdf and .png")
# 

###############################################################
## GWIF genome wide integration frequency
###############################################################

### any difference between LV and Random?
# allIS_lv <- sample_n(IS_gdf_molten[which(IS_gdf_molten$SampleID == "Lenti"), ], subsampling_n_is)
allIS_lv <- IS_gdf_molten[which(IS_gdf_molten$SampleID == "Lenti"), ]
# allIS_cr1 <- sample_n(IS_gdf_molten[which(IS_gdf_molten$SampleID == "CR1"), ], subsampling_n_is)
allIS_cr1 <- IS_gdf_molten[which(IS_gdf_molten$SampleID == "CR1"), ]
# allIS_rte3 <- sample_n(IS_gdf_molten[which(IS_gdf_molten$SampleID == "RTE3"), ], subsampling_n_is)
allIS_rte3 <- IS_gdf_molten[which(IS_gdf_molten$SampleID == "RTE3"), ]
# allIS_rte25 <- sample_n(IS_gdf_molten[which(IS_gdf_molten$SampleID == "RTE25"), ], subsampling_n_is)
allIS_rte25 <- IS_gdf_molten[which(IS_gdf_molten$SampleID == "RTE25"), ]
# allIS_vingi <- sample_n(IS_gdf_molten[which(IS_gdf_molten$SampleID == "Vingi"), ], subsampling_n_is)
allIS_vingi <- IS_gdf_molten[which(IS_gdf_molten$SampleID == "Vingi"), ]
# allIS_random already exists (see above)

id_cols_min <- c("chr", "integration_locus", "integration_strand", "GeneName", "GeneStrand")

allvectors_fisherFreq <- NULL

# even with subsampling vs Random
# RTE3
rte3_vs_random <- compareGeneFrequency_Fisher(df_g1 = allIS_random_all[colnames(allIS_rte3)], 
                                              df_g2 = allIS_rte3, 
                                              min_is_per_gene = 1)
rte3_vs_random <- rte3_vs_random %>%
  mutate(Prevalence = case_when(
    FC_G1G2 < 1 ~ "RTE3",
    FC_G1G2 == 1 ~ "Neutral",
    FC_G1G2 > 1 ~ "Random"
  ))
rte3_vs_random <- merge(x = rte3_vs_random, y = msk_db, by = c("GeneName"), all.x = T)
rte3_vs_random$KnownGeneClass <- ifelse(!is.na(rte3_vs_random$`Gene Type`), 
                                        (ifelse((rte3_vs_random$`Gene Type`) %in% c("ONCOGENE", "ONCOGENE_AND_TSG"), 
                                                "OncoGene", "TumSuppressor")), 
                                        "Other")
rte3_vs_random[is.na(rte3_vs_random)] <- NA
rte3_vs_random$positive_outlier_and_significant <- ifelse((!is.na(rte3_vs_random$FisherTest_pvalue_fdr) & rte3_vs_random$FisherTest_pvalue_fdr < 0.05), TRUE, FALSE)
rte3_vs_random$Vector <- "RTE3"

# LV
lv_vs_random <- compareGeneFrequency_Fisher(df_g1 = allIS_random_all[colnames(allIS_lv)], 
                                            df_g2 = allIS_lv, 
                                            min_is_per_gene = 1)
lv_vs_random <- lv_vs_random %>%
  mutate(Prevalence = case_when(
    FC_G1G2 < 1 ~ "Lenti",
    FC_G1G2 == 1 ~ "Neutral",
    FC_G1G2 > 1 ~ "Random"
  ))
lv_vs_random <- merge(x = lv_vs_random, y = msk_db, by = c("GeneName"), all.x = T)
lv_vs_random$KnownGeneClass <- ifelse(!is.na(lv_vs_random$`Gene Type`), 
                                      (ifelse((lv_vs_random$`Gene Type`) %in% c("ONCOGENE", "ONCOGENE_AND_TSG"), 
                                              "OncoGene", "TumSuppressor")), 
                                      "Other")
lv_vs_random[is.na(lv_vs_random)] <- NA
lv_vs_random$positive_outlier_and_significant <- ifelse((!is.na(lv_vs_random$FisherTest_pvalue_fdr) & lv_vs_random$FisherTest_pvalue_fdr < 0.05), TRUE, FALSE)
lv_vs_random$Vector <- "Lenti"

# RTE25
rte25_vs_random <- compareGeneFrequency_Fisher(df_g1 = allIS_random_all[colnames(allIS_rte25)], 
                                               df_g2 = allIS_rte25,
                                               min_is_per_gene = 1)
rte25_vs_random <- rte25_vs_random %>%
  mutate(Prevalence = case_when(
    FC_G1G2 < 1 ~ "RTE25",
    FC_G1G2 == 1 ~ "Neutral",
    FC_G1G2 > 1 ~ "Random"
  ))
rte25_vs_random <- merge(x = rte25_vs_random, y = msk_db, by = c("GeneName"), all.x = T)
rte25_vs_random$KnownGeneClass <- ifelse(!is.na(rte25_vs_random$`Gene Type`), 
                                         (ifelse((rte25_vs_random$`Gene Type`) %in% c("ONCOGENE", "ONCOGENE_AND_TSG"), 
                                                 "OncoGene", "TumSuppressor")), 
                                         "Other")
rte25_vs_random[is.na(rte25_vs_random)] <- NA
rte25_vs_random$positive_outlier_and_significant <- ifelse((!is.na(rte25_vs_random$FisherTest_pvalue_fdr) & rte25_vs_random$FisherTest_pvalue_fdr < 0.05), TRUE, FALSE)
rte25_vs_random$Vector <- "RTE25"

# Vingi
vingi_vs_random <- compareGeneFrequency_Fisher(df_g1 = allIS_random_all[colnames(allIS_vingi)], df_g2 = allIS_vingi,
                                               min_is_per_gene = 1)
vingi_vs_random <- vingi_vs_random %>%
  mutate(Prevalence = case_when(
    FC_G1G2 < 1 ~ "Vingi",
    FC_G1G2 == 1 ~ "Neutral",
    FC_G1G2 > 1 ~ "Random"
  ))
vingi_vs_random <- merge(x = vingi_vs_random, y = msk_db, by = c("GeneName"), all.x = T)
vingi_vs_random$KnownGeneClass <- ifelse(!is.na(vingi_vs_random$`Gene Type`), 
                                         (ifelse((vingi_vs_random$`Gene Type`) %in% c("ONCOGENE", "ONCOGENE_AND_TSG"), 
                                                 "OncoGene", "TumSuppressor")), 
                                         "Other")
vingi_vs_random[is.na(vingi_vs_random)] <- NA
vingi_vs_random$positive_outlier_and_significant <- ifelse((!is.na(vingi_vs_random$FisherTest_pvalue_fdr) & vingi_vs_random$FisherTest_pvalue_fdr < 0.05), TRUE, FALSE)
vingi_vs_random$Vector <- "Vingi"

# CR1
cr1_vs_random <- compareGeneFrequency_Fisher(df_g1 = allIS_random_all[colnames(allIS_cr1)], df_g2 = allIS_cr1, 
                                             min_is_per_gene = 1)
cr1_vs_random <- cr1_vs_random %>%
  mutate(Prevalence = case_when(
    FC_G1G2 < 1 ~ "CR1",
    FC_G1G2 == 1 ~ "Neutral",
    FC_G1G2 > 1 ~ "Random"
  ))
cr1_vs_random <- merge(x = cr1_vs_random, y = msk_db, by = c("GeneName"), all.x = T)
cr1_vs_random$KnownGeneClass <- ifelse(!is.na(cr1_vs_random$`Gene Type`), 
                                       (ifelse((cr1_vs_random$`Gene Type`) %in% c("ONCOGENE", "ONCOGENE_AND_TSG"), 
                                               "OncoGene", "TumSuppressor")), 
                                       "Other")
cr1_vs_random[is.na(cr1_vs_random)] <- NA
cr1_vs_random$positive_outlier_and_significant <- ifelse((!is.na(cr1_vs_random$FisherTest_pvalue_fdr) & cr1_vs_random$FisherTest_pvalue_fdr < 0.05), TRUE, FALSE)
cr1_vs_random$Vector <- "CR1"

# now combin all data
significance_threshold_minus_log_p <- -log(0.05, base = 10)
allvectors_fisherFreq <- do.call(rbind, list(rte25_vs_random, rte3_vs_random, cr1_vs_random, vingi_vs_random, lv_vs_random))
# allvectors_fisherFreq <- do.call(rbind, list(vingi_vs_random, lv_vs_random))
allvectors_fisherFreq$minlog_FisherTest_pvalue <- -log(x = allvectors_fisherFreq$FisherTest_pvalue, base = 10)

allvectors_fisherFreq <- merge(x = allvectors_fisherFreq, y = vector_rename, by = "Vector", all.x=T)
faceting_order_vector_extended <- c("Random", "RTE-25_Lmi", "RTE-3_BF", "CR1-1_PH", "Vingi-1_Acar", "Lenti")
allvectors_fisherFreq$Vector_Order_Extended <- factor(allvectors_fisherFreq$ExtendedVector, levels = faceting_order_vector_extended)

plot_isprofile_all_vs_random <- 
  ggplot(data = allvectors_fisherFreq, aes( # y = log_pvalue_fdr, 
    x = log_FC,
    y = log_pvalue_fdr,
    # y = minlog_FisherTest_pvalue,
    color = Prevalence, fill = Prevalence), na.rm = T, se = TRUE) +
    # color = FisherTest_pvalue_significant, fill = FisherTest_pvalue_significant), na.rm = T, se = TRUE) +
  # ), na.rm = T, se = TRUE) +
  geom_point(size = 3, alpha = .5) +

  scale_color_manual(values = c("orange", "red4", "gray80", "green4", "blue", "purple")) +
  # facet_grid( . ~ Vector) +
  facet_grid( . ~ Vector_Order_Extended) +
  geom_label_repel(
    data = subset(allvectors_fisherFreq, FisherTest_pvalue < 0.05 & n_IS_perGene_G2 >= 3),
    aes(label = GeneName, fill = KnownGeneClass),
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    color = "black",
    segment.color = "black"
    # show.legend = FALSE
  ) + 
  scale_fill_manual(
    values = c(
      OncoGene = "firebrick2",
      TumSuppressor   = "orange",
      Other    = "white"
    )
  ) + 
  theme_bw() +
  theme(strip.text.y = element_text(size = 16, colour = "gray50", angle = 0), strip.text.x = element_text(size = 16, colour = "black", angle = 0)) +
  theme(strip.text = element_text(face="bold", size=16)) +
  # theme(strip.background = element_rect(fill="darkblue", colour="white", size=1)) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(legend.position="none") +
  labs(title = paste("Gene-based Integration frequency - All vectors vs Random"), 
       y = "P-value Fisher test (-log10(p))", 
       x = "Integration frequency (log10)", 
       color = "Vector", fill = "Oncogene/TSG",
       subtitle = paste0("Genes labeled if P-value < 0.05 post FDR correction.")) 

png(file = paste(source_folder, analysis_folder_date, ".GeneFreq.AllVectors_vs_Random.png", sep = ""), height=6, width=14, units = "in", res = 300)
plot_isprofile_all_vs_random
dev.off()
pdf(file = paste(source_folder, analysis_folder_date, ".GeneFreq.AllVectors_vs_Random.pdf", sep = ""), height=6, width=14)
plot_isprofile_all_vs_random
dev.off()



###############################################################
## Frequency of integration against LV
###############################################################

# RTE3
rte3_vs_lv <- compareGeneFrequency_Fisher(df_g1 = allIS_lv, 
                                          df_g2 = allIS_rte3, 
                                          min_is_per_gene = 1)
rte3_vs_lv <- rte3_vs_lv %>%
  mutate(Prevalence = case_when(
    FC_G1G2 < 1 ~ "RTE3",
    FC_G1G2 == 1 ~ "Neutral",
    FC_G1G2 > 1 ~ "LV"
  ))
rte3_vs_lv <- merge(x = rte3_vs_lv, y = msk_db, by = c("GeneName"), all.x = T)
rte3_vs_lv$KnownGeneClass <- ifelse(!is.na(rte3_vs_lv$`Gene Type`), 
                                        (ifelse((rte3_vs_lv$`Gene Type`) %in% c("ONCOGENE", "ONCOGENE_AND_TSG"), 
                                                "OncoGene", "TumSuppressor")), 
                                        "Other")
rte3_vs_lv[is.na(rte3_vs_lv)] <- NA
rte3_vs_lv$positive_outlier_and_significant <- ifelse((!is.na(rte3_vs_lv$FisherTest_pvalue_fdr) & rte3_vs_lv$FisherTest_pvalue_fdr < 0.05), TRUE, FALSE)
rte3_vs_lv$Vector <- "RTE3"

# RTE25
rte25_vs_lv <- compareGeneFrequency_Fisher(df_g1 = allIS_lv, 
                                           df_g2 = allIS_rte25,
                                           min_is_per_gene = 1)
rte25_vs_lv <- rte25_vs_lv %>%
  mutate(Prevalence = case_when(
    FC_G1G2 < 1 ~ "RTE25",
    FC_G1G2 == 1 ~ "Neutral",
    FC_G1G2 > 1 ~ "LV"
  ))
rte25_vs_lv <- merge(x = rte25_vs_lv, y = msk_db, by = c("GeneName"), all.x = T)
rte25_vs_lv$KnownGeneClass <- ifelse(!is.na(rte25_vs_lv$`Gene Type`), 
                                         (ifelse((rte25_vs_lv$`Gene Type`) %in% c("ONCOGENE", "ONCOGENE_AND_TSG"), 
                                                 "OncoGene", "TumSuppressor")), 
                                         "Other")
rte25_vs_lv[is.na(rte25_vs_lv)] <- NA
rte25_vs_lv$positive_outlier_and_significant <- ifelse((!is.na(rte25_vs_lv$FisherTest_pvalue_fdr) & rte25_vs_lv$FisherTest_pvalue_fdr < 0.05), TRUE, FALSE)
rte25_vs_lv$Vector <- "RTE25"

# Vingi
vingi_vs_lv <- compareGeneFrequency_Fisher(df_g1 = allIS_lv, 
                                           df_g2 = allIS_vingi,
                                           min_is_per_gene = 1)
vingi_vs_lv <- vingi_vs_lv %>%
  mutate(Prevalence = case_when(
    FC_G1G2 < 1 ~ "Vingi",
    FC_G1G2 == 1 ~ "Neutral",
    FC_G1G2 > 1 ~ "LV"
  ))
vingi_vs_lv <- merge(x = vingi_vs_lv, y = msk_db, by = c("GeneName"), all.x = T)
vingi_vs_lv$KnownGeneClass <- ifelse(!is.na(vingi_vs_lv$`Gene Type`), 
                                         (ifelse((vingi_vs_lv$`Gene Type`) %in% c("ONCOGENE", "ONCOGENE_AND_TSG"), 
                                                 "OncoGene", "TumSuppressor")), 
                                         "Other")
vingi_vs_lv[is.na(vingi_vs_lv)] <- NA
vingi_vs_lv$positive_outlier_and_significant <- ifelse((!is.na(vingi_vs_lv$FisherTest_pvalue_fdr) & vingi_vs_lv$FisherTest_pvalue_fdr < 0.05), TRUE, FALSE)
vingi_vs_lv$Vector <- "Vingi"

# CR1
cr1_vs_lv <- compareGeneFrequency_Fisher(df_g1 = allIS_lv, 
                                         df_g2 = allIS_cr1, 
                                         min_is_per_gene = 1)
cr1_vs_lv <- cr1_vs_lv %>%
  mutate(Prevalence = case_when(
    FC_G1G2 < 1 ~ "CR1",
    FC_G1G2 == 1 ~ "Neutral",
    FC_G1G2 > 1 ~ "LV"
  ))
cr1_vs_lv <- merge(x = cr1_vs_lv, y = msk_db, by = c("GeneName"), all.x = T)
cr1_vs_lv$KnownGeneClass <- ifelse(!is.na(cr1_vs_lv$`Gene Type`), 
                                       (ifelse((cr1_vs_lv$`Gene Type`) %in% c("ONCOGENE", "ONCOGENE_AND_TSG"), 
                                               "OncoGene", "TumSuppressor")), 
                                       "Other")
cr1_vs_lv[is.na(cr1_vs_lv)] <- NA
cr1_vs_lv$positive_outlier_and_significant <- ifelse((!is.na(cr1_vs_lv$FisherTest_pvalue_fdr) & cr1_vs_lv$FisherTest_pvalue_fdr < 0.05), TRUE, FALSE)
cr1_vs_lv$Vector <- "CR1"

# now combin all data
significance_threshold_minus_log_p <- -log(0.05, base = 10)
allvectors_fisherFreq_vs_lv <- do.call(rbind, list(rte25_vs_lv, rte3_vs_lv, cr1_vs_lv, vingi_vs_lv))
# allvectors_fisherFreq_vs_lv <- do.call(rbind, list(vingi_vs_lv, lv_vs_lv))
allvectors_fisherFreq_vs_lv$minlog_FisherTest_pvalue <- -log(x = allvectors_fisherFreq_vs_lv$FisherTest_pvalue, base = 10)

allvectors_fisherFreq_vs_lv <- merge(x = allvectors_fisherFreq_vs_lv, y = vector_rename, by = "Vector", all.x=T)
faceting_order_vector_extended <- c("Random", "RTE-25_Lmi", "RTE-3_BF", "CR1-1_PH", "Vingi-1_Acar", "Lenti")
allvectors_fisherFreq_vs_lv$Vector_Order_Extended <- factor(allvectors_fisherFreq_vs_lv$ExtendedVector, levels = faceting_order_vector_extended)

plot_isprofile_all_vs_lv <- 
  ggplot(data = allvectors_fisherFreq_vs_lv, aes( # y = log_pvalue_fdr, 
    x = log_FC,
    y = log_pvalue_fdr,
    # y = minlog_FisherTest_pvalue,
    color = Prevalence, fill = Prevalence), na.rm = T, se = TRUE) +
  # color = FisherTest_pvalue_significant, fill = FisherTest_pvalue_significant), na.rm = T, se = TRUE) +
  # ), na.rm = T, se = TRUE) +
  geom_point(size = 3, alpha = .5) +
  scale_color_manual(values = c(t(tessera_color_schema[1,]))) +
  # scale_fill_manual(values = c(t(tessera_color_schema[1,]))) +
  # scale_color_manual(values = c("orange", "red4", "gray80", "green4", "blue", "purple")) +
  # facet_grid( . ~ Vector) +
  facet_grid( . ~ Vector_Order_Extended) +
  geom_label_repel(
    data = subset(allvectors_fisherFreq_vs_lv, FisherTest_pvalue < 0.05),
    aes(label = GeneName, fill = KnownGeneClass),
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    color = "black",
    segment.color = "black"
    # show.legend = FALSE
  ) + 
  scale_fill_manual(
    values = c(
      OncoGene = "firebrick2",
      TumSuppressor   = "orange",
      Other    = "white"
    )
  ) + 
  theme_bw() +
  theme(strip.text.y = element_text(size = 16, colour = "gray50", angle = 0), strip.text.x = element_text(size = 16, colour = "black", angle = 0)) +
  theme(strip.text = element_text(face="bold", size=16)) +
  # theme(strip.background = element_rect(fill="darkblue", colour="white", size=1)) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(legend.position="none") +
  labs(title = paste("Gene-based Integration frequency - All vectors vs LV"), 
       y = "P-value Fisher test (-log10(p))", 
       x = "Integration frequency (log10)", 
       color = "Vector", fill = "Oncogene/TSG",
       subtitle = paste0("Genes labeled if P-value < 0.05 post FDR correction.")) 

png(file = paste(source_folder, analysis_folder_date, ".GeneFreq.AllVectors_vs_lv.png", sep = ""), height=6, width=14, units = "in", res = 300)
plot_isprofile_all_vs_lv
dev.off()
pdf(file = paste(source_folder, analysis_folder_date, ".GeneFreq.AllVectors_vs_lv.pdf", sep = ""), height=6, width=14)
plot_isprofile_all_vs_lv
dev.off()




###############################################################
## CIS
###############################################################

# -------- CIS --------- hotspots -----
cis_df <- NULL
basic_id_cols <- c("chr", "integration_locus", "integration_strand", "GeneName", "GeneStrand")
subsampling_n_is_forCIS <- 3500
for (vec in setdiff(colnames(IS_gdf), uni_id_cols)) {
  message(paste("[AP]\tProcessing:", vec))
  vec_IS_mat <- IS_gdf[c(basic_id_cols, vec)]
  vec_IS_mat <- compactDfByRows(df = vec_IS_mat, data_columns = vec, annotation_columns = basic_id_cols)
  vec_IS_mat <- sample_n(vec_IS_mat, subsampling_n_is_forCIS)
  vec_IS_mat$chr <- gsub("chr", "", vec_IS_mat$chr)
  names(vec_IS_mat) <- c("chr", "integration_locus", "strand", "GeneName", "GeneStrand", vec)
  # compute CIS
  # NOTE: Provide path to gene-based annotation file if available
  vec_cis <- CISGrubbs(df = vec_IS_mat,
                       annotation_cols = c("chr", "integration_locus", "strand", "GeneName", "GeneStrand", "GeneDistance"),
                       genomic_annotation_genebased_file = NULL)
  vec_cis$Assay <- "Uniseq"
  vec_cis$Vector <- vec
  write.table(x = vec_cis, 
              file = gzfile(paste(source_folder, analysis_folder_date, ".CIS.", vec, ".tsv.gz", sep = "")), 
              sep = "\t", quote = FALSE, row.names = F, col.names = T, na = '')
  # subset only the potive cis genes
  vec_cis_topgenes <- vec_cis[which(!is.na(vec_cis$tdist_positive_and_corrected)), c("chr", "integration_locus", "integration_locus", "GeneName", "n_IS_perGene")]
  if (nrow(vec_cis_topgenes) > 0) {
  write.xlsx(x = vec_cis_topgenes, 
             file = paste(source_folder, analysis_folder_date, ".CIS.", vec, ".positivegenes.xlsx", sep = ""), 
             rowNames = F)
  }
  if (!is.null(dim(vec_cis))) {
    # create a full data df
    if (!is.null(dim(cis_df))) {
      cis_df <- rbind(cis_df, vec_cis)
    } else {
      cis_df <- vec_cis
    } # if (length(cis_df) > 0)
  } # if (!(dim(mat_df) == NULL))
}

# # the random
# allIS_random_all$strand <- allIS_random_all$integration_strand
# random_cis <- CISGrubbs(df = allIS_random_all, 
#                      annotation_cols = c("chr", "integration_locus", "strand", "GeneName", "GeneStrand", "GeneDistance"), 
#                      genomic_annotation_genebased_file = "/Users/andreacalabria/Library/CloudStorage/OneDrive-VL58/Utils/GitHub/rna-writer_prototype/metadata/hg38.refGene.oracle.tsv.gz")
# vec_cis$Assay <- "Uniseq"
# vec_cis$Vector <- "Random"
# write.table(x = vec_cis, 
#             file = gzfile(paste(source_folder, analysis_folder_date, ".CIS.", "Random", ".tsv.gz", sep = "")), 
#             sep = "\t", quote = FALSE, row.names = F, col.names = T, na = '')
# # subset only the potive cis genes
# vec_cis_topgenes <- vec_cis[which(!is.na(vec_cis$tdist_positive_and_corrected)), c("chr", "integration_locus", "integration_locus", "GeneName", "n_IS_perGene")]

# add annotations
cis_df$mlog10_tdist_fdr <- -log(cis_df$tdist_fdr, base = 10)
cis_df <- merge(x = cis_df, y = msk_db, by = c("GeneName"), all.x = T)
cis_df$KnownGeneClass <- ifelse(!is.na(cis_df$`Gene Type`), 
                                (ifelse((cis_df$`Gene Type`) %in% c("ONCOGENE", "ONCOGENE_AND_TSG"), 
                                        "OncoGene", "TumSuppressor")), 
                                "Other")
cis_df[is.na(cis_df)] <- NA
cis_df$positive_outlier_and_significant <- ifelse((!is.na(cis_df$tdist_fdr) & cis_df$tdist_fdr < 0.05), TRUE, FALSE)

cis_df$Vector <- factor(cis_df$Vector, levels = faceting_order_vector)
cis_df <- merge(x = cis_df, y = vector_rename, by = "Vector", all.x=T)
faceting_order_vector_extended <- c("Random", "RTE-25_Lmi", "RTE-3_BF", "CR1-1_PH", "Vingi-1_Acar", "Lenti")
cis_df$Vector_Order_Extended <- factor(cis_df$ExtendedVector, levels = faceting_order_vector_extended)

write.table(x = cis_df,
            file = gzfile(paste(source_folder, analysis_folder_date, ".CIS.", "AllVectorsMerged", ".tsv.gz", sep = "")),
            sep = "\t", quote = FALSE, row.names = F, col.names = T, na = '')

# cis_df <- cis_df[which(cis_df$average_TxLen>300 & cis_df$n_IS_perGene>=3),] # remove annotations of genomic elements too short (not real genes, like MIR), and min IS per gene = 2
cis_df <- cis_df[which(cis_df$average_TxLen>300),] # remove annotations of genomic elements too short (not real genes, like MIR), and min IS per gene = 2


# get top hit genes
cis_df_topgenes <- cis_df[order(cis_df$geneIS_frequency_byHitIS, decreasing = T), 
                          c("chr", "integration_locus", "integration_locus", "GeneName", "geneIS_frequency_byHitIS")]
cis_df_topgenes_slice <- unique(cis_df_topgenes[1:30, c("GeneName")])

# cis_df_topgenes <- cis_df %>% 
#   group_by(chr, GeneName) %>%
#   arrange(desc(geneIS_frequency_byHitIS), .by_group = TRUE) %>%
#   filter(row_number()==1)
# 
# cis_df_topgenes <- cis_df[which(!is.na(cis_df$tdist_positive_and_corrected)), c("chr", "integration_locus", "integration_locus", "GeneName", "n_IS_perGene")]
# bed_rte1_d_cis_topgenes <- rte1_d_cis_topgenes
# bed_rte1_d_cis_topgenes$strand <- "+"
# names(bed_rte1_d_cis_topgenes) <- bed_colnames
# write.xlsx(x = rte1_d_cis_topgenes, 
#            file = paste(source_folder, analysis_folder_date, ".CIS.RTE1.d.positivegenes.xlsx", sep = ""), 
#            rowNames = F)



# plot CIS volcano facet vector
plot_cis_byvector <- 
  ggplot(data = cis_df, 
         aes(x = neg_zscore_minus_log2_integration_freq_withtolerance,
             y =  mlog10_tdist_fdr
             # color = KnownGeneClass, fill = KnownGeneClass
         ), na.rm = T, se = TRUE) +
  geom_point(size = 3, alpha = .5) +
  # scale_color_manual(values = c("red", "green", "blue")) +
  # scale_size(name = "Stdev P-value") + 
  facet_wrap( . ~ Vector_Order_Extended, ncol = 3) +
  geom_hline(yintercept = significance_threshold_minus_log_p, color='firebrick4', size=1, show.legend = T, linetype="dashed") +
  # scale_x_continuous(breaks = seq(-4, 4,2)) +
  geom_label_repel(
    # data = subset(cis_df, ( tdist_fdr < 0.05 & n_IS_perGene >= 3 & average_TxLen > 300) ),
    data = subset(cis_df, ( (tdist_fdr < 0.05 | mlog10_tdist_fdr > 0.5) & n_IS_perGene >= 3 & average_TxLen > 300) ),
    aes(label = GeneName, fill = KnownGeneClass),
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    color = 'black',
    segment.color = 'black',
    # max.overlaps = Inf
  ) + 
  scale_fill_manual(
    values = c(
      OncoGene = "firebrick2",
      TumSuppressor   = "orange",
      Other    = "white"
    )
  ) + 
  theme(strip.text.y = element_text(size = 16, colour = "black", angle = 0), 
        strip.text.x = element_text(size = 16, colour = "blue", angle = 0)) +
  theme(strip.text = element_text(face="bold", size=16)) +
  # theme(strip.background = element_rect(fill="darkblue", colour="white", size=1)) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste("CIS genes"), 
       y = "P-value Grubbs test (-log10(p))", 
       x = "Integration frequency (log2)", 
       # color = "Gene class", 
       subtitle = paste0("P-value < 0.05 (FDR adjusted, -log(p)=", round(significance_threshold_minus_log_p, 3), ").\nAbove the dashed red line, CIS genes passing FDR correction.")) 

png(file = paste(source_folder, analysis_folder_date, ".CISgenes.png", sep = ""), height=6, width=9, units = "in", res = 300)
plot_cis_byvector
dev.off()
pdf(file = paste(source_folder, analysis_folder_date, ".CISgenes.pdf", sep = ""), height=6, width=9)
plot_cis_byvector
dev.off()


# plot CIS volcano, all together
plot_cis_byvector_allinone_v1 <- 
  ggplot(data = cis_df, 
         aes(x = neg_zscore_minus_log2_integration_freq_withtolerance,
             y =  mlog10_tdist_fdr,
             color = KnownGeneClass, fill = KnownGeneClass
         ), na.rm = T, se = TRUE) +
  geom_point(size = 3, alpha = .5) +
  # scale_color_manual(values = c("red", "green", "blue")) +
  # scale_size(name = "Stdev P-value") + 
  facet_wrap( . ~ Vector_Order_Extended, ncol = 3) +
  # scale_x_continuous(breaks = seq(-4, 4,2)) +
  # geom_hline(yintercept = significance_threshold_minus_log_p, color='navy', size=1, show.legend = T, linetype="dashed") +
  geom_label_repel(
    # data = subset(cis_df, (tdist_positive_and_corrected >0 & n_IS_perGene > 3 & average_TxLen > 300) ),
    data = subset(cis_df, (tdist_fdr <0.05 & n_IS_perGene >= 3 & average_TxLen > 300 ) ),
    aes(label = GeneName), 
    # aes(label = GeneName), 
    # size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    color = 'white',
    segment.color = 'black',
    max.overlaps = nrow(cis_df)
  ) + 
  # geom_label_repel(
  #   data = subset(cis_df, (tdist_fdr <0.1 & !is.na(`Gene Type`) & n_IS_perGene >= 3 & average_TxLen > 300 )),
  #   aes(label = GeneName, fill = KnownGeneClass),
  #   # aes(label = GeneName),
  #   # size = 5,
  #   box.padding = unit(0.35, "lines"),
  #   point.padding = unit(0.3, "lines"),
  #   color = 'white',
  #   segment.color = 'black'
  # ) +
  theme_bw() +
  theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0), strip.text.x = element_text(size = 16, colour = "blue", angle = 0)) +
  theme(strip.text = element_text(face="bold", size=16)) +
  # theme(strip.background = element_rect(fill="darkblue", colour="white", size=1)) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste("CIS genes"), 
       y = "P-value Grubbs test (-log10(p))", 
       x = "Integration frequency (log2)", 
       # color = "Gene class", 
       subtitle = paste0("P-value < 0.05 (FDR adjusted, -log(p)=", round(significance_threshold_minus_log_p, 3), ").")) 

png(file = paste(source_folder, analysis_folder_date, ".CISgenes.PaperFigTmp.v1.png", sep = ""), height=10, width=16, units = "in", res = 300)
plot_cis_byvector_allinone_v1
dev.off()
pdf(file = paste(source_folder, analysis_folder_date, ".CISgenes.PaperFigTmp.v1.pdf", sep = ""), height=10, width=16)
plot_cis_byvector_allinone_v1
dev.off()

# plot CIS volcano, all together
plot_cis_byvector_allinone_v2 <- 
  ggplot(data = cis_df, 
         aes(x = neg_zscore_minus_log2_integration_freq_withtolerance,
             y =  mlog10_tdist_fdr,
             color = KnownGeneClass, fill = KnownGeneClass
         ), na.rm = T, se = TRUE) +
  geom_point(size = 3, alpha = .5) +
  # scale_color_manual(values = c("red", "green", "blue")) +
  # scale_size(name = "Stdev P-value") + 
  facet_wrap( . ~ Vector_Order_Extended, ncol = 3) +
  # scale_x_continuous(breaks = seq(-4, 4,2)) +
  # geom_hline(yintercept = significance_threshold_minus_log_p, color='navy', size=1, show.legend = T, linetype="dashed") +
  geom_label_repel(
    # data = subset(cis_df, (tdist_positive_and_corrected >0 & n_IS_perGene > 3 & average_TxLen > 300) ),
    data = subset(cis_df, ( (tdist_fdr < 0.05 | mlog10_tdist_fdr > 0.5) & n_IS_perGene >= 3 & average_TxLen > 300) ),
    aes(label = GeneName), 
    # aes(label = GeneName), 
    # size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    color = 'white',
    segment.color = 'black',
    max.overlaps = nrow(cis_df)
  ) + 
  # geom_label_repel(
  #   data = subset(cis_df, (tdist_fdr <0.1 & !is.na(`Gene Type`) )),
  #   aes(label = GeneName, fill = KnownGeneClass),
  #   # aes(label = GeneName),
  #   # size = 5,
  #   box.padding = unit(0.35, "lines"),
  #   point.padding = unit(0.3, "lines"),
  #   color = 'white',
  #   segment.color = 'black'
  # ) +
  theme_bw() +
  theme(strip.text.y = element_text(size = 16, colour = "black", angle = 0), strip.text.x = element_text(size = 16, colour = "black", angle = 0)) +
  theme(strip.text = element_text(face="bold", size=16)) +
  # theme(strip.background = element_rect(fill="darkblue", colour="white", size=1)) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste("CIS genes"), 
       y = "P-value Grubbs test (-log10(p))", 
       x = "Integration frequency (log2)", 
       # color = "Gene class", 
       subtitle = paste0("P-value < 0.05 (FDR adjusted, -log(p)=", round(significance_threshold_minus_log_p, 3), ").")) 

png(file = paste(source_folder, analysis_folder_date, ".CISgenes.PaperFigTmp.v2.png", sep = ""), height=10, width=16, units = "in", res = 300)
plot_cis_byvector_allinone_v2
dev.off()
pdf(file = paste(source_folder, analysis_folder_date, ".CISgenes.PaperFigTmp.v2.pdf", sep = ""), height=10, width=16)
plot_cis_byvector_allinone_v2
dev.off()


