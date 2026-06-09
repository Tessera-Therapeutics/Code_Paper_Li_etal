###### Init UNI-seq analysis ########

uni_id_cols <- c("chr", "integration_locus", "integration_strand", "GeneName", "GeneStrand", "GeneDistance")
old_id_cols <- c("refIS_chr", "refIS_locus", "integration_strand", "feature_geneName", "feature_strand", "feature_distance")

tessera_color_schema <- data.frame("txblue" = "#9999FF",
                                   ## "CD34" = "green4",
                                   "txred" = "#FF0066",
                                   "txpurple" = "#333399",
                                   "txgray" = "#999999",
                                   "orange" = "#FF9933",
                                   "red" = "#FF3333",
                                   "azure"  ="#3399FF")

paper_barplot_colors <- data.frame("ppgray" = "#a4a4a4",
                                   "ppblue" = "#3c3c93",
                                   "ppwater" = "#44acc3",
                                   "ppolivegreen" = "#9bbb56"
                                   )

vector_rename <- data.frame("Vector" = c("CR1", "RTE1", "RTE3", "RTE25", "Vingi", "Lenti", "Random"),
                            "ExtendedVector" = c("CR1-1_PH", "RTE-1_MD", "RTE-3_BF", "RTE-25_Lmi", "Vingi-1_Acar", "Lenti", "Random") )

# source("/Users/andreacalabria/Library/CloudStorage/OneDrive-VL58/Utils/GitHub/isatk/script/R/isa_utils_functions.R")
source("/Users/andreacalabria/Library/CloudStorage/OneDrive-VL58/Utils/GitHub/rna-writer_prototype/R/functions/isa_utils.R")

# load Master file
master_file <- "/Users/andreacalabria/Library/CloudStorage/OneDrive-VL58/Utils/GitHub/rna-writer_prototype/results/2023.MasterFileRun.Uniseq.xlsx"
uniseq_master_df <- read.xlsx(xlsxFile = master_file, sheet = "uniseq_sample_list")
vector_df <- read.xlsx(xlsxFile = master_file, sheet = "vector")
# uniseq_master_df <- merge(x = uniseq_master_df, y = vector_df, by = c("VectorName_BAM"), all.x = T)
# target_genome_GTF_file <- "/Users/andreacalabria/Library/CloudStorage/OneDrive-VL58/Project DISTALseq/metadata/hg38/hg38.refGene.sorted.transcript.gtf.gz"
target_genome_fa_file <- "/Users/andreacalabria/Library/CloudStorage/OneDrive-VL58/Genomes/hg38/index/hg38.p13.fa"
target_genome_GTF_file <- "/Users/andreacalabria/Library/CloudStorage/OneDrive-VL58/Project DISTALseq/metadata/hg38/hg38.ncbiRefSeq.sorted.transcript.gtf.gz" # hg38.ncbiRefSeq.sorted.transcript.gtf.gz
target_genome_TSS_file <- "/Users/andreacalabria/Library/CloudStorage/OneDrive-VL58/Project DISTALseq/metadata/hg38/hg38.ncbiRefSeq.sorted.TSS2.bed.gz" # hg38.ncbiRefSeq.sorted.TSS.bed.gz
target_genome_exons_file <- "/Users/andreacalabria/Library/CloudStorage/OneDrive-VL58/Project DISTALseq/metadata/hg38/hg38.ncbiRefSeq.sorted_slice_exons.gtf.gz" # hg38.ncbiRefSeq.sorted_slice_exons.gtf.gz
target_genome_blacklist_mappability_GTF <- "/Users/andreacalabria/Library/CloudStorage/OneDrive-VL58/Utils/GitHub/rna-writer_prototype/metadata/hg38.ucsc.encode_grc.lowmappability.2024.sort.gtf"
# target_genome_GTF_file <- "/Users/andreacalabria/Library/CloudStorage/OneDrive-VL58/Project DISTALseq/metadata/hg38/ucsc.hg38.refGene.2309.gtf"
# target_genome_GTF_file <- "/Users/andreacalabria/Library/CloudStorage/OneDrive-VL58/Project DISTALseq/metadata/hg38/hg38.ncbiRefSeq.sorted.slice_exon.gtf"
# faceting_order_ddpcr <- c("ddPCR_UTR5", "ddPCR_Transgene5p", "ddPCR_Transgenemid", "ddPCR_Transgene3p", "ddPCR_Promoter", "ddPCR_UTR3", "ddPCR_NoPrimer")

threshold_IS_span <- 10


# cancer genes
msk_db <- read.csv("/Users/andreacalabria/Library/CloudStorage/OneDrive-VL58/Utils/GitHub/isatk/script/R/publicdb/2601.msk.cancerGeneList.tsv", header=TRUE, fill=T, sep='\t', check.names = FALSE)
msk_db$GeneName <- msk_db$`Hugo Symbol` # just add a name that I am used to merge
msk_db <- msk_db[which((msk_db$`Gene Type` %in% c("ONCOGENE", "ONCOGENE_AND_TSG", "TSG")) & 
                         msk_db$`# of occurrence within resources (Column J-P)` > 1),] # avoid unknown genes and false positive / weak genes

