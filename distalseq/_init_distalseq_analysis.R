###### Init DISTAL-seq analysis ########

# token replacement file for RG tag update
dict_replacement_file <- "metadata/Vingi.Clones.knownIS.xlsx"
dict_clones <- read.xlsx(xlsxFile = dict_replacement_file)

ProjectID <- "DS"
# id_cols <- c("chr", "integration_locus", "strand", "GeneName", "GeneStrand")
id_cols <- c("chr", "integration_locus", "integration_strand", "GeneName", "GeneStrand", "GeneDistance")
old_id_cols <- c("refIS_chr", "refIS_locus", "integration_strand", "feature_geneName", "feature_strand", "feature_distance")

tessera_color_schema <- data.frame("txblue" = "#9999FF",
                                   ## "CD34" = "green4",
                                   "txred" = "#FF0066",
                                   "txpurple" = "#333399",
                                   "txgray" = "#999999",
                                   "orange" = "#FF9933",
                                   "red" = "#FF3333",
                                   "azure"  ="#3399FF")

# vector_cytoband_file <- "metadata/genomes/PLV5194.cytoband"
# vector_chr <- "PLV5194"
# zoomed_chr_index = 23
# species = "hg38"
# bp_res= 300
# color_vector_fwd="orange"
# color_vector_rev="green"
# color_target_rev="violet"
# color_target_fwd="blue"

template_string_retag <- "#!/bin/bash
cd \'TOKEN_FILEFOLDER\'

# generate sub BAM and index
picard FilterSamReads I=TOKEN_INPUTRAWBAM O=TOKEN_OUTPUT_SINGLE_IS_BAM READ_LIST_FILE=TOKEN_READLIST FILTER=includeReadList
samtools index TOKEN_OUTPUT_SINGLE_IS_BAM

# retag this file BAM
picard AddOrReplaceReadGroups I=TOKEN_OUTPUT_SINGLE_IS_BAM O=TOKEN_FINAL_SINGLE_IS RGID='TOKEN_KNOWNIS_ID' RGLB='TOKEN_CLONEID' RGPL='TOKEN_SAMPLEID' RGPU='TOKEN_CLONEID' RGSM='TOKEN_VECTORID'
samtools index TOKEN_FINAL_SINGLE_IS

mv TOKEN_OUTPUT_SINGLE_IS_BAM TOKEN_OUTPUT_SINGLE_IS_BAM.bai TOKEN_READLIST TOKEN_FINAL_SINGLE_IS TOKEN_FINAL_SINGLE_IS.bai TOKEN_OUTPUTFOLDER

"

# source("/Users/andreacalabria/Library/CloudStorage/OneDrive-VL58/Utils/GitHub/isatk/script/R/isa_utils_functions.R")
# source("/Users/andreacalabria/Library/CloudStorage/OneDrive-VL58/Utils/GitHub/rna-writer_prototype/R/functions/distalseq_utils.R")

threshold_IS_span <- 30
# load Master file
master_file <- "/Users/andreacalabria/Library/CloudStorage/OneDrive-VL58/Utils/GitHub/rna-writer_prototype/results/2023.MasterFileRun.DISTALseq.xlsx"
master_df <- read.xlsx(xlsxFile = master_file, sheet = "sample_list")
vector_df <- read.xlsx(xlsxFile = master_file, sheet = "vector", startRow = 2)
flankseq_df <- read.xlsx(xlsxFile = master_file, sheet = "IS_Sequence", startRow = 1)
master_df <- merge(x = master_df, y = vector_df, by = c("VectorName_BAM"), all.x = T)
master_df <- merge(x = master_df, y = flankseq_df, by = c("CommonPrefixAllFiles"), all.x = T)
# target_genome_GTF_file <- "/Users/andreacalabria/Library/CloudStorage/OneDrive-VL58/Project DISTALseq/metadata/hg38/hg38.ncbiRefSeq.sorted.slice_exon.gtf"
# target_genome_GTF_file <- "/Users/andreacalabria/Library/CloudStorage/OneDrive-VL58/Project DISTALseq/metadata/hg38/hg38.refGene.sorted.transcript.gtf.gz"

target_genome_GTF_file <- "/Users/andreacalabria/Library/CloudStorage/OneDrive-VL58/Project DISTALseq/metadata/hg38/hg38.ncbiRefSeq.sorted.transcript.gtf.gz" # hg38.ncbiRefSeq.sorted.transcript.gtf.gz
target_genome_TSS_file <- "/Users/andreacalabria/Library/CloudStorage/OneDrive-VL58/Project DISTALseq/metadata/hg38/hg38.ncbiRefSeq.sorted.TSS2.bed.gz" # hg38.ncbiRefSeq.sorted.TSS.bed.gz
target_genome_exons_file <- "/Users/andreacalabria/Library/CloudStorage/OneDrive-VL58/Project DISTALseq/metadata/hg38/hg38.ncbiRefSeq.sorted_slice_exons.gtf.gz" # hg38.ncbiRefSeq.sorted_slice_exons.gtf.gz
target_genome_blacklist_mappability_GTF <- "/Users/andreacalabria/Library/CloudStorage/OneDrive-VL58/Utils/GitHub/rna-writer_prototype/metadata/hg38.ucsc.encode_grc.lowmappability.2024.sort.gtf"

# faceting_order_ddpcr <- c("ddPCR_UTR5", "ddPCR_Transgene5p", "ddPCR_Transgenemid", "ddPCR_Transgene3p", "ddPCR_Promoter", "ddPCR_UTR3", "ddPCR_NoPrimer")
faceting_order_ddpcr <- c("ddPCR_Set1", "ddPCR_Set2", "ddPCR_NoPrimer")
chromosomes_to_use <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")

## for PAF ingest
map_pafcols_to_mycols <- data.frame("mycols" = c("chr", "start", "end", "ratio_bp_aligned_on_raw", "query_start", "query_end", "strand", "name", "read_len"),
                                    "pafcols" = c("tname", "tstart", "tend", "mapq", "qstart", "qend", "strand", "qname", "qlen") )



# cancer genes
msk_db <- read.csv("/Users/andreacalabria/Library/CloudStorage/OneDrive-VL58/Utils/GitHub/isatk/script/R/publicdb/2601.msk.cancerGeneList.tsv", header=TRUE, fill=T, sep='\t', check.names = FALSE)
msk_db$GeneName <- msk_db$`Hugo Symbol` # just add a name that I am used to merge
msk_db <- msk_db[which((msk_db$`Gene Type` %in% c("ONCOGENE", "ONCOGENE_AND_TSG", "TSG")) & 
                         msk_db$`# of occurrence within resources (Column J-P)` > 1),] # avoid unknown genes and false positive / weak genes

