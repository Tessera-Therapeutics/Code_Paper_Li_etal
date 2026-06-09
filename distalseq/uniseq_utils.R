### my functions, for Uniseq


###############################################################
#' @title Collect BED files as output of the pipeline Uniseq, post merge.
#' 
#' @author Andrea Calabria
#' @details version 0.1 (23-12-2)
#'
#' @rdname readUniseqBed
#' @docType methods
#' @aliases readUniseqBed
#'
#' @param df the input dataframe, output of DISTAL-seq pipeline
#' @param df_column_name the column name of the Clone ID to look for in the input df. Default: "KnownInts"
#' @param dict_tokens Dictionary of clones, df from an Excel file showing a map between token and value to replace in the bash command.
#' @param template_string template bash command with explicit TOKEN string. The dictionay will replace the tokens.
#'
#' @return void
#' @usage TODO
#' @note : TODO
###############################################################
readUniseqBed <- function(rootfolder_id, 
                          rootfolder_id_colname = "RootFolder",
                          extra_path = "exact_integration_sites",
                          # molten = TRUE,
                          uniseq_master_df,
                          common_suffix_re = "*.integration-sites.sorted.merged-5.bed",
                          suffix_toremove = ".integration-sites.sorted.merged-5.bed",
                          bed_colnames = c("chr", "integration_locus", "integration_locus_end", "string_score", "SC", "integration_strand", "string_strand", "string_loci"),
                          common_prefix_string = "" # look for this in the master file
                          ){
  ## ---- preliminar operations and checks ------ ##
  # check that THIS obj is not duplicated in the master file
  file_row_index <- NULL
  if (nrow(uniseq_master_df[which(uniseq_master_df[,rootfolder_id_colname] == rootfolder_id & uniseq_master_df$ToProcess == T),]) >1) {
    stop(paste0("[AP]\tError: Duplicated filename (", filename, ") in the metadata file. "))
  } else {
    file_row_index <- rownames(uniseq_master_df[which(uniseq_master_df[,rootfolder_id_colname] == rootfolder_id & uniseq_master_df$ToProcess == T),])
  }
  # look for those files and get data
  files_in_folder <- list.files(path = paste0(uniseq_master_df[file_row_index, rootfolder_id_colname], extra_path), 
                                pattern = common_suffix_re, full.names = T)
  common_prefix_string <- uniseq_master_df[file_row_index, "CommonPrefixAllFiles"]
  
  # for each file, collect BED info and sample ID from the name
  is_data <- NULL
  for (bedfile in files_in_folder) {
    # open, read, acquire data
    is_in_bed <- read.csv(file = bedfile, 
                          header=F, fill=T, sep='\t', check.names = FALSE)
    names(is_in_bed) <- bed_colnames
    rownames(is_in_bed) <- apply(is_in_bed[c("chr", "integration_locus", "integration_strand")], 1, function(x) {
      paste0(x[1], "_", as.character(as.numeric(x[2])), "_", x[3])
      }
    )
    # get metadata (assuming from filename) and add metadata to IS data
    sample_id <- gsub(suffix_toremove, "", basename(bedfile))
    is_in_bed$SampleID <- sample_id
    is_in_bed$SampleName <- uniseq_master_df[file_row_index, "SampleName"]
    is_in_bed$Vector <- uniseq_master_df[file_row_index, "Vector"]
    is_in_bed$PoolID <- uniseq_master_df[file_row_index, "PoolID"]
    # global add
    if (length(is_data)>0){
      is_data <- rbind(is_data, is_in_bed)
    } else {
      is_data <- is_in_bed
    } # if (length(is_data)>0){
  } # for bed files
  
  return (is_data)
}
