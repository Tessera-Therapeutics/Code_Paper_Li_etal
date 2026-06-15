### my functions, starting from DISTALseq


###############################################################
#' @title Retag a BAM file for clones
#' 
#' @author Andrea Calabria
#' @details version 0.1 (23-7-1)
#'
#' @rdname reTagGroupByKnownClone
#' @docType methods
#' @aliases reTagGroupByKnownClone
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
reTagGroupByKnownClone <- function(df, 
                                   df_column_name = "KnownInts", 
                                   dict_tokens = dict_clones, 
                                   template_string = template_string_retag) {
  # other checks to do: folder exists, etc.
  for (row_to_process in rownames(dict_tokens[which(dict_tokens$ToProcess == T),])) {
    
    message(paste("[AP]\tProcessing Row:", row_to_process))
    
    # init 
    this_template_string <- gsub("TOKEN_FILEFOLDER", as.character(dict_tokens[row_to_process, "TOKEN_FILEFOLDER"]), template_string)
    this_template_string <- gsub("TOKEN_INPUTRAWBAM", as.character(dict_tokens[row_to_process, "TOKEN_INPUTRAWBAM"]), this_template_string)
    this_template_string <- gsub("TOKEN_OUTPUT_SINGLE_IS_BAM", as.character(dict_tokens[row_to_process, "TOKEN_OUTPUT_SINGLE_IS_BAM"]), this_template_string)
    this_template_string <- gsub("TOKEN_READLIST", as.character(dict_tokens[row_to_process, "TOKEN_READLIST"]), this_template_string)
    this_template_string <- gsub("TOKEN_FINAL_SINGLE_IS", as.character(dict_tokens[row_to_process, "TOKEN_FINAL_SINGLE_IS"]), this_template_string)
    this_template_string <- gsub("TOKEN_KNOWNIS_ID", as.character(dict_tokens[row_to_process, "TOKEN_KNOWNIS_ID"]), this_template_string)
    this_template_string <- gsub("TOKEN_CLONEID", as.character(dict_tokens[row_to_process, "TOKEN_CLONEID"]), this_template_string)
    this_template_string <- gsub("TOKEN_VECTORID", as.character(dict_tokens[row_to_process, "TOKEN_VECTORID"]), this_template_string)
    this_template_string <- gsub("TOKEN_OUTPUTFOLDER", as.character(dict_tokens[row_to_process, "TOKEN_OUTPUTFOLDER"]), this_template_string)
    this_template_string <- gsub("TOKEN_SAMPLEID", as.character(dict_tokens[row_to_process, "TOKEN_SAMPLEID"]), this_template_string)
    
    # slice df
    isid <- as.character(dict_tokens[row_to_process, "TOKEN_KNOWNIS_ID"])
    df_slice <- df[which(df[,df_column_name] == isid & df$SampleID == as.character(dict_tokens[row_to_process, "SampleID"]) ), ]
    # write only those reads
    readlist <- paste0(as.character(dict_tokens[row_to_process, "TOKEN_FILEFOLDER"]), as.character(dict_tokens[row_to_process, "TOKEN_READLIST"]))
    write.table(x = rownames(df_slice), file = readlist, sep = "\t", quote = FALSE, row.names = FALSE, col.names = F, na = '')
    
    # run it 
    system(this_template_string, intern = TRUE,
           ignore.stdout = FALSE, ignore.stderr = FALSE,
           show.output.on.console = TRUE)
  }
  # then run picard MergeSamFiles \
  # I=Vingi.CloneC1C3_chr16.rgIS.bam \
  # I=Vingi.CloneC1C3_chr3.rgIS.bam \
  # I=Vingi.CloneC1C3_chr4.rgIS.bam \
  # I=Vingi.CloneC1C3_chr5.rgIS.bam \
  # I=Vingi.CloneC1C3_chr6.rgIS.bam \
  # I=Vingi.CloneC1C3_chr9.rgIS.bam \
  # I=Vingi.CloneC2F2_chr1.rgIS.bam \
  # I=Vingi.CloneC2F2_chr18.rgIS.bam \
  # I=Vingi.CloneC2F2_chr19.rgIS.bam \
  # I=Vingi.CloneC2F2_chr2.rgIS.bam \
  # O=AllClones.AllIS.rg.bam
  
}


###############################################################
#' @title Retag a BAM file for buckets
#' 
#' @author Andrea Calabria
#' @details version 0.1 (24-6-10)
#'
#' @rdname reTagGroupByBucket
#' @docType methods
#' @aliases reTagGroupByBucket
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
reTagGroupByKnownClone <- function(df, 
                                   df_column_name = "KnownInts", 
                                   dict_tokens = dict_clones, 
                                   template_string = template_string_retag) {
  # other checks to do: folder exists, etc.
  for (row_to_process in rownames(dict_tokens[which(dict_tokens$ToProcess == T),])) {
    
    message(paste("[AP]\tProcessing Row:", row_to_process))
    
    # init 
    this_template_string <- gsub("TOKEN_FILEFOLDER", as.character(dict_tokens[row_to_process, "TOKEN_FILEFOLDER"]), template_string)
    this_template_string <- gsub("TOKEN_INPUTRAWBAM", as.character(dict_tokens[row_to_process, "TOKEN_INPUTRAWBAM"]), this_template_string)
    this_template_string <- gsub("TOKEN_OUTPUT_SINGLE_IS_BAM", as.character(dict_tokens[row_to_process, "TOKEN_OUTPUT_SINGLE_IS_BAM"]), this_template_string)
    this_template_string <- gsub("TOKEN_READLIST", as.character(dict_tokens[row_to_process, "TOKEN_READLIST"]), this_template_string)
    this_template_string <- gsub("TOKEN_FINAL_SINGLE_IS", as.character(dict_tokens[row_to_process, "TOKEN_FINAL_SINGLE_IS"]), this_template_string)
    this_template_string <- gsub("TOKEN_KNOWNIS_ID", as.character(dict_tokens[row_to_process, "TOKEN_KNOWNIS_ID"]), this_template_string)
    this_template_string <- gsub("TOKEN_CLONEID", as.character(dict_tokens[row_to_process, "TOKEN_CLONEID"]), this_template_string)
    this_template_string <- gsub("TOKEN_VECTORID", as.character(dict_tokens[row_to_process, "TOKEN_VECTORID"]), this_template_string)
    this_template_string <- gsub("TOKEN_OUTPUTFOLDER", as.character(dict_tokens[row_to_process, "TOKEN_OUTPUTFOLDER"]), this_template_string)
    this_template_string <- gsub("TOKEN_SAMPLEID", as.character(dict_tokens[row_to_process, "TOKEN_SAMPLEID"]), this_template_string)
    
    # slice df
    isid <- as.character(dict_tokens[row_to_process, "TOKEN_KNOWNIS_ID"])
    df_slice <- df[which(df[,df_column_name] == isid & df$SampleID == as.character(dict_tokens[row_to_process, "SampleID"]) ), ]
    # write only those reads
    readlist <- paste0(as.character(dict_tokens[row_to_process, "TOKEN_FILEFOLDER"]), as.character(dict_tokens[row_to_process, "TOKEN_READLIST"]))
    write.table(x = rownames(df_slice), file = readlist, sep = "\t", quote = FALSE, row.names = FALSE, col.names = F, na = '')
    
    # run it 
    system(this_template_string, intern = TRUE,
           ignore.stdout = FALSE, ignore.stderr = FALSE,
           show.output.on.console = TRUE)
  }
  # then run picard MergeSamFiles \
  # I=Vingi.CloneC1C3_chr16.rgIS.bam \
  # I=Vingi.CloneC1C3_chr3.rgIS.bam \
  # I=Vingi.CloneC1C3_chr4.rgIS.bam \
  # I=Vingi.CloneC1C3_chr5.rgIS.bam \
  # I=Vingi.CloneC1C3_chr6.rgIS.bam \
  # I=Vingi.CloneC1C3_chr9.rgIS.bam \
  # I=Vingi.CloneC2F2_chr1.rgIS.bam \
  # I=Vingi.CloneC2F2_chr18.rgIS.bam \
  # I=Vingi.CloneC2F2_chr19.rgIS.bam \
  # I=Vingi.CloneC2F2_chr2.rgIS.bam \
  # O=AllClones.AllIS.rg.bam
  
}

###############################################################
#' @title Find genomic intervals in a given list of locations
#' 
#' @author Andrea Calabria
#' @details version 0.1 (23-7-10)
#'
#' @rdname findGenomicIntervals
#' @docType methods
#' @aliases findGenomicIntervals
#'
#' @param df the input dataframe, output of DISTAL-seq pipeline
#'
#' @return df with interval IDs
#' @usage TODO
#' @note : TODO is_id_col_sorted <- findGenomicIntervals(df = ispileup, show_status_bar = T)
#' @todo Strand specific evaluation
###############################################################
findGenomicIntervals <- function(df, 
                                 threshold_IS_span = 10, 
                                 show_status_bar = T,
                                 in_chr_colname = "targetRegion_chr",
                                 in_locus_colname = "integration_locus",
                                 in_strand_colname = "integration_strand",
                                 out_colnames = c("RefID", "Ref_chr", "RefID_start", "RefID_end", "RefID_span", "RefID_strand"),
                                 quiet = T) {
  message(paste0("[AP]\tFind genomic intervals"))
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = nrow(df), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 100,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  # df[,in_locus_colname] <- as.numeric(df[,in_locus_colname]) # grant a numeric field
  is_id <- 1
  first_element <- TRUE
  id_start <- NULL
  id_end <- NULL
  id_span <- NULL
  id_df <- NULL
  for (i in seq(1, nrow(df))) {
    if (!quiet) {message(paste0("[AP]\t ================\n\t\ti=",i))}
    if (!first_element) {
      if (!quiet) {message(paste0("[AP]\tnot first element - ",i))}
      if ( ((df[i, in_chr_colname] == df[i-1, in_chr_colname]) & 
            ( abs( df[i, in_locus_colname] - df[i-1, in_locus_colname] ) > threshold_IS_span ) ) |
           (df[i, in_chr_colname] != df[i-1, in_chr_colname])
      ) { # init a NEW ref interval
        if (!quiet) {message(paste0("[AP]\tidentical chr and loci delta above thr:: df[i,in_chr_colname]:",df[i, in_chr_colname],"---","df[i-1, in_chr_colname]:",df[i-1, in_chr_colname]))}
        is_id <- is_id + 1
        ## ---- update df ----
        id_start <- as.numeric(df[i, in_locus_colname])
        id_end <- as.numeric(df[i, in_locus_colname])
        id_span <- id_end-id_start
        id_df <- rbind( id_df, data.frame("RefID" = paste0(df[i, in_chr_colname], "_", is_id), df[i, in_chr_colname], id_start, id_end, id_span, df[i, in_strand_colname]))
        if (!quiet) {message(paste0("[AP]\tis_id = ",is_id))}
      } else { # this is belongs to the same (previous) interval
        if (!quiet) {message(paste0("[AP]\telse branch = ",is_id, " -> df[i,in_chr_colname]:",df[i, in_chr_colname],"---","df[i-1, in_chr_colname]:",df[i-1, in_chr_colname]))}
        ## ---- update df ----
        id_end <- as.numeric(df[i, in_locus_colname])
        id_span <- id_end-id_start
        id_df <- rbind( id_df, data.frame("RefID" = paste0(df[i, in_chr_colname], "_", is_id), df[i, in_chr_colname], id_start, id_end, id_span, df[i, in_strand_colname]))
      } # if ( ((df[i, in_chr_colname] == df[i-1, in_chr_colname])
    } else {
      if (!quiet) {message(paste0("[AP]\tfirst elem, is_id = ",is_id))}
      first_element <- F
      ## ---- init df ----
      id_start <- as.numeric(df[i, in_locus_colname])
      id_end <- id_start
      id_span <- 0
      id_df <- data.frame("RefID" = paste0(df[i, in_chr_colname], "_", is_id), df[i, in_chr_colname], id_start, id_end, id_span, df[i, in_strand_colname])
      # id_df <- data.frame("RefID" = is_id,
      #                     "RefID_start" = id_start,
      #                     "RefID_end" = id_end,
      #                     "RefID_span" = id_span)
    } # if (!first_element)
    # Sets the progress bar to the current state
    if (show_status_bar) {setTxtProgressBar(pb, i)}
  } # for (i in seq(1, nrow(df))
  names(id_df) <- out_colnames
  return (id_df)
}

###############################################################
#' @title Split df by chr
#' 
#' @author Andrea Calabria
#' @details version 0.1 (23-9-13)
#'
#' @rdname splitDfByChr
#' @docType methods
#' @aliases splitDfByChr
#'
#' @param df the input dataframe of IS
#'
#' @return df with interval IDs
#' @usage TODO
#' @note : 
#' @todo 
###############################################################
splitDfByChr <- function(df, chr_colname = "chr") {
  message(paste0("[AP]\tSplit in chrs"))
  if (chr_colname %in% colnames(df)) {
    outlist <- split(df, f = df[,chr_colname])
  } else {
    message(paste0("[AP]\tError: input chr column not found: ", chr_colname))
  }
  return (outlist)
}

###############################################################
#' @title Classify single IS as complete or not
#' 
#' @author Andrea Calabria
#' @details version 0.1 (23-11-23 - Thanksgiving)
#'
#' @rdname classifyCompleteIntegrationByIS
#' @docType methods
#' @aliases classifyCompleteIntegrationByIS
#'
#' @param df (tidy) of reads (with reference IS) having the following columns: "RefID", "SampleID", "Complete"
#'
#' @return input df with corrected classification (unique per IS)
#' @usage TODO
#' @logic Groupd reads by the Complete tag. If complete is unique (True or False) by IS, then report it, else: if True is > False -> Complete = T.
#' @note : 
#' @todo 
###############################################################
classifyCompleteIntegrationByIS <- function(df, 
                                            group_by_cols = c("RefID", "SampleID", "Complete"),
                                            locus_refID_colname = "RefID",
                                            complete_colname = "Complete", 
                                            complete_positive_val = "True", 
                                            complete_negative_val = "False",
                                            output_plot_file_prefix = "QC.Complete_vs_Truncated.MixedISCases",
                                            locus_sample_col = "SampleID",
                                            height=7, width=9, plot_qc = F
                                            ) {
  message(paste0("[AP]\tClassify IS as Complete or not through reads (correcting PCR errors)."))
  return_df <- NULL
  # check colnames in df
  if (FALSE %in% (group_by_cols %in% colnames(df))) {
    stop(paste0("[AP]\tError: input column/s not found. Your input is: ", paste0(group_by_cols, collapse = " + ")))
  } else {
    # check vals
    if (!(complete_positive_val %in% levels(factor(df[,complete_colname])) | complete_negative_val %in% levels(factor(df[,complete_colname])))) {
      stop(paste0("[AP]\tError: Wrong values in DF. Check complete_positive_val or complete_negative_val"))
    }
    # create df of refid (IS) with N of reads complete or not by sample
    df_group <- df %>% 
      group_by_at(group_by_cols) %>% 
      dplyr::summarise(countCompleteN = n())
    # cast data on "complete" to compare the counts
    cast_cols_string <- setdiff(group_by_cols, complete_colname)
    df_group_cast <- dcast(data = df_group, 
                       paste(paste0(cast_cols_string, collapse = " + "), "~", complete_colname),
                       value.var = "countCompleteN", fun.aggregate = sum)
    # complete_colname to be padded!
    if ( FALSE %in% (c(complete_positive_val, complete_negative_val) %in% colnames(df_group_cast)) ){
      if (complete_negative_val %in% colnames(df_group_cast)) {
        df_group_cast[,complete_positive_val] <- 0
      } else {
        df_group_cast[,complete_negative_val] <- 0
      }
    } 
    # evaluate complete
    df_group_cast$CompleteIntegration <- apply(df_group_cast[c(complete_positive_val, complete_negative_val)], 1, function(x) {
      ifelse(x[1]>=x[2], "Full", "Partial")
      }
    )
    # merge data
    return_df <- merge(x = df, 
                       y = df_group_cast, 
                       by = c(cast_cols_string),
                       all.x = T)
    
    if (plot_qc) {
      # before lcosing everythig, do QC plots
      set_c <- unique(df[which(df[,complete_colname] == complete_positive_val), locus_refID_colname])
      set_t <- unique(df[which(df[,complete_colname] == complete_negative_val), locus_refID_colname])
      set_m <- intersect(set_c, set_t)
      df_group_m <- df[which(df[,locus_refID_colname] %in% set_m),]
      # assumption: min is THE real complete IS, max is the not complete (if group by without the field complete).
      df_group_m_summary <- df_group_m %>% 
        group_by_at(c(locus_refID_colname, locus_sample_col)) %>% 
        dplyr::summarise(MaxDelta = (max(as.numeric(vec_seg1_start)) - min(as.numeric(vec_seg1_start))), 
                         nReads = n()) %>%
        filter(MaxDelta>0)
      
      plot_histplot_notcomplete <- 
        ggplot(df_group_m_summary, aes(x = MaxDelta ) ) +
        geom_density(
          aes(
            # colour = factor(cyl),
            # fill = after_scale(alpha(colour, 0.3)),
            y = after_stat(count / sum(count))
          ) ) +
        # scale_y_continuous(labels=scales::percent) +
        scale_y_sqrt(labels=scales::percent) +
        # scale_x_log10() +
        theme_bw() +
        theme(axis.text.x = element_text(size=14, angle = 0), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=22)) +
        labs(title = paste0("Mixed IS (complete and not), alignment distribution\nDistance from the 5'UTR end."), 
             x = "Distance (bp) from the end inward the vector", y = "Density", color = "SampleID", fill = "SampleID",
             subtitle = paste0("N. IS = ", length(unique(df_group_m[,locus_refID_colname])), "; N. reads = ", nrow(df_group_m), "; Percentage of IS on total input IS = ", round((length(unique(df_group_m[,locus_refID_colname]))/length(unique( pull(df_group, locus_refID_colname) )))*100, digits=2) ,"%.") ) 
      
      pdf(file = paste(output_plot_file_prefix, ".density.pdf", sep = ""), height=height, width=width)
      plot(plot_histplot_notcomplete)
      dev.off()
      png(file = paste(output_plot_file_prefix, ".density.png", sep = ""), height=height, width=width, units = "in", res = 300)
      plot(plot_histplot_notcomplete)
      dev.off()
      
      # now with the ridge plot
      plot_ridges_facet <- 
        ggplot(df_group_m_summary, aes(x = MaxDelta, y = SampleID) ) +
        # ggplot(full_df, aes(x = Read_Len, y = facet_on_col, fill = facet_on_col)) +
        geom_density_ridges() +
        theme_ridges() + 
        theme(legend.position = "none") +
        theme(strip.text = element_text(face="plain", size=14)) +
        theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0),
              strip.text.y = element_text(size = 14, colour = "darkblue", angle = 270)) +
        theme(axis.text.x = element_text(size=14, angle = 30, hjust = 1, vjust = 1),
              axis.text.y = element_text(size=16),
              axis.title = element_text(size=16),
              plot.title = element_text(size=22)) +
        labs(title = paste0("Mixed IS (complete and not), alignment distribution\nDistance from the 5'UTR end."),
             x = "Base-pairs from the end of the vector", y = "Density", color = "Sample", fill = "Sample")
      
      pdf(file = paste(output_plot_file_prefix, ".density.ridges.pdf", sep = ""), height=height, width=width)
      plot(plot_ridges_facet)
      dev.off()
      png(file = paste(output_plot_file_prefix, ".density.ridges.png", sep = ""), height=height, width=width, units = "in", res = 300)
      plot(plot_ridges_facet)
      dev.off()
    } # if plot qc
    
  }
  return (return_df)
}

# add genomics info of the vector hits, col "Vector_Hits"
overlaps <- function(start1, end1, start2, end2) {
  return(!(end1 < start2 || end2 < start1))
}

# function to extract info from a single summarized string of alignment, format: <chr>:<start>-<end><strand>
getAlmInfo <- function(this_is) {
  # this_is <- paste0("chr", strsplit(as.character(x), ';', fixed = T)[[1]][1])
  strand <- substr(this_is, start = nchar(this_is), stop = nchar(this_is))
  sp1 <- strsplit(this_is, ':', fixed = T)
  chr <- sp1[[1]][1]
  coord <- strsplit(sp1[[1]][2], '-', fixed = T)
  start <- as.numeric(coord[[1]][1])
  end <- as.numeric(strsplit(coord[[1]][2], split = "[+-]")[[1]][1])
  return(c(chr, start, end, strand))
}

###############################################################
#' @title Acquire the TSV output file of the distalseq pipeline as df
#' 
#' @author Andrea Calabria
#' @details version 0.1 (23-9-13)
#'
#' @rdname parseDistalseqPipeTSVOutput
#' @docType methods
#' @aliases parseDistalseqPipeTSVOutput
#'
#' @param filename to ingest of the TSV output (usually with the suffix "combined.ToHGandPLV.ReadInfo.additions.tsv.gz")
#'
#' @return df with interval IDs
#' @usage t <- parseDistalseqPipeTSVOutput(filename = "2023_07_21_RTE1_gDNA77261_B1.combined.ToHGandPLV.ReadInfo.additions.tsv.gz")
#' @note : 
#' @todo : check all required cols -> or generalize as input par
###############################################################
parseDistalseqPipeTSVOutput <- function(filename,
                                        master_df,
                                        required_metadata_cols = c("SampleID", "SampleName", "Replica", "Vector", "PoolID", "GroupName", "Transgene", "SampleType", "Separator", "Promoter_start", "Promoter_end", "pA_start", "pA_end", "Functional_min", "Functional_max", "VectorOrientation", "MinReadLen"),
                                        import_only_proper = TRUE,
                                        rownames_as_id = FALSE,
                                        nastring = c("NONE", "NA", "NULL", "NaN", "ND", ""),
                                        threshold_IS_span = 30,
                                        analyze_ddPCR_primer_set = FALSE,
                                        chromosomes_to_use = c(seq(1, 23), "X", "Y") ) {
  ## ---- preliminar operations and checks ------ ##
  # # read metadata
  # if (file.exists(master_metadata_file)) {
  #   master_df <- read.xlsx(xlsxFile = master_metadata_file)  
  # } else {
  #   stop(paste0("[AP]\tError: metadata file does not exist. ", master_metadata_file))
  # }
  # check all required cols in the metadata file
  if (FALSE %in% (required_metadata_cols %in% colnames(master_df))) {
    stop(paste0("[AP]\tError: Missing colnames in metadata df"))
  }
  # check that THIS FILENAME is not duplicated in the master file
  file_row_index <- NULL
  if (nrow(master_df[which(master_df$FileName == filename & master_df$ToProcess == T),]) >1) {
    stop(paste0("[AP]\tError: Duplicated filename (", filename, ") in the metadata file. "))
  } else {
    file_row_index <- rownames(master_df[which(master_df$FileName == filename & master_df$ToProcess == T),])
  }
  ## ---- read input files ------ ##
  file_name <- paste0(as.character(master_df[file_row_index, "FileFolder"]), as.character(master_df[file_row_index, "FileName"]))
  message(paste("[AP]\tProcessing File:", file_name))
  # message(paste("[AP]\t\tGet data"))
  if (as.character(master_df[file_row_index,"Separator"]) == ",") {
    mat_df <- read.csv(file = file_name, header=TRUE, fill=T, sep=',', check.names = FALSE, na.strings = nastring)
  } else {
    mat_df <- read.csv(file = file_name, header=TRUE, fill=T, sep='\t', check.names = FALSE, na.strings = nastring)
  }
  # if (FALSE %in% (colnames(mat_df) %in% pivolta_colnames)) {
  #   stop(paste0("[AP]\tError: some pivotal col names are missing in the input file. Required cols: ", paste(pivolta_colnames, collapse = ' - ')))
  # }
  # TMP: to prevent some issues with duplicated rows, try to import only proper reads
  if (import_only_proper) {
    mat_df <- mat_df[which(mat_df$Proper == "True"),]
  }
  if (nrow(mat_df) == 0) {
    message(paste("[AP]\t\tSKIPPING THIS FILE SINCE EMPTY:", file_name))
    return (NULL) 
  } else {
    # apply rownames if true
    if (rownames_as_id) {
      # remove doublets, produced by error (splitting reads in the pipe nextflow)
      allID <- as.data.frame(table(mat_df$Read_Name))
      rownames(allID) <- allID$Var1
      uniqueID <- rownames(allID[which(allID$Freq == 1),])
      message(paste("[AP]\t\tUnique ID: Loosing", (nrow(allID) - length(uniqueID)), "but keeping", length(uniqueID), 
                    "[loss =", round((nrow(allID) - length(uniqueID))/(nrow(allID))*100, digits = 3), "%]"))
      # now use ONLY the unique IDs
      mat_df <- mat_df[which(mat_df$Read_Name %in% uniqueID),]
      rownames(mat_df) <- paste(mat_df$Read_Name, as.character(master_df[file_row_index,"SampleID"]), sep = "_")
    }
    
    # adjust colnames to latest names (Jason, June 5)
    names(mat_df) <- gsub("vec_FirstSegment_", "vec_seg1_", colnames(mat_df))
    names(mat_df) <- gsub("vec_SecondSegment_", "vec_seg2_", colnames(mat_df))
    names(mat_df) <- gsub("vec_ThirdSegment_", "vec_seg3_", colnames(mat_df))
    
    # add annotation col by sample
    mat_df$SampleName <- as.character(master_df[file_row_index,"SampleName"])
    mat_df$SampleID <- as.character(master_df[file_row_index,"SampleID"])
    mat_df$SampleType <- as.character(master_df[file_row_index,"SampleType"])
    mat_df$PoolID <- as.character(master_df[file_row_index,"PoolID"])
    mat_df$Replica <- as.character(master_df[file_row_index,"Replica"])
    mat_df$Vector <- as.character(master_df[file_row_index,"Vector"])
    mat_df$Transgene <- as.character(master_df[file_row_index,"Transgene"])
    mat_df$gDNA <- as.character(master_df[file_row_index,"DNA_ID"])
    mat_df$GroupName <- as.character(master_df[file_row_index,"GroupName"])
    
    # min read len?
    mat_df$AboveMinReadLen <- ifelse( mat_df$Read_Len >= master_df[file_row_index,"MinReadLen"],
                                      TRUE, 
                                      FALSE
    )
    mat_df$TargetGenomeAlmSize <- ifelse(!is.na(mat_df$Human_AlnStart), abs(mat_df$Human_AlnEnd - mat_df$Human_AlnStart), NA)
    mat_df$TargetGenomeAlmSize_AboveMinLen <- ifelse( is.na(mat_df$Human_AlnStart), 
                                                      FALSE,
                                                      ifelse( mat_df$TargetGenomeAlmSize >= master_df[file_row_index,"MinTargetAlmLen"], TRUE, FALSE ) )
    
    # redefine "proper" label
    mat_df$Chimeric <- ifelse(mat_df$Num_HumanAln >= 1 & mat_df$Num_VectorAln > 0, TRUE, FALSE)
    
    # some ad hoc cols
    # functional is present from pA to promoter
    # OLD:: mat_df$Functionality <- ifelse( (mat_df$vec_seg1_start <= master_df[file_row_index,"pA_start"]) & (mat_df$vec_seg1_end >= master_df[file_row_index,"Promoter_end"]-1),
    #                                 "Functional", "Defective"
    # ) # ifelse
    mat_df$Functionality <- ifelse( 
      (mat_df$vec_seg1_start <= master_df[file_row_index,"Functional_min"]) & 
        (mat_df$vec_seg1_end >= master_df[file_row_index,"Functional_max"]),
      "Functional", 
      "Defective"
    ) # ifelse
    
    # do we have other features?
    mat_df$Promoter <- apply(mat_df[c("vec_seg1_start", "vec_seg1_end")], 1, function(x) {
      ifelse(!is.na(master_df[file_row_index,"Promoter_start"]) ,
             ifelse( ((x[1] <= master_df[file_row_index,"Promoter_start"]) & 
                        (x[2] >= master_df[file_row_index,"Promoter_end"])), 
                     "Full", 
                     "Truncated"),
             NA) # ifelse
    })
    mat_df$RearrPromoter <- apply(mat_df[c("Num_VectorAln", "vec_seg2_start", "vec_seg2_end", "vec_seg3_start", "vec_seg3_end")], 1, function(x) {
      ifelse((!is.na(master_df[file_row_index,"Promoter_start"]) & (x[1] > 1)),
             ifelse( ((x[2] <= master_df[file_row_index,"Promoter_start"]) & 
                        (x[3] >= master_df[file_row_index,"Promoter_end"])), 
                     "Full", 
                     "Truncated"),
             NA) # ifelse
    })
    mat_df$IntegrationType <- ifelse( mat_df$Num_VectorAln == 1, "Single", "Rearranged") # old:: "Structural Change" | Linear vs Rearranged
    # mat_df$InvertedFragment <- ifelse( mat_df$Inversion == "True", "Inverted", "Not inverted")
    mat_df$InvertedFragment <- apply(mat_df[c("Num_VectorAln", "Vector_AlnStrands")], 1, function(x){
      strand_vect <- unique(as.character(strsplit(as.character(x[2]), ',', fixed = T)[[1]]))
      ifelse( length(strand_vect) > 1, "Inverted", "Not inverted" )
      } 
    )
    
    # add genomics info of the IS, from col "Human_Hits"
    mat_df_targetIS_info <- as.data.frame(
      t(
        apply(mat_df[c("Human_Hits", "Read_Name")], 1, function(x) {
          this_is <- paste0("chr", strsplit(as.character(x[1]), ';', fixed = T)[[1]][1])
          getAlmInfo(this_is)
        }
        )
      )
    )
    names(mat_df_targetIS_info) <- paste("targetRegion", c("chr", "start", "end", "strand"), sep = "_")
    mat_df <- cbind(mat_df, mat_df_targetIS_info)
    mat_df$human_IS_min <- mat_df$human_IS - threshold_IS_span
    mat_df$human_IS_max <- mat_df$human_IS + threshold_IS_span
    
    # this test data is for the next apply on vector hits
    test_f <- data.frame("Vector_Hits" = c("PLV5194:5581-6639-;PLV5194:4741-6000-", "PLV5194:5568-5639-;PLV5194:4449-4550-", "PLV5194:1482-5638-;PLV5194:4626-4877-", "PLV5194:5581-6639-"), 
                         "Read_Name" = c(1, 2, 3, 4))
    mat_df_vechits_info <- 
      data.frame( "OverlappingVectorHits" =
      # t(
        apply(mat_df[c("Vector_Hits", "Read_Name")], 1, function(x) {
        # apply(test_f[c("Vector_Hits", "Read_Name")], 1, function(x) {
          vhits <- strsplit(as.character(x[1]), ';', fixed = T)[[1]]
          df_hits <- NULL
          for (vh in vhits) {
            this_is <- vh
            df_hits <- rbind( df_hits, c(x[2], getAlmInfo(this_is) ) )
          } # for
          overlap <- FALSE
          if (nrow(df_hits)>1) {
            overlap <- overlaps(df_hits[1,3], df_hits[1,4], df_hits[2,3], df_hits[2,4])
          } 
          # overlap
          ifelse(overlap, "Overlapping", "Gapped")
          } # apply function
        ) # apply
      # ) # t
    ) # df
    mat_df <- cbind(mat_df, mat_df_vechits_info)
    
    # our IS locus is needed now
    mat_df$integration_locus <- ifelse( ((mat_df$vec_seg1_strand == "+" & mat_df$targetRegion_strand == "+") | ((mat_df$vec_seg1_strand == "-" & mat_df$targetRegion_strand == "-"))),
                                        mat_df$targetRegion_end,
                                        ifelse( ((mat_df$vec_seg1_strand == "+" & mat_df$targetRegion_strand == "-") | (mat_df$vec_seg1_strand == "-" & mat_df$targetRegion_strand == "+") ),
                                                mat_df$targetRegion_start,
                                                NA)
    )
    # IS strand, logics: if vector and human strand are concordant -> negative, else positive
    mat_df$integration_strand <- ifelse( 
      ((mat_df$vec_seg1_strand == "+" & mat_df$targetRegion_strand == "+") | ((mat_df$vec_seg1_strand == "-" & mat_df$targetRegion_strand == "-"))),
      "-",
      ifelse( ((mat_df$vec_seg1_strand == "+" & mat_df$targetRegion_strand == "-") | (mat_df$vec_seg1_strand == "-" & mat_df$targetRegion_strand == "+") ),
              "+",
              NA)
    )
    
    # deal with Human hits, in particular with the cases of multiple hits (not captured as proper chimera in the pipe 1.9)
    # hh <- "12:103946744-103947202-;2:169347472-169349176-"
    mat_df_humanhits_info <- 
      data.frame( t(
                    apply(mat_df[c("TargetAln_Order", "Vector_Hits", "Vector_AlnStrands", "Human_Hits", "Read_Name", "Num_VectorAln", "Num_HumanAln")], 1, function(x) {
                      if (x[6] > 0 & x[7] > 1) {
                        # get df of human hits
                        thits <- strsplit(as.character(x[4]), ';', fixed = T)[[1]]
                        df_hits <- NULL
                        for (th in thits) {
                          df_hits <- rbind( df_hits, c(x[5], getAlmInfo(th)) )
                        } # for
                        df_hits <- as.data.frame(df_hits)
                        names(df_hits) <- c("Read_Name", "chr", "start", "end", "strand")
                        # get order of almn
                        alm_order <- strsplit(as.character(x[1]), '-', fixed = T)[[1]]
                        vector_position <- which((alm_order %in% as.character(gsub("chr", "", chromosomes_to_use))) == FALSE)[1] # index of the occurrence of the vector
                        if (vector_position > 1) { # case > 1 => take IS at the LEFT of the vector -> strand + !. this means take the flanking human hit (!). i.e. "10-9-PL11181" take chr 9
                          chr_flankingIS <- alm_order[(vector_position-1)] # get the left flanking hit
                          if (nrow(df_hits[which(df_hits[2] == chr_flankingIS),]) > 1) { # if the same chr has been targeted >1, capture the first/last one
                            ref_thit <- df_hits[which(df_hits[2] == chr_flankingIS),][2,]
                          } else {
                            ref_thit <- df_hits[which(df_hits[2] == chr_flankingIS),]
                          } # if > 1 hit on the same chr
                        } else { # case == 1 -> - orientation, get RIGHT alm
                          chr_flankingIS <- alm_order[(vector_position+1)] # get the right flanking hit
                          ref_thit <- df_hits[1,]
                        }
                        # identify the locus using the vector strand and the hg strand alm
                        putative_locus <- NULL
                        if (x[3] == '+') { # + vector strand -> left flanking IS -> IS alm end has the IS
                          if (ref_thit[5] == "+") {
                            # the end is the IS
                            putative_locus <- ref_thit$end
                          } else {
                            # the start is the IS
                            putative_locus <- ref_thit$start
                          }
                        } else { # vector is -, IS in + is the start, in - the end
                          if (ref_thit[5] == "+") {
                            # the start is the IS
                            putative_locus <- ref_thit$start
                          } else {
                            # the end is the IS
                            putative_locus <- ref_thit$end
                          }
                        }
                        ref_thit$putative_locus <- putative_locus
                        # now the vector... get its features...
                        tvh <- getAlmInfo( strsplit(as.character(x[2]), ';', fixed = T)[[1]][1] )
                        ref_thit$flankingTargetRegion_functionality <- ifelse( 
                          (tvh[2] <= master_df[file_row_index,"Functional_min"]) & 
                            (tvh[3] >= master_df[file_row_index,"Functional_max"]),
                          "Functional", 
                          "Defective"
                        ) # ifelse
                        ref_thit$flankingTargetRegion_complete <- ifelse( 
                          (tvh[2] <= master_df[file_row_index,"UTR5_start"]),
                          "True", 
                          "False"
                        ) # ifelse
                        ref_thit$fromalm_vec_seg1_start <- tvh[2]
                        ref_thit$fromalm_vec_seg1_end <- tvh[3]
                        ref_thit$fromalm_vec_seg1_strand <- tvh[4]
                        # return all info
                        unlist(ref_thit)
                      } else {
                        c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
                      } # if x6 x7
                    } # apply function

                  ) # apply
                ) # t
      ) # df
    names(mat_df_humanhits_info) <- c("rname", "flankingTargetRegion_chr", "flankingTargetRegion_start", "flankingTargetRegion_end", "flankingTargetRegion_strand", "flankingTargetRegion_locus", "flankingTargetRegion_functionality", "flankingTargetRegion_complete", "fromalm_vec_seg1_start", "fromalm_vec_seg1_end", "fromalm_vec_seg1_strand")
    mat_df <- cbind(mat_df, mat_df_humanhits_info)
    
    # only for the case of 2 human hits and 1 vector hit, correct the main fields
    mat_df$Complete <- ifelse((is.na(mat_df$Complete) & !is.na(mat_df$flankingTargetRegion_complete)), 
                              mat_df$flankingTargetRegion_complete,
                              mat_df$Complete)
    mat_df$Functionality <- ifelse((is.na(mat_df$Functionality) & !is.na(mat_df$flankingTargetRegion_functionality)), 
                              mat_df$flankingTargetRegion_functionality,
                              mat_df$Functionality)
    mat_df$integration_locus <- ifelse((is.na(mat_df$integration_locus) & !is.na(mat_df$flankingTargetRegion_locus)), 
                              mat_df$flankingTargetRegion_locus,
                              mat_df$integration_locus)
    mat_df$human_IS <- ifelse((is.na(mat_df$human_IS) & !is.na(mat_df$flankingTargetRegion_locus)), 
                              mat_df$flankingTargetRegion_locus,
                              mat_df$human_IS)
    mat_df$integration_strand <- ifelse((is.na(mat_df$integration_strand) & !is.na(mat_df$flankingTargetRegion_strand)),
                                        ifelse(
                                          ((mat_df$Vector_AlnStrands == "+" & mat_df$flankingTargetRegion_strand == "+") | ((mat_df$Vector_AlnStrands == "-" & mat_df$flankingTargetRegion_strand == "-"))),
                                          "-",
                                          ifelse( ((mat_df$Vector_AlnStrands == "+" & mat_df$flankingTargetRegion_strand == "-") | (mat_df$Vector_AlnStrands == "-" & mat_df$flankingTargetRegion_strand == "+") ),
                                                  "+",
                                                  NA)
                                        ),
                                        # mat_df$flankingTargetRegion_strand,
                                        mat_df$integration_strand)
    mat_df$Human_AlnStart <- ifelse((is.na(mat_df$Human_AlnStart) & !is.na(mat_df$flankingTargetRegion_start)), 
                              mat_df$flankingTargetRegion_start,
                              mat_df$Human_AlnStart)
    mat_df$Human_AlnEnd <- ifelse((is.na(mat_df$Human_AlnEnd) & !is.na(mat_df$flankingTargetRegion_end)), 
                              mat_df$flankingTargetRegion_end,
                              mat_df$Human_AlnEnd)
    mat_df$Human_Strand <- ifelse((is.na(mat_df$Human_Strand) & !is.na(mat_df$flankingTargetRegion_strand)), 
                              mat_df$flankingTargetRegion_strand,
                              mat_df$Human_Strand)
    mat_df$Human_Chr <- ifelse((is.na(mat_df$Human_Chr) & !is.na(mat_df$flankingTargetRegion_chr)), 
                              mat_df$flankingTargetRegion_chr,
                              mat_df$Human_Chr)
    mat_df$targetRegion_chr <- ifelse(!is.na(mat_df$flankingTargetRegion_chr), 
                              paste0("chr", as.character(mat_df$flankingTargetRegion_chr)),
                              mat_df$targetRegion_chr)
    mat_df$targetRegion_start <- ifelse(!is.na(mat_df$flankingTargetRegion_start), 
                              mat_df$flankingTargetRegion_start,
                              mat_df$targetRegion_start)
    mat_df$targetRegion_end <- ifelse(!is.na(mat_df$flankingTargetRegion_end), 
                              mat_df$flankingTargetRegion_end,
                              mat_df$targetRegion_end)
    mat_df$targetRegion_strand <- ifelse(!is.na(mat_df$flankingTargetRegion_strand), 
                              mat_df$flankingTargetRegion_strand,
                              mat_df$targetRegion_strand)
    mat_df$vec_seg1_start <- ifelse((is.na(mat_df$vec_seg1_start) | mat_df$vec_seg1_start == "NA") , 
                              mat_df$fromalm_vec_seg1_start,
                              mat_df$vec_seg1_start)
    mat_df$vec_seg1_end <- ifelse((is.na(mat_df$vec_seg1_end) | mat_df$vec_seg1_end == "NA") , 
                              mat_df$fromalm_vec_seg1_end,
                              mat_df$vec_seg1_end)
    mat_df$vec_seg1_strand <- ifelse((is.na(mat_df$vec_seg1_strand) | mat_df$vec_seg1_strand == "NA") , 
                              mat_df$fromalm_vec_seg1_strand,
                              mat_df$vec_seg1_strand)
    
    # does it contain UTR primers?
    if (analyze_ddPCR_primer_set) {
      
      # mat_df$UTR5_Primer <- ifelse(!is.na(master_df[file_row_index,"UTR5_Primer_start"]) ,
      #                              ifelse( (mat_df$vec_seg1_start <= master_df[file_row_index,"UTR5_Primer_start"]+1) & (mat_df$vec_seg1_end >= master_df[file_row_index,"UTR5_Primer_end"]-1), "UTR5Primer", "Not UTR5Primer"),
      #                              NA) # ifelse
      # mat_df$UTR3_Primer <- ifelse(!is.na(master_df[file_row_index,"UTR3_Primer_start"]) ,
      #                              ifelse( (mat_df$vec_seg1_start <= master_df[file_row_index,"UTR3_Primer_start"]+1) & (mat_df$vec_seg1_end >= master_df[file_row_index,"UTR3_Primer_end"]-1), "UTR3Primer", "Not UTR3Primer"),
      #                              NA) # ifelse
      
      # ddPCR probes
      mat_df <- mat_df %>%
        mutate(ddPCR_probe = case_when(
          vec_seg1_start <= master_df[file_row_index,"ddPCR_Set1_start"] ~"ddPCR_Set1",
          vec_seg1_start <= master_df[file_row_index,"ddPCR_Set2_start"] ~"ddPCR_Set2",
          # vec_seg1_start <= master_df[file_row_index,"ddPCR_Transgenemid_start"]+1 ~"ddPCR_Transgenemid",
          # vec_seg1_start <= master_df[file_row_index,"ddPCR_Transgene3p_start"]+1 ~"ddPCR_Transgene3p",
          # vec_seg1_start <= master_df[file_row_index,"ddPCR_Promoter_start"]+1 ~"ddPCR_Promoter",
          # vec_seg1_start <= master_df[file_row_index,"ddPCR_UTR3_start"]+1 ~"ddPCR_UTR3",
          # vec_seg1_start > master_df[file_row_index,"ddPCR_UTR3_start"]+1 ~"ddPCR_NoPrimer"
          vec_seg1_start > master_df[file_row_index,"ddPCR_Set2_start"] ~"ddPCR_NoPrimer"
        ))
    } # if (analyze_ddPCR_primer_set)
    
    return (mat_df) 
  } # if (nrow(mat_df) == 0)
}



###############################################################
#' @title Acquire the TSV output file of vector only reads from distalseq pipeline
#' 
#' @author Andrea Calabria
#' @details version 0.1 (23-11-14)
#'
#' @rdname parseDistalseqPipeTSVOutput_VectorOnly
#' @docType methods
#' @aliases parseDistalseqPipeTSVOutput_VectorOnly
#'
#' @param filename to ingest of the TSV output (usually with the suffix "combined.ToHGandPLV.ReadInfo.additions.tsv.gz")
#'
#' @return df with interval IDs
#' @usage t <- parseDistalseqPipeTSVOutput(filename = "2023_07_21_RTE1_gDNA77261_B1.combined.ToHGandPLV.ReadInfo.additions.tsv.gz")
#' @note : 
#' @todo : check all required cols -> or generalize as input par
###############################################################
parseDistalseqPipeTSVOutput_VectorOnly <- function(filename,
                                        master_df,
                                        required_metadata_cols = c("SampleID", "SampleName", "Replica", "Vector", "PoolID", "GroupName", "Transgene", "SampleType", "Separator", "Promoter_start", "Promoter_end", "pA_start", "pA_end", "Functional_min", "Functional_max", "VectorOrientation", "MinReadLen"),
                                        # import_only_proper = TRUE,
                                        rownames_as_id = FALSE,
                                        prefix_filename_col = "CommonPrefixAllFiles",
                                        suffix_vectoronly_reads = ".combined.VectorOnly.AlnInfo.tsv.gz",
                                        nastring = c("NONE", "NA", "NULL", "NaN", "ND", "")) {
  
  # check all required cols in the metadata file
  if (FALSE %in% (required_metadata_cols %in% colnames(master_df))) {
    stop(paste0("[AP]\tError: Missing colnames in metadata df"))
  }
  # check that THIS FILENAME is not duplicated in the master file
  file_row_index <- NULL
  if (nrow(master_df[which(master_df$FileName == filename & master_df$ToProcess == T),]) >1) {
    stop(paste0("[AP]\tError: Duplicated filename (", filename, ") in the metadata file. "))
  } else {
    file_row_index <- rownames(master_df[which(master_df$FileName == filename & master_df$ToProcess == T),])
  }
  ## ---- read input files ------ ##
  file_name <- paste0(as.character(master_df[file_row_index, "FileFolder"]), 
                      as.character(master_df[file_row_index, prefix_filename_col]), 
                      suffix_vectoronly_reads)
  message(paste("[AP]\tProcessing File:", file_name))
  # message(paste("[AP]\t\tGet data"))
  if (as.character(master_df[file_row_index,"Separator"]) == ",") {
    mat_df <- read.csv(file = file_name, header=TRUE, fill=T, sep=',', check.names = FALSE, na.strings = nastring)
  } else {
    mat_df <- read.csv(file = file_name, header=TRUE, fill=T, sep='\t', check.names = FALSE, na.strings = nastring)
  }
  if (nrow(mat_df) == 0) {
    message(paste("[AP]\t\tSKIPPING THIS FILE SINCE EMPTY:", file_name))
    return (NULL)
  } else {
    # apply rownames if true
    if (rownames_as_id) {
      rownames(mat_df) <- paste(mat_df$Query_Name, as.character(master_df[file_row_index,"SampleID"]), sep = "_")
    }
    # min read len?
    mat_df$AboveMinReadLen <- ifelse( mat_df$Query_Len >= master_df[file_row_index,"MinReadLen"],
                                      TRUE,
                                      FALSE
    )
    # redefine "proper" label
    mat_df$Alignment <- paste0(mat_df$Target_Name, ":", mat_df$Target_Start, "-", mat_df$Target_End, mat_df$Strand )

    # group by read name and process it
    mat_df_rn <- mat_df %>% 
      group_by(Query_Name, .drop = FALSE) %>%
      dplyr::summarise(VectorAlmHits = paste0(Alignment, collapse = ";"), 
                VectorAlmStrand = paste0(Strand, collapse = ";"), 
                Num_VectorAln = n(),
                # Read_Name = Query_Name,	
                Read_Len = mean(Query_Len)
                )
      # mutate(VectorAlmStrand = paste0(Strand, collapse = ";"))  %>%
      # dplyr::count() 
    
    # add annotation col by sample
    mat_df_rn$SampleName <- as.character(master_df[file_row_index,"SampleName"])
    mat_df_rn$SampleID <- as.character(master_df[file_row_index,"SampleID"])
    mat_df_rn$SampleType <- as.character(master_df[file_row_index,"SampleType"])
    mat_df_rn$PoolID <- as.character(master_df[file_row_index,"PoolID"])
    mat_df_rn$Replica <- as.character(master_df[file_row_index,"Replica"])
    mat_df_rn$Vector <- as.character(master_df[file_row_index,"Vector"])
    mat_df_rn$Transgene <- as.character(master_df[file_row_index,"Transgene"])
    mat_df$gDNA <- as.character(master_df[file_row_index,"DNA_ID"])
    
    mat_df_rn$AboveMinReadLen <- ifelse( mat_df_rn$Read_Len >= master_df[file_row_index,"MinReadLen"],
                                      TRUE,
                                      FALSE
    )
    
    mat_df_rn$IntegrationType <- ifelse( mat_df_rn$Num_VectorAln == 1, "Single", "Rearranged") # old:: "Structural Change" | Linear vs Rearranged
    # count by split
    mat_df_rn$InvertedFragment <- apply(mat_df_rn[c("Num_VectorAln", "VectorAlmStrand")], 1, function(x){
      strand_vect <- unique(as.character(strsplit(as.character(x[2]), ';', fixed = T)[[1]]))
      ifelse( length(strand_vect) > 1, "Inverted", "Not inverted" )
      } 
    )

    return (mat_df_rn)
  } # if (nrow(mat_df) == 0)
}

###############################################################
#' @title Parse the output stats files
#' 
#' @author Andrea Calabria
#' @details version 0.1 (23-9-22)
#'
#' @rdname parseDistalseqPipeStatsFiles
#' @docType methods
#' @aliases parseDistalseqPipeStatsFiles
#'
#' @param df source DF
#' @param field_name columnd to look for the IS
#' @param dest_file
#'
#' @return df with stats
#' @usage TODO
#' @note : !!!!WARNING!!!! Vector alignments have the SAME prefix search of hybrid genome. here we are using the SECOND file, but it could change or be wrong on next updates!!!
#' @todo : check it
###############################################################
parseDistalseqPipeStatsFiles <- function(filename, 
                                         molten = TRUE,
                                         master_df,
                                         common_suffix_re = "*.NumReads.txt",
                                         common_prefix_string = "", # look for this in the master file
                                         Total_Reads_suffixFile = "All.NumReads.txt",
                                         Vector_Reads_suffixFile = "PriAlnToVector.NumReads.txt",
                                         TargetGenome_Reads_suffixFile = "PriAlnToHG.PriAlnToVector.NumReads.txt",
                                         PassingP5filter_Reads_suffixFile = "PriAlnToHG.PriAlnToVector.P5UMI_PASS.NumReads.txt",
                                         Vector_PassingP5filter_Reads_suffixFile = "VectorOnly.P5PASS.NumReads.txt",
                                         perc_scale100 = TRUE){
  ## ---- preliminar operations and checks ------ ##
  # check that THIS FILENAME is not duplicated in the master file
  file_row_index <- NULL
  if (nrow(master_df[which(master_df$FileName == filename & master_df$ToProcess == T),]) >1) {
    stop(paste0("[AP]\tError: Duplicated filename (", filename, ") in the metadata file. "))
  } else {
    file_row_index <- rownames(master_df[which(master_df$FileName == filename & master_df$ToProcess == T),])
  }
  # look for those files and get data
  files_in_folder <- list.files(path = master_df[file_row_index, "FileFolder"], 
                                pattern = common_suffix_re, full.names = T)
  common_prefix_string <- master_df[file_row_index, "CommonPrefixAllFiles"]
  
  # the numbers
  total_reads <- ifelse( file.exists(grep(Total_Reads_suffixFile, files_in_folder, value = T)), 
                         as.numeric(read.csv(file = grep(Total_Reads_suffixFile, files_in_folder, value = T), header=F, fill=T, sep='\t', check.names = FALSE)), 
                         NA )
  vector_reads <- ifelse( length(grep(paste0(common_prefix_string, ".", Vector_Reads_suffixFile), files_in_folder, value = T))>0 , 
          # as.numeric(read.csv(file = grep(Vector_Reads_suffixFile, files_in_folder, value = T)[-1], header=F, fill=T, sep='\t', check.names = FALSE)), 
          as.numeric(read.csv(file = grep(paste0(common_prefix_string, ".", Vector_Reads_suffixFile), files_in_folder, value = T), header=F, fill=T, sep='\t', check.names = FALSE)), 
          NA )
  targetgenome_reads <- ifelse( file.exists(grep(TargetGenome_Reads_suffixFile, files_in_folder, value = T)), 
                                as.numeric(read.csv(file = grep(TargetGenome_Reads_suffixFile, files_in_folder, value = T), header=F, fill=T, sep='\t', check.names = FALSE)), 
                                NA )
  chimera_passingP5_reads <- ifelse( file.exists(grep(PassingP5filter_Reads_suffixFile, files_in_folder, value = T)), 
                             as.numeric(read.csv(file = grep(PassingP5filter_Reads_suffixFile, files_in_folder, value = T), header=F, fill=T, sep='\t', check.names = FALSE)), 
                             NA )
  vector_passingP5_reads <- ifelse( length(grep(Vector_PassingP5filter_Reads_suffixFile, files_in_folder, value = T))>0, 
                                    as.numeric(read.csv(file = grep(Vector_PassingP5filter_Reads_suffixFile, files_in_folder, value = T), header=F, fill=T, sep='\t', check.names = FALSE)), 
                                    NA )
  
  stats_df <- NULL
  if (!molten) {
    message("[AP]\t\tNon molten")
    if (length(files_in_folder)>0 | !(NA %in% c(total_reads, vector_reads, targetgenome_reads, chimera_passingP5_reads)) ) {
      stats_df <- data.frame(
        "FileName" = filename,
        "FileFolder" = master_df[file_row_index, "FileFolder"],
        "SampleName" = as.character(master_df[file_row_index,"SampleName"]),
        "SampleID" = as.character(master_df[file_row_index,"SampleID"]),
        "SampleType" = as.character(master_df[file_row_index,"SampleType"]),
        "PoolID" = as.character(master_df[file_row_index,"PoolID"]),
        "Replica" = as.character(master_df[file_row_index,"Replica"]),
        "Vector" = as.character(master_df[file_row_index,"Vector"]),
        "Transgene" = as.character(master_df[file_row_index,"Transgene"]),
        "Total_Reads" = total_reads,
        "WithVector_Reads" = vector_reads,
        "WithoutVector_Reads" = total_reads - vector_reads,
        "WithVector_Reads_perc" = vector_reads / total_reads,
        "WithoutVector_Reads_perc" = (total_reads - vector_reads) / total_reads,
        "Chimeric_Reads" = targetgenome_reads,
        "VectorOnly_Reads" = vector_reads - targetgenome_reads,
        "Chimeric_Reads_perc" = targetgenome_reads / vector_reads,
        "VectorOnly_Reads_perc" = (vector_reads - targetgenome_reads) / vector_reads,
        "ChimeraPassingP5filter_Reads" = chimera_passingP5_reads,
        "VectorOnlyPassingP5filter_Reads" = vector_passingP5_reads,
        "ChimeraPassingP5filter_Reads_perc" = chimera_passingP5_reads / targetgenome_reads,
        "VectorOnlyPassingP5filter_Reads_perc" = vector_passingP5_reads / targetgenome_reads
      )  
    } else {
      stats_df <- data.frame(
        "FileName" = filename,
        "FileFolder" = master_df[file_row_index, "FileFolder"],
        "SampleName" = as.character(master_df[file_row_index,"SampleName"]),
        "SampleID" = as.character(master_df[file_row_index,"SampleID"]),
        "SampleType" = as.character(master_df[file_row_index,"SampleType"]),
        "PoolID" = as.character(master_df[file_row_index,"PoolID"]),
        "Replica" = as.character(master_df[file_row_index,"Replica"]),
        "Vector" = as.character(master_df[file_row_index,"Vector"]),
        "Transgene" = as.character(master_df[file_row_index,"Transgene"]),
        "Total_Reads" = NA,
        "WithVector_Reads" = NA,
        "WithoutVector_Reads" = NA,
        "WithVector_Reads_perc" = NA,
        "WithoutVector_Reads_perc" = NA,
        "Chimeric_Reads" = NA,
        "VectorOnly_Reads" = NA,
        "Chimeric_Reads_perc" = NA,
        "VectorOnly_Reads_perc" = NA,
        "ChimeraPassingP5filter_Reads" = NA,
        "VectorOnlyPassingP5filter_Reads" = NA,
        "ChimeraPassingP5filter_Reads_perc" = NA,
        "VectorOnlyPassingP5filter_Reads_perc" = NA
      )
    } # if else
  } else {
    first_fields <- c(filename,
                      master_df[file_row_index, "FileFolder"],
                      as.character(master_df[file_row_index,"SampleName"]),
                      as.character(master_df[file_row_index,"SampleID"]),
                      as.character(master_df[file_row_index,"SampleType"]),
                      as.character(master_df[file_row_index,"PoolID"]),
                      as.character(master_df[file_row_index,"Replica"]),
                      as.character(master_df[file_row_index,"Vector"]),
                      as.character(master_df[file_row_index,"Transgene"]))
    if (length(files_in_folder)>0 | !(NA %in% c(total_reads, vector_reads, targetgenome_reads, chimera_passingP5_reads)) ) {
      stats_df <- rbind(stats_df, as.data.frame(t(c(first_fields, total_reads, "WithVector", vector_reads, vector_reads / total_reads))))
      stats_df <- rbind(stats_df, as.data.frame(t(c(first_fields, total_reads, "WithoutVector", total_reads - vector_reads, (total_reads - vector_reads) / total_reads ))))
      stats_df <- rbind(stats_df, as.data.frame(t(c(first_fields, total_reads, "Chimeric", targetgenome_reads, targetgenome_reads / vector_reads))))
      stats_df <- rbind(stats_df, as.data.frame(t(c(first_fields, total_reads, "VectorOnly", vector_reads - targetgenome_reads, (vector_reads - targetgenome_reads) / vector_reads) )))
      stats_df <- rbind(stats_df, as.data.frame(t(c(first_fields, total_reads, "ChimeraPassingP5filter", chimera_passingP5_reads, chimera_passingP5_reads / targetgenome_reads))))
      stats_df <- rbind(stats_df, as.data.frame(t(c(first_fields, total_reads, "VectorOnlyPassingP5filter", vector_passingP5_reads, vector_passingP5_reads / targetgenome_reads))))
      names(stats_df) <- c("FileName", "FileFolder", "SampleName", "SampleID", "SampleType", "PoolID", "Replica", "Vector", "Transgene", "Total_Reads", "PipelineStep", "N_Reads", "Perc_Reads")
      # stats_df$Perc_Reads <- as.numeric(stats_df$Perc_Reads)
      stats_df$Perc_Reads <- ifelse(perc_scale100 & !is.na(stats_df$Perc_Reads), as.numeric(stats_df$Perc_Reads)*100, as.numeric(stats_df$Perc_Reads))
      stats_df$N_Reads <- as.numeric(stats_df$N_Reads)
      stats_df$Total_Reads <- as.numeric(stats_df$Total_Reads)
    } # if (length(files_in_folder)>0
  }# if molten
  return (stats_df)
}


###############################################################
#' @title Write BED file from DF
#' 
#' @author Andrea Calabria
#' @details version 0.1 (23-6-10)
#'
#' @rdname write_bed
#' @docType methods
#' @aliases write_bed
#'
#' @param df source DF
#' @param field_name columnd to look for the IS
#' @param dest_file
#'
#' @return df with interval IDs
#' @usage TODO
#' @note : 
#' @todo : check it
###############################################################
write_bed <- function(df, field_name, dest_file){
  if (field_name %in% colnames(df)) {
    write.table(x = df[which(df$Proper == "True"), c(field_name)], file = dest_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = F, na = '')
  } else {
    message(paste("[AP]\tNo column named", field_name))
  }
}


###############################################################
#' @title Get number of elements in ggplot
#' 
#' @author Andrea Calabria, from web
#' @details version 0.x
#'
#' @rdname n_fun
#' @docType methods
#' @aliases n_fun
#'
#' @param x element to explore and count
#' @param k position in the plot (percentage, scale 0-100)
#'
#' @return void
#' @usage TODO
#' @note : TODO
###############################################################
n_fun <- function(x, k){
  return(data.frame(y = c(k),
                    label = paste0("N=", length(x)) ))
}


###############################################################
#' @title Get number of reads from BAM file
#'
#' @author Andrea Calabria
#' @details version 0.1 (23-7-1)
#'
#' @rdname getStatsFromBAM
#' @docType methods
#' @aliases getStatsFromBAM
#'
#' @param stats_df the stats df output of the function parseDistalseqPipeStatsFiles
#' @param output_plot_file_prefix outptput file name (will create both png and pdf, so no extensions needed; full path is better)
#'
#' @return void
#' @usage TODO
#' @note : TODO
###############################################################
plotPoolStats <- function(stats_df,
                          output_plot_file_prefix, # with abs path
                          height=7, 
                          width=10,
                          perc_scale100 = FALSE,
                          color_schema = c(t(tessera_color_schema[1,])),
                          ... ) {
  # # melt data to have categories
  # stats_df_molten <- melt(data = stats_df, 
  #                          id.vars = setdiff(colnames(stats_df), grep("_perc", colnames(stats_df), value = T)),
  #                          variable.name = "PipelineStep", 
  #                          na.rm = T, 
  #                          value.name = "Percentage")
  # stats_df_molten$PipelineStep <- gsub("_Reads_perc", "", stats_df_molten$PipelineStep)
  # stats_df_molten_slice <- stats_df_molten[which(stats_df_molten$PipelineStep %in% c("Chimeric", "WithVector")),]
  # stats_df_molten_slice$Percentage <- ifelse(stats_df_molten_slice$PipelineStep == "WithVector", 
  #                                            stats_df_molten_slice$Percentage * -1, 
  #                                            stats_df_molten_slice$Percentage)
  # stats_df_molten_slice <- melt(data = stats_df_molten_slice, 
  #                         id.vars = setdiff(colnames(stats_df_molten_slice), grep("_Reads", colnames(stats_df_molten_slice), value = T)),
  #                         variable.name = "PipelineFeature", 
  #                         na.rm = T, 
  #                         value.name = "Reads")
  # stats_df_molten_slice$PipelineFeature <- gsub("_Reads", "", stats_df_molten_slice$PipelineFeature)
  # stats_df_molten_slice <- stats_df_molten_slice[which(stats_df_molten_slice$PipelineFeature %in% c("Total", "WithVector")),]
  stats_df$Perc_Reads <- ifelse(perc_scale100 == TRUE & !is.na(stats_df$Perc_Reads), stats_df$Perc_Reads*100, stats_df$Perc_Reads)
  stats_df$Perc_Reads <- ifelse(stats_df$PipelineStep == "WithVector",
                                stats_df$Perc_Reads * -1,
                                stats_df$Perc_Reads)
  stats_df$SampleName_forPlot <- gsub("\\.|_", " ", stats_df$SampleName)
  # breaks_values <- pretty(stats_df$Perc_Reads, n = 10)
  breaks_values <- seq(-100, 100, 20)
  
  # barplot divergent
  plot_barplot_facet <- ggplot(stats_df[which(stats_df$PipelineStep %in% c("Chimeric", "WithVector")),], 
                               aes(x = SampleID, y = Perc_Reads, color = PipelineStep, fill = PipelineStep)) +
    scale_color_manual(values = color_schema) +
    scale_fill_manual(values = color_schema) +
    geom_bar(stat = "identity") +
    # geom_text(aes(label = paste0(N_Reads, "/", Total_Reads) ), 
    geom_text(aes(label = N_Reads), 
              hjust=1.0, fontface = "bold",
              color="black"
              ) +
    # geom_text(aes(label = Total_Reads), 
    #           vjust = 0, hjust = 0,
    #           fontface = "bold",
    #           color="white"
    # ) +
    coord_flip()+
    scale_y_continuous(breaks = breaks_values,
                      labels = abs(breaks_values),
                      limits = c(-100, 100)) +
    theme_bw() +
    # facet_wrap(. ~ SampleName, nrow = 1, scales = "free") +
    # facet_wrap(. ~ Transgene, nrow = 1, scales = "free") +
    facet_grid(SampleName_forPlot ~ ., scales = "free", space = "free_y", labeller = label_wrap_gen(width=6)) +
    # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
    # theme(strip.text = element_text(face="bold", size=14)) +
    theme(strip.text = element_text(face="plain", size=14)) +
    # theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0), strip.text.y = element_text(size = 16, colour = "darkblue", angle = 270)) +
    theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0), 
          strip.text.y = element_text(size = 14, colour = "darkblue", angle = 270)) +
    theme(axis.text.x = element_text(size=14, angle = 0), 
          axis.text.y = element_text(size=16), 
          axis.title = element_text(size=16), 
          plot.title = element_text(size=22)) +
    theme(legend.position="bottom", 
          legend.text = element_text(size=14)) + 
    labs(title = paste0("DISTAL-seq Pipeline steps"), 
         x = "Samples", y = "Perc. of reads", color = "Steps", fill = "Steps") 
  
  pdf(file = paste(output_plot_file_prefix, ".Barplot.pdf", sep = ""), height=height, width=width)
  plot(plot_barplot_facet)
  dev.off()
  png(file = paste(output_plot_file_prefix, ".Barplot.png", sep = ""), height=height, width=width, units = "in", res = 300)
  plot(plot_barplot_facet)
  dev.off()
  
}


###############################################################
#' @title Plot number of reads in the filtering steps
#'
#' @author Andrea Calabria
#' @details version 0.1 (23-10-09)
#'
#' @rdname plotPoolStats_FilteringSteps
#' @docType methods
#' @aliases plotPoolStats_FilteringSteps
#'
#' @param stats_df the stats df output of the function getPoolStats_FilteringSteps + integrated with the previous stats
#' @param output_plot_file_prefix outptput file name (will create both png and pdf, so no extensions needed; full path is better)
#'
#' @return void
#' @usage TODO
#' @note : TODO
###############################################################
plotPoolStats_FilteringSteps <- function(stats_df,
                          output_plot_file_prefix, # with abs path
                          height=7, 
                          width=10,
                          perc_scale100 = FALSE,
                          pipeline_steps_order = c("WithoutVector", "WithVector", "VectorOnly", "VectorOnlyPassingP5filter", "Chimeric", "ChimeraPassingP5filter", "ChimericPassingMinLenFilter", "ChimericPassingProperFilter"),
                          cols_to_plot = c("ChimeraPassingP5filter", "ChimericPassingMinLenFilter", "ChimericPassingProperFilter"),
                          color_schema = c(t(tessera_color_schema[1,])),
                          ... ) {
  # do checks ---- TODO
  # adjust perc
  stats_df$Perc_Reads <- ifelse(perc_scale100 == TRUE & !is.na(stats_df$Perc_Reads), stats_df$Perc_Reads*100, stats_df$Perc_Reads)
  stats_df$Perc_Reads <- ifelse(stats_df$PipelineStep == "ChimeraPassingP5filter",
                                stats_df$Perc_Reads * -1,
                                stats_df$Perc_Reads)
  stats_df$PipelineStep <- factor(stats_df$PipelineStep, levels = pipeline_steps_order)
  stats_df$SampleName_forPlot <- gsub("\\.", " ", stats_df$SampleName)
  # breaks_values <- pretty(stats_df$Perc_Reads, n = 10)
  breaks_values <- seq(-100, 100, 20)
  
  # barplot divergent
  plot_barplot_facet <- ggplot(stats_df[which(stats_df$PipelineStep %in% cols_to_plot),], 
                               aes(x = SampleID, y = Perc_Reads, color = PipelineStep, fill = PipelineStep)) +
    scale_color_manual(values = color_schema) +
    scale_fill_manual(values = color_schema) +
    geom_bar(stat = "identity") +
    # geom_text(aes(label = N_Reads), 
    #           # hjust=-1.0, 
    #           fontface = "plain",
    #           color="black"
    # ) +
    coord_flip()+
    scale_y_continuous(breaks = breaks_values,
                       labels = abs(breaks_values),
                       limits = c(-100, 100)) +
    theme_bw() +
    # facet_wrap(. ~ SampleName, nrow = 1, scales = "free") +
    # facet_wrap(. ~ Transgene, nrow = 1, scales = "free") +
    facet_grid(SampleName_forPlot ~ ., scales = "free", space = "free_y", labeller = label_wrap_gen(width=6)) +
    # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
    # theme(strip.text = element_text(face="bold", size=14)) +
    theme(strip.text = element_text(face="plain", size=14)) +
    # theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0), strip.text.y = element_text(size = 16, colour = "darkblue", angle = 270)) +
    theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0), 
          strip.text.y = element_text(size = 14, colour = "darkblue", angle = 270)) +
    theme(axis.text.x = element_text(size=14, angle = 0), 
          axis.text.y = element_text(size=16), 
          axis.title = element_text(size=16), 
          plot.title = element_text(size=22)) +
    theme(legend.position="bottom", 
          legend.text = element_text(size=14)) + 
    labs(title = paste0("DISTAL-seq Pipeline steps"), 
         x = "Samples", y = "Perc. of reads", color = "Steps", fill = "Steps") 
  
  pdf(file = paste(output_plot_file_prefix, ".Barplot.pdf", sep = ""), height=height, width=width)
  plot(plot_barplot_facet)
  dev.off()
  png(file = paste(output_plot_file_prefix, ".Barplot.png", sep = ""), height=height, width=width, units = "in", res = 300)
  plot(plot_barplot_facet)
  dev.off()
  
}


###############################################################
#' @title Get pipe stats from reads (after P5 filtering)
#'
#' @author Andrea Calabria
#' @details version 0.1 (23-10-09)
#'
#' @rdname getPoolStats_FilteringSteps
#' @docType methods
#' @aliases getPoolStats_FilteringSteps
#'
#' @param df the input dataframe of the parsed reads post pipe
#'
#' @return void
#' @usage filtering_steps_stats <- getPoolStats_FilteringSteps(df = full_df, perc_scale100 = T, incremental_delta_perc = T)
#' @note : TODO
###############################################################
getPoolStats_FilteringSteps <- function(df,
                          perc_scale100 = FALSE,
                          out_colnames = c("SampleID", "SampleName", "PipelineStep", "N_Reads", "Perc_Reads"),
                          incremental_delta_perc = TRUE,
                          ... ) {
  # get some stats
  classification_read_counts <- df %>% 
    # group_by(targetRegion_chr, integration_locus, targetRegion_strand, .drop = FALSE) %>%
    group_by(SampleID, SampleName, AboveMinReadLen, Proper, .drop = FALSE) %>%
    dplyr::count() 
  # do percentages
  classification_read_counts <- classification_read_counts %>% 
    group_by(SampleID, SampleName) %>% 
    dplyr::mutate(PassingP5Filter_Reads = sum(n))
  classification_read_counts <- classification_read_counts %>% 
    group_by(SampleID, SampleName, AboveMinReadLen) %>% 
    dplyr::mutate(MinLenFilter_Reads = sum(n))
  classification_read_counts$ByLenProper_Perc_Reads <- ifelse(perc_scale100 & !is.na(classification_read_counts$MinLenFilter_Reads),
                                                              100 * (classification_read_counts$n / classification_read_counts$PassingP5Filter_Reads),
                                                              (classification_read_counts$n / classification_read_counts$PassingP5Filter_Reads) )
  if (incremental_delta_perc) {
    classification_read_counts$ByLen_Perc_Reads <- ifelse(perc_scale100 & !is.na(classification_read_counts$MinLenFilter_Reads), 
                                                          100 * (classification_read_counts$MinLenFilter_Reads / classification_read_counts$PassingP5Filter_Reads) - classification_read_counts$ByLenProper_Perc_Reads,
                                                          (classification_read_counts$MinLenFilter_Reads / classification_read_counts$PassingP5Filter_Reads)  - classification_read_counts$ByLenProper_Perc_Reads )
  } else {
    classification_read_counts$ByLen_Perc_Reads <- ifelse(perc_scale100 & !is.na(classification_read_counts$MinLenFilter_Reads), 
                                                          100 * (classification_read_counts$MinLenFilter_Reads / classification_read_counts$PassingP5Filter_Reads),
                                                          (classification_read_counts$MinLenFilter_Reads / classification_read_counts$PassingP5Filter_Reads)  )
  }
  
  # filter only TRUE vals (=passing filters reads) and create a melt structure
  classification_read_counts_slice <- classification_read_counts[which(classification_read_counts$AboveMinReadLen == TRUE & classification_read_counts$Proper == "True"),]
  classification_read_counts_slice_len <- classification_read_counts_slice[c("SampleID", "SampleName", "MinLenFilter_Reads", "ByLen_Perc_Reads")]
  names(classification_read_counts_slice_len) <- c("SampleID", "SampleName", "N_Reads", "Perc_Reads")
  classification_read_counts_slice_len$PipelineStep <- "ChimericPassingMinLenFilter"
  classification_read_counts_slice_proper <- classification_read_counts_slice[c("SampleID", "SampleName", "n", "ByLenProper_Perc_Reads")]
  names(classification_read_counts_slice_proper) <- c("SampleID", "SampleName", "N_Reads", "Perc_Reads")
  classification_read_counts_slice_proper$PipelineStep <- "ChimericPassingProperFilter"
  classification_read_counts_slice_bind <- rbind(classification_read_counts_slice_len[out_colnames], 
                                                 classification_read_counts_slice_proper[out_colnames])
  return( classification_read_counts_slice_bind )
}

###############################################################
#' @title Get pipe stats from reads (after P5 filtering), vector only reads
#'
#' @author Andrea Calabria
#' @details version 0.1 (23-11-13)
#'
#' @rdname getPoolStats_FilteringSteps_VectorOnly
#' @docType methods
#' @aliases getPoolStats_FilteringSteps_VectorOnly
#'
#' @param df the input dataframe of the parsed reads post pipe
#'
#' @return void
#' @usage filtering_steps_stats <- getPoolStats_FilteringSteps_VectorOnly(df = full_df, perc_scale100 = T, incremental_delta_perc = T)
#' @note : TODO
###############################################################
getPoolStats_FilteringSteps_VectorOnly <- function(df,
                                        perc_scale100 = FALSE,
                                        out_colnames = c("SampleID", "SampleName", "PipelineStep", "N_Reads", "Perc_Reads"),
                                        # incremental_delta_perc = TRUE,
                                        ... ) {
  # get some stats
  classification_read_counts <- df %>% 
    # group_by(targetRegion_chr, integration_locus, targetRegion_strand, .drop = FALSE) %>%
    group_by(SampleID, SampleName, AboveMinReadLen, .drop = FALSE) %>%
    dplyr::count() 
  # do percentages
  classification_read_counts <- classification_read_counts %>% 
    group_by(SampleID, SampleName) %>% 
    dplyr::mutate(PassingP5Filter_Reads = sum(n))
  classification_read_counts <- classification_read_counts %>% 
    group_by(SampleID, SampleName, AboveMinReadLen) %>% 
    dplyr::mutate(MinLenFilter_Reads = sum(n))
  classification_read_counts$ByLen_Perc_Reads <- ifelse(perc_scale100 & !is.na(classification_read_counts$MinLenFilter_Reads), 
                                                        100 * (classification_read_counts$MinLenFilter_Reads / classification_read_counts$PassingP5Filter_Reads),
                                                        (classification_read_counts$MinLenFilter_Reads / classification_read_counts$PassingP5Filter_Reads)  )
  
  # filter only TRUE vals (=passing filters reads) and create a melt structure
  classification_read_counts_slice <- classification_read_counts[which(classification_read_counts$AboveMinReadLen == TRUE),]
  classification_read_counts_slice_len <- classification_read_counts_slice[c("SampleID", "SampleName", "MinLenFilter_Reads", "ByLen_Perc_Reads")]
  names(classification_read_counts_slice_len) <- c("SampleID", "SampleName", "N_Reads", "Perc_Reads")
  classification_read_counts_slice_len$PipelineStep <- "VectorOnlyPassingMinLenFilter"
  
  return( classification_read_counts_slice_len[out_colnames] )
}


###############################################################
#' @title Reshape stats by changing the denominator
#'
#' @author Andrea Calabria
#' @details version 0.1 (23-10-23)
#'
#' @rdname rescaleStatsPerStep
#' @docType methods
#' @aliases rescaleStatsPerStep
#'
#' @param df stats_extended for all steps with N_Reads and Perc_Reads cols
#' @param colnames colnames of the actual df stats. default: c("SampleID", "SampleName", "PipelineStep", "N_Reads", "Perc_Reads")
#'
#' @return void
#' @usage 
#' @note : TODO
###############################################################
rescaleStatsPerStep <- function(df_steps, 
                                sorted_cols_to_rescale = c("ChimericPassingProperFilter", "ChimericPassingMinLenFilter", "ChimeraPassingP5filter"),
                                ref_col_denominator = c("Chimeric"),
                                annotation_cols = c("SampleID", "SampleName"),
                                perc_scale100 = FALSE,
                                exclude_missing_values = TRUE,
                                exclude_missing_values_refcol = "ChimericPassingMinLenFilter",
                                ... ) {
  # cast data
  steps_cast <- dcast(data = df_steps, 
                      SampleID + SampleName ~ PipelineStep, 
                      value.var = "N_Reads", fun.aggregate = sum)
  rownames(steps_cast) <- steps_cast$SampleID
  # stats_df <- rbind(stats_df, as.data.frame(t(c(first_fields, total_reads, "WithVector", vector_reads, vector_reads / total_reads))))
  steps_cast_incremental <- steps_cast[c(sorted_cols_to_rescale, ref_col_denominator)]
  if (length(sorted_cols_to_rescale) > 1) {
    for (k in seq(length(sorted_cols_to_rescale), 2)){
      # update values
      steps_cast_incremental[k] <- steps_cast_incremental[k] - steps_cast_incremental[k-1]
    }
  }
  # now in percentage
  steps_cast_incremental_perc <- steps_cast_incremental[sorted_cols_to_rescale]
  for (k in seq(1, length(sorted_cols_to_rescale))){
    # update values
    steps_cast_incremental_perc[k] <- steps_cast_incremental_perc[k] / steps_cast[ref_col_denominator]
  }
  if (perc_scale100){
    steps_cast_incremental_perc <- steps_cast_incremental_perc*100
  }
  # add annotations
  steps_cast_incremental <- cbind(steps_cast[annotation_cols], steps_cast_incremental)
  steps_cast_incremental_perc <- cbind(steps_cast[annotation_cols], steps_cast_incremental_perc)
  if (exclude_missing_values) {
    missing_samples  <- rownames(steps_cast[which(steps_cast[,exclude_missing_values_refcol] %in% c(0, NA)),])
    steps_cast_incremental <- steps_cast_incremental[setdiff(rownames(steps_cast_incremental), missing_samples),]
    steps_cast_incremental_perc <- steps_cast_incremental_perc[setdiff(rownames(steps_cast_incremental_perc), missing_samples),]
  }
  # return data after molten
  steps_cast_incremental_molten <- 
    melt(data = steps_cast_incremental, 
         id.vars = annotation_cols,
         variable.name = "PipelineStep", 
         na.rm = F, 
         value.name = "N_Reads")
  steps_cast_incremental_perc_molten <- 
    melt(data = steps_cast_incremental_perc, 
         id.vars = annotation_cols,
         variable.name = "PipelineStep", 
         na.rm = F, 
         value.name = "Perc_Reads")
  steps_cast_incremental_full <- merge(x = steps_cast_incremental_molten, 
                                       y = steps_cast_incremental_perc_molten,
                                       by = c(annotation_cols, "PipelineStep"))
  return( steps_cast_incremental_full )
}

###############################################################
#' @title Read sequences from Distalseq plugin (DA)
#'
#' @author Andrea Calabria
#' @details version 0.1 (24-04-24, the day 0)
#'
#' @rdname importFlankingSeqs
#' @docType methods
#' @aliases importFlankingSeqs
#'
#' @param file to read.
#' @param ... ...
#'
#' @return df
#' @usage 
#' @note : TODO
###############################################################
importFlankingSeqs <- function(seq_file, ... ) {
  message(paste0("[AP]\tImport flanking sequences from ", seq_file))
  if (file.exists(seq_file)) {
    seq_df <- read.csv(file = seq_file, 
                       header=TRUE, fill=T, 
                       sep='\t', check.names = FALSE, 
                       na.strings = c("NONE", "NA", "NULL", "NaN", "ND", ""))
    rownames(seq_df) <- seq_df$ReadID
    return (seq_df)
  } else {
    message(paste0("[AP]\tERROR: ", seq_file, " NOT FOUND!"))
  }
}



