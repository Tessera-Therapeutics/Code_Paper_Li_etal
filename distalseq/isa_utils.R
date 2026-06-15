#### ISA Utils ####

###############################################################
#' @title Write BED file with a window up/downdtream for sequence extraction
#' 
#' @author Andrea Calabria
#' @details version 0.1 (23-12-19)
#'
#' @rdname getFaFromExtendBedWindow
#' @docType methods
#' @aliases getFaFromExtendBedWindow
#'
#' @param df the input dataframe of a BED file! not a matrix.
#' @param size +/- window size. Default: 20: output will have X-20, X+20 bp.
#'
#' @return a new df. NB: it runs bedtools with system call!
#' @usage TODO
#' @note : TODO
###############################################################
getFaFromExtendBedWindow <- function(df, 
                                     size = 20,
                                     bedfile_output, 
                                     target_genome_fasta_file, 
                                     tab_outformat = T) {
  # TODO: need to check strand specific options. Currently the window is centered
  df_extendedbed <- as.data.frame( t(
    apply(df, 1, function(x) {
      c(x[1], (as.numeric(x[2])-size), (as.numeric(x[3])+size), paste0(x[4], ":", x[6]), x[5], x[6])
    }
    ) )
  )
  write.table(x = df_extendedbed, 
              file = bedfile_output, 
              sep = "\t", col.names = F, row.names = F, na = '', quote = F)
  bed_bname <- basename(bedfile_output)
  bed_dname <- dirname(bedfile_output)
  bed_suffix <- gsub(".bed", "", bed_bname)
  bed_sorted <- paste(bed_dname, "/", bed_suffix, ".sorted.bed", sep = "")
  fa_stranded <- paste(bed_dname, "/", bed_suffix, ".sorted.fa", sep = "")
  system(command = paste0("bedtools sort -i '", 
                          bedfile_output, "' > '",
                          bed_sorted, "'"))
  if (tab_outformat) {
    system(command = paste0("bedtools getfasta -name -tab -s -fi '", 
                            target_genome_fasta_file, 
                            "' -bed '",
                            bed_sorted, "' -fo '", 
                            fa_stranded, "'"))
  } else {
    system(command = paste0("bedtools getfasta -name -s -fi '", 
                            target_genome_fasta_file, 
                            "' -bed '",
                            bed_sorted, "' -fo '", 
                            fa_stranded, "'"))
  }
  return (list("bed_sorted_file" = bed_sorted, 
               "fa_stranded" = fa_stranded,
               "bed_sorted_df" = df_extendedbed) )
}


###############################################################
#' @title Write rownames in DF
#' 
#' @author Andrea Calabria
#' @details version 0.1 (23-6-10)
#'
#' @rdname rownamesAsIS
#' @docType methods
#' @aliases rownamesAsIS
#'
#' @param df source DF
#' @param id_cols annotation cols to combine
#'
#' @return df with annotated rows
#' @usage TODO
#' @note : 
#' @todo : check it
###############################################################
rownamesAsIS <- function(df, 
                         id_cols = c("chr", "integration_locus", "integration_strand", "GeneName", "GeneStrand", "GeneDistance"),
                         only_coordinates = FALSE
                        ){
  if (!(FALSE %in% (id_cols %in% colnames(df)) )) {
    if (only_coordinates) {
      rownames(df) <- apply(df[id_cols], 1, function(x) {
        paste0(x[1], "_", as.character(as.numeric(x[2])), "_", x[3])
          }
        )
    } else {
      rownames(df) <- apply(df[id_cols], 1, function(x) {
        paste0(x[1], "_", as.character(as.numeric(x[2])), "_", x[3], "_", x[4], "_", x[5])
          }
        )
    } # if (only_coordinates) {
  } else {
    stop(paste0("[AP]\tError: one or more columns in id_cols are not present in your input df. "))
  }
  return (df)
}



###############################################################
#' @title Compute pairwise Fisher test on gene frequencies from CIS results
#' 
#' @author Andrea Calabria
#' @details version 0.1, 2020-07-08
#'
#' @rdname compareGeneFrequency_Fisher
#' @docType methods
#' @aliases compareGeneFrequency_Fisher
#'
#' @param g1_cis_df CIS results df group 1.
#' @param g2_cis_df CIS results df group 2.
#' @param min_is_per_gene Min number of IS per gene to be reported as CIS gene.
#'
#' @return a df of gene frequences 
#' @usage TODO
#' @description This function is aimed at comparing the groups (G1 and G2) in terms of gene frequencies
#' by applying Fisher exact test on the confusion matrix 2x2 by gene (FDR corrected).
#' @note : 
#'
###############################################################
compareGeneFrequency_Fisher <- function(df_g1, df_g2, 
                                        # annotation_cols = c("chr", "integration_locus", "strand", "GeneName", "GeneStrand"),
                                        # gene_name_col = "GeneName",
                                        # gene_strand_col = "GeneStrand",
                                        # chr_name_col = "chr",
                                        group_by_cols = c("chr", "GeneName", "GeneStrand"),
                                        # g1_cis_df, g2_cis_df, 
                                        min_is_per_gene = 3, 
                                        gene_set_method = "INTERSECTION",
                                        # scary_genes_toannotate = c("Lmo2", "Smg6", "Mecom", "Mds", "Ccnd2"),
                                        remove_only_unbalanced_0 = TRUE ){
  message(paste("[AP]\tComparing G1 and G2 gene frequencies with Fisher test."))
  # filename_infix <- ".Intersection" # ".Intersection" ".Union"
  # gene_set_method <- "INTERSECTION" # "INTERSECTION" "UNION"
  # min_is_per_gene <- 3
  
  g1_cis_df <- df_g1 %>% 
    dplyr::group_by_at(group_by_cols) %>% 
    dplyr::summarise(n_IS_perGene = n())
    # mutate(nISperGene = n())
  g1_cis_df$TotalIS <- sum(g1_cis_df$n_IS_perGene)
  
  g2_cis_df <- df_g2 %>% 
    dplyr::group_by_at(group_by_cols) %>% 
    dplyr::summarise(n_IS_perGene = n())
    # mutate(nISperGene = n())
  g2_cis_df$TotalIS <- sum(g2_cis_df$n_IS_perGene)
  
  g1_cis_df_forfreqplot <- g1_cis_df[which(g1_cis_df$n_IS_perGene >= min_is_per_gene),]
  g2_cis_df_forfreqplot <- g2_cis_df[which(g2_cis_df$n_IS_perGene >= min_is_per_gene),]
  
  study_mergegenename <- NULL
  if (gene_set_method == "UNION") {
    ### -------- version 1: Group 1 OR Group 2 genes ----------- equivalent to UNION ###
    union_of_genes <- union(as.character(g1_cis_df_forfreqplot[,c("GeneName")]), as.character(g2_cis_df_forfreqplot[,c("GeneName")]))
    new_overall_genet_cast_merge <- merge(x = g1_cis_df[which(g1_cis_df$GeneName %in% union_of_genes),], 
                                          y = g2_cis_df[which(g2_cis_df$GeneName %in% union_of_genes),], 
                                          by = group_by_cols)
    names(new_overall_genet_cast_merge) <- gsub("\\.x", "_G1", colnames(new_overall_genet_cast_merge))
    names(new_overall_genet_cast_merge) <- gsub("\\.y", "_G2", colnames(new_overall_genet_cast_merge))
    study_mergegenename <- new_overall_genet_cast_merge  
  } else if (gene_set_method == "INTERSECTION") {
    ### -------- version 2: Group 1 AND Group 2 genes ----------- equivalent to INTERSECTION ###
    # gene intersection
    # overall_genet_cast_merge <- merge(x = g1_cis_df_forfreqplot, y = g2_cis_df_forfreqplot, by = annotation_cols_to_get)
    overall_genet_cast_merge <- merge(x = g1_cis_df_forfreqplot, 
                                      y = g2_cis_df_forfreqplot, 
                                      by = group_by_cols)
    names(overall_genet_cast_merge) <- gsub("\\.x", "_G1", colnames(overall_genet_cast_merge))
    names(overall_genet_cast_merge) <- gsub("\\.y", "_G2", colnames(overall_genet_cast_merge))
    # use a backup
    study_mergegenename <- overall_genet_cast_merge
  } else {
    message(paste0("[AP]\tERROR: wrong method selected for gene_set_method."))
  }
  
  # check len and return results
  if (nrow(study_mergegenename) == 0) {
    message(paste0("[AP]\tWARNING: No genes in comomn resulted in G1 and G2 for frequency comparison. Rerutning NULL."))
    return (NULL)
  } else {
    ### --- now move ahead wuth the rest --- 
    # now do comparison by gene frequency
    # study_mergegenename$scary <- ifelse(study_mergegenename$GeneName %in% scary_genes_toannotate, TRUE, FALSE)
    selected_col_for_fisher <- c("n_IS_perGene_G1", "TotalIS_G1", "n_IS_perGene_G2", "TotalIS_G2")
    # selected_col_for_fisher <- c("rd027_NISperGene", "rd027_TotalIS", "tm035_NISperGene", "tm035_TotalIS")
    ft_pval <- data.frame("FisherTest_pvalue" = apply(study_mergegenename[selected_col_for_fisher], 1, function(x) {
      .m <- matrix(c(x[1], x[2]-x[1], x[3], x[4]-x[3]), nrow = 2, dimnames = list(RD027 = c("ISofGene", "TotalIS"), TM035 = c("ISofGene", "TotalIS")))
      .ft <- fisher.test(.m)
      .ft$p.value
    })
    )
    study_mergegenename <- cbind(study_mergegenename, ft_pval)
    study_mergegenename$FisherTest_pvalue_significant <- ifelse(study_mergegenename$FisherTest_pvalue < 0.05, TRUE, FALSE)
    
    ### ---- to correct the pvalue, first try including or not the 0
    # get the min number of IS for 0 genes in the other sample
    mean_nis_pergene_a <- ceiling(mean(study_mergegenename$n_IS_perGene_G1[study_mergegenename$n_IS_perGene_G1>0]))
    mean_nis_pergene_c <- ceiling(mean(study_mergegenename$n_IS_perGene_G2[study_mergegenename$n_IS_perGene_G2>0]))
    study_mergegenename$to_exclude_from_test <- apply(study_mergegenename[selected_col_for_fisher], 1, function(x) {
      ifelse( (x[1] == 0 | x[3] == 0),
              ifelse( (((x[1] < mean_nis_pergene_a) & (x[3] == 0)) | ((x[1] == 0) & (x[3] < mean_nis_pergene_c))),
                      TRUE,
                      FALSE),
              FALSE
      )
    })
    
    if (remove_only_unbalanced_0) {
      study_mergegenename <- study_mergegenename[which(study_mergegenename$to_exclude_from_test == FALSE),] 
    }
    
    # corrected pvalue
    # study_mergegenename$FisherTest_pvalue_bonferroni <- p.adjust(study_mergegenename$FisherTest_pvalue, method = "bonferroni", n = length(study_mergegenename$FisherTest_pvalue))
    study_mergegenename$FisherTest_pvalue_fdr <- p.adjust(study_mergegenename$FisherTest_pvalue, method = "fdr", n = length(study_mergegenename$FisherTest_pvalue))
    study_mergegenename$FisherTest_pvalue_benjamini <- p.adjust(study_mergegenename$FisherTest_pvalue, method = "BY", n = length(study_mergegenename$FisherTest_pvalue))
    # add vars for volcano plot
    study_mergegenename$geneIS_frequency_byHitIS_G1 <- study_mergegenename$n_IS_perGene_G1 / study_mergegenename$TotalIS_G1
    study_mergegenename$geneIS_frequency_byHitIS_G2 <- study_mergegenename$n_IS_perGene_G2 / study_mergegenename$TotalIS_G2
    study_mergegenename$log_FC <- log(x = (study_mergegenename$geneIS_frequency_byHitIS_G1 / study_mergegenename$geneIS_frequency_byHitIS_G2), base = 10)
    study_mergegenename$FC_G1G2 <- study_mergegenename$geneIS_frequency_byHitIS_G1 / study_mergegenename$geneIS_frequency_byHitIS_G2
    study_mergegenename$log_pvalue <- -log(x = study_mergegenename$FisherTest_pvalue, base = 10)
    study_mergegenename$log_pvalue_fdr <- -log(x = study_mergegenename$FisherTest_pvalue_fdr, base = 10)
    # study_mergegenename$log_pvalue <- -log(x = study_mergegenename$FisherTest_pvalue, base = 10)
    study_mergegenename$log2_FC <- log(x = (study_mergegenename$geneIS_frequency_byHitIS_G1 / study_mergegenename$geneIS_frequency_byHitIS_G2), base = 2)
    # study_mergegenename$log2_FC_TPM <- log(x = (study_mergegenename$IS_per_kbGeneLen_perMDepth_TPM_G1 / study_mergegenename$IS_per_kbGeneLen_perMDepth_TPM_G2), base = 2)
    
    return(study_mergegenename)
  } # if (nrow(study_mergegenename) == 0)
  
}


###############################################################
#' @title Plot chimera in CIRCOS
#' 
#' @author Andrea Calabria
#' @details version 0.1, 2020-11-19; version 0.2 Tessera, 2023-12-07
#'
#' @rdname plot_circos_rearrangement_imp
#' @docType methods
#' @aliases plot_circos_rearrangement_imp
#'
#' @param df_alm input df with alignment results
#' @param min_cigar_alm_width minimum alignment subread widt for each CIGAR string. default = 5. This value will be applied to ALL tags (indels included)
#'
#' @return plot
#' @description parse CIGAR string and get data for AAV studies (and general chimera). Definition of chimera (as subread): the ONLY portion of the read (identified by CIGAR) that flanks target genome alignment.
#' @note : 
#' @usage 
#' plot_circos_rearrangement_imp(df = alm_reads_for_rearrangements_withstatsbyread, sample_read = "m64047_200620_005306/19269857/ccs", 
#' outfile_png = paste0(dest_dir, "results/", ProjectID, ".plot_alm_reads_stats.circos.InVivo.V-V-V-V-X.01.png", sep = ""), 
#' cols_to_search = c("chr", "start", "end", "ratio_bp_aligned_on_raw"), vector_chr = "chrV", zoomed_chr_index = 23, vector_cytoband_file = "source/metadata/CAG_Tomato_withBackbone/AAV-CAG-tdTomato.withBackbone.cytoband")
###############################################################
plot_circos_rearrangement_imp <- function(df, sample_read, outfile_png, outfile_pdf,
                                          cols_to_search = c("chr", "start", "end", "ratio_bp_aligned_on_raw"), 
                                          vector_chr = "chrV", 
                                          zoomed_chr_index = 23, 
                                          vector_cytoband_file = "source/metadata/CAG_Tomato_withBackbone/AAV-CAG-tdTomato.withBackbone.cytoband",
                                          species = "hg38",
                                          bp_res= 300,
                                          color_vector_fwd="orange",
                                          color_vector_rev="green",
                                          color_target_rev="violet",
                                          color_target_fwd="blue",
                                          Elem5p_name = "UTR5",
                                          Elem3p_name = "UTR3",
                                          force_order_from_vector = T, 
                                          title = '',
                                          subtitle_vectorgenomename = '',
                                          subtitle_targetgenomename = 'Hg38') {
  if (FALSE %in% (cols_to_search %in% colnames(df)) ) {
    message(paste0("[AP]\tERROR: Columns in df are different from expected: ", paste(cols_to_search, collapse = ',')))
    return (NULL)
  } else {
    # slice data
    t <- sample_read
    slice_t <- df[which(df$name == t),]
    slice_t <- slice_t[order(slice_t$query_start, decreasing = F),]
    # case specific to DISTALseq
    if (force_order_from_vector) {
      # force order to have the first alignment in the vector (starting with chrV or NOT chr)
      if (startsWith(slice_t[1,"chr"], "chr")) {
        slice_t <- slice_t[order(slice_t$query_start, decreasing = T),]
        if (slice_t[1,"strand"] == "+") {
          slice_t$strand <- ifelse(slice_t$strand == "+", "-", "+")
        }
      }
    }
    
    bed_1 <- slice_t[(1:(nrow(slice_t)-1)), cols_to_search]
    names(bed_1) <- c("chr", "start", "end", "value1")
    bed_2 <- slice_t[(2:(nrow(slice_t))), cols_to_search]
    names(bed_2) <- c("chr", "start", "end", "value1")
    
    # orientations
    bed_all_stranded <- slice_t[c("chr", "start", "end", "strand")]
    bed_all_stranded$strand <- ifelse(bed_all_stranded$strand == '+', 1, 0)
    names(bed_all_stranded) <- c("chr", "start", "end", "value1")
    
    # get start-end of each read
    bed_all_long <- slice_t[c("chr", "start", "end", "strand", "query_start", "query_end")]
    outbed_line_sorted_start <- NULL
    if (nrow(bed_all_long) > 1) {
      for (i in seq(1, (nrow(bed_all_long)-1))) {
        if (bed_all_long[i,"strand"] == '+') {
          from_chr <- bed_all_long[i,c("chr")]
          from_pos <- bed_all_long[i,c("end")]
          if (bed_all_long[i+1,"strand"] == '+') {
            to_chr <- bed_all_long[i+1,c("chr")]
            to_pos <- bed_all_long[i+1,c("start")]
          } else {
            to_chr <- bed_all_long[i+1,c("chr")]
            to_pos <- bed_all_long[i+1,c("end")]
          } # if (bed_all_long[i+1,"strand"] == '+')
        } else {
          from_chr <- bed_all_long[i,c("chr")]
          from_pos <- bed_all_long[i,c("start")]
          if (bed_all_long[i+1,"strand"] == '+') {
            to_chr <- bed_all_long[i+1,c("chr")]
            to_pos <- bed_all_long[i+1,c("start")]
          } else {
            to_chr <- bed_all_long[i+1,c("chr")]
            to_pos <- bed_all_long[i+1,c("end")]
          } # if (bed_all_long[i+1,"strand"] == '+')
        } # if (bed_all_long[i,"strand"] == '+')
        if (length(outbed_line_sorted_start) > 0) {
          outbed_line_sorted_start <- rbind(outbed_line_sorted_start, data.frame("from_chr" = from_chr,
                                                                                 "from_pos" = from_pos,
                                                                                 "to_chr" = to_chr,
                                                                                 "to_pos" = to_pos))  
        } else {
          outbed_line_sorted_start <- data.frame("from_chr" = from_chr,
                                                 "from_pos" = from_pos,
                                                 "to_chr" = to_chr,
                                                 "to_pos" = to_pos)
        } # if (length(outbed_line_sorted_start) > 0)
        
      } # for (i in seq(1, nrow(bed_all_long)))
    } # if (nrow(bed_all_long) > 1)
    outbed_line_sorted_start$value1 <- 1
    bed_1_line <- outbed_line_sorted_start[grep("from", colnames(outbed_line_sorted_start))]
    names(bed_1_line) <- c("chr", "start")
    bed_1_line$end <- bed_1_line$start
    bed_1_line$value1 <- 1
    bed_2_line <- outbed_line_sorted_start[grep("to", colnames(outbed_line_sorted_start))]
    names(bed_2_line) <- c("chr", "start")
    bed_2_line$end <- bed_2_line$start
    bed_2_line$value1 <- 1
    
    target_chr <- gsub("chr", "", setdiff( levels(factor(slice_t$chr)), vector_chr))
    # if 0 -> only AAV (chrV), else target too
    if (length(target_chr) > 0) {
      # go into genomics
      human_cytoband <- read.cytoband(species = species)$df
      vector_cytoband <- read.csv(file = vector_cytoband_file, header=F, fill=T, check.names = FALSE, sep = '\t')
      cytoband_rbind <- rbind(human_cytoband, vector_cytoband)
      cytoband <- read.cytoband(cytoband_rbind)
      
      cytoband_df = cytoband$df
      chromosome = cytoband$chromosome
      # acquire target chr
      if (target_chr %in% c("X", "Y")) {
        if (target_chr == "X") {target_chr <- grep("X", chromosome)}
        if (target_chr == "Y") {target_chr <- grep("Y", chromosome)}
      } else {
        target_chr <- as.numeric(target_chr)
      }
      # determine range of zoom
      xrange = c(cytoband$chr.len, cytoband$chr.len[vector_chr])
      normal_chr_index = target_chr
      zoomed_chr_index = zoomed_chr_index
      sector.width = c(xrange[normal_chr_index] / sum(xrange[normal_chr_index]), 
                       xrange[zoomed_chr_index] / sum(xrange[zoomed_chr_index])) 
      # PNG
      png(file = outfile_png, height=5, width=5, units = "in", res = 600)
      extended <-extend_chromosomes(cytoband_df, vector_chr)
      chr_target_ff<- bed_all_stranded[bed_all_stranded$chr!=vector_chr,]$chr
      s <- bed_all_stranded[bed_all_stranded$chr==chr_target_ff,]$start
      e <- bed_all_stranded[bed_all_stranded$chr==chr_target_ff,]$end
      m <- (s+e)/2
      r <- extended[(extended$V1==chr_target_ff) & (extended$V2<m) & (extended$V3>m) ,]
      increase<-slice_t$read_len[1]+1000
      r$V2 <- m-(increase/2)
      r$V3 <- m+(increase/2)
      extended[(extended$V1==chr_target_ff) & (extended$V2<m) & (extended$V3>m) ,] <- r
      itr5 <- extended[extended$V4==Elem5p_name & extended$V1==vector_chr,]$V2
      itr3 <- extended[extended$V4==Elem3p_name & extended$V1==vector_chr,]$V3
      extended2 <- extended[((extended$V1==chr_target_ff) & (extended$V2<m) & (extended$V3>m)) | (extended$V1!=chr_target_ff & extended$V1!=vector_chr)
                            | ((extended$V1==vector_chr & (extended$V2>=itr5 & extended$V3<=itr3))),]
      extended2$V6 <-''
      extended2[extended2$V1==vector_chr,]$V6 <- extended2[extended2$V1==vector_chr,]$V4
      
      circos.initializeWithIdeogram(extended2, 
                                    chromosome.index = c(chromosome[target_chr], vector_chr),
                                    sector.width = sector.width, tickLabelsStartFromZero = FALSE, major.by = bp_res)
      
      f = colorRamp2(breaks = c(0, 1,2,3), colors = c(color_vector_rev, color_vector_fwd,color_target_rev,color_target_fwd))
      l <- seq(m-5000,m+5000,600)
      for (subread in seq(1, nrow(bed_all_stranded))) {
        circos.genomicTrackPlotRegion(bed_all_stranded[subread,], stack = TRUE, 
                                      track.height = 0.05, bg.col = NA, bg.border = "gray90",
                                      panel.fun = function(region, value, ...) {
                                        x1=region$start
                                        y1=region$end
                                        if(value==1){
                                          v<-value
                                          if (bed_all_stranded[subread,]$chr!=vector_chr)
                                            v<-v+2
                                          circos.arrow(x1,y1, col = f(v), 
                                                       border = 1, arrow.head.length = cm_x(0.2) ) 
                                        }
                                        else{
                                          v<-value
                                          if (bed_all_stranded[subread,]$chr!=vector_chr)
                                            v<-v+2
                                          circos.arrow(x1, y1, col = f(v), 
                                                       border = 1, arrow.position = "start", arrow.head.length = cm_x(0.2))
                                        }
                                        
                                        
                                        i = getI(...)
                                        cell.xlim = get.cell.meta.data("cell.xlim")
                                        # circos.lines(cell.xlim, c(i, i), lty = 2, col = "#FFFFFF")
                                      })
        
        
      }
      # circos.genomicLink(bed_1, bed_2, col = rand_color(nrow(bed_1), transparency = 0.5))
      circos.genomicLink(bed_1_line, bed_2_line, col = "black", directional = 1)
      title(main = title)
      text(-0.9, -0.8, paste0("Target\ngenome\n", subtitle_targetgenomename))
      text(0.9, 0.8, paste0("Vector\ngenome\n", subtitle_vectorgenomename))
      dev.off()
      circos.clear()
      # PDF
      pdf(file = outfile_pdf, height=5, width=5)
      circos.initializeWithIdeogram(extended2, 
                                    chromosome.index = c(chromosome[target_chr], vector_chr),
                                    sector.width = sector.width, tickLabelsStartFromZero = FALSE, major.by = bp_res)
      
      f = colorRamp2(breaks = c(0, 1,2,3), colors = c(color_vector_rev, color_vector_fwd,color_target_rev,color_target_fwd))
      for (subread in seq(1, nrow(bed_all_stranded))) {
        circos.genomicTrackPlotRegion(bed_all_stranded[subread,], stack = TRUE, 
                                      track.height = 0.05, bg.col = NA, bg.border = "gray90",
                                      panel.fun = function(region, value, ...) {
                                        x1=region$start
                                        y1=region$end
                                        if(value==1){
                                          v<-value
                                          if (bed_all_stranded[subread,]$chr!=vector_chr)
                                            v<-v+2
                                          circos.arrow(x1,y1, col = f(v), 
                                                       border = 1, arrow.head.length = cm_x(0.2) ) 
                                        }
                                        else{
                                          v<-value
                                          if (bed_all_stranded[subread,]$chr!=vector_chr)
                                            v<-v+2
                                          circos.arrow(x1, y1, col = f(v), 
                                                       border = 1, arrow.position = "start", arrow.head.length = cm_x(0.2))
                                        }
                                        
                                        
                                        i = getI(...)
                                        cell.xlim = get.cell.meta.data("cell.xlim")
                                        # circos.lines(cell.xlim, c(i, i), lty = 2, col = "#FFFFFF")
                                      })
      }
      # circos.genomicLink(bed_1, bed_2, col = rand_color(nrow(bed_1), transparency = 0.5))
      circos.genomicLink(bed_1_line, bed_2_line, col = "black", directional = 1)
      title(main = title)
      text(-0.9, -0.8, paste0("Target\ngenome\n", subtitle_targetgenomename))
      text(0.9, 0.8, paste0("Vector\ngenome\n", subtitle_vectorgenomename))
      dev.off()
      circos.clear()
      
    } else {
      vector_cytoband <- read.csv(file = vector_cytoband_file, header=F, fill=T, check.names = FALSE, sep = '\t')
      # cytoband_rbind <- rbind(human_cytoband, vector_cytoband)
      cytoband <- read.cytoband(vector_cytoband)
      
      cytoband_df = cytoband$df
      chromosome = cytoband$chromosome
      # PNG
      png(file = outfile_png, height=5, width=5, units = "in", res = 600)
      circos.par("start.degree" = 90)
      itr5 <- extended[extended$V4==Elem5p_name & extended$V1==vector_chr,]$V2
      itr3 <- extended[extended$V4==Elem3p_name & extended$V1==vector_chr,]$V3
      cytoband_df <- cytoband_df[(cytoband_df$V2>=itr5 & cytoband_df$V3<=itr3),]
      
      circos.initializeWithIdeogram(cytoband = cytoband_df,
                                    chromosome.index = c(vector_chr), 
                                    tickLabelsStartFromZero = FALSE, major.by = bp_res)
      f = colorRamp2(breaks = c(0, 1,2,3), colors = c(color_vector_rev, color_vector_fwd,color_target_rev,color_target_fwd))
      for (subread in seq(1, nrow(bed_all_stranded))) {
        circos.genomicTrackPlotRegion(bed_all_stranded[subread,], stack = TRUE, 
                                      track.height = 0.05, bg.col = NA, bg.border = "gray90",
                                      panel.fun = function(region, value, ...) {
                                        x1=region$start
                                        y1=region$end
                                        if(value==1){
                                          v<-value
                                          if (bed_all_stranded[subread,]$chr!=vector_chr)
                                            v<-v+2
                                          circos.arrow(x1,y1, col = f(v), 
                                                       border = 1, arrow.head.length = cm_x(0.2) ) 
                                        }
                                        else{
                                          v<-value
                                          if (bed_all_stranded[subread,]$chr!=vector_chr)
                                            v<-v+2
                                          circos.arrow(x1, y1, col = f(v), 
                                                       border = 1, arrow.position = "start", arrow.head.length = cm_x(0.2))
                                        }
                                        
                                        
                                        i = getI(...)
                                        cell.xlim = get.cell.meta.data("cell.xlim")
                                        # circos.lines(cell.xlim, c(i, i), lty = 2, col = "#FFFFFF")
                                      })
      }
      # circos.genomicLink(bed_1, bed_2, col = rand_color(nrow(bed_1), transparency = 0.5))
      circos.genomicLink(bed_1_line, bed_2_line, col = "black", directional = 1)
      title(main = title)
      text(-0.9, -0.8, paste0("Target\ngenome\n", subtitle_targetgenomename))
      text(0.9, 0.8, paste0("Vector\ngenome\n", subtitle_vectorgenomename))
      dev.off()
      circos.clear()
      # PDF
      pdf(file = outfile_pdf, height=5, width=5)
      circos.par("start.degree" = 90)
      circos.initializeWithIdeogram(cytoband = cytoband_df, chromosome.index = c(vector_chr),
                                    tickLabelsStartFromZero = FALSE, major.by = bp_res)
      f = colorRamp2(breaks = c(0, 1,2,3), colors = c(color_vector_rev, color_vector_fwd,color_target_rev,color_target_fwd))
      for (subread in seq(1, nrow(bed_all_stranded))) {
        circos.genomicTrackPlotRegion(bed_all_stranded[subread,], stack = TRUE, 
                                      track.height = 0.05, bg.col = NA, bg.border = "gray90",
                                      panel.fun = function(region, value, ...) {
                                        x1=region$start
                                        y1=region$end
                                        if(value==1){
                                          v<-value
                                          if (bed_all_stranded[subread,]$chr!=vector_chr)
                                            v<-v+2
                                          circos.arrow(x1,y1, col = f(v), 
                                                       border = 1, arrow.head.length = cm_x(0.2) ) 
                                        }
                                        else{
                                          v<-value
                                          if (bed_all_stranded[subread,]$chr!=vector_chr)
                                            v<-v+2
                                          circos.arrow(x1, y1, col = f(v), 
                                                       border = 1, arrow.position = "start", arrow.head.length = cm_x(0.2))
                                        }
                                        
                                        
                                        i = getI(...)
                                        cell.xlim = get.cell.meta.data("cell.xlim")
                                        # circos.lines(cell.xlim, c(i, i), lty = 2, col = "#FFFFFF")
                                      })
      }
      # circos.genomicLink(bed_1, bed_2, col = rand_color(nrow(bed_1), transparency = 0.5))
      circos.genomicLink(bed_1_line, bed_2_line, col = "black", directional = 1)
      title(main = title)
      text(-0.9, -0.8, paste0("Target\ngenome\n", subtitle_targetgenomename))
      text(0.9, 0.8, paste0("Vector\ngenome\n", subtitle_vectorgenomename))
      dev.off()
      circos.clear()
    }  # if (length(target_chr) > 0) 
  } # if (FALSE %in% (cols_to_search %in% colnames(df)) )
  
} # function

extend_chromosomes = function(bed, chromosome, prefix = "zoom_") {
  zoom_bed = bed[bed[[1]] %in% chromosome, , drop = FALSE]
  zoom_bed[[1]] = paste0(prefix, zoom_bed[[1]])
  rbind(bed, zoom_bed)
}


###############################################################
#' @title Plot IS CIRCOS
#' 
#' @author Andrea Calabria
#' @details version 0.1, 2024-01-24
#'
#' @rdname plot_circos_IS
#' @docType methods
#' @aliases plot_circos_IS
#'
#' @param df_alm input df with alignment results
#' @param min_cigar_alm_width minimum alignment subread widt for each CIGAR string. default = 5. This value will be applied to ALL tags (indels included)
#'
#' @return plot
#' @description TODO
#' @note : 
#' @usage 
#' plot_circos_IS(df = alm_reads_for_rearrangements_withstatsbyread, sample_read = "m64047_200620_005306/19269857/ccs", 
#' outfile_png = paste0(dest_dir, "results/", ProjectID, ".plot_alm_reads_stats.circos.InVivo.V-V-V-V-X.01.png", sep = ""), 
#' cols_to_search = c("chr", "start", "end", "ratio_bp_aligned_on_raw"), vector_chr = "chrV", zoomed_chr_index = 23, vector_cytoband_file = "source/metadata/CAG_Tomato_withBackbone/AAV-CAG-tdTomato.withBackbone.cytoband")
###############################################################
plot_circos_IS <- function(df, outfile_png, outfile_pdf,
                           vector_cols = c("Vector", "vec_seg1_start"), 
                           targetgenome_cols = c("targetRegion_chr", "integration_locus"), 
                           vector_chr = "chrV", 
                           junction_strand_col = c("vec_seg1_strand"),
                           integration_strand_col = c("integration_strand"),
                           bed_colnames = c("chr","start","end","value1", "score", "strand"),
                           basenumber_for_full_integration = 5,
                           # zoomed_chr_index = 23, 
                           vector_cytoband_file = "BEDFILE",
                           species = "hg38",
                           # bp_res= 300,
                           # color_vector_fwd="orange",
                           # color_vector_rev="green",
                           # color_target_rev="violet",
                           # color_target_fwd="blue",
                           Elem5p_name = "UTR5",
                           Elem3p_name = "UTR3",
                           exclude_backbone = TRUE,
                           force_order_from_vector = T, 
                           title = '',
                           subtitle_vectorgenomename = '',
                           subtitle_targetgenomename = 'Hg38') {
  if (FALSE %in% ( c(vector_cols, targetgenome_cols) %in% colnames(df)) ) {
    message(paste0("[AP]\tERROR: Columns in df are different from expected: ", paste(cols_to_search, collapse = ',')))
    return (NULL)
  } else {
    # ========== STABLE  ========================= ##
    human_cytoband <- read.cytoband(species = species)$df
    vector_cytoband <- read.csv(file = vector_cytoband_file, header=F, fill=T, check.names = FALSE, sep = '\t')
    cytoband_rbind <- rbind(human_cytoband, vector_cytoband)
    cytoband <- read.cytoband(cytoband_rbind)
    
    cytoband_df = cytoband$df
    chromosome = cytoband$chromosome
    
    xrange = c(cytoband$chr.len[grep("chr", names(cytoband$chr.len) )], cytoband$chr.len[grep(vector_chr, names(cytoband$chr.len) )])
    normal_chr_index = grep("chr", names(xrange))
    zoomed_chr_index = grep(vector_chr, names(xrange))
    
    # normalize in normal chromsomes and zoomed chromosomes separately
    sector.width = c(xrange[normal_chr_index] / sum(xrange[normal_chr_index]), 
                     xrange[zoomed_chr_index] / sum(xrange[zoomed_chr_index])) 
    
    # links IS - junction
    vector_junction_df <- df[which(df$targetRegion_chr %in% chromosome), vector_cols]
    names(vector_junction_df) <- c("chr", "pos")
    vector_junction_df$chr <- vector_chr
    
    target_IS_df <- df[which(df$targetRegion_chr %in% chromosome), targetgenome_cols]
    target_IS_df$integration_locus <- as.numeric(target_IS_df$integration_locus)
    
    bed_vector_start <- df[which(df$targetRegion_chr %in% chromosome), c(vector_cols, vector_cols[2], "Read_Name", "Read_Len", "vec_seg1_strand")]
    bed_vector_start$Vector <- vector_chr
    names(bed_vector_start) <- bed_colnames
    bed_vector_start$start <- bed_vector_start$start + 1
    bed_vector_start$end <- bed_vector_start$end + 1
    
    bed_is <- df[which(df$targetRegion_chr %in% chromosome), c(targetgenome_cols, targetgenome_cols[2], "Read_Name", "Read_Len", "integration_strand")]
    names(bed_is) <- bed_colnames
    bed_is$start <- as.numeric(bed_is$start)
    bed_is$end <- as.numeric(bed_is$end)
    bed_df <- rbind(bed_is, bed_vector_start[1:10,])
    
    # arrange view on vector, visible ONLY within UTRs
    if (exclude_backbone) {
      itr5 <- cytoband_df[cytoband_df$V4==Elem5p_name & cytoband_df$V1==vector_chr,]$V2
      itr3 <- cytoband_df[cytoband_df$V4==Elem3p_name & cytoband_df$V1==vector_chr,]$V3
      # cytoband_df <- cytoband_df[(cytoband_df$V2>=itr5 & cytoband_df$V3<=itr3),]
      cytoband_df <- cytoband_df[(!(cytoband_df$V1==vector_chr & (cytoband_df$V2<itr5 | cytoband_df$V3>itr3))),]
    }
    # circos
    png(file = outfile_png, height=5, width=5, units = "in", res = 600)
    
    circos.par(start.degree = 5)
    circos.initializeWithIdeogram(cytoband_df, 
                                  chromosome.index = c(chromosome[grep("chr", chromosome)], vector_chr),
                                  sector.width = sector.width,
                                  tickLabelsStartFromZero = FALSE)
    # circos.genomicDensity(bed_df, col = c("navyblue"), track.height = theight)
    # circos.genomicDensity(bed_vector_start[1:100,], col = c("green"), track.height = theight)
    rcols <- scales::alpha(ifelse(vector_junction_df$pos < basenumber_for_full_integration, "#f46d43", "#66c2a5"), alpha=0.2)
    circos.genomicLink(region1 = vector_junction_df, 
                       region2 = target_IS_df, 
                       col = rcols, border = NA)
    title(main = title)
    text(-0.9, -0.8, paste0("Target\ngenome\n", subtitle_targetgenomename))
    text(0.9, 0.8, paste0("Vector\ngenome\n", subtitle_vectorgenomename))
    dev.off()
    
    circos.clear()
    # PDF
    pdf(file = outfile_pdf, height=5, width=5)
    circos.par(start.degree = 5)
    circos.initializeWithIdeogram(cytoband_df, 
                                  chromosome.index = c(chromosome[grep("chr", chromosome)], vector_chr),
                                  sector.width = sector.width,
                                  tickLabelsStartFromZero = FALSE)
    # circos.genomicDensity(bed_df, col = c("navyblue"), track.height = theight)
    # circos.genomicDensity(bed_vector_start[1:100,], col = c("green"), track.height = theight)
    rcols <- scales::alpha(ifelse(vector_junction_df$pos < basenumber_for_full_integration, "#f46d43", "#66c2a5"), alpha=0.2)
    circos.genomicLink(region1 = vector_junction_df, 
                       region2 = target_IS_df, 
                       col = rcols, border = NA)
    title(main = title)
    text(-0.9, -0.8, paste0("Target\ngenome\n", subtitle_targetgenomename))
    text(0.9, 0.8, paste0("Vector\ngenome\n", subtitle_vectorgenomename))
    dev.off()
    circos.clear()
    ## =================================== ##
  }
} # function





# valore di riferimento: numero di IS condivise tra campioni -> fissato un paziente, per ogni tempo|campione fai rapporto di quante sono le IS condivise prma e dopo (diventa un rapporto N_condivise/N_osservate)
getSharedISratio <- function(df, compact = TRUE, left_to_rigth_reading_output = TRUE, prune_selected_rows = FALSE, pruning_rows_labels = c()) {
  # in questa funzione, a partire dal df (RxC), per ogni colonna c  in C fissi le is di c (ovvero quelle >0 not NA) e conti quante sono condivise restituendo un vettore
  message("[AP]\tConverting input df by adding 0 to NA")
  df[is.na(df)] <- 0 # avoid here inside NA
  if (compact) {
    df <- compactDfByColumns(df)
  }
  list_container <- list() # init the list of resulting object
  message(paste("[AP]\tNow looping over columns to comput sharing results"))
  c_index <- 1
  for (c in colnames(df)) {
    message(paste("[AP]\t-> processing the colum\t", c, "\tposition", as.character(c_index), "of", as.character(ncol(df)), "\t[", as.character(round(c_index*100/ncol(df),2)), "%]"))
    slice_df <- df[which(df[c]>0),] # slice df
    list_container <- c( list_container, data.frame(c = apply(slice_df, 2, function(x) {length(x[x>0])/nrow(slice_df)})) ) # the relativ eprecentage of contaminations
    c_index <- c_index + 1
  }
  names(list_container) <- colnames(df) # rename list objects
  r <- as.data.frame(list_container)
  rownames(r) <- colnames(slice_df)
  ### NB: in questo momento la matrice (df r) ha lettura alto basso (!!!!), e non sinistra destra. in base all'opzione invertila nel return
  if (left_to_rigth_reading_output) {
    message(paste("[AP]\tTranspose output matrix: Left->Right orientation"))
    # if need to prune_selected_rows
    if (prune_selected_rows & length(pruning_rows_labels)>0) {
      message(paste("[AP]\tRemoving selected rows"))
      r <- r[,!colnames(r) %in% pruning_rows_labels]
      return (t(r))
    } else {
      return (t(r))
    }
  } 
  else {
    message(paste("[AP]\tOutput matrix orientation: Top-Down"))
    if (prune_selected_rows & length(pruning_rows_labels)>0) {
      message(paste("[AP]\tRemoving selected rows"))
      r <- r[!rownames(r) %in% pruning_rows_labels,]
      return (r)
    } else {
      return (r)
    }
  }
}

# come la gemella getSharedISratio, cambia il tipo di calcolo
getSharedISnumber <- function(df, compact = TRUE, left_to_rigth_reading_output = TRUE, prune_selected_rows = FALSE, pruning_rows_labels = c()) {
  # in questa funzione, a partire dal df (RxC), per ogni colonna c  in C fissi le is di c (ovvero quelle >0 not NA) e conti quante sono condivise restituendo un vettore
  message("[AP]\tConverting input df by adding 0 to NA")
  df[is.na(df)] <- 0 # avoid here inside NA
  if (compact) {
    df <- compactDfByColumns(df)
  }
  list_container <- list() # init the list of resulting object
  message(paste("[AP]\tNow looping over columns to comput sharing results"))
  c_index <- 1
  for (c in colnames(df)) {
    message(paste("[AP]\t-> processing the colum\t", c, "\tposition", as.character(c_index), "of", as.character(ncol(df)), "\t[", as.character(round(c_index*100/ncol(df),2)), "%]"))
    slice_df <- df[which(df[c]>0),] # slice df
    list_container <- c( list_container, data.frame(c = apply(slice_df, 2, function(x) {length(x[x>0])})) ) # the relativ eprecentage of contaminations
    c_index <- c_index + 1
  }
  names(list_container) <- colnames(df) # rename list objects
  r <- as.data.frame(list_container)
  rownames(r) <- colnames(slice_df)
  ### NB: in questo momento la matrice (df r) ha lettura alto basso (!!!!), e non sinistra destra. in base all'opzione invertila nel return
  if (left_to_rigth_reading_output) {
    message(paste("[AP]\tTranspose output matrix: Left->Right orientation"))
    # if need to prune_selected_rows
    if (prune_selected_rows & length(pruning_rows_labels)>0) {
      message(paste("[AP]\tRemoving selected rows"))
      r <- r[,!colnames(r) %in% pruning_rows_labels]
      return (t(r))
    } else {
      return (t(r))
    }
  } 
  else {
    message(paste("[AP]\tOutput matrix orientation: Top-Down"))
    if (prune_selected_rows & length(pruning_rows_labels)>0) {
      message(paste("[AP]\tRemoving selected rows"))
      r <- r[!rownames(r) %in% pruning_rows_labels,]
      return (r)
    } else {
      return (r)
    }
  }
}


compactDfByColumns <- function(df) {
  # compct df by pruning columns without any value (all rows of these columns are NA or 0)
  # 1. identify null columns
  # 2. return the same df witout these columns
  colnames(df) 
  notnull_cols <- apply(df, 2, function(x) {nrow(df) - length(x[ x == 0 | is.na(x) ])})
  null_col_names <- as.vector(names(notnull_cols[notnull_cols==0]))
  null_col_index <- colnames(df) %in% null_col_names
  return (df[!null_col_index])
}

compactDfByRows <- function(df, data_columns, annotation_columns = NULL) {
  # compact df by pruning rows without any value (meaning that no value for any columns, all NA or 0)
  # 1. identify null rows
  # 2. return the same df witout these rows
  # Example usage: t <- compactDfByRows(gdf, data_columns = grep("BThal", colnames(gdf)))
  message("[AP]\tCompacting by rows (pruning rows without values)")
  if (length(data_columns) == 0) {
    stop("[AP]\tERROR: Arguments data_columns must be specified. Try with grep, example: compactDfByRows(gdf, data_columns = grep('BThal', colnames(gdf)))")
  } 
  col_sum <- data.frame("SUM" = apply(df[data_columns], 1, function(x) {abs(sum(x, na.rm = TRUE))}), 
                        "ID" = rownames(df[data_columns])
  ) # use absolute values so that only positive values are the rows to keep
  rows_to_keep <- rownames(col_sum[which(col_sum$SUM > 0),])
  if (length(annotation_columns) > 0) {
    return (df[rows_to_keep, c(annotation_columns, data_columns)])
  } else {
    return (df[rows_to_keep, data_columns]) 
  }
}



aggregateDfColumnsByName <- function(df, metadata_df, key_field, starting_data_col_index = 2, number_of_last_cols_to_remove = 0) {
  # given the input df, look up the column  names and aggregate by name
  # so far,  this is a by-hand procedure... waiting for better ideas
  # for each col name, slice df and apply function
  colID_df <- data.frame("colID" = colnames(df)[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)]) # get column names of df into a new df
  rownames(colID_df) <- colnames(df)[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)]
  colmetadata_df <- merge(colID_df, metadata_df, by=0, all.x = TRUE) # get the output merge df of annotations with all x metadata
  # get names from the key field
  colmetadata_df_outnames <- colmetadata_df[,c(key_field)]
  
  # produce a warning IF two columns share the same pattern matching (prefix or suffix)
  if (length(grep(key_field, colnames(colmetadata_df)))>1) {
    message(paste("[AP]\tWARNING: in your metadata file you have two columns with a very similar name, this is not allowed. I.e.: Group and TestGroup."))
  }
  
  output_cols <- as.character(levels(factor(colmetadata_df_outnames)))
  #message(paste("[AP]\tThese will be the output/aggregated columns:", output_cols))
  out_df <- NULL
  for (k in output_cols) {
    # get col IDs from metadata
    #metadata_df[which(metadata_df[grep(key_field, colnames(metadata_df))] == k), c("colID")]
    message(paste("[AP]\tProcessing column(s)", k))
    slice_cols <- as.character(colmetadata_df[which(colmetadata_df[grep(key_field, colnames(colmetadata_df))] == k), c("colID")])
    if (length(out_df) == 0) { # if this is the first loop
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        out_df <- df[slice_cols]
        out_df[is.na(out_df)] <- 0
        names(out_df) <- k
      } else {
        out_df <- as.data.frame(apply(df[, slice_cols], 1, function(x) {sum(x, na.rm = TRUE)}))
        names(out_df) <- k
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    } else { # if (length(out_df) == 0) { # if this is the first loop
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        actual_colnames <- colnames(out_df)
        this_slice <- df[, slice_cols]
        this_slice[is.na(this_slice)] <- 0
        out_df <- cbind(out_df, this_slice )
        names(out_df) <- c(actual_colnames, k)
      } else {
        actual_colnames <- colnames(out_df)
        out_df <- cbind(out_df, as.data.frame(apply(df[, slice_cols], 1, function(x) {sum(x, na.rm = TRUE)})) )
        names(out_df) <- c(actual_colnames, k)
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    }
  } # else of if (length(out_df) == 0) { # if this is the first loop
  # return df
  return (out_df)
}


computeAbundancePercentage <- function(df = df, metadata_df = NULL, key_field = "QuantificationSum", starting_data_col_index = 1, last_data_col_index = NA){
  # given a matrix of data only, compute abundance
  message(paste("[AP]\tCalculating abundance from input matrix"))
  # evaluate column length
  if (is.na(last_data_col_index)) {
    last_data_col_index <- length(colnames(df))
  }
  source_df <- df[starting_data_col_index:last_data_col_index]
  # compute internally the sum df, and call the column of sum as key_field
  # fill NA with 0
  source_df[is.na(source_df)] <- 0
  metadata_df <- data.frame( 
    "QuantificationSum" = apply(source_df, 2, function(x) {sum(x, na.rm=T)} ),
    "NumIS" = apply(source_df, 2, function(x) {length(x[x>0])} )
  )
  names(metadata_df) <- c(key_field, "NumIS")
  # first check if you have any emopty column (no rows with values >0)
  overall_sum <- apply(source_df, 2, function(x) {sum(x, na.rm = TRUE)})
  if (min(overall_sum) == 0) {
    message(paste("[AP]\t\tWARNING: You have empty cols (cols with no values)"))
  }
  # do abundance
  abundance_df <- data.frame(sapply(colnames(source_df), function(x) {
    # standard relative %:
    .quantification_sum <- metadata_df[x, key_field]
    if (.quantification_sum > 0) {
      ( source_df[,x]*100 / (metadata_df[x, key_field]) )
    } else {
      ( source_df[,x]*0 )
    }
  } )
  ) # do the abundance
  rownames(abundance_df) <- rownames(df) # change names
  colnames(abundance_df) <- colnames(source_df)
  # do warnings
  if (min(abundance_df)<0) {
    message(paste("[AP]\t\tThe abundance is producing values <0 (!!!)"))
  }
  # return out df
  return (abundance_df)
}

getTopHitGenes <- function(df, data_column_to_use, group_by_colums, maxGenesToReturn = 10, column_suffix = c("Gene Symbol", "Annotated IS")) {
  # input: df, data columns to extract results
  # logics: get the top HIT genes (not by abundance)
  # todo: control on max number of returning elements:: if dat aelements are < than max -> check it
  # output: df with values
  result_df <- NULL 
  for (k in data_column_to_use) {
    message(paste("[AP]\tPocessing column(s)", k))
    aggregated_slice <- aggregate(df[which(df[k] > 0), c(k)], by=list(df[which(df[k] > 0), group_by_colums]), FUN=function(x) {length(x[x>0])})
    aggregated_slice_sorted <- aggregated_slice[order(-c(aggregated_slice$x)),]
    # get output with values
    if (length(result_df) == 0) {
      # init out df
      result_df <- as.data.frame(aggregated_slice_sorted[1:maxGenesToReturn, ])
      names(result_df) <- c(paste(k, column_suffix[1]), paste(k, column_suffix[2]))
    }
    else {
      actual_colnames <- colnames(result_df)
      result_df <- cbind(result_df, as.data.frame(aggregated_slice_sorted[1:maxGenesToReturn, ]) )
      names(result_df) <- c(actual_colnames, paste(k, column_suffix[1]), paste(k, column_suffix[2]))
    }
  }
  return (result_df)
}

getTopCIS_withIScount <- function(df, data_column_to_use, group_by_colum, maxGenesToReturn = 10, column_suffix = c("Gene Symbol", "Annotated IS")) {
  # input: df of melted data (!!!!), data values to extract results
  # logics: get first the top CIS and their MAX score, then combine their number of annotated landed IS (IS count)
  # todo: control on max number of returning elements:: if dat aelements are < than max -> check it
  # output: df with values
  result_df <- NULL 
  for (k in data_column_to_use) {
    message(paste("[AP]\tPocessing column(s)", k))
    aggregated_slice <- aggregate(df[which(df[k] > 0), c(k)], by=list(df[which(df[k] > 0), group_by_colums]), FUN=function(x) {length(x[x>0])})
    aggregated_slice_sorted <- aggregated_slice[order(-c(aggregated_slice$x)),]
    # get output with values
    if (length(result_df) == 0) {
      # init out df
      result_df <- as.data.frame(aggregated_slice_sorted[1:maxGenesToReturn, ])
      names(result_df) <- c(paste(k, column_suffix[1]), paste(k, column_suffix[2]))
    }
    else {
      actual_colnames <- colnames(result_df)
      result_df <- cbind(result_df, as.data.frame(aggregated_slice_sorted[1:maxGenesToReturn, ]) )
      names(result_df) <- c(actual_colnames, paste(k, column_suffix[1]), paste(k, column_suffix[2]))
    }
  }
  return (result_df)
}



###############################################################
#' @title CIS with Grubbs
#' 
#' @author Andrea Calabria
#' @details version 0.1 (20 September 2019) 
#'
#' @rdname CISGrubbs
#' @docType methods
#' @aliases CISGrubbs
#'
#' @param df an input dataframe of IS matrix.
#' @param annotation_cols columns of annotation
#' @param genomic_annotation_genebased_file Annotation file: gene based genomic annotation file. It must be formatted UCSC based, see description for details.
#' @param grubbs_flankingene_bp Extra bp for flanking genes. Default 100000.
#' @param threshold_alpha P-value for filtering significant genes. Default 0.05
#' @param gene_name_col Column name in input df for column gene. Default = "GeneName",
#' @param chr_name_col Column name in input df for column chr. Default = "chr",
#' @param add_standard_padjust Do you want to compute standard p.adjust methods? default = TRUE
#' @param compactDfByRows If the matrix requires be compacted, set this to TRUE (Default = F)
#'
#' @return DF of CIS results
#' @usage todo
#' @description Do CIS analysis with Grubbs test for outliers (as EM and DC designed).
#' 
#' Gene name File formatting:
#'      name2 chrom strand min_txStart max_txEnd minmax_TxLen average_TxLen      name min_cdsStart max_cdsEnd minmax_CdsLen average_CdsLen
# 1     A1BG chr19      -    58858171  58864865         6694          6694 NM_130786     58858387   58864803          6416        6416.00
# 2 A1BG-AS1 chr19      +    58863335  58866549         3214          3214 NR_015380     58866549   58866549             0           0.00
# 3     A1CF chr10      -    52559168  52645435        86267         86267 NM_138933     52566488   52619700         53212       48635.50
# 4      A2M chr12      -     9220303   9268825        48522         48522 NM_000014      9220418    9268445         48027       46276.75
# 5  A2M-AS1 chr12      +     9217772   9220651         2879          2879 NR_026971      9220651    9220651             0           0.00
# 6    A2ML1 chr12      +     8975067   9029377        54310         43044 NM_144670      8975247    9027607         52360       41100.00
#' 
###############################################################
CISGrubbs <- function(df,
                      annotation_cols = c("chr", "integration_locus", "strand", "GeneName", "GeneStrand"),
                      gene_name_col = "GeneName",
                      chr_name_col = "chr",
                      genomic_annotation_genebased_file, 
                      grubbs_flankingene_bp = 100000, 
                      threshold_alpha = 0.05,
                      add_standard_padjust = TRUE, 
                      compactDfByRows = FALSE
) {
  
  library(sqldf)
  
  ok_checks_passed <- TRUE
  if (!file.exists(genomic_annotation_genebased_file)) {
    ok_checks_passed <- FALSE
  }
  
  if (ok_checks_passed) {
    if (compactDfByRows) {
      message(paste("[AP]\tCompacting input dataframe (NOTE: be sure that all data columns are real data columns, no 'all' cols accepted."))
      df <- compactDfByRows(df = df, annotation_columns = annotation_cols, data_columns = setdif(colnames(df), annotation_cols))
    }
    
    # annotations
    message(paste("[AP]\tCIS analysis started, files are ok. Importing data."))
    refgenes <- read.csv(file = genomic_annotation_genebased_file, 
                         header=TRUE, fill=T, sep='\t', 
                         check.names = FALSE, na.strings = c("NONE", "NA", "NULL", "NaN", ""))
    
    df_bygene <- sqldf("select chr, integration_locus, GeneName, GeneStrand, count(distinct integration_locus) as n_IS_perGene, 
                          min(integration_locus) as min_bp_integration_locus, max(integration_locus) as max_bp_integration_locus, (max(integration_locus) - min(integration_locus)) as IS_span_bp, avg(integration_locus) as avg_bp_integration_locus, median(integration_locus) as median_bp_integration_locus, count(distinct strand) as distinct_orientations
                       from df 
                       where 1 
                       group by GeneName, chr")
    df_bygene$chr <- gsub("^", "chr", df_bygene$chr)
    # df_bygene$TotIS_asSumByGene <- sum(df_bygene$n_IS_perGene) # sum of IS
    
    df_bygene_withannotation <- sqldf("select df_bygene.chr, df_bygene.integration_locus, df_bygene.GeneName, df_bygene.n_IS_perGene, refgenes.average_TxLen, df_bygene.min_bp_integration_locus, df_bygene.max_bp_integration_locus, df_bygene.avg_bp_integration_locus, df_bygene.median_bp_integration_locus, df_bygene.distinct_orientations, df_bygene.IS_span_bp
                                      from df_bygene, refgenes
                                      where df_bygene.chr = refgenes.chrom and df_bygene.GeneStrand = refgenes.strand and df_bygene.GeneName = refgenes.name2 ")
    df_bygene_withannotation$TotIS_asDfRow <- nrow(df)
    df_bygene_withannotation$geneIS_frequency_byHitIS <- df_bygene_withannotation$n_IS_perGene / df_bygene_withannotation$TotIS_asDfRow
    
    # add columns as grubbs test requires
    message(paste("[AP]\tComputing tests."))
    extrabp <- grubbs_flankingene_bp
    n_elements <- nrow(df_bygene_withannotation)
    df_bygene_withannotation$raw_gene_integration_frequency <- df_bygene_withannotation$n_IS_perGene / df_bygene_withannotation$average_TxLen
    df_bygene_withannotation$integration_frequency_withtolerance <- (df_bygene_withannotation$n_IS_perGene / (df_bygene_withannotation$average_TxLen + extrabp)) * 1000
    df_bygene_withannotation$minus_log2_integration_freq_withtolerance <- -log(x = df_bygene_withannotation$integration_frequency_withtolerance, base = 2)
    # average_minus_log2_integration_freq <- mean(df_bygene_withannotation$minus_log2_integration_freq_withtolerance)
    # stdev_minus_log2_integration_freq <- sd(df_bygene_withannotation$minus_log2_integration_freq_withtolerance)
    # df_bygene_withannotation$ratioZEM_minus_log2_integration_freq_withtolerance <- (average_minus_log2_integration_freq - df_bygene_withannotation$minus_log2_integration_freq_withtolerance) / stdev_minus_log2_integration_freq
    
    df_bygene_withannotation$zscore_minus_log2_integration_freq_withtolerance <- scale(-log(x = df_bygene_withannotation$integration_frequency_withtolerance, base = 2))
    df_bygene_withannotation$neg_zscore_minus_log2_integration_freq_withtolerance <- -scale(-log(x = df_bygene_withannotation$integration_frequency_withtolerance, base = 2))
    df_bygene_withannotation$t_z_mlif <- sqrt( (n_elements * (n_elements - 2) * (df_bygene_withannotation$neg_zscore_minus_log2_integration_freq_withtolerance)^2) /
                                                 ( ((n_elements -1)^2) - (n_elements * (df_bygene_withannotation$neg_zscore_minus_log2_integration_freq_withtolerance)^2))
    )
    
    # df_bygene_withannotation$minus_log2_integration_freq <- -log(x = ((df_bygene_withannotation$raw_gene_integration_frequency)*1000), base = 2)
    # df_bygene_withannotation$z_minus_log_integration_freq <- scale(df_bygene_withannotation$minus_log2_integration_freq)
    # df_bygene_withannotation$t_z_mlif <- sqrt( (nrow(df_bygene_withannotation) * (nrow(df_bygene_withannotation) - 2) * (df_bygene_withannotation$z_minus_log_integration_freq)^2) /
    #                                              ( ((nrow(df_bygene_withannotation) -1)^2) - (nrow(df_bygene_withannotation) * (df_bygene_withannotation$z_minus_log_integration_freq)^2))
    #                                            )
    
    df_bygene_withannotation$tdist2t <- T.DIST.2T(df_bygene_withannotation$t_z_mlif, df = (nrow(df_bygene_withannotation) - 2))
    # df_bygene_withannotation$tdist <- dt(x = df_bygene_withannotation$t_z_mlif, df = nrow(df_bygene_withannotation) - 2) * 2
    # df_bygene_withannotation$tdist <- (1 - dt(x = df_bygene_withannotation$t_z_mlif, df = nrow(df_bygene_withannotation) - 2))
    df_bygene_withannotation$tdist_pt <- pt(q = df_bygene_withannotation$t_z_mlif, df = nrow(df_bygene_withannotation) - 2)
    
    df_bygene_withannotation$tdist_bonferroni_AC <- ifelse(df_bygene_withannotation$tdist2t * nrow(df_bygene_withannotation) > 1, 1, df_bygene_withannotation$tdist2t * nrow(df_bygene_withannotation))
    if (add_standard_padjust) {
      df_bygene_withannotation$tdist_bonferroni <- p.adjust(df_bygene_withannotation$tdist2t, method = "bonferroni", n = length(df_bygene_withannotation$tdist2t))
      df_bygene_withannotation$tdist_fdr <- p.adjust(df_bygene_withannotation$tdist2t, method = "fdr", n = length(df_bygene_withannotation$tdist2t))
      df_bygene_withannotation$tdist_benjamini <- p.adjust(df_bygene_withannotation$tdist2t, method = "BY", n = length(df_bygene_withannotation$tdist2t))
    }
    
    df_bygene_withannotation$tdist_positive_and_corrected <- ifelse((df_bygene_withannotation$tdist_bonferroni_AC < threshold_alpha & df_bygene_withannotation$neg_zscore_minus_log2_integration_freq_withtolerance > 0), 
                                                                    df_bygene_withannotation$tdist_bonferroni_AC, 
                                                                    NA)
    df_bygene_withannotation$tdist_positive <- ifelse((df_bygene_withannotation$tdist2t < threshold_alpha & df_bygene_withannotation$neg_zscore_minus_log2_integration_freq_withtolerance > 0), 
                                                      df_bygene_withannotation$tdist2t, 
                                                      NA)
    EM_correction_N <- length(df_bygene_withannotation$tdist_positive[!is.na(df_bygene_withannotation$tdist_positive)])
    df_bygene_withannotation$tdist_positive_and_correctedEM <- ifelse((df_bygene_withannotation$tdist2t * EM_correction_N < threshold_alpha & df_bygene_withannotation$neg_zscore_minus_log2_integration_freq_withtolerance > 0), 
                                                                      df_bygene_withannotation$tdist2t * EM_correction_N, 
                                                                      NA)
    
    return (df_bygene_withannotation)
    
  } else {
    message(paste("[AP]\tERROR: Some problems with input files."))
  }  # if (ok_checks_passed)
  
}

###############################################################
#' @title MS Excel like T.DIST.2T
#' 
#' @author Andrea Calabria
#' @details version 0.1 (20 September 2019) 
#'
#' @rdname T.DIST.2T
#' @docType methods
#' @aliases T.DIST.2T
#'
#' @param x an input vector.
#' @param df degrees of freedom
#'
#' @return T.DIST.2T
#' @usage todo
#' @description ToDo
#' 
###############################################################

T.DIST.2T <- function(x, df) {
  return ( (1 - pt(x, df))*2 )
}


###############################################################
#' @title Annotate IS matrix
#' 
#' @author Andrea Calabria
#' @details version 0.1 (24-4-9)
#'
#' @rdname annotateISMatrix
#' @docType methods
#' @aliases annotateISMatrix
#'
#' @param df IS matrix
#' @param id_cols Annotation columns (che start strand) to use for annotation
#'
#' @return a new df. NB: it runs bedtools with system call!
#' @usage TODO
#' @note : TODO
#' @protype annotateISMatrix()
###############################################################
annotateISMatrix <- function(df, 
                             cols_chr_start_end_strand = c("chr", "integration_locus", "integration_locus", "integration_strand"),
                             id_cols = c("chr", "integration_locus", "integration_strand"),
                             features_to_annotate_gtf_bed,
                             feature_output_names = c("GeneName", "GeneStrand", "GeneDistance"),
                             df_bedfile_towrite = "ISmatrix.bed",
                             df_bedfile_annotated = "ISmatrix.annotated.bed",
                             bedtools_options = " -D ref -t first ",
                             ref_bed_format = F,
                             sort_bed = FALSE
                             ) {
  message(paste("[AP]\tAnnotate IS matrix using the annotation file\n\t", features_to_annotate_gtf_bed)) 
  # add rownames for next join
  df <- rownamesAsIS(df, id_cols = c("chr", "integration_locus", "integration_strand"), only_coordinates = T)
  
  # write tmp input file
  write.table(x = df[cols_chr_start_end_strand], 
              file = df_bedfile_towrite, 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = F, na = '')
  
  # sort and annotate 
  if (sort_bed) {
    system(command = paste0("bedtools sort -i '", df_bedfile_towrite, "' > '", paste0(df_bedfile_towrite, ".sorted.bed"), "'"))
    # run bedtools such this:
    system(command = paste0("bedtools closest -b '", 
                            features_to_annotate_gtf_bed, 
                            "' -a '", 
                            paste0(df_bedfile_towrite, ".sorted.bed"), 
                            "' ", bedtools_options, " > '", 
                            df_bedfile_annotated, "'"))
  } else {
    # run bedtools such this:
    system(command = paste0("bedtools closest -b '", 
                            features_to_annotate_gtf_bed, 
                            "' -a '", 
                            df_bedfile_towrite, 
                            "' ", bedtools_options, " > '", 
                            df_bedfile_annotated, "'")) 
  }
  
  if (ref_bed_format) {
    # if the annotation file ref is in BED file format, acquire a different number of fields
    # read the annotation of all the reads
    df_annotated <- read.csv(file = df_bedfile_annotated, 
                             header=F, fill=T, sep='\t', check.names = FALSE, 
                             na.strings = c("NONE", "NA", "NULL", "NaN", "ND", ""))
    names(df_annotated) <- c("chr", "integration_locus", "integration_locus_copy", "integration_locus_copy2", "integration_locus_copy3",
                             "integration_strand", "feature_chr", "feature_start", "feature_end", "FeatureName", "feature_score", "FeatureStrand", "FeatureDistance")
    rownames(df_annotated) <- apply(df_annotated[c("chr", "integration_locus", "integration_strand")], 1, function(x) {
      paste0(x[1], "_", as.character(as.numeric(x[2])), "_", x[3])})
  } else {
    # read the annotation of all the reads
    df_annotated <- read.csv(file = df_bedfile_annotated, 
                             header=F, fill=T, sep='\t', check.names = FALSE, 
                             na.strings = c("NONE", "NA", "NULL", "NaN", "ND", ""))
    names(df_annotated) <- c("chr", "integration_locus", "integration_locus_copy", "integration_locus_copy2", "integration_locus_copy3",
                             "integration_strand", "feature_chr", "feature_source", "feature_type", 
                             "feature_start", "feature_end", "feature_dot", "FeatureStrand", "feature_dot2", "feature_string", "FeatureDistance")
    rownames(df_annotated) <- apply(df_annotated[c("chr", "integration_locus", "integration_strand")], 1, function(x) {
      paste0(x[1], "_", as.character(as.numeric(x[2])), "_", x[3])})
    df_annotated$FeatureName <- apply(df_annotated[c("feature_string", "FeatureDistance")], 1, function(x) {
      strsplit( strsplit(x[1], ';', fixed = T)[[1]][1] , ' ', fixed = T)[[1]][2]
    } )
  }
  
  # return df annotated features
  if (sort_bed) {
    df_annotated_full <- merge(x = df, y = df_annotated[c(id_cols, "FeatureName", "FeatureStrand", "FeatureDistance")], by = id_cols, all.x = T)
    names(df_annotated_full)[(length(df_annotated_full)-2):(length(df_annotated_full))] <- feature_output_names
  } else{ 
    df_annotated_full <- cbind(df, df_annotated[c("FeatureName", "FeatureStrand", "FeatureDistance")])
    df_annotated_full <- df_annotated_full[c(id_cols, "FeatureName", "FeatureStrand", "FeatureDistance", setdiff(colnames(df), id_cols))]
    names(df_annotated_full) <- c(id_cols, feature_output_names, setdiff(colnames(df), id_cols) )  
  }
  return(df_annotated_full)
}


###############################################################
#' @title Assign an IS to a specific annotated feature
#' 
#' @author Andrea Calabria
#' @details version 0.1 (24-6-11)
#'
#' @rdname uniqueAnnotation
#' @docType methods
#' @aliases uniqueAnnotation
#'
#' @param elements vector of elements
#'
#' @return 
#' @usage TODO
#' @note : TODO
#' @protype uniqueAnnotation
###############################################################
uniqueAnnotation <- function(df, 
                             # annotation_prefix = c("Exon", "Gene", "TSS"),
                             annotation_cols_name_and_order = c("GeneName", "GeneStrand", "GeneDistance", "ExonName", "ExonStrand", "ExonDistance", "TSSName", "TSSStrand", "TSSDistance"),
                             id_cols = c("chr", "integration_locus", "integration_strand"),
                             annotation_suffix = c("Distance"),
                             tss_upstream_kb = 10000,
                             ...) {
  message(paste0("[AP]\tUnivocally label each IS with a feature.\n\n\t!!WARNING!! THE ORDER OF THE COLUMNS MUST BE: GeneName, GeneStrand, GeneDistance, ExonName, ExonStrand, ExonDistance, TSSName, TSSStrand, TSSDistance!!!\n"))
  df_sorted <- df[annotation_cols_name_and_order]
  # parse annotations
  # df_sorted$AnnotationClass <- apply(df[grep(paste(paste0(annotation_prefix, annotation_suffix), collapse = "|"), colnames(df), value = T)], 1, 
  # df_sorted$AnnotationClass <- apply(df[grep(paste(paste0(annotation_prefix, annotation_suffix), collapse = "|"), colnames(df), value = T)], 1, 
  df$AnnotationClass <- apply(df_sorted[grep(annotation_suffix, colnames(df_sorted), value = T)], 1, 
                              function(x) {
    if(x[2]==0) {
      "Exon"
    } else {
      # if((abs(x[2])>0) & (x[1]==0)) {
      if(x[1]==0) {
        "Intron"
      } else {
        # if( (abs(x[2])>0) & (abs(x[1])>0) & (x[3]>=(-tss_upstream_kb) & x[3]<=0) ) {
        if( x[3]>=(-tss_upstream_kb) & x[3]<=0 ) {
          paste0("TSS_", tss_upstream_kb)
        } else {
        "Intergenic"
        }
      }
    }
  }
  )
  return(df)
}


###############################################################
#' @title Create the list of comparisons
#' 
#' @author Andrea Calabria
#' @details version 0.1 (24-5-6)
#'
#' @rdname generateListOfPairs
#' @docType methods
#' @aliases generateListOfPairs
#'
#' @param elements vector of elements
#'
#' @return a new df. NB: it runs bedtools with system call!
#' @usage TODO
#' @note : TODO
#' @protype annotateISMatrix()
###############################################################
generateListOfPairs <- function(elements) {
  message(paste0("[AP]\tGenerate the list of pairs from: ", paste(elements, collapse = "~")))
  # outlist <- NULL
  outlist <- list()
  for (i in seq(1, (length(elements)-1))) {
    for (k in seq((i+1), length(elements))){
      # message(paste(i, k))
      # outlist <- list(outlist, c(elements[i], elements[i+k]))
      if (length(outlist) == 0) {
        # outlist <- data.frame("p1" = elements[i], "p2" = elements[k])
        outlist <- list(c(elements[i], elements[k]))
      } else {
        # outlist <- rbind(outlist, data.frame("p1" = elements[i], "p2" = elements[k]))
        outlist <- c(outlist, list(c(elements[i], elements[k])))
      }
    }
  }
  # return(list(outlist))
  return(outlist)
}

###############################################################
#' @title Load gene annotated as onco or tumor suppressors
#' 
#' @author Andrea Calabria
#' @details version 0.1 (27 September 2019) 
#'
#' @rdname loadOncoTSgenes
#' @docType methods
#' @aliases loadOncoTSgenes
#'
#' @param onco_db_file Path to UniProt oncogene database file (from config)
#' @param tumsup_db_file Path to UniProt tumor suppressor database file (from config)
#' @param species Default = "all" # alternatives: human, mouse
#' 
#'
#' @return dataframe
#' @usage todo
#' @description ToDo
#' 
###############################################################

loadOncoTSgenes <- function(onco_db_file = UNIPROT_ONCOGENE,
                            tumsup_db_file = UNIPROT_TUMOR_SUPPRESSOR,
                            species = "all" # alternatives: human, mouse
) {
  require(Hmisc)
  error <- ""
  ok_checks_passed <- TRUE
  if (!file.exists(onco_db_file) | !file.exists(tumsup_db_file)) {
    ok_checks_passed <- FALSE
    error <- paste(error, "\n[AP]\t\t\t-> Input files are not existing or wrong path.")
  }
  if (!(species %in% c("all", "human", "mouse"))) {
    ok_checks_passed <- FALSE
    error <- paste(error, "\n[AP]\t\t\t-> Specie not existing, choose one among: all, human, mouse")
  }
  
  if (ok_checks_passed) {
    message(paste("[AP]\tLoading annotated genes"))
    
    # acquire DB
    onco_db <- read.csv(file = onco_db_file, header=TRUE, fill=T, sep='\t', check.names = FALSE)
    tumsup_db <- read.csv(file = tumsup_db_file, header=TRUE, fill=T, sep='\t', check.names = FALSE)
    # onco_db <- read.csv("/Users/calabria.andrea/Dropbox (FONDAZIONE TELETHON)/Project Tumour Development/source/publicdb/uniprot-oncogene_Mouse.tab", header=TRUE, fill=T, sep='\t', check.names = FALSE)
    # tumsup_db <- read.csv("/Users/calabria.andrea/Dropbox (FONDAZIONE TELETHON)/Project Tumour Development/source/publicdb/uniprot-tumorsuppressor_Mouse.tab", header=TRUE, fill=T, sep='\t', check.names = FALSE)
    oncots_df_touse <- NULL # output df
    if (species == "mouse") {
      specie <- "Mus musculus (Mouse)"
      message(paste("[AP]\t-> Specie selected:", species, " -> output gene names will be Title case"))
      # subset data for this specie
      mouse_onco_db <- onco_db[which(onco_db$Organism == specie),]
      # mouse_onco_db <- onco_db
      # rownames(mouse_onco_db) <- mouse_onco_db$`Gene names  (primary )`
      mouse_tumsup_db <- tumsup_db[which(tumsup_db$Organism == specie),]
      # mouse_tumsup_db <- tumsup_db
      # rownames(mouse_tumsup_db) <- mouse_tumsup_db$`Gene names  (primary )`
      
      # get gene list: only reviewd genes and all gene aliases (from all species, thus for mice you must do lowe case and capitalize the gene name)
      mouse_onco_db_genes <- unique(as.character(mouse_onco_db$`Gene names  (primary )`))
      mouse_onco_db_genes_allnames <- unique(unlist(strsplit(as.character(unique(as.character(mouse_onco_db[which(mouse_onco_db$Status == "reviewed" & !is.na(mouse_onco_db$`Gene names`)), c("Gene names")]))), " ", fixed = TRUE), function(x) {c(x)}))
      mouse_onco_db_genes_allnames <- gsub(";", "", unique(capitalize(tolower(mouse_onco_db_genes_allnames)))) # merge datasets
      mouse_onco_db_genes_allnames_df <- data.frame("OncoGene" = mouse_onco_db_genes_allnames)
      rownames(mouse_onco_db_genes_allnames_df) <- mouse_onco_db_genes_allnames_df$OncoGene
      
      mouse_tumsup_db_genes <- unique(as.character(mouse_tumsup_db$`Gene names  (primary )`))
      mouse_tumsup_db_genes_allnames <- unique(unlist(strsplit(as.character(unique(as.character(mouse_tumsup_db[which(mouse_tumsup_db$Status == "reviewed" & !is.na(mouse_tumsup_db$`Gene names`)), c("Gene names")]))), " ", fixed = TRUE), function(x) {c(x)}))
      mouse_tumsup_db_genes_allnames <- gsub(";", "", unique(capitalize(tolower(mouse_tumsup_db_genes_allnames))) ) # merge datasets
      mouse_tumsup_db_genes_allnames_df <- data.frame("TumorSuppressor" = mouse_tumsup_db_genes_allnames)
      rownames(mouse_tumsup_db_genes_allnames_df) <- mouse_tumsup_db_genes_allnames_df$TumorSuppressor
      
      # merge df
      mouse_oncotumsup_db_genes_allnames_df <- merge(x = mouse_onco_db_genes_allnames_df, y = mouse_tumsup_db_genes_allnames_df, by = 0, all = T)
      rownames(mouse_oncotumsup_db_genes_allnames_df) <- mouse_oncotumsup_db_genes_allnames_df$Row.names
      names(mouse_oncotumsup_db_genes_allnames_df) <- c("GeneName", "OncoGene", "TumorSuppressor")
      mouse_oncotumsup_db_genes_allnames_df <- cbind(mouse_oncotumsup_db_genes_allnames_df, data.frame(
        "Onco1_TS2" = apply(mouse_oncotumsup_db_genes_allnames_df, 1, function(x) {
          .onco <- ifelse(!(is.na(x[2])), 1, 0)
          .tums <- ifelse(!(is.na(x[3])), 1, 0)
          ifelse( (.onco > 0 & .tums > 0 ), 3,
                  ifelse(.onco > 0, 1,
                         ifelse(.tums > 0, 2, NA)
                  )
          )
        })
      )
      )
      oncots_df_touse <- mouse_oncotumsup_db_genes_allnames_df
      return (oncots_df_touse)
    }
    
    if (species == "human") {
      specie <- "Homo sapiens (Human)"
      message(paste("[AP]\t-> Specie selected:", species, " -> output gene names will be Upper case."))
      # subset data for this specie
      human_onco_db <- onco_db[which(onco_db$Organism == specie),]
      # human_onco_db <- onco_db
      # rownames(human_onco_db) <- human_onco_db$`Gene names  (primary )`
      human_tumsup_db <- tumsup_db[which(tumsup_db$Organism == specie),]
      # human_tumsup_db <- tumsup_db
      # rownames(human_tumsup_db) <- human_tumsup_db$`Gene names  (primary )`
      
      # get gene list: only reviewd genes and all gene aliases (from all species, thus for mice you must do lowe case and capitalize the gene name)
      human_onco_db_genes <- unique(as.character(human_onco_db$`Gene names  (primary )`))
      human_onco_db_genes_allnames <- unique(unlist(strsplit(as.character(unique(as.character(human_onco_db[which(human_onco_db$Status == "reviewed" & !is.na(human_onco_db$`Gene names`)), c("Gene names")]))), " ", fixed = TRUE), function(x) {c(x)}))
      # human_onco_db_genes_allnames <- gsub(";", "", unique(capitalize(tolower(human_onco_db_genes_allnames)))) # merge datasets
      human_onco_db_genes_allnames <- gsub(";", "", unique(human_onco_db_genes_allnames)) # merge datasets
      human_onco_db_genes_allnames_df <- data.frame("OncoGene" = human_onco_db_genes_allnames)
      rownames(human_onco_db_genes_allnames_df) <- human_onco_db_genes_allnames_df$OncoGene
      
      human_tumsup_db_genes <- unique(as.character(human_tumsup_db$`Gene names  (primary )`))
      human_tumsup_db_genes_allnames <- unique(unlist(strsplit(as.character(unique(as.character(human_tumsup_db[which(human_tumsup_db$Status == "reviewed" & !is.na(human_tumsup_db$`Gene names`)), c("Gene names")]))), " ", fixed = TRUE), function(x) {c(x)}))
      # human_tumsup_db_genes_allnames <- gsub(";", "", unique(capitalize(tolower(human_tumsup_db_genes_allnames))) ) # merge datasets
      human_tumsup_db_genes_allnames <- gsub(";", "", unique(human_tumsup_db_genes_allnames)) # merge datasets
      human_tumsup_db_genes_allnames_df <- data.frame("TumorSuppressor" = human_tumsup_db_genes_allnames)
      rownames(human_tumsup_db_genes_allnames_df) <- human_tumsup_db_genes_allnames_df$TumorSuppressor
      
      # merge df
      human_oncotumsup_db_genes_allnames_df <- merge(x = human_onco_db_genes_allnames_df, y = human_tumsup_db_genes_allnames_df, by = 0, all = T)
      rownames(human_oncotumsup_db_genes_allnames_df) <- human_oncotumsup_db_genes_allnames_df$Row.names
      names(human_oncotumsup_db_genes_allnames_df) <- c("GeneName", "OncoGene", "TumorSuppressor")
      human_oncotumsup_db_genes_allnames_df <- cbind(human_oncotumsup_db_genes_allnames_df, data.frame(
        "Onco1_TS2" = apply(human_oncotumsup_db_genes_allnames_df, 1, function(x) {
          .onco <- ifelse(!(is.na(x[2])), 1, 0)
          .tums <- ifelse(!(is.na(x[3])), 1, 0)
          ifelse( (.onco > 0 & .tums > 0 ), 3, 
                  ifelse(.onco > 0, 1, 
                         ifelse(.tums > 0, 2, NA)
                  )
          )
        })
      )
      )
      oncots_df_touse <- human_oncotumsup_db_genes_allnames_df
      return (oncots_df_touse)
    }
    
    if (species == "all") {
      message(paste("[AP]\t-> Specie selected:", species, " -> output gene names will be Title case"))
      # subset data for this specie
      # allspecies_onco_db <- onco_db[which(onco_db$Organism == specie),]
      allspecies_onco_db <- onco_db
      # rownames(allspecies_onco_db) <- allspecies_onco_db$`Gene names  (primary )`
      # allspecies_tumsup_db <- tumsup_db[which(tumsup_db$Organism == specie),]
      allspecies_tumsup_db <- tumsup_db
      # rownames(allspecies_tumsup_db) <- allspecies_tumsup_db$`Gene names  (primary )`
      
      # get gene list: only reviewd genes and all gene aliases (from all species, thus for mice you must do lowe case and capitalize the gene name)
      allspecies_onco_db_genes <- unique(as.character(allspecies_onco_db$`Gene names  (primary )`))
      allspecies_onco_db_genes_allnames <- unique(unlist(strsplit(as.character(unique(as.character(allspecies_onco_db[which(allspecies_onco_db$Status == "reviewed" & !is.na(allspecies_onco_db$`Gene names`)), c("Gene names")]))), " ", fixed = TRUE), function(x) {c(x)}))
      allspecies_onco_db_genes_allnames <- gsub(";", "", unique(capitalize(tolower(allspecies_onco_db_genes_allnames)))) # merge datasets
      allspecies_onco_db_genes_allnames_df <- data.frame("OncoGene" = allspecies_onco_db_genes_allnames)
      rownames(allspecies_onco_db_genes_allnames_df) <- allspecies_onco_db_genes_allnames_df$OncoGene
      
      allspecies_tumsup_db_genes <- unique(as.character(allspecies_tumsup_db$`Gene names  (primary )`))
      allspecies_tumsup_db_genes_allnames <- unique(unlist(strsplit(as.character(unique(as.character(allspecies_tumsup_db[which(allspecies_tumsup_db$Status == "reviewed" & !is.na(allspecies_tumsup_db$`Gene names`)), c("Gene names")]))), " ", fixed = TRUE), function(x) {c(x)}))
      allspecies_tumsup_db_genes_allnames <- gsub(";", "", unique(capitalize(tolower(allspecies_tumsup_db_genes_allnames))) ) # merge datasets
      allspecies_tumsup_db_genes_allnames_df <- data.frame("TumorSuppressor" = allspecies_tumsup_db_genes_allnames)
      rownames(allspecies_tumsup_db_genes_allnames_df) <- allspecies_tumsup_db_genes_allnames_df$TumorSuppressor
      
      # merge df
      allspecies_oncotumsup_db_genes_allnames_df <- merge(x = allspecies_onco_db_genes_allnames_df, y = allspecies_tumsup_db_genes_allnames_df, by = 0, all = T)
      rownames(allspecies_oncotumsup_db_genes_allnames_df) <- allspecies_oncotumsup_db_genes_allnames_df$Row.names
      names(allspecies_oncotumsup_db_genes_allnames_df) <- c("GeneName", "OncoGene", "TumorSuppressor")
      allspecies_oncotumsup_db_genes_allnames_df <- cbind(allspecies_oncotumsup_db_genes_allnames_df, data.frame(
        "Onco1_TS2" = apply(allspecies_oncotumsup_db_genes_allnames_df, 1, function(x) {
          .onco <- ifelse(!(is.na(x[2])), 1, 0)
          .tums <- ifelse(!(is.na(x[3])), 1, 0)
          ifelse( (.onco > 0 & .tums > 0 ), 3, 
                  ifelse(.onco > 0, 1, 
                         ifelse(.tums > 0, 2, NA)
                  )
          )
        })
      )
      )
      
      oncots_df_touse <- allspecies_oncotumsup_db_genes_allnames_df
      return (oncots_df_touse)
      
    } # if all
    
  } else {
    message(paste("[AP]\tERROR: Some problems with input files:", error))
    return (NULL)
  }  # if (ok_checks_passed)
  
}


###############################################################
#' @title Load gene annotated as onco or tumor suppressors
#' 
#' @author Andrea Calabria
#' @details version 0.1 (14 October 2019) 
#'
#' @rdname loadOncoTSgenes_fromCancerMine
#' @docType methods
#' @aliases loadOncoTSgenes_fromCancerMine
#'
#' @param cancermine_collated_file Path to CancerMine database file (from config) 
#'
#' @return dataframe
#' @usage todo
#' @description ToDo
#' 
###############################################################

loadOncoTSgenes_fromCancerMine <- function(cancermine_collated_file = CANCERMINE_DB) {
  require(Hmisc)
  error <- ""
  ok_checks_passed <- TRUE
  if (!file.exists(cancermine_collated_file)) {
    ok_checks_passed <- FALSE
    error <- paste(error, "\n[AP]\t\t\t-> Input files are not existing or wrong path.")
  }
  
  if (ok_checks_passed) {
    message(paste("[AP]\tLoading annotated genes"))
    
    # acquire DB
    cancermine_db <- read.csv(file = cancermine_collated_file, header=TRUE, fill=T, sep='\t', check.names = FALSE)
    
    cancermine_db$GeneName <- cancermine_db$gene_normalized
    cancermine_db$OncoTS <- cancermine_db$role
    cancermine_db$OncoTS <- gsub("Oncogene", "OncoGene", gsub("_", "", cancermine_db$OncoTS))
    
    return (cancermine_db)
    
    
  } else {
    message(paste("[AP]\tERROR: Some problems with input files:", error))
    return (NULL)
  }  # if (ok_checks_passed)
  
}

