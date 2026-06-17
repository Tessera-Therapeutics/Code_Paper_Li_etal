#!/usr/bin/env python3
import pandas as pd 
import numpy as np
import scanpy as sc
import infercnvpy as cnv
from math import floor

"""
# Calculating Rates of Aneuploid Cells Per Treatment
In the Replogle et al paper (Cell, July 2022), the authors did a Perturb-seq experiment in which they characterized 
aneuploidy at a per-cell level in response to a pooled genome-wide CRISPRi screening library. They used infercnvpy 
as we did above, with the same parameters. To call cells with "karyotpic events," they used the following standard:
For each single-cell, infercnv generated a vector of CIN values across the genome. To label cells with likely 
karyotpic abnormalities, unstable karyotypic cells were heuristically defined as having >= 1 chromosome with evidence 
of changes in chromosomal copy number (nonzero CIN values) for >80% of the chromosomal length

Below, we apply the same standard to each timepoint, and summarize the rates of aneuploidy per treatment.
A key piece of information to be able to summarize across columns is knowing how the windows in the result matrix
map to their corresponding chromosomes. The method authors have stored the positions of where the chromosomes 
start relative to the windows represented in the result matrix (`anndata_obj.obsm['X_cnv']`) in `anndata_obj.uns["cnv"]`. 
"""

def determine_aneuploidy_for_chr_width(cnv_mat, 
                                        start, 
                                        stop, 
                                        abs_cnv_threshold = 0, 
                                        pct_aberrant_chr_threshold = 0.8):
    eighty_pct_chr_width = floor((stop - start) * pct_aberrant_chr_threshold) # threshold of num aberrand windows for this chr (based on how many windows this chr has)
    cnv_mat_subset = cnv_mat[:, start:stop] # extract windows for this chr
    cnv_events = abs(cnv_mat_subset) > abs_cnv_threshold # calculate num windows above CNV threshold [CIN value]
    event_sums = cnv_events.sum(axis = 1).tolist() # count how many windows in this chr are abnormal (one number per cell)
    flattened_event_sums = [event_sum[0] for event_sum in event_sums] # flatten to a list
    aneuploid_for_chr = [num_events >= eighty_pct_chr_width for num_events in flattened_event_sums] # call aneuploidy for each cell if threshold is met
    return(list(map(int, aneuploid_for_chr))) # convert booleans to 0/1

def determine_aneuploidy_for_chr_width_exclude(cnv_mat,
                                                start,
                                                stop,
                                                chr_name,
                                                exclude_chromosomes,
                                                abs_cnv_threshold=0,
                                                pct_aberrant_chr_threshold=0.8):
    if chr_name in exclude_chromosomes:
        return [0] * cnv_mat.shape[0]
    return determine_aneuploidy_for_chr_width(cnv_mat, start, stop, abs_cnv_threshold, pct_aberrant_chr_threshold)


def get_abnormal_cells(anndata_obj,
                        abs_cnv_threshold = 0,
                        pct_aberrant_chr_threshold = 0.8,
                        exclude_chromosomes = None):
    '''
    get chr positions from cnv data
    chr_pos = {
    'chr1': 0,
    'chr2': 120,
    'chr3': 210,
    ...
    }
    chr1 starts at col 0, chr2 at 120 etc

    exclude_chromosomes: optional list of chromosome names (e.g. ['chr6', 'chr15'])
                         to exclude from the aneuploidy calculation.
    '''
    if exclude_chromosomes is None:
        exclude_chromosomes = []

    chr_pos = anndata_obj.uns["cnv"]['chr_pos']
    chr_names = list(chr_pos.keys())
    chr_widths = {}
    num_infercnv_res_columns = anndata_obj.obsm['X_cnv'].shape[1] # total number of genomic windows
    # Build chromosome start-end table
    for idx in range(0, len(chr_pos)):
        current_chr = chr_names[idx]
        current_chr_width = {}
        if idx == len(chr_names) - 1:
            current_chr_width = {'start': chr_pos[current_chr], 'end': num_infercnv_res_columns}
        else:
            next_chr = chr_names[idx + 1]
            current_chr_width = {'start': chr_pos[current_chr], 'end': chr_pos[next_chr]}
        chr_widths[current_chr] = current_chr_width
    # convert to dataframe and compute widths
    chr_widths = pd.DataFrame.from_dict(chr_widths).transpose().reset_index()
    chr_widths.rename(columns={'index': 'chr_name'}, inplace=True)
    chr_widths['width'] = chr_widths['end'] - chr_widths['start']
    chr_widths = chr_widths[chr_widths['width'] > 1]  # In some datasets, infercnvpy only calculates one window for a chr. Ignore these

    cnv_mat = anndata_obj.obsm['X_cnv'] # rows = cells; columns = genomic windows
    # run aneuploidy calculation for each chromosome in chr_widths (it will do all cells automatically)
    if exclude_chromosomes:
        per_chr_aneuploid_cells = chr_widths.apply(
            lambda row: determine_aneuploidy_for_chr_width_exclude(
                cnv_mat, row['start'], row['end'], row['chr_name'],
                exclude_chromosomes, abs_cnv_threshold, pct_aberrant_chr_threshold), axis=1)
    else:
        per_chr_aneuploid_cells = chr_widths.apply(
            lambda row: determine_aneuploidy_for_chr_width(
                cnv_mat, row['start'], row['end'], abs_cnv_threshold, pct_aberrant_chr_threshold), axis=1)
    # put aneuploidy calculation in a dataframe with rows = cells, cols = chromosomes
    per_chr_aneuploidy_df = pd.DataFrame({name: arr for name, arr in zip(chr_widths['chr_name'], per_chr_aneuploid_cells)})
    # count num abnormal chr per cell
    abnormal_cell_df = pd.DataFrame(per_chr_aneuploidy_df.sum(axis=1),
                                    index=list(range(0, per_chr_aneuploidy_df.shape[0], 1)),
                                    columns=['num_abnormal_chr'])
    # prepare cell metadata as anndata obj
    anndata_abnormal = anndata_obj.obs[['cell_barcode', 'Timepoint', 'sample_id', 'Treatment','Run']].copy()
    # append num abnormal chr to anndata object as cell property
    anndata_abnormal['num_abnormal_chr'] = abnormal_cell_df['num_abnormal_chr'].tolist()
    # call each cell normal or abnormal
    anndata_abnormal['has_karyotypic_abnormalities'] = anndata_abnormal['num_abnormal_chr'] > 0
    # aggregate data by timepoint, treatment, AND sample_id (replicate)
    summary_df = anndata_abnormal.groupby(['Timepoint', 'Treatment', 'Run', 'sample_id'], observed=True).agg(
        num_karyotypically_abnormal=('has_karyotypic_abnormalities', 'sum'),
        num_cells_in_treatment=('Treatment', 'count')  # you could use any column here, it is just counting the non-null vals
    )
    # compute pct abnormal
    summary_df = summary_df.assign(pct_karyotypically_abnormal_cells=summary_df['num_karyotypically_abnormal'] / summary_df['num_cells_in_treatment'] * 100)

    if exclude_chromosomes:
        print(f"  NOTE: Excluded chromosomes from aneuploidy calculation: {exclude_chromosomes}")

    return summary_df, anndata_abnormal


