#!/usr/bin/env python3

import pandas as pd
import numpy as np
import scanpy as sc
import infercnvpy as cnv
import pickle
import os
from aneuploid_rate_calculation_byrep import get_abnormal_cells

"""
Run infercnv per timepoint group so that each group's cases are compared against
matched controls. This avoids confounding from timepoint-specific batch effects.

Each entry in TIMEPOINT_GROUPS defines one infercnv run:
  - case_timepoints: which timepoints' case cells to include
  - ref_timepoints:  which timepoints' reference cells to use as the normal baseline

day10 has no Activated-Tcells controls, so it is grouped with day8 and they
share the day8 controls in a single infercnv run.
"""

# ========================= USER SETTINGS =========================
# UPDATE: Set these paths to your output directories

ANNDATA_FILE = "full_run_preprocessed.pickle"
indir = "results/02_prep_infercnv"  # From step 01
outdir = "results/03_infercnv_timepoints/"

# Reference treatment used as the normal baseline
REFERENCE_TREATMENT = "Activated-Tcells"

# Case treatments to evaluate (set to None to include everything that isn't reference)
CASE_TREATMENTS = [
    "Etop 500 nM",
    "Etop 250nM",
    "Nucleofection-only",
    "Vingi EN- RT- 100ng",
    "Vingi WT 100ng",
    "Cas9-B2M-TRAC-DKO",
]

# Each entry defines one infercnv run:
#   case_timepoints  – timepoints whose case cells go into the run
#   ref_timepoints   – timepoints from which to draw REFERENCE_TREATMENT cells
# Day10 has no controls, so it shares a run with day8 using day8 controls.

TIMEPOINT_GROUPS = [
    {"case_timepoints": ["day3"],          "ref_timepoints": ["day3"]},
    {"case_timepoints": ["day4"],          "ref_timepoints": ["day4"]},
    {"case_timepoints": ["day8","day10"],  "ref_timepoints": ["day8"]},
]

# infercnv parameters
INFERCNV_PARAMS = dict(
    window_size=100,
    step=1,
    chunksize=5000,
    n_jobs=64,
)

# Aneuploidy calling threshold (passed to get_abnormal_cells)
ABS_CNV_THRESHOLD = 0

# ========================= END SETTINGS ==========================

os.makedirs(outdir, exist_ok=True)

# ---- Load data ----
with open(os.path.join(indir, ANNDATA_FILE), "rb") as f:
    full_run = pickle.load(f)
print(f"Loaded anndata: {full_run.shape[0]} cells x {full_run.shape[1]} genes")

available_timepoints = full_run.obs["Timepoint"].unique().tolist()
available_treatments = full_run.obs["Treatment"].unique().tolist()
print(f"Available timepoints: {sorted(available_timepoints)}")
print(f"Available treatments: {sorted(available_treatments)}")

# ---- Run per-timepoint group ----
all_summaries = []

for group in TIMEPOINT_GROUPS:
    case_tps = group["case_timepoints"]
    ref_tps = group["ref_timepoints"]
    group_label = "+".join(case_tps)

    print(f"\n{'='*60}")
    print(f"Processing case timepoints: {case_tps}  (reference controls from: {ref_tps})")
    print(f"{'='*60}")

    # Select reference cells from the specified timepoints
    ref_mask = (
        (full_run.obs["Treatment"] == REFERENCE_TREATMENT)
        & (full_run.obs["Timepoint"].isin(ref_tps))
    )

    # Select case cells at the case timepoints
    if CASE_TREATMENTS is not None:
        case_mask = (
            (full_run.obs["Treatment"].isin(CASE_TREATMENTS))
            & (full_run.obs["Timepoint"].isin(case_tps))
        )
    else:
        case_mask = (
            (full_run.obs["Treatment"] != REFERENCE_TREATMENT)
            & (full_run.obs["Timepoint"].isin(case_tps))
        )

    n_ref = ref_mask.sum()
    n_case = case_mask.sum()
    print(f"  Reference cells: {n_ref}  |  Case cells: {n_case}")

    if n_ref == 0:
        print(f"  WARNING: no reference cells found — skipping {group_label}")
        continue
    if n_case == 0:
        print(f"  WARNING: no case cells found — skipping {group_label}")
        continue

    subset = full_run[ref_mask | case_mask].copy()
    print(f"  Total cells in subset: {subset.shape[0]}")
    print(f"  Treatments: {sorted(subset.obs['Treatment'].unique())}")
    print(f"  Sample IDs: {sorted(subset.obs['sample_id'].unique())}")

    cnv.tl.infercnv(
        subset,
        reference_key="Treatment",
        reference_cat=REFERENCE_TREATMENT,
        **INFERCNV_PARAMS,
    )

    aneuploidy_summary, anndata_abnormal = get_abnormal_cells(
        subset, abs_cnv_threshold=ABS_CNV_THRESHOLD
    )
    print(f"\n  Aneuploidy summary for {group_label}:")
    print(aneuploidy_summary)
    aneuploidy_summary = aneuploidy_summary.reset_index()
    aneuploidy_summary.to_csv(os.path.join(outdir, f"cases_vs_controls_{group_label}.csv"), index=False)
    
    if "Cas9-B2M-TRAC-DKO" in CASE_TREATMENTS:
        dko = subset[subset.obs["Treatment"] == "Cas9-B2M-TRAC-DKO"]
        aneuploidy_summary_b2m_ko, anndata_abnormal_b2m_ko = get_abnormal_cells(dko, abs_cnv_threshold=ABS_CNV_THRESHOLD, exclude_chromosomes=["chr14"])
        aneuploidy_summary_b2m_ko = aneuploidy_summary_b2m_ko.reset_index()
        aneuploidy_summary_b2m_ko["Treatment"] = aneuploidy_summary_b2m_ko["Treatment"].str.replace("Cas9-B2M-TRAC-DKO", "Cas9-B2M-KO")
        aneuploidy_summary_b2m_ko.to_csv(os.path.join(outdir, f"cases_vs_controls_{group_label}_Cas9_B2M_KO.csv"), index=False)

        aneuploidy_summary_trac_ko, anndata_abnormal_trac_ko = get_abnormal_cells(dko, abs_cnv_threshold=ABS_CNV_THRESHOLD, exclude_chromosomes=["chr15"])
        aneuploidy_summary_trac_ko = aneuploidy_summary_trac_ko.reset_index()
        aneuploidy_summary_trac_ko["Treatment"] = aneuploidy_summary_trac_ko["Treatment"].str.replace("Cas9-B2M-TRAC-DKO", "Cas9-TRAC-KO")
        aneuploidy_summary_trac_ko.to_csv(os.path.join(outdir, f"cases_vs_controls_{group_label}_Cas9_TRAC_KO.csv"), index=False)

        aneuploidy_summary_no_target, anndata_abnormal_no_target = get_abnormal_cells(dko, abs_cnv_threshold=ABS_CNV_THRESHOLD, exclude_chromosomes=["chr14", "chr15"])
        aneuploidy_summary_no_target = aneuploidy_summary_no_target.reset_index()
        aneuploidy_summary_no_target["Treatment"] = aneuploidy_summary_no_target["Treatment"].str.replace("Cas9-B2M-TRAC-DKO", "Cas9-NO-TARGET")
        aneuploidy_summary_no_target.to_csv(os.path.join(outdir, f"cases_vs_controls_{group_label}_Cas9_NO-TARGET.csv"), index=False)
    
    all_summaries.append(aneuploidy_summary)
    all_summaries.append(aneuploidy_summary_b2m_ko)
    all_summaries.append(aneuploidy_summary_trac_ko)
    all_summaries.append(aneuploidy_summary_no_target)

# ---- Combine all timepoints into one CSV ----
if all_summaries:
    combined = pd.concat(all_summaries, ignore_index=True)
    combined.to_csv(os.path.join(outdir, "cases_vs_untreated.csv"), index=False)
    print(f"\nCombined results: {len(combined)} rows")
    print(f"Saved to: {outdir}/cases_vs_untreated.csv")

print("\nDone!")
