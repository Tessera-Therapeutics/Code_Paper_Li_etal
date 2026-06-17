# Single-cell RNA-seq Analysis

Analysis pipeline for detecting karyotypic abnormalities (aneuploidy) in single-cell RNA-seq data using InferCNV.

## Overview

This pipeline processes 10X Genomics single-cell RNA-seq data through three main steps:

1. **Data preparation and QC** - Load cellranger outputs, filter cells/genes, prepare for CNV analysis
2. **InferCNV analysis by timepoint** - Detect copy number variations using infercnvpy
3. **Statistical modeling** - Non-linear mixed models to analyze aneuploidy trajectories over time

## Pipeline Steps

### Step 1: Prepare InferCNV from Counts

**Script**: `01_prep_infercnv_from_counts.ipynb` (Jupyter notebook)

Loads individual cellranger count outputs, filters cells and genes by quality metrics, and prepares a single AnnData object for CNV inference.

**Inputs:**
- Sample metadata: `../metadata/sample_info/scRNAseq.Samples.csv`
- Cellranger outputs: `filtered_feature_bc_matrix.h5` files (paths in metadata CSV)
- Gene annotations: GTF file (hg38)

**QC Filters:**
- Min 200 genes per cell
- Min 1000 counts per cell
- <10% mitochondrial content
- Genes must be in ≥3% of cells
- Canonical chromosomes only (chr1-22)

**Output:**
- `results/02_prep_infercnv/full_run_preprocessed.pickle` - Filtered and normalized AnnData object

### Step 2: Run InferCNV by Timepoint

**Script**: `02_run_infercnv_by_timepoint.py`

Runs infercnvpy separately for each timepoint group, comparing treated cells against matched control timepoints to avoid batch effects.

**Key features:**
- Groups timepoints with matched controls (day3, day4, day8+day10)
- Reference treatment: Activated-Tcells
- Case treatments: Etop, Nucleofection, Vingi, Cas9-B2M-TRAC-DKO
- Calls aneuploidy using Replogle et al. (Cell, 2022) criteria:
  - ≥1 chromosome with CNV changes in >80% of chromosomal length
- Special handling for Cas9-B2M-TRAC-DKO: reports with/without target chromosomes excluded

**InferCNV parameters:**
```python
window_size = 100
step = 1
chunksize = 5000
n_jobs = 64
```

**Output:**
- `results/03_infercnv_timepoints/cases_vs_controls_*.csv` - Per-timepoint aneuploidy rates
- `results/03_infercnv_timepoints/cases_vs_untreated.csv` - Combined results

### Step 3: Non-Linear Mixed Model Analysis

**Script**: `03_nlmm_aneuploidy_baseline_corrected.py`

Fits non-linear mixed models to baseline-corrected aneuploidy trajectories.

**Model structure:**
```
Y_ij = β₀ + β₁·day + β₂·exp(-day)
     + β₃·treatment + β₄·day×treatment + β₅·exp(-day)×treatment
     + u_i + ε_ij
```

**Baseline correction:**
- Etop, Cas9DKO, Nucleofection → subtract mean(Activated-Tcells) per day
- Vingi WT → subtract mean(Vingi EN- RT-) per day

**Statistical tests:**
- Trajectory comparison (linear + exponential terms)
- Slope vs. zero tests
- Treatment comparisons

**Output:**
- `results/05_nlmm_timepoints_baseline_corrected/` - Model results, plots, statistics

### Utility Module

**Script**: `aneuploid_rate_calculation_byrep.py`

Reusable functions for calculating aneuploidy rates per treatment/replicate:

- `get_abnormal_cells()` - Applies Replogle criteria to call aneuploid cells
- `determine_aneuploidy_for_chr_width()` - Per-chromosome aneuploidy detection
- Supports excluding specific chromosomes (for Cas9 target analysis)

## Requirements

### Python packages

```bash
pip install scanpy infercnvpy pandas numpy scipy matplotlib seaborn statsmodels plotnine kneed
```

### Key dependencies

- **scanpy** (≥1.9) - Single-cell analysis
- **infercnvpy** (≥0.4) - Copy number variation inference
- **statsmodels** - Mixed linear models
- **Jupyter** - For running notebooks

## Setup

1. **Prepare sample metadata**

   Edit `../metadata/sample_info/scRNAseq.Samples.csv`:
   - Set `Cellranger_counts_h5` paths to your cellranger output locations
   - Paths should point to `molecule_info.h5` files (code auto-replaces with `filtered_feature_bc_matrix.h5`)

2. **Update paths in notebook**

   In `01_prep_infercnv_from_counts.ipynb`:
   ```python
   GENOME_GTF_PATH = Path("path/to/your/genes.gtf.gz")
   ```

3. **Run pipeline in order**

   ```bash
   # Step 1: Prepare data (Jupyter notebook)
   jupyter notebook 01_prep_infercnv_from_counts.ipynb
   
   # Step 2: Run InferCNV
   python 02_run_infercnv_by_timepoint.py
   
   # Step 3: Statistical analysis
   python 03_nlmm_aneuploidy_baseline_corrected.py
   ```

## Methods Reference

### Aneuploidy Calling Criteria

Following Replogle et al. (Cell, 2022):

> Unstable karyotypic cells were heuristically defined as having ≥1 chromosome with evidence of changes in chromosomal copy number (nonzero CIN values) for >80% of the chromosomal length.

For each cell:
1. InferCNV generates a vector of CIN values across the genome
2. For each chromosome, count windows with |CNV| > threshold (default: 0)
3. If ≥80% of chromosome windows are aberrant → chromosome is aneuploid
4. If ≥1 chromosome is aneuploid → cell has karyotypic abnormalities

## Configuration

Key parameters can be adjusted in each script:

**02_run_infercnv_by_timepoint.py:**
- `REFERENCE_TREATMENT` - Control treatment (default: "Activated-Tcells")
- `CASE_TREATMENTS` - Treatments to evaluate
- `TIMEPOINT_GROUPS` - How to group timepoints for analysis
- `INFERCNV_PARAMS` - Window size, step, parallelization
- `ABS_CNV_THRESHOLD` - Minimum CNV value to call aberrant (default: 0)

**03_nlmm_aneuploidy_baseline_corrected.py:**
- `BASELINE_MAP` - Which control to use for each treatment
- `TRAJECTORY_TEST_PARAMS` - Which model terms to test
- `COMPARISONS` - Treatment pairs to compare

## Output Structure

```
results/
├── 02_prep_infercnv/
│   ├── full_run_preprocessed.pickle
│   └── preprocessing_filter_stats.csv
├── 03_infercnv_timepoints/
│   ├── cases_vs_controls_day3.csv
│   ├── cases_vs_controls_day4.csv
│   ├── cases_vs_controls_day8+day10.csv
│   └── cases_vs_untreated.csv
└── 05_nlmm_timepoints_baseline_corrected/
    ├── model_results.csv
    ├── trajectory_plots.pdf
    └── statistics.txt
```

## Notes

- **Day 10 handling**: Day10 has no Activated-Tcells controls, so it's grouped with day8 and shares day8 controls
- **Cas9-B2M-TRAC-DKO**: Special handling excludes chr14 (B2M), chr15 (TRAC), or both to assess off-target aneuploidy
- **Computational requirements**: InferCNV can be memory-intensive; recommend ≥64GB RAM for large datasets
- **Runtime**: Full pipeline takes ~2-4 hours depending on dataset size and CPU cores

## Citation

If you use this code, please cite:

```
Li et al. (2026). "In vivo Gene Writing with engineered retrotransposons."
[Journal]. DOI: [to be added]
```

And the InferCNV method:

```
Replogle et al. (2022). "Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq."
Cell 185(14): 2559-2575.e28. DOI: 10.1016/j.cell.2022.05.013
```