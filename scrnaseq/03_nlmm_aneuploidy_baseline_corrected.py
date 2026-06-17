#!/usr/bin/env python3
"""
Non-Linear Mixed Model Analysis for Aneuploidy Trajectories
(Baseline-Corrected)
------------------------------------------------------------
Each treatment's replicates are corrected by subtracting the mean of their
respective control at each timepoint:
  - Etop 500 nM, Etop 250nM, Cas9DKO, Nucleofection  → subtract mean(Activated-Tcells) per day
  - Vingi WT 100ng           → subtract mean(Vingi EN- RT- 100ng) per day

Model structure on corrected data:
    Y_ij = β₀ + β₁·day + β₂·exp(-day)
         + β₃·treatment + β₄·day×treatment + β₅·exp(-day)×treatment
         + u_i + ε_ij

"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.regression.mixed_linear_model import MixedLM
import os


# ========================= USER SETTINGS =========================
# UPDATE: Set these paths to your input/output directories

TIMEPOINT = "_timepoints"
indir = f"results/03_infercnv{TIMEPOINT}/"  # From step 02
input_csv = os.path.join(indir, "cases_vs_untreated.csv")
outdir = f"results/05_nlmm{TIMEPOINT}_baseline_corrected"
os.makedirs(outdir, exist_ok=True)
LEGEND = True

# ========================= END SETTINGS ==========================

# ============================================================================
# CONFIGURATION CAS9 PLOT WITH NUCLEOFECTION 
# ============================================================================

TRAJECTORY_TEST_PARAMS = "both"  # "linear_only", "exponential_only", or "both"
SLOPE_VS_ZERO_PARAMS = {"Nucleofection-only": "both"}
SLOPE_VS_ZERO_DEFAULT = "both"
SLOPE_VS_ZERO_SKIP = ["Cas9-TRAC-KO", "Cas9-B2M-KO", "Cas9-B2M-TRAC-DKO"]
SLOPE_VS_ZERO_ANNOTATE = ["Nucleofection-only"]
COMPARISONS = []

BASELINE_MAP = {
    "Cas9-B2M-TRAC-DKO": "Activated-Tcells", # we keep activated T cells here instead of nucleofection; we want to see the nucleofection effect and compare with Cas9 samples after removing target chromosomes
    "Nucleofection-only": "Activated-Tcells",
    "Cas9-B2M-KO": "Activated-Tcells",
    "Cas9-TRAC-KO": "Activated-Tcells",
    "Cas9-NO-TARGET": "Activated-Tcells",}

SAMPLE_FILTERS = {}
NAMEPIC = "cas9"

###
# ============================================================================
# CONFIGURATION VINGI PLOT
# ============================================================================

# TRAJECTORY_TEST_PARAMS = "both"  # "linear_only", "exponential_only", or "both"
# SLOPE_VS_ZERO_PARAMS = {"Vingi WT 100ng": "both",}
# SLOPE_VS_ZERO_DEFAULT = "both"
# SLOPE_VS_ZERO_SKIP = []
# SLOPE_VS_ZERO_ANNOTATE = ["Vingi WT 100ng"]

# BASELINE_MAP = {
#     "Etop 500 nM": "Activated-Tcells",
#     "Etop 250nM": "Activated-Tcells",
#     "Vingi WT 100ng": "Vingi EN- RT- 100ng",
# }

# SAMPLE_FILTERS = {}

# COMPARISONS = [
#     ("Vingi WT 100ng", "Etop 250nM"),
#     ("Vingi WT 100ng", "Etop 500 nM"),
#     # ("Vingi WT 100ng", "Cas9-B2M-TRAC-DKO"),
# ]
# NAMEPIC = "vingi"

# ============================================================================
# STEP 1: Load and prepare the data
# ============================================================================

print("=" * 70)
print("DATA LOADING AND PREPARATION")
print("=" * 70)

df_all = pd.read_csv(input_csv)
df_all = df_all.drop_duplicates(
    subset=["sample_id", "Timepoint", "Treatment", "pct_karyotypically_abnormal_cells"]
)

all_treatments = list(BASELINE_MAP.keys()) + list(set(BASELINE_MAP.values()))
df = df_all[df_all["Treatment"].isin(all_treatments)].copy()

for trt, filt_fn in SAMPLE_FILTERS.items():
    mask = df["Treatment"] == trt
    df = pd.concat([df[~mask], filt_fn(df[mask])], ignore_index=True)

response_col = "pct_karyotypically_abnormal_cells"

# ============================================================================
# STEP 2: Extract and create variables
# ============================================================================

df["day_num"] = (
    df["Timepoint"]
    .astype(str)
    .str.extract(r"(\d+)", expand=False)
    .astype(int)
)

def extract_replicate(row):
    import re
    sid = str(row["sample_id"])
    m = re.search(r"(rep\d+)", sid)
    if m:
        return m.group(1)
    m = re.search(r"D\d+_(D\d+)_", sid)
    if m:
        return m.group(1)
    return sid

df["replicate"] = df.apply(extract_replicate, axis=1)
df["replicate_id"] = df["Treatment"].astype(str) + "_" + df["replicate"].astype(str)
df["Treatment"] = df["Treatment"].str.strip()

print(f"\nLoaded {len(df)} observations from {input_csv}")
print(f"Treatments: {df['Treatment'].unique().tolist()}")
print(f"Timepoints (days): {sorted(df['day_num'].unique())}")

# ============================================================================
# STEP 3: Baseline correction
# ============================================================================

print("\n" + "=" * 70)
print("BASELINE CORRECTION")
print("=" * 70)

act_day8 = df[(df["Treatment"] == "Activated-Tcells") & (df["day_num"] == 8)].copy()
act_day8["day_num"] = 10
act_day8["Timepoint"] = "day10"
df = pd.concat([df, act_day8], ignore_index=True)

baseline_means = {}
baseline_sems = {}
for baseline_trt in set(BASELINE_MAP.values()):
    bl_data = df[df["Treatment"] == baseline_trt]
    bl_means = bl_data.groupby("day_num")[response_col].mean()
    bl_sems = bl_data.groupby("day_num")[response_col].sem()
    if baseline_trt == "Vingi EN- RT- 100ng" and 8 not in bl_means.index and 10 in bl_means.index:
        bl_means[8] = bl_means[10]
        bl_sems[8] = bl_sems[10]
        print(f"\n  NOTE: Using day10 as proxy for day8 in {baseline_trt}")
    baseline_means[baseline_trt] = bl_means
    baseline_sems[baseline_trt] = bl_sems
    print(f"\n  Baseline means for {baseline_trt}:")
    for day, val in sorted(bl_means.items()):
        print(f"    day {day}: {val:.3f}% (SEM: {bl_sems[day]:.3f})")

corrected_rows = []
for trt, baseline_trt in BASELINE_MAP.items():
    trt_data = df[df["Treatment"] == trt].copy()
    bl_means = baseline_means[baseline_trt]
    for _, row in trt_data.iterrows():
        day = row["day_num"]
        if day in bl_means.index:
            corrected_val = row[response_col] - bl_means[day]
            new_row = row.copy()
            new_row[response_col] = corrected_val
            corrected_rows.append(new_row)

df_corr = pd.DataFrame(corrected_rows)

treatments = list(BASELINE_MAP.keys())
df_corr["Treatment"] = pd.Categorical(df_corr["Treatment"], categories=treatments, ordered=True)

print(f"\nAfter correction: {len(df_corr)} observations")
print(f"Treatments in corrected data: {df_corr['Treatment'].cat.categories.tolist()}")
print(f"\nCorrected data summary per treatment and day:")
print(df_corr.groupby(["Treatment", "day_num"])[response_col].agg(["mean", "std", "count"]))

# ============================================================================
# STEP 4: Create design matrix for fixed effects
# ============================================================================

print("\n" + "=" * 70)
print("BUILDING DESIGN MATRIX")
print("=" * 70)

REFERENCE_TREATMENT = treatments[0]
treatment_dummies = pd.get_dummies(df_corr["Treatment"], drop_first=True, prefix="Treatment")

print(f"\nReference treatment (first in model): {REFERENCE_TREATMENT}")

X_fixed = pd.DataFrame({
    'Intercept': np.ones(len(df_corr), dtype=float),
    'day_num': df_corr['day_num'].values.astype(float),
    'day_num_exp': np.exp(-df_corr['day_num'].values).astype(float),
})

for col in treatment_dummies.columns:
    X_fixed[col] = treatment_dummies[col].values.astype(float)

for col in treatment_dummies.columns:
    interaction_name = f"day_num:{col}"
    X_fixed[interaction_name] = (df_corr['day_num'].values * treatment_dummies[col].values).astype(float)

for col in treatment_dummies.columns:
    interaction_name = f"day_num_exp:{col}"
    X_fixed[interaction_name] = (np.exp(-df_corr['day_num'].values) * treatment_dummies[col].values).astype(float)

print(f"\nFinal design matrix shape: {X_fixed.shape}")
print(f"Predictors in the model:")
for i, col in enumerate(X_fixed.columns):
    print(f"  {i+1}. {col}")

# ============================================================================
# STEP 5: Fit the non-linear mixed model
# ============================================================================

print("\n" + "=" * 70)
print("FITTING NON-LINEAR MIXED MODEL")
print("=" * 70)

y = df_corr[response_col].values
groups = df_corr["replicate_id"]

print(f"\nModel specification:")
print(f"  Outcome: {response_col} (baseline-corrected)")
print(f"  Fixed effects: {', '.join(X_fixed.columns)}")
print(f"  Random effects: Random intercept for each replicate_id")
print(f"  Number of groups (replicates): {groups.nunique()}")
print(f"  Number of observations: {len(y)}")

model = MixedLM(endog=y, exog=X_fixed, groups=groups)

try:
    result = model.fit(method="powell", reml=True, maxiter=500)
except Exception as e:
    print(f"Powell method failed: {e}")
    print("Trying with BFGS method instead...")
    result = model.fit(method="bfgs", reml=True, maxiter=500)

print("\n" + "=" * 70)
print("MODEL RESULTS")
print("=" * 70)
print(result.summary())

# ============================================================================
# STEP 6: Extract and interpret parameter estimates
# ============================================================================

print("\n" + "=" * 70)
print("PARAMETER INTERPRETATION")
print("=" * 70)

fe_params = result.fe_params
fe_stderr = result.bse_fe
fe_zscores = fe_params / fe_stderr
fe_pvalues = 2 * (1 - stats.norm.cdf(np.abs(fe_zscores)))

params_table = pd.DataFrame({
    'Parameter': list(fe_params.index),
    'Estimate': fe_params.values,
    'Std_Error': fe_stderr.values,
    'z_statistic': fe_zscores.values,
    'P_value': fe_pvalues,
})

print("\nFixed effects parameters:")
print(params_table.to_string(index=False))

random_effects_var = result.cov_re.values[0, 0]
residual_var = result.scale
icc = random_effects_var / (random_effects_var + residual_var)

print(f"\n--- Variance components ---")
print(f"  Random intercept variance (between replicates): {random_effects_var:.3f}")
print(f"  Residual variance (within replicates): {residual_var:.3f}")
print(f"  Intraclass correlation (ICC): {icc:.3f}")
print(f"    (i.e., {icc*100:.1f}% of variance is between replicates)")

# ============================================================================
# STEP 7: Covariance matrix for predictions
# ============================================================================

cov_matrix = result.cov_params().iloc[:len(fe_params), :len(fe_params)]

# ============================================================================
# STEP 8: Trajectory shape comparisons (Wald test)
# ============================================================================

print("\n" + "=" * 70)
print("TRAJECTORY SHAPE COMPARISONS BETWEEN TREATMENTS")
print("=" * 70)
if TRAJECTORY_TEST_PARAMS == "linear_only":
    print("\nTesting if the LINEAR component of trajectory differs between treatments")
elif TRAJECTORY_TEST_PARAMS == "exponential_only":
    print("\nTesting if the EXPONENTIAL component of trajectory differs between treatments")
else:
    print("\nTesting if the TRAJECTORY SHAPE differs between treatments")
    print("by jointly testing trajectory-related parameters:")
    print("  - Treatment effect on linear component (day_num:Treatment)")
    print("  - Treatment effect on exponential component (day_num_exp:Treatment)")
print("\nNote: We exclude intercept differences to focus on the DYNAMICS of change.")


def compare_trajectories_wald(treatment1, treatment2, fe_params, cov_matrix,
                              reference_treatment, test_params):
    param_names = list(fe_params.index)
    n_params = len(param_names)
    contrasts = []
    param_labels = []

    if treatment1 == reference_treatment and treatment2 != reference_treatment:
        target = treatment2
        sign = 1.0
    elif treatment2 == reference_treatment and treatment1 != reference_treatment:
        target = treatment1
        sign = 1.0
    else:
        target = None

    if target is not None:
        if test_params in ["linear_only", "both"]:
            name = f"day_num:Treatment_{target}"
            if name in param_names:
                L = np.zeros(n_params)
                L[param_names.index(name)] = sign
                contrasts.append(L)
                param_labels.append(name)
        if test_params in ["exponential_only", "both"]:
            name = f"day_num_exp:Treatment_{target}"
            if name in param_names:
                L = np.zeros(n_params)
                L[param_names.index(name)] = sign
                contrasts.append(L)
                param_labels.append(name)
    else:
        if test_params in ["linear_only", "both"]:
            n1 = f"day_num:Treatment_{treatment1}"
            n2 = f"day_num:Treatment_{treatment2}"
            if n1 in param_names and n2 in param_names:
                L = np.zeros(n_params)
                L[param_names.index(n1)] = 1.0
                L[param_names.index(n2)] = -1.0
                contrasts.append(L)
                param_labels.append(f"day_num:{treatment1} - day_num:{treatment2}")
        if test_params in ["exponential_only", "both"]:
            n1 = f"day_num_exp:Treatment_{treatment1}"
            n2 = f"day_num_exp:Treatment_{treatment2}"
            if n1 in param_names and n2 in param_names:
                L = np.zeros(n_params)
                L[param_names.index(n1)] = 1.0
                L[param_names.index(n2)] = -1.0
                contrasts.append(L)
                param_labels.append(f"day_num_exp:{treatment1} - day_num_exp:{treatment2}")

    if len(contrasts) == 0:
        return None, None, 0, []

    L_matrix = np.array(contrasts)
    L_beta = np.dot(L_matrix, fe_params.values)
    L_cov_L = np.dot(L_matrix, np.dot(cov_matrix.values, L_matrix.T))
    wald_stat = np.dot(L_beta, np.dot(np.linalg.inv(L_cov_L), L_beta))
    dof = len(contrasts)
    p_value = 1 - stats.chi2.cdf(wald_stat, dof)
    return wald_stat, p_value, dof, param_labels


trajectory_comparisons = []

print("\n" + "=" * 70)
print("PAIRWISE TRAJECTORY SHAPE COMPARISONS (Wald Test)")
print("=" * 70)

for treatment_a, treatment_b in COMPARISONS:
    wald_stat, p_value, dof, param_labels = compare_trajectories_wald(
        treatment_a, treatment_b, fe_params, cov_matrix,
        REFERENCE_TREATMENT, TRAJECTORY_TEST_PARAMS
    )
    if wald_stat is not None:
        trajectory_comparisons.append({
            'Comparison': f"{treatment_a} vs {treatment_b}",
            'Wald_statistic': wald_stat,
            'df': dof,
            'P_value': p_value,
            'Parameters_tested': ', '.join(param_labels),
        })
        sig = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "n.s."
        print(f"\n{treatment_a} vs {treatment_b}:")
        print(f"  Wald statistic: {wald_stat:.3f} (df={dof})")
        print(f"  P-value: {p_value:.4f} {sig}")
        print(f"  Parameters tested: {', '.join(param_labels)}")

trajectory_table = pd.DataFrame(trajectory_comparisons)
trajectory_table.to_csv(os.path.join(outdir, "trajectory_comparisons.csv"), index=False)

# ============================================================================
# STEP 8b: Test whether each treatment's slope differs from zero
# ============================================================================

print("\n" + "=" * 70)
print("SLOPE-VS-ZERO TESTS FOR EACH TREATMENT")
print("=" * 70)
print("\nFor each treatment, testing whether the trajectory slope")
print("(linear and/or exponential component) is significantly different from 0.")
print("This answers: 'Is aneuploidy actually changing over time for this treatment?'")

param_names_8b = list(fe_params.index)
n_params_8b = len(param_names_8b)

slope_vs_zero_rows = []

for treatment in treatments:
    if treatment in SLOPE_VS_ZERO_SKIP:
        continue
    test_mode = SLOPE_VS_ZERO_PARAMS.get(treatment, SLOPE_VS_ZERO_DEFAULT)
    contrasts = []
    contrast_labels = []

    if test_mode in ["linear_only", "both"]:
        L_lin = np.zeros(n_params_8b)
        L_lin[param_names_8b.index('day_num')] = 1.0
        if treatment != REFERENCE_TREATMENT:
            name = f"day_num:Treatment_{treatment}"
            if name in param_names_8b:
                L_lin[param_names_8b.index(name)] = 1.0
        contrasts.append(L_lin)
        contrast_labels.append("linear")

    if test_mode in ["exponential_only", "both"]:
        L_exp = np.zeros(n_params_8b)
        L_exp[param_names_8b.index('day_num_exp')] = 1.0
        if treatment != REFERENCE_TREATMENT:
            name = f"day_num_exp:Treatment_{treatment}"
            if name in param_names_8b:
                L_exp[param_names_8b.index(name)] = 1.0
        contrasts.append(L_exp)
        contrast_labels.append("exponential")

    if len(contrasts) == 0:
        continue

    print(f"\n  {treatment}:")

    # Joint Wald test (only meaningful when df > 1)
    L_matrix = np.array(contrasts)
    L_beta = np.dot(L_matrix, fe_params.values)
    L_cov_L = np.dot(L_matrix, np.dot(cov_matrix.values, L_matrix.T))
    wald_stat = np.dot(L_beta, np.dot(np.linalg.inv(L_cov_L), L_beta))
    dof = len(contrasts)
    p_joint = 1 - stats.chi2.cdf(wald_stat, dof)
    sig_joint = "***" if p_joint < 0.001 else "**" if p_joint < 0.01 else "*" if p_joint < 0.05 else "n.s."

    if dof > 1:
        print(f"    Joint test (slope ≠ 0): Wald={wald_stat:.3f}, df={dof}, p={p_joint:.4f} {sig_joint}")
        slope_vs_zero_rows.append({
            'Treatment': treatment,
            'Component': 'joint',
            'Estimate': np.nan,
            'SE': np.nan,
            'Wald_statistic': round(wald_stat, 3),
            'df': dof,
            'P_value': round(p_joint, 5),
            'Significance': sig_joint,
        })

    # Individual component tests (z-tests)
    for L, label in zip(contrasts, contrast_labels):
        estimate = np.dot(L, fe_params.values)
        se = np.sqrt(np.dot(L, np.dot(cov_matrix.values, L)))
        z = estimate / se
        p = 2 * (1 - stats.norm.cdf(np.abs(z)))
        sig_i = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."
        print(f"    {label:12s} slope: estimate={estimate:+.4f} ± {se:.4f}, z={z:+.2f}, p={p:.4f} {sig_i}")
        slope_vs_zero_rows.append({
            'Treatment': treatment,
            'Component': label,
            'Estimate': round(estimate, 5),
            'SE': round(se, 5),
            'Wald_statistic': round(z ** 2, 3),
            'df': 1,
            'P_value': round(p, 5),
            'Significance': sig_i,
        })

slope_vs_zero_df = pd.DataFrame(slope_vs_zero_rows)
slope_vs_zero_df.to_csv(os.path.join(outdir, "slope_vs_zero_tests.csv"), index=False)
print(f"\n  Saved: {os.path.join(outdir, 'slope_vs_zero_tests.csv')}")

# ============================================================================
# STEP 9: Pointwise contrasts at each observed timepoint
# ============================================================================

print("\n" + "=" * 70)
print("POINTWISE CONTRASTS AT EACH TIMEPOINT")
print("=" * 70)

param_names = list(fe_params.index)
n_params = len(param_names)


def build_prediction_contrast(treatment, reference, day, fe_params, cov_matrix):
    param_names = list(fe_params.index)
    n_params = len(param_names)
    L = np.zeros(n_params)

    ref_cat = REFERENCE_TREATMENT

    if treatment != ref_cat:
        trt_dummy = f"Treatment_{treatment}"
        if trt_dummy in param_names:
            L[param_names.index(trt_dummy)] = 1.0
        trt_day = f"day_num:Treatment_{treatment}"
        if trt_day in param_names:
            L[param_names.index(trt_day)] = day
        trt_exp = f"day_num_exp:Treatment_{treatment}"
        if trt_exp in param_names:
            L[param_names.index(trt_exp)] = np.exp(-day)

    if reference != ref_cat:
        ref_dummy = f"Treatment_{reference}"
        if ref_dummy in param_names:
            L[param_names.index(ref_dummy)] = -1.0
        ref_day = f"day_num:Treatment_{reference}"
        if ref_day in param_names:
            L[param_names.index(ref_day)] = -day
        ref_exp = f"day_num_exp:Treatment_{reference}"
        if ref_exp in param_names:
            L[param_names.index(ref_exp)] = -np.exp(-day)

    diff = np.dot(L, fe_params.values)
    se = np.sqrt(np.dot(L, np.dot(cov_matrix.values, L)))
    z = diff / se
    p = 2 * (1 - stats.norm.cdf(np.abs(z)))
    return diff, se, z, p


pointwise_rows = []
for treatment_a, treatment_b in COMPARISONS:
    trt_days = sorted(df_corr.loc[df_corr["Treatment"] == treatment_a, "day_num"].unique())
    for day in trt_days:
        diff, se, z, p = build_prediction_contrast(
            treatment_a, treatment_b, day, fe_params, cov_matrix
        )
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."
        pointwise_rows.append({
            'Comparison': f"{treatment_a} vs {treatment_b}",
            'Day': day,
            'Difference': round(diff, 3),
            'SE': round(se, 3),
            'z_statistic': round(z, 3),
            'P_value': round(p, 5),
            'Significance': sig,
        })
        print(f"  {treatment_a} vs {treatment_b} | day{day:>2d}: "
              f"diff={diff:+7.3f} ± {se:.3f}  z={z:+6.2f}  p={p:.4f} {sig}")

pointwise_df = pd.DataFrame(pointwise_rows)
pointwise_df.to_csv(os.path.join(outdir, "pointwise_contrasts.csv"), index=False)

# ============================================================================
# STEP 10: Generate model predictions for plotting
# ============================================================================

print("\n" + "=" * 70)
print("GENERATING MODEL PREDICTIONS")
print("=" * 70)

day_min = df_corr["day_num"].min()
day_max = df_corr["day_num"].max()
time_grid = np.linspace(day_min, day_max, 200)

predictions = []
for treatment in df_corr['Treatment'].cat.categories:
    for day in time_grid:
        row_data = {'Intercept': 1.0, 'day_num': day, 'day_num_exp': np.exp(-day)}
        for other_treatment in df_corr['Treatment'].cat.categories[1:]:
            dummy_name = f"Treatment_{other_treatment}"
            row_data[dummy_name] = 1.0 if treatment == other_treatment else 0.0
        for other_treatment in df_corr['Treatment'].cat.categories[1:]:
            interaction_name = f"day_num:Treatment_{other_treatment}"
            row_data[interaction_name] = day if treatment == other_treatment else 0.0
        for other_treatment in df_corr['Treatment'].cat.categories[1:]:
            interaction_name = f"day_num_exp:Treatment_{other_treatment}"
            row_data[interaction_name] = np.exp(-day) if treatment == other_treatment else 0.0
        predictions.append({'Treatment': treatment, 'day_num': day, **row_data})

pred_df = pd.DataFrame(predictions)
X_pred = pred_df[X_fixed.columns]
pred_df['fitted'] = np.dot(X_pred.values, fe_params.values)
pred_variances = np.sum((np.dot(X_pred.values, cov_matrix.values)) * X_pred.values, axis=1)
pred_df['se_fitted'] = np.sqrt(pred_variances)
pred_df['ci_lower'] = pred_df['fitted'] - 1.96 * pred_df['se_fitted']
pred_df['ci_upper'] = pred_df['fitted'] + 1.96 * pred_df['se_fitted']

print(f"\nGenerated {len(pred_df)} prediction points")
print(f"  {len(time_grid)} time points x {len(df_corr['Treatment'].cat.categories)} treatments")

# ============================================================================
# STEP 11: Create visualization
# ============================================================================

print("\n" + "=" * 70)
print("CREATING VISUALIZATION")
print("=" * 70)

sns.set_theme(style="whitegrid", context="talk")
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


def format_pval_stars(pval):
    if pval < 0.001:
        return "***"
    elif pval < 0.01:
        return "**"
    elif pval < 0.05:
        return "*"
    else:
        return "n.s."


palette = {
    "Etop 250nM": "#6498c6ff",
    "Etop 500 nM": "#3d4ca3ff",
    "Vingi WT 100ng": "#f9a72dff",
    "Cas9-B2M-TRAC-DKO": "#a71a1aff",
    "Cas9-B2M-KO": "#E97451",
    "Cas9-TRAC-KO": "#F35682",
    "Cas9-NO-TARGET": "#6495ED",
    "Nucleofection-only": "#54a24bff",
    "Vingi EN- RT- 100ng": "#4a4a4aff",
}


fig, ax1 = plt.subplots(1, 1, figsize=(6, 4))

for treatment in df_corr['Treatment'].cat.categories:
    treatment_data = df_corr[df_corr['Treatment'] == treatment]
    for rep_id in treatment_data['replicate_id'].unique():
        rep_data = treatment_data[treatment_data['replicate_id'] == rep_id].sort_values('day_num')
        ax1.plot(
            rep_data['day_num'], rep_data[response_col],
            marker='o', markersize=4, linestyle='-', linewidth=0.8, alpha=0.35,
            color=palette.get(treatment, None),
        )

for treatment in df_corr['Treatment'].cat.categories:
    pred_subset = pred_df[pred_df['Treatment'] == treatment]
    ax1.fill_between(
        pred_subset['day_num'], pred_subset['ci_lower'], pred_subset['ci_upper'],
        alpha=0.18, color=palette.get(treatment, None),
    )
    ax1.plot(
        pred_subset['day_num'], pred_subset['fitted'],
        linewidth=2, color=palette.get(treatment, None), label=treatment,
    )

observed_means = (
    df_corr.groupby(['Treatment', 'day_num'], observed=False)[response_col]
    .mean().reset_index()
)
for treatment in df_corr['Treatment'].cat.categories:
    mean_subset = observed_means[observed_means['Treatment'] == treatment]
    ax1.scatter(
        mean_subset['day_num'], mean_subset[response_col],
        s=40, color=palette.get(treatment, None),
        edgecolor='black', linewidth=1, zorder=5,
    )

ax1.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.5)

bl_avg_sems = []
for baseline_trt in set(BASELINE_MAP.values()):
    per_day_sems = baseline_sems[baseline_trt].dropna()
    bl_avg_sems.append(per_day_sems.mean())
bl_sem = max(bl_avg_sems)
bl_band = 1.96 * bl_sem
ax1.axhspan(-bl_band, bl_band, color='black', alpha=0.08,
            label='Baseline variability (95% CI)')

day3_value = float(df_corr['day_num'].min())
all_annotations = []

for annotate_trt in SLOPE_VS_ZERO_ANNOTATE:
    test_mode = SLOPE_VS_ZERO_PARAMS.get(annotate_trt, SLOPE_VS_ZERO_DEFAULT)
    if test_mode == "both":
        component = 'joint'
    elif test_mode == "exponential_only":
        component = 'exponential'
    else:
        component = 'linear'
    slope_row = slope_vs_zero_df[
        (slope_vs_zero_df['Treatment'] == annotate_trt) &
        (slope_vs_zero_df['Component'] == component)
    ]
    if not slope_row.empty:
        trt_p = slope_row.iloc[0]['P_value']
        trt_pred = pred_df[(pred_df['Treatment'] == annotate_trt) &
                           (np.abs(pred_df['day_num'] - day3_value) < 0.1)]
        if not trt_pred.empty:
            all_annotations.append({
                'type': 'slope_vs_zero',
                'y1': 0,
                'y2': trt_pred.iloc[0]['fitted'],
                'p_value': trt_p,
            })

for t1, t2 in COMPARISONS:
    comp_name = f"{t1} vs {t2}"
    pval_row = trajectory_table[trajectory_table['Comparison'] == comp_name]
    if not pval_row.empty:
        pred1 = pred_df[(pred_df['Treatment'] == t1) &
                        (np.abs(pred_df['day_num'] - day3_value) < 0.1)]
        pred2 = pred_df[(pred_df['Treatment'] == t2) &
                        (np.abs(pred_df['day_num'] - day3_value) < 0.1)]
        if not pred1.empty and not pred2.empty:
            all_annotations.append({
                'type': 'between_treatment',
                'y1': pred1.iloc[0]['fitted'],
                'y2': pred2.iloc[0]['fitted'],
                'p_value': pval_row.iloc[0]['P_value'],
            })

n_comparisons_plotted = len(all_annotations)
for ann in all_annotations:
    ann['p_value_bonferroni'] = min(ann['p_value'] * n_comparisons_plotted, 1.0)

for i, ann in enumerate(all_annotations):
    x_pos = df_corr['day_num'].min() - 0.3 - i * 0.3
    stars = format_pval_stars(ann['p_value_bonferroni'])

    ax1.plot([x_pos, x_pos], [ann['y1'], ann['y2']],
             color='black', linewidth=1, linestyle='-', zorder=10)
    tick_length = 0.1
    ax1.plot([x_pos - tick_length, x_pos + tick_length], [ann['y1'], ann['y1']],
             color='black', linewidth=1, linestyle='-', zorder=10)
    ax1.plot([x_pos - tick_length, x_pos + tick_length], [ann['y2'], ann['y2']],
             color='black', linewidth=1, linestyle='-', zorder=10)

    y_middle = (ann['y1'] + ann['y2']) / 2
    ax1.text(x_pos - 0.1, y_middle, stars,
             ha='center', va='center', fontsize=12, fontweight='bold',
             rotation=90, zorder=11)

    ann['x_position'] = x_pos

x_left = min(a['x_position'] for a in all_annotations) - 0.5 if all_annotations else df_corr['day_num'].min() - 0.5
ax1.set_xlabel('Day', weight='bold', fontsize=7)
ax1.set_ylabel('% Karyotypically Abnormal Cells', weight='bold', fontsize=7)
ax1.set_xticks(sorted(df_corr['day_num'].unique()))
ax1.set_xlim(x_left, 10 + 0.5)
ax1.set_ylim(-1, 11)
if LEGEND:
    ax1.legend(title='Treatment', frameon=True, loc='best', fontsize=6, title_fontsize=7)
ax1.grid(True, alpha=0.3)

ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_color('black')
ax1.spines['bottom'].set_color('black')
ax1.spines['left'].set_linewidth(2)
ax1.spines['bottom'].set_linewidth(2)
ax1.tick_params(axis='both', which='major', direction='out', length=8, width=2,
                colors='black', labelsize=6)
ax1.tick_params(axis='both', which='minor', direction='out', length=5, width=1.5,
                colors='black')
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')

plt.tight_layout()
plot_path = os.path.join(outdir, f'nlmm_aneuploidy_trajectories_baseline_corrected_{NAMEPIC}.pdf')
plt.savefig(plot_path, bbox_inches='tight')
plot_path = os.path.join(outdir, f'nlmm_aneuploidy_trajectories_baseline_corrected_{NAMEPIC}.png')
plt.savefig(plot_path, bbox_inches='tight', dpi=600)
plt.close()
print(f"\nPlot saved to: {plot_path}")

# ============================================================================
# STEP 12: Save results
# ============================================================================

print("\n" + "=" * 70)
print("SAVING RESULTS")
print("=" * 70)

with open(os.path.join(outdir, f'nlmm_full_model_summary_{NAMEPIC}.txt'), 'w') as f:
    f.write("=" * 70 + "\n")
    f.write("NON-LINEAR MIXED MODEL SUMMARY (BASELINE-CORRECTED)\n")
    f.write("=" * 70 + "\n\n")
    f.write("Baseline correction:\n")
    for trt, bl in BASELINE_MAP.items():
        f.write(f"  {trt} -= mean({bl}) per timepoint\n")
    f.write("\n")
    f.write(str(result.summary()))
    f.write("\n\n" + "=" * 70 + "\n")
    f.write("TRAJECTORY SHAPE COMPARISONS (Wald Test)\n")
    f.write("=" * 70 + "\n\n")
    if len(trajectory_table) > 0:
        f.write(trajectory_table[['Comparison', 'Wald_statistic', 'df', 'P_value']].to_string(index=False))
        f.write("\n\n")
        for _, row in trajectory_table.iterrows():
            sig = "***" if row['P_value'] < 0.001 else "**" if row['P_value'] < 0.01 else "*" if row['P_value'] < 0.05 else "n.s."
            f.write(f"  {row['Comparison']}: {sig}\n")
    f.write("\n\n" + "=" * 70 + "\n")
    f.write("SLOPE-VS-ZERO TESTS\n")
    f.write("=" * 70 + "\n\n")
    f.write(slope_vs_zero_df.to_string(index=False))
    f.write("\n\n" + "=" * 70 + "\n")
    f.write("POINTWISE CONTRASTS\n")
    f.write("=" * 70 + "\n\n")
    f.write(pointwise_df.to_string(index=False))
    f.write("\n\n" + "=" * 70 + "\n")
    f.write("MODEL NOTES\n")
    f.write("=" * 70 + "\n\n")
    f.write("Model equation (on baseline-corrected data):\n")
    f.write("Y = b0 + b1*day + b2*exp(-day) + treatment effects + interactions\n\n")
    f.write(f"Trajectory test mode: {TRAJECTORY_TEST_PARAMS}\n")

print(f"  Saved: {os.path.join(outdir, 'nlmm_full_model_summary.txt')}")
print(f"  Saved: {os.path.join(outdir, 'trajectory_comparisons.csv')}")
print(f"  Saved: {os.path.join(outdir, 'slope_vs_zero_tests.csv')}")
print(f"  Saved: {os.path.join(outdir, 'pointwise_contrasts.csv')}")

# ---- Final summary ----
print("\n" + "=" * 70)
print("TRAJECTORY COMPARISON SUMMARY")
print("=" * 70)
if len(trajectory_table) > 0:
    print(trajectory_table[['Comparison', 'Wald_statistic', 'df', 'P_value']].to_string(index=False))
    print()
    for _, row in trajectory_table.iterrows():
        sig = "***" if row['P_value'] < 0.001 else "**" if row['P_value'] < 0.01 else "*" if row['P_value'] < 0.05 else "n.s."
        print(f"  {row['Comparison']}: {sig}")

print("\n" + "=" * 70)
print("SLOPE-VS-ZERO SUMMARY")
print("=" * 70)
print(slope_vs_zero_df.to_string(index=False))

print("\n" + "=" * 70)
print("POINTWISE CONTRAST SUMMARY")
print("=" * 70)
print(pointwise_df.to_string(index=False))

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
