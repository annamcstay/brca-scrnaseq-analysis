# Leave-one-out robustness check — what went wrong and how to fix it

**File:** `Notebooks/22_phase4_leaveoneout_robustness.ipynb`, cell 6 (the "Cell 3 (corrected)" cell)
**Affects:** the GSE114725 macrophage result only. The HER2+ memory T-cell half of the notebook is correctly built and its result stands.
**Thesis sections to rewrite once re-run:** Results 3.4.6, Discussion 4.3.1, Limitations 4.6.1, Abstract.

---

## 1. The problem, in plain English

Each of the eight patients in GSE114725 donated cells from four different places: the tumour, nearby healthy tissue, blood, and a lymph node. The main analysis keeps those separate, which is the whole point — you are comparing tumour to healthy tissue, so a patient's tumour cells and their healthy cells have to be two separate entries. That gives **twelve entries to compare: eight tumour, four healthy.**

The robustness check instead pools **all** of a patient's cells into a single lump — tumour, healthy, blood and lymph node all mixed together — and then labels that one lump with whichever tissue happened to appear first for that patient in the file.

An analogy. Imagine testing whether a drug lowers blood pressure, where every participant gave four readings: morning, afternoon, evening and night. The real analysis compares morning readings against evening readings. The robustness check instead averages each person's four readings into a single number, then labels that number "morning" or "evening" depending on which of their readings happened to be typed in first — and quietly discards anyone whose first-typed reading was "afternoon" or "night".

Three things follow from that, and all three are visible in the notebook's own printed output:

1. **Three of the eight patients silently disappear.** BC1, BC2 and BC4 had blood or lymph node listed first, so they were never in the comparison at all. That is why "excluding BC1", "excluding BC2" and "excluding BC4" produce *byte-identical* numbers (FN1 log2FC 3.075921, padj 0.077732 in all three rows). You cannot remove someone who was never there. Those three rows are the tell.

2. **Every surviving entry is a mixture.** A single number blending tumour, healthy, blood and lymph-node cells, filed under one arbitrary label. So the test is comparing mixtures against mixtures, not tumour against healthy.

3. **Twelve entries become five.** Every iteration prints "3 TUMOR / 2 NORMAL" or "2 TUMOR / 2 NORMAL", against the true 8 and 4. With so few, the statistical test would not run at all — which is why the minimum-evidence threshold had to be relaxed from 10 to 3. That relaxation was the second warning sign; it was treated as the bug rather than as a symptom of it.

**The one-line summary:** a robustness check is supposed to re-run *the same analysis* on *slightly less* data. This re-ran a *different* analysis on *much less* data. Its verdict is about that different analysis, so it tells us nothing about the macrophage finding either way.

For scale: FN1 comes out at log2FC 3.08, padj 0.078 in the broken version, against **6.80, padj 0.0093** in the real analysis.

**What this does not mean.** It does not mean the macrophage finding is confirmed. It means it is currently untested. Once re-run properly it may still look fragile — four healthy-tissue samples is four healthy-tissue samples, and that limitation is real and worth stating. But it will at least be a limitation of the finding rather than of the harness.

**Where the correct code already lives.** Notebook 05 does this properly and carries an inline comment explaining exactly why:

> `# patients can have multiple tissues — need per-(patient,tissue) pseudobulk,`
> `# not per-patient alone, since a patient contributes separately to each tissue`

The fix is to reuse that, not to write anything new.

**And to be fair about how this happened.** Two lines of abandoned first-attempt code are still sitting in notebook 05, immediately *above* that warning comment:

```python
11 | for ct in viable_cell_types_1:
12 |     counts_df, meta_df = build_pseudobulk(          # per-patient — result discarded at line 35
13 |         adata1_raw, sample_col="patient", cell_type=ct, min_cells=10
14 |     )
15 |     # restrict metadata to tissue == TUMOR or NORMAL for this comparison
16 |     tissue_lookup = adata1_raw.obs.drop_duplicates("patient").set_index("patient")   # DEAD, never read
17 |     # patients can have multiple tissues — need per-(patient,tissue) pseudobulk,
18 |     # not per-patient alone, since a patient contributes separately to each tissue
19 |     subset = adata1_raw[adata1_raw.obs["cell_type"] == ct]
```

Line 12's result is overwritten at line 35 and `tissue_lookup` is never read, so notebook 05's own output is unaffected — it correctly reports 8 TUMOR / 4 NORMAL and 12 significant genes. But those two orphaned lines are exactly what notebook 22 picked up. So this was not a case of ignoring the warning; it was copying the dead code sitting directly above it, which still looked like working code. The root cause is that the correct logic was written inline inside notebook 05's loop rather than factored into a shared function, so notebook 22 reached for the *other* function in scope — `build_pseudobulk(..., "patient", ...)` — which does something subtly different.

---

## 2. Current code (notebook 22, cell 6)

The two defects are marked.

```python
patients = adata1_raw.obs["patient"].cat.categories.tolist()
target_genes = ["FN1", "HSPA1A", "HSPA1B"]
loo_results = []

for excluded_patient in patients:
    adata_loo = adata1_raw[adata1_raw.obs["patient"] != excluded_patient].copy()

    # ▼ DEFECT 1 — one pseudobulk per PATIENT, pooling all four tissues into one lump.
    #   The main analysis builds one per (patient, tissue).
    counts_df, meta_df = build_pseudobulk(adata_loo, "patient", "Macrophages")

    # ▼ DEFECT 2 — labels each patient-lump with whichever tissue appeared FIRST.
    #   Patients whose first row was BLOOD or LYMPHNODE get filtered out on the
    #   next line and vanish from the comparison without warning.
    tissue_lookup = adata_loo.obs.drop_duplicates("patient").set_index("patient")["tissue"]
    meta_df["tissue"] = meta_df["patient"].map(tissue_lookup)

    meta_sub = meta_df[meta_df["tissue"].isin(["TUMOR", "NORMAL"])].copy()
    meta_sub["tissue"] = meta_sub["tissue"].astype(str)
    counts_sub = counts_df[meta_sub.index]

    n_tumor = (meta_sub["tissue"] == "TUMOR").sum()
    n_normal = (meta_sub["tissue"] == "NORMAL").sum()
    print(f"Excluding {excluded_patient}: {n_tumor} TUMOR / {n_normal} NORMAL samples remaining")

    if n_tumor < 2 or n_normal < 2:
        print(f"  SKIPPED — insufficient samples after exclusion")
        continue

    # ▼ CONSEQUENCE of defects 1-2 — only ~5 samples survive, so the project's
    #   standard >=10 threshold could not run and was lowered to >=3.
    #   The threshold was never the problem.
    gene_filter = (counts_sub > 0).sum(axis=1) >= 3
    counts_sub_filtered = counts_sub[gene_filter]

    try:
        dds = DeseqDataSet(counts=counts_sub_filtered.T.astype(int), metadata=meta_sub,
                            design_factors="tissue", refit_cooks=True, quiet=True)
        dds.deseq2()
        ds = DeseqStats(dds, contrast=["tissue", "TUMOR", "NORMAL"], quiet=True)
        ds.summary()
        results = ds.results_df.copy()
        ...
```

---

## 3. Proposed fix

Two changes. First add a shared helper that reproduces notebook 05's per-(patient, tissue) construction; then rewrite the loop to use it together with notebook 05's own `run_pydeseq2` runner.

### 3.1 New helper cell — insert before the leave-one-out loop

Note `run_pydeseq2` is defined in notebook 05 but **not** in notebook 22, so it has to be copied across verbatim as well — reproduced in full in §3.2 below. Copying it unchanged is deliberate: the defence is that the robustness check calls the identical function as the primary analysis.

```python
# ----------------------------
# Per-(patient, tissue) pseudobulk — lifted verbatim from notebook 05's
# GSE114725 loop, where it was written inline. A patient contributes
# SEPARATELY to each tissue, so the sample key must be patient + tissue,
# never patient alone.
#
# Factored out here so the robustness check and the primary analysis
# cannot drift apart again. Ideally notebook 05 is refactored to import
# this same function rather than keeping its own inline copy.
# ----------------------------
def build_pseudobulk_by_patient_tissue(adata_raw, cell_type,
                                       tissues=("TUMOR", "NORMAL"),
                                       min_cells=10):
    """
    Returns (counts_df, meta_df) for one cell type.
      counts_df : genes x samples, raw summed counts
      meta_df   : one row per sample, with a 'tissue' column
    One sample == one (patient, tissue) pair, matching notebook 05.
    """
    subset = adata_raw[adata_raw.obs["cell_type"] == cell_type]
    subset = subset[subset.obs["tissue"].isin(list(tissues))]
    subset.obs["sample_id"] = (subset.obs["patient"].astype(str)
                               + "_" + subset.obs["tissue"].astype(str))

    pseudobulk_samples, sample_ids, tissue_labels, n_cells_list = [], [], [], []
    for sid in subset.obs["sample_id"].unique():
        mask = (subset.obs["sample_id"] == sid).values
        n_cells = mask.sum()
        if n_cells < min_cells:
            continue
        pseudobulk_samples.append(np.asarray(subset.X[mask].sum(axis=0)).flatten())
        sample_ids.append(sid)
        tissue_labels.append(subset.obs.loc[mask, "tissue"].iloc[0])
        n_cells_list.append(n_cells)

    counts_df = pd.DataFrame(pseudobulk_samples, index=sample_ids,
                             columns=subset.var_names).T
    meta_df = pd.DataFrame({"tissue": tissue_labels, "n_cells": n_cells_list},
                           index=sample_ids)
    return counts_df, meta_df
```

### 3.2 Replacement for cell 6

```python
# ----------------------------
# Leave-one-out robustness check: Macrophages, Tumour vs Normal.
#
# Each iteration drops one patient and re-runs the IDENTICAL analysis
# reported in notebook 05 — same per-(patient, tissue) pseudobulk, same
# >=10-sample gene filter, same PyDESeq2 runner. A robustness check must
# differ from the primary analysis in the excluded patient and in
# NOTHING else, or it is testing a different question.
# ----------------------------
patients = adata1_raw.obs["patient"].cat.categories.tolist()
target_genes = ["FN1", "HSPA1A", "HSPA1B"]
loo_results = []

for excluded_patient in patients:
    adata_loo = adata1_raw[adata1_raw.obs["patient"] != excluded_patient].copy()

    counts_df, meta_df = build_pseudobulk_by_patient_tissue(
        adata_loo, cell_type="Macrophages", tissues=("TUMOR", "NORMAL"), min_cells=10
    )

    n_tumor  = (meta_df["tissue"] == "TUMOR").sum()
    n_normal = (meta_df["tissue"] == "NORMAL").sum()

    # Sanity gate — the full analysis has 8 TUMOR / 4 NORMAL, so removing
    # one patient must leave 7 tumour samples and 3 or 4 normal. Anything
    # far below that means samples are being dropped silently, which is
    # exactly the failure this rewrite exists to prevent.
    print(f"Excluding {excluded_patient}: {n_tumor} TUMOR / {n_normal} NORMAL "
          f"({len(meta_df)} samples total)")
    assert n_tumor >= 6, f"Only {n_tumor} tumour samples — samples are being lost"

    results = run_pydeseq2(                      # notebook 05's runner, unchanged
        counts_df, meta_df,
        group_col="tissue", group_a="TUMOR", group_b="NORMAL",
        cell_type="Macrophages",
        comparison_name=f"Tumor_vs_Normal_excl_{excluded_patient}",
        dataset_name="GSE114725",
        min_genes_expressed=10,                  # matches the primary analysis
    )

    if results is None:
        print(f"  SKIPPED — model did not fit")
        del adata_loo; gc.collect()
        continue

    row = {"excluded_patient": excluded_patient,
           "n_tumor": n_tumor, "n_normal": n_normal,
           "n_sig_total": int((results["padj"] < 0.05).sum())}
    for gene in target_genes:
        if gene in results.index:
            row[f"{gene}_log2FC"]      = results.loc[gene, "log2FoldChange"]
            row[f"{gene}_padj"]        = results.loc[gene, "padj"]
            row[f"{gene}_significant"] = bool(results.loc[gene, "padj"] < 0.05)
        else:
            row[f"{gene}_log2FC"] = None
            row[f"{gene}_padj"] = None
            row[f"{gene}_significant"] = False
    loo_results.append(row)

    del adata_loo, results
    gc.collect()

loo_df = pd.DataFrame(loo_results)
print("\n=== Leave-one-out results (per-(patient,tissue) pseudobulk) ===")
print(loo_df.to_string(index=False))
loo_df.to_csv(RESULTS_DIR / "GSE114725_macrophage_leaveoneout_FIXED.csv", index=False)

# Guard against the failure mode that produced the original result:
# identical rows mean a patient was never in the comparison to begin with.
dupes = loo_df.duplicated(subset=[f"{g}_log2FC" for g in target_genes], keep=False)
if dupes.any():
    print("\nWARNING — identical result rows detected for: "
          f"{loo_df.loc[dupes, 'excluded_patient'].tolist()}. "
          "Those patients are probably not entering the comparison.")
```

### 3.3 What to expect when it runs

- **Eight iterations, not six.** With per-(patient, tissue) pseudobulk, dropping any single patient still leaves at least 7 tumour and 3 healthy samples, so every exclusion is viable. The thesis's current "6 valid exclusions given sample-size constraints on the remaining 2" is itself an artefact of the broken version and needs correcting too.
- **The printed sample counts are the check.** They should read 7 TUMOR / 3 NORMAL when dropping BC1, BC2, BC3 or BC7 (the four who contributed healthy tissue) and 7 TUMOR / 4 NORMAL when dropping BC4, BC5, BC6 or BC8. If any iteration prints something much smaller, stop — the pooling bug is back. The `assert` is there to catch that.
- Some (patient, tissue) pairs may still legitimately drop out under the `min_cells=10` macrophage minimum — BC3's healthy sample is small. That is a real, principled exclusion, and the printed counts will show it.
- **Do not lower the `>=10` gene filter.** With 10–11 samples per iteration it is once again workable, and keeping it identical to notebook 05 is the whole point: the defence in a viva is "the robustness check is literally the same function as the primary analysis". If you want to be strictly proportional, `round(10/12 * n_samples)` is the equivalent — but run the plain `>=10` version as the headline.
- **No prediction about the outcome.** It may still show fragility. Report whatever it gives.

### 3.4 While you are in there

Two things worth doing in the same pass, since they touch the same code:

- **Factor notebook 05 onto the shared helper too**, so there is exactly one implementation of per-(patient, tissue) pseudobulk in the project. The duplication is the root cause here, and leaving two copies invites the same regression a third time.
- **Consider a paired design.** Four patients (BC1, BC2, BC3, BC7) contributed both tumour and healthy tissue, so patient is partly confounded with tissue in the current unpaired 8-vs-4 comparison. Running `~patient + tissue`, or a paired analysis restricted to those four patients, is a stronger test and a good sensitivity check to report alongside. Separate from this fix, but the same cell.
