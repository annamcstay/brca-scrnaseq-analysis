# Methods: Differential Expression, Pathways, and Cell-Cell Interactions (Phase 3)

Written as analysis progresses, per the same practice used for `methods_clustering.md`.

---

## 1. Differential expression design: comparison group selection

Pseudobulk DE requires enough independent patient-level samples per group to estimate
biological variance reliably (a practical floor of ≥3 patients per group, each contributing
enough cells to form a representative pseudobulk profile — used ≥10 cells per patient per
cell type as the minimum for a usable pseudobulk sample).

### GSE114725

Initial plan was a tumour-vs-blood comparison (matching the tissue conditions in the
original dataset design). Checking patient-level sample sizes before committing to this
comparison revealed only 2 of 8 patients had any blood samples at all (BC1, BC4) — too few
to estimate variance in a blood group, regardless of cell counts within those two patients.
**Comparison switched to tumour vs normal-adjacent tissue**, which has better patient
coverage (8 tumour patients, 4 normal patients).

A per-patient, per-cell-type cell count check was then run to determine which cell types
have adequate data in both groups:

| Cell type | Tumour patients viable (≥10 cells) | Normal patients viable (≥10 cells) | Included in DE? |
|---|---|---|---|
| T cells | 8/8 | 3/4 | Yes |
| CD8/Effector T cells | 8/8 | 4/4 | Yes |
| NK/Cytotoxic T cells | 8/8 | 3/4 | Yes |
| B cells | 8/8 | 3/4 | Yes |
| Macrophages | 8/8 | 4/4 | Yes |
| Monocytes/DC | 3/8 | 2/4 | No — insufficient patients in both groups |
| Mast cells | 6/8 | 1/4 | No — essentially absent from normal tissue in most patients |
| pDC | 5/8 | 0/4 | No — zero patients clear the floor in normal tissue |
| Mixed/stromal-contaminated | 3/8 | 1/4 | No — excluded on both power and identity grounds (known contamination cluster, not a genuine cell type; see `methods_clustering.md` §5) |

**GSE114725 pseudobulk tumour-vs-normal DE run on 5 of 9 cell types**: T cells,
CD8/Effector T cells, NK/Cytotoxic T cells, B cells, Macrophages. The remaining 4 cell
types are either too rare or too sparsely distributed across patients in normal tissue to
support reliable pseudobulk aggregation, and were excluded from DE by design rather than
run and reported with an unreliable result.

### GSE176078

Comparison groups are the three breast cancer subtypes present in the dataset (ER+, n=11
samples; HER2+, n=5 samples; TNBC, n=10 samples). The same per-sample, per-cell-type
viability check was applied (≥10 cells per sample, ≥3 viable samples per subtype required).

| Cell type | ER+ viable | HER2+ viable | TNBC viable | Included in DE? |
|---|---|---|---|---|
| Endothelial cells | 11/11 | 5/5 | 9/10 | Yes |
| CAFs | 10/11 | 5/5 | 10/10 | Yes |
| PVL | 10/11 | 5/5 | 9/10 | Yes |
| B cells | 8/11 | 5/5 | 5/10 | Yes |
| CD8 T cells | 10/11 | 5/5 | 8/10 | Yes |
| NK cells | 10/11 | 5/5 | 8/10 | Yes |
| Memory T cells | 10/11 | 5/5 | 8/10 | Yes |
| T cells | 10/11 | 5/5 | 8/10 | Yes |
| Cycling epithelial | 7/11 | 3/5 | 8/10 | Yes |
| Macrophages | 11/11 | 5/5 | 10/10 | Yes |
| Epithelial (ambiguous) | 6/11 | 4/5 | 8/10 | Yes |
| Luminal epithelial | 9/11 | 4/5 | 9/10 | Yes |
| Basal epithelial | 2/11 | 2/5 | 3/10 | No — rare cell type dataset-wide (1,073 cells total), insufficient in all three subtypes |
| pDC | 2/11 | 2/5 | 3/10 | No — rare cell type dataset-wide (315 cells total), insufficient in all three subtypes |
| Plasma cells | 7/11 | 2/5 | 5/10 | No — fails specifically in HER2+ (smallest subtype group) |
| Cycling T cells | 2/11 | 5/5 | 6/10 | No — fails specifically in ER+; consistent with ER+ having the lowest cycling T cell proportion in composition analysis (Phase 2) |
| Unassigned (n=28, artefact) | 0/11 | 0/5 | 1/10 | No — known contamination cluster, excluded on identity grounds as well as power |

**GSE176078 pseudobulk subtype DE run on 12 of 17 cell types.**

**Comparison structure: pairwise, not a single 3-group test.** ER+ vs TNBC, ER+ vs HER2+,
and TNBC vs HER2+ run separately per cell type, rather than one omnibus 3-group test. This
gives more directly interpretable results for the write-up (which specific pair differs,
and in which direction) at the cost of more individual tests — mitigated by FDR correction
applied within each comparison.

## 2. Scope: broad cell types first, sub-cluster DE as a possible extension

Pseudobulk DE was run on the broad Phase 2 `cell_type` labels (5 for GSE114725, 12 for
GSE176078 — see above), not the finer T cell/macrophage sub-cluster labels from Phase 2's
stretch sub-clustering (`tcell_subtype`, `mac_subtype`). This matches Adrien's core brief,
which asks for within-cell-type DE at the level already established, not automatically at
the sub-cluster level.

Sub-cluster-level DE (e.g. "how do LAM-like macrophages specifically differ tumour vs
normal") is a natural follow-on question, but was deliberately deferred rather than
attempted alongside the core analysis: sub-clusters are smaller slices of already
sample-size-limited cell types (e.g. macrophages had only 4 viable normal-tissue patients
at the broad level), so a per-sub-type viability check is needed before committing to it,
the same way it was needed at the broad level. Revisited as a stretch task once core DE,
pathway enrichment, and cell-cell interaction analysis are complete.

## 3. Pseudobulk DE results: run outcome and key caveats

All 5 GSE114725 (Tumour vs Normal) and 36 GSE176078 (12 cell types × 3 pairwise subtype
comparisons) pseudobulk DE tests completed successfully. Full results and FDR-significant
(padj<0.05) gene lists saved per cell type/comparison in `results/phase3_pseudobulk_de/`;
summary table in `phase3_DE_summary_all_comparisons.csv`.

**Technical note — categorical dtype bug.** Initial runs of the GSE176078 pairwise
comparisons failed uniformly with `LinAlgError: Singular matrix`. Cause: `subtype` was
stored as a pandas `Categorical` retaining all 3 original levels even after filtering to
2 groups for a pairwise comparison, causing PyDESeq2 to build a design-matrix column for
the absent third category (all zeros), which cannot be inverted. Fixed by casting the
group column to plain string after filtering, immediately before fitting.

**Key caveat — HER2+ vs ER+ comparisons are underpowered.** HER2+_vs_ER+ showed 0
significant genes in 8 of 12 GSE176078 cell types (Endothelial, CAFs, PVL, B cells,
CD8 T cells, NK cells, Memory T cells, Macrophages), a uniform pattern across nearly the
entire immune/stromal compartment. Given HER2+ has only 5 samples (vs 10-11 for ER+/TNBC),
this most likely reflects reduced statistical power rather than genuinely absent biological
difference, and should not be interpreted as evidence that HER2+ and ER+ microenvironments
are transcriptionally identical. Where HER2+_vs_ER+ *did* show substantial hits — Luminal
epithelial (565 genes) and Cycling epithelial (28 genes) — these are cases where HER2
status drives a large enough intrinsic transcriptional effect to remain detectable despite
low power, consistent with HER2 signalling's well-established direct effect on epithelial
cell state.

**TNBC comparisons dominate significant findings across nearly all cell types**, consistent
with TNBC's established position as the most transcriptionally distinct breast cancer
subtype. The largest single effect was Luminal epithelial cells, TNBC vs ER+ (1,602 of
15,800 genes tested, ~10%) — biologically plausible given the two groups represent
maximally different differentiation programs, but also the best-powered comparison (9 vs 9
samples), so effect size should not be directly compared to less-powered comparisons
without accounting for this.

**Minor caveat — B cells.** All three B cell comparisons used DESeq2's mean-based
dispersion trend fallback (parametric fitting did not converge), on the smallest gene
universe of any cell type (3,910-5,879 genes tested). B cell TNBC-related results (82 and
159 significant genes) should be treated with somewhat more caution than better-powered,
better-converged cell types until corroborated by pathway-level analysis.

## 4. Deferred / optional refinement — flagged for later, not blocking

**GSE114725 Tumour-vs-Normal: partial patient overlap not accounted for in the design.**
All 4 patients contributing Normal samples (BC1, BC2, BC3, BC7) also contribute Tumour
samples — the comparison is not built from fully independent samples per group, even
though it's run as an unpaired test (`~tissue`). This could make p-values mildly
anti-conservative, since some tumour/normal variance is actually stable within-patient
signal a blocked design (`~patient + tissue`) would separate out properly. Not fixed now:
with only 12 total samples and partial (not full) overlap, a blocked model may not have
enough degrees of freedom to add real value, and the current results are otherwise sound.
Worth revisiting only if reviewing for maximum statistical rigour with time to spare —
not required for a complete, defensible analysis.

**Pathway gene set coverage — future work item.** Macrophage Tumour-vs-Normal DE genes
showed biological coherence at the individual gene level (FN1: matrix remodeling;
HSPA1A/B: stress response) but did not reach pathway-level significance — only 3 of 12
significant genes were annotated within the MSigDB Hallmark gene set collection (50 broad,
curated pathways), insufficient for Fisher's exact test to detect enrichment after
correction.

This is not an isolated case. Both GSE176078 B cell comparisons involving TNBC
(TNBC vs ER+ and TNBC vs HER2+) showed volcano plots visually dominated by a heat-shock/
chaperone gene cluster (HSPA1A, HSPA1B, DNAJB1, DNAJB4 among the top hits in both), yet
neither comparison's pathway dot plot contains any heat-shock-specific term — Hallmark has
no dedicated heat-shock response gene set, so this signal is invisible at the pathway
level despite being one of the clearest patterns in the gene-level results. TNF-alpha
Signaling via NF-kB tops both B cell comparisons instead, likely reflecting whichever
genes in the significant list do have Hallmark annotations, not necessarily the dominant
underlying biology.

**Future work:** re-run comparisons showing this gene-level/pathway-level disconnect
against a broader gene set library (GO Biological Process or Reactome — both substantially
larger than Hallmark's 50 pathways, and GO in particular has dedicated heat-shock-response
terms) to check whether currently-unannotated significant genes map to real pathways those
libraries do cover.

## 5. Cell-cell communication (LIANA)

Run per condition rather than once per dataset, matching the DE/pathway comparison
structure: GSE114725 Tumour vs Normal (2 networks), GSE176078 ER+/HER2+/TNBC (3 networks).
Input: `adata.raw` (log-normalised, unscaled, full gene set) — the correct input for
LIANA, consistent with the convention already established for marker validation in Phase
2. Method: `li.mt.rank_aggregate`, LIANA's consensus method combining multiple ligand-
receptor scoring approaches (including CellPhoneDB-style analysis). Cell types with <50
cells in a given condition, and the known contamination artefact cluster, excluded from
that condition's network.

**Ranking metric — specificity_rank, not magnitude_rank.** Initial results sorted by
`magnitude_rank` were dominated almost entirely by a single interaction (B2M → KLRD1,
MHC-I sensing by NK cells) repeated across nearly every row with only the source cell type
varying. This is genuine, correct biology, but not distinctive — B2M is a near-universally
expressed MHC-I component, so magnitude-based ranking (driven by absolute expression level)
mathematically favours it over more cell-type-specific signalling regardless of biological
interest. Switched to `specificity_rank` (how unusually strong an interaction is for a
*particular* cell-type pair, not just how highly expressed the genes are overall) and
de-duplicated by ligand-receptor pair, surfacing a substantially more varied and
distinctive set of interactions per network.

**Key findings:**
- **GSE114725 Tumour**: broad CXCL8 (IL-8) inflammatory signalling hub radiating from
  Monocytes/DC and Macrophages to nearly every other retained cell type (B cells, Mast
  cells, pDC, the stromal-contaminated cluster) via CXCR1/CXCR2/ACKR1/CD79A. Tight,
  specific S100A8/S100A9 → CD36/CD68 signalling between Macrophages and Mast cells.
- **GSE114725 Normal**: same core myeloid axis present but substantially narrower —
  signalling concentrated almost entirely within Monocyte/DC ↔ Macrophage crosstalk
  (S100A8/9-CD68/CD36, MMP9-LRP1), without Tumour's broader cast of receiving cell types.
  Consistent with — and independently corroborating — the pseudobulk DE finding that
  Macrophages showed the largest tumour-vs-normal transcriptional shift of any tested cell
  type (FN1, HSPA1A/B; see §3), now shown as an expanded signalling network as well as an
  expression change.
- **GSE176078 (all three subtypes)**: a conserved stromal-vascular signalling programme —
  collagen (COL1A1/COL1A2 → CD93/FLT4) and complement (C1QA/C1QB → CD93/LRP1) signalling
  from CAFs/Macrophages into Endothelial cells — present with largely the same partners and
  mechanism across ER+, HER2+, and TNBC. This is a legitimate finding in its own right:
  what is *conserved* across subtypes, rather than only what differs.
- **pDC → Endothelial cells via SCT-RAMP2/3** appeared consistently across all three
  GSE176078 subtypes. Checked for robustness given pDC's small size (315 cells total)
  before reporting: SCT is expressed in 60-69% of pDC cells across all three subtypes
  (not driven by a handful of outlier cells), with consistent mean expression (~1.0-1.2)
  in each. Confirmed as a genuine, if under-studied in this context, signalling axis rather
  than an artefact — secretin/RAMP signalling is a recognised vasoactive pathway in other
  physiological contexts, worth noting as a candidate for further investigation rather than
  a definitive functional claim.

Full results and top-10 (by specificity_rank) tables per network saved in
`results/phase3_liana/`; dot plots per network in `figures/phase3_liana/`.

## 6. Compositional analysis (scCODA)

Standard non-parametric tests (used as a first pass in Phase 2) don't formally account
for the compositional nature of cell-type proportions (they sum to 100% per sample, so
categories aren't statistically independent). scCODA — a Bayesian Dirichlet-Multinomial
model built specifically for compositional data — used here as the proper replacement.

**Reference cell type chosen deliberately per comparison**, not defaulted: T cells for
GSE114725 (large, stable, and not the a priori candidate for driving tumour/normal
difference — that's macrophages, per existing DE/LIANA findings), PVL for GSE176078
(structural/stromal, unlikely to be a primary driver of subtype-specific differences
compared to CAFs or Macrophages, which already showed real signal in DE/pathway analysis).
A reference plausibly central to the biology under test would risk biasing results toward
finding everything else "changed" relative to it.

**Technical note — same categorical dtype bug as the DE notebook.** Initial GSE114725 run
produced an Effects table listing three tissue levels (including a phantom
`tissue[T.LYMPHNODE]` with zero actual samples in this comparison) and every single effect
came back exactly 0.0 across all 8 cell types — a strong sign of a degenerate fit rather
than a genuine uniform result. Cause: identical to the PyDESeq2 issue (§3) — the tissue
column retained unused categorical levels after filtering. Fixed the same way, casting to
plain string before model fitting.

**GSE114725 Tumour vs Normal — result.** No cell type showed a credible compositional
change (`credible_effects()` all `False`). Importantly, this should not be read as
"confidently no difference" — inclusion probabilities (the model's posterior confidence
that each effect is real) clustered around 0.47-0.57 for every cell type, i.e. near chance
level, rather than near zero (which would indicate confident exclusion). This reflects
insufficient statistical power from the small sample (4 tumour, 4 normal patients) rather
than genuine compositional equivalence.

**Notable detail:** Macrophages showed the highest inclusion probability of any cell type
(0.572), still short of credibility but directionally consistent with — and worth reading
alongside — the macrophage-specific findings already established: 12 significant
pseudobulk DE genes (§3, FN1/HSPA1A/B) and an expanded cell-cell signalling network in
tumour vs normal tissue (§5). Taken together, this suggests macrophage *abundance* does
not credibly shift between tumour and normal tissue, while macrophage *transcriptional
state and signalling activity* clearly do — a reprogramming-without-recruitment pattern,
consistent across three independent analytical approaches (DE, LIANA, and now
compositional analysis), rather than a simple increase in macrophage numbers.

**GSE176078 subtype comparisons — result: Memory T cells is a credible, directional, dual-
confirmed finding.** Across the three pairwise comparisons (TNBC vs ER+, HER2+ vs ER+,
TNBC vs HER2+), Memory T cells was the only cell type of 16 to reach credible status, and
it did so consistently: `credible_effects()` = `True` for both **HER2+ vs ER+**
(log2-fold change = +1.42, i.e. higher in HER2+) and **TNBC vs HER2+** (log2-fold change =
-1.35, i.e. lower in TNBC relative to HER2+). Taken together, these two independent
credible results describe the same underlying pattern: **Memory T cells are specifically
elevated in HER2+ tumours relative to both other subtypes**, not merely "different between
some pair." This corroborates the raw composition numbers already observed descriptively
in Phase 2 (HER2+ 20.5% Memory T cells vs ER+ 5.7%, TNBC 8.7% — the largest subtype
difference in that table), now confirmed with a proper compositional statistical model
rather than left as a descriptive observation. TNBC vs ER+ showed no credible effects,
consistent with the HER2+-specific nature of this pattern rather than a general
TNBC-vs-others difference.

**Technical note.** All three GSE176078 runs reported *"Zero counts encountered in data!
Added a pseudocount of 0.5"* — standard, expected handling for compositional models when
a rare cell type (e.g. pDC, Basal epithelial) has zero cells in at least one sample; not
indicative of an error.
