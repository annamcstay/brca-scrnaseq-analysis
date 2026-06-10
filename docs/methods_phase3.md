Phase 3 — Cell-Cell Communication

Cell-Cell Communication Inference with LIANA

Cell-cell communication analysis was performed using LIANA (Ligand-receptor Interaction ANAlysis) version 1.7.3, which aggregates multiple ligand-receptor inference methods into a single consensus ranking. LIANA was run using the rank\_aggregate function, which combines five independent methods: CellPhoneDB, Connectome, log2FC, NATMI, and SingleCellSignalR. Using multiple methods reduces the false positive rate associated with any individual approach and produces a more robust ranking of interactions.

The consensus ligand-receptor resource was used, which is a curated combination of multiple publicly available databases including CellPhoneDB, CellChatDB, CellTalkDB, and others. Interactions were filtered to ligand-receptor pairs where the ligand was expressed in at least 10% of cells in the source cell type (expr\_prop=0.1) and where a minimum of 10 cells were present per cell type (min\_cells=10). Analysis was performed on log-normalised expression values in adata.X.

For GSE114725, LIANA was run across all six annotated immune cell types. For GSE176078, analysis was restricted to immune cell populations (T cells, CD8 T cells, NK cells, Naive/memory T cells, B cells, Plasma cells, Macrophages, Monocytes/DC, pDC) as the full dataset including epithelial and stromal populations exceeded available RAM for this computation step.

Significant interactions were defined as those with a specificity rank below 0.05, reflecting interactions that are specific to particular cell type pairs rather than broadly expressed across all pairs.

Results Summary

GSE114725 — 4,644 total ligand-receptor pairs were tested, of which 458 (9.9%) were statistically specific (specificity rank < 0.05). Key interactions included S100A8/S100A9 signalling to CD36 and CD68 receptors from Monocytes/DC, reflecting alarmin-mediated myeloid activation in the tumour microenvironment. CXCL8 signalling via CXCR2 and ACKR1 was identified as a prominent inflammatory chemokine interaction. IL1B→IL1R2 signalling between Activated T cells reflected inflammatory cytokine crosstalk. LTF→IL1RL1 interactions were identified between Activated T cells and Monocytes/DC, suggesting a regulatory axis between these populations.

GSE176078 — 25,434 total ligand-receptor pairs were tested across nine immune cell types, of which 3,086 (12.1%) were significant. The most specific interactions were dominated by CCL19→CCR7 signalling from Monocytes/DC to multiple immune populations including T cells, B cells, CD8 T cells and NK cells, consistent with a key role for Monocytes/DC in lymphocyte recruitment and homing within the breast tumour microenvironment. CXCL9→FCGR2A macrophage-to-macrophage signalling was identified, reflecting autocrine macrophage activation. CD14→TLR9 interactions between Macrophages and pDC indicated innate immune crosstalk. PLAU→LRP1 plasminogen activator signalling between Macrophages and Monocytes/DC was also identified, consistent with remodelling activity in the TME.

Limitations

LIANA predicts interactions based on co-expression of ligand-receptor pairs and does not confirm physical binding or downstream signalling activity. The consensus resource covers approximately 28% of human ligand-receptor pairs, meaning interactions outside the resource cannot be detected. For GSE176078, restricting analysis to immune cells means stromal-immune interactions (e.g. CAF→T cell or Endothelial→immune) were not captured in this analysis and represent a direction for future work.

