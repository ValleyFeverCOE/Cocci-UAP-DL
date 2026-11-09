# cocci-n-of-1-pathways-zscore

Directional N-of-1 MixEnrich on cocci log2 fold-change data (7,122 genes x 51 patients),
producing lnOR, variance, and z-score matrices for 5,179 GO terms x 2 directions (`_up`/`_down`).

- `data/cocci_in.rds` — list of `GO_anno` (annotation) and `log2fcs` (genes x subjects); gene
  symbols repaired from Excel date-mangling (e.g. `2-Sep` -> `SEPT2`). Same fix applied in
  `data/cocci_log2fc_7122_genes_51_patients.csv`.
- Pipeline code taken unchanged from inter-nof1 (`MixEnrich.R`, `runMixEnrich.R`,
  `differential_expression_direction.R`, `ef_method_pipeline.R`).
- Parameters: posterior probability cutoff 0.99, |log2FC| > log2(1.2), two-sided.

Run:

```sh
Rscript run_all.R   # writes output/nof1_{lnOR,var,zscore}_matrix_7122.csv
```

Requires R with `mixtools`, `doParallel`, `data.table`.

