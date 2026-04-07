# Mini DESeq2 Demo Fixture

This synthetic dataset is meant for quick manual runs of `q2-deseq2` and for
visual inspection of the generated visualization.

## Design

- Six samples total, with duplicate replicates for each condition.
- `control`: `CTRL_1`, `CTRL_2`
- `treated_a`: `TRTA_1`, `TRTA_2`
- `treated_b`: `TRTB_1`, `TRTB_2`
- The metadata also includes a `batch` column so the same fixture can be used
  with `fixed_effects_formula='batch + condition'`.

## Expected Signal

- `gene_alpha` is strongly up in `treated_a` only.
- `gene_beta` is strongly down in `treated_a` only.
- `gene_gamma` is strongly up in `treated_b` only.
- `gene_delta` is strongly down in `treated_b` only.
- `gene_shared` is up in both treated groups relative to `control`.
- `gene_housekeeping`, `gene_balanced`, and `gene_batch_shift` are mostly
  stable, with `gene_batch_shift` showing a mild batch effect.
- `gene_low_count` has a total count of `5`, so it should be removed whenever
  `min_total_count >= 10`.

## Files

- `feature-table.tsv` is the human-readable count matrix.
- `feature-table.biom` is the import-ready BIOM version of the same matrix.
- `sample-metadata.tsv` contains the `condition` and `batch` columns.
- `reference.gff3` adds gene names and products to the visualization.

## Suggested Use

- Import `feature-table.biom` as `FeatureTable[Frequency]`.
- Use the `condition` column from `sample-metadata.tsv` with
  `reference_level=control` for the simple two-step or pipeline workflows.
- Use `fixed_effects_formula='batch + condition'` and
  `reference_levels=['condition::control']` to exercise `estimate-model`.
- Use `reference.gff3` if you want gene names and product annotations in the
  visualization tables and plots.
