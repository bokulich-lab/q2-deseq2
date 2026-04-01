# Mini Microbiome Demo Fixture

This synthetic dataset simulates a sparse microbiome ASV table for testing
`q2-deseq2` with data that stresses the zero-handling and `poscounts` size
factor estimation.

## Design

- Six samples, three per body site.
- `gut`: `GUT_1`, `GUT_2`, `GUT_3`
- `skin`: `SKIN_1`, `SKIN_2`, `SKIN_3`
- The metadata also includes a `subject` column (paired design).

## Expected Signal (skin vs gut)

- `asv_gut_dominant` is abundant in gut, near-zero in skin (strong down).
- `asv_skin_dominant` is abundant in skin, near-zero in gut (strong up).
- `asv_gut_enriched` is moderately enriched in gut (moderate down).
- `asv_skin_enriched` is moderately enriched in skin (moderate up).
- `asv_shared_high` is present in both at similar levels (no signal).
- `asv_ubiquitous` is abundant and stable everywhere (no signal).
- `asv_rare_gut` and `asv_rare_skin` are very sparse, site-specific features.
- `asv_sparse_shared` and `asv_zero_heavy` are extremely sparse with no clear
  differential signal.
- `asv_low_count` has a total count of 2, so it should be filtered out
  whenever `min_total_count >= 5`.
- `asv_variable` has high within-group variance and no clear signal.

## Sparsity

The table is ~30% zeros (22 out of 72 cells), which is moderate for a
pre-filtered ASV table. Several features have entire rows or near-entire rows
of zeros.

## Files

- `feature-table.tsv` is the human-readable count matrix.
- `feature-table.biom` is the import-ready BIOM version.
- `sample-metadata.tsv` contains `condition` and `subject` columns.
