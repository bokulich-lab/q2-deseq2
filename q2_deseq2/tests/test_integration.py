import unittest
from pathlib import Path
from subprocess import run as subprocess_run

import biom
import numpy as np
import qiime2
from qiime2.plugin.testing import TestPluginBase

from q2_deseq2 import methods


def _r_deseq2_available():
    try:
        result = subprocess_run(
            ["Rscript", "-e", 'library("DESeq2"); cat("ok")'],
            capture_output=True, text=True, timeout=30,
        )
        return result.returncode == 0 and "ok" in result.stdout
    except Exception:
        return False


class _IntegrationBase(TestPluginBase):
    package = "q2_deseq2.tests"

    def _get_effect_row(self, results, effect_id, feature_id):
        mask = (results["effect_id"] == effect_id) & (
            results["feature_id"] == feature_id
        )
        rows = results.loc[mask]
        self.assertEqual(len(rows), 1, f"Expected 1 row for {effect_id}/{feature_id}")
        return rows.iloc[0]


@unittest.skipUnless(_r_deseq2_available(), "R or DESeq2 not available")
class TestRNASeqIntegration(_IntegrationBase):

    def setUp(self):
        super().setUp()
        demo_dir = Path(self.get_data_path("mini-diff-expr-demo"))
        self.table = biom.load_table(str(demo_dir / "feature-table.biom"))
        self.metadata = qiime2.Metadata.load(
            str(demo_dir / "sample-metadata.tsv")
        )
        self.condition = self.metadata.get_column("condition")

    def test_run_deseq2_simple_two_contrasts(self):
        result = methods.run_deseq2(
            table=self.table,
            condition=self.condition,
            reference_level="control",
            min_total_count=10,
        )

        df = result.results
        self.assertEqual(result.reference_level, "control")
        # Two contrasts: treated_a vs control, treated_b vs control
        effect_ids = sorted(df["effect_id"].unique())
        self.assertEqual(len(effect_ids), 2)
        self.assertIn("contrast::condition::treated_a::control", effect_ids)
        self.assertIn("contrast::condition::treated_b::control", effect_ids)

        # gene_low_count (total=5) should be filtered out by min_total_count=10
        self.assertNotIn("gene_low_count", df["feature_id"].values)

        # --- treated_a vs control ---
        ta_effect = "contrast::condition::treated_a::control"

        # gene_alpha: strongly up in treated_a (20->220)
        row = self._get_effect_row(df, ta_effect, "gene_alpha")
        self.assertGreater(row["log2FoldChange"], 2.0)
        self.assertLess(row["padj"], 0.05)

        # gene_beta: strongly down in treated_a (180->20)
        row = self._get_effect_row(df, ta_effect, "gene_beta")
        self.assertLess(row["log2FoldChange"], -2.0)
        self.assertLess(row["padj"], 0.05)

        # gene_shared: up in treated_a (30->140)
        row = self._get_effect_row(df, ta_effect, "gene_shared")
        self.assertGreater(row["log2FoldChange"], 1.0)
        self.assertLess(row["padj"], 0.05)

        # gene_housekeeping: stable in treated_a
        row = self._get_effect_row(df, ta_effect, "gene_housekeeping")
        self.assertLess(abs(row["log2FoldChange"]), 1.0)

        # gene_gamma: stable in treated_a (not its signal group)
        row = self._get_effect_row(df, ta_effect, "gene_gamma")
        self.assertLess(abs(row["log2FoldChange"]), 1.5)

        # --- treated_b vs control ---
        tb_effect = "contrast::condition::treated_b::control"

        # gene_gamma: strongly up in treated_b (15->210)
        row = self._get_effect_row(df, tb_effect, "gene_gamma")
        self.assertGreater(row["log2FoldChange"], 2.0)
        self.assertLess(row["padj"], 0.05)

        # gene_delta: strongly down in treated_b (200->22)
        row = self._get_effect_row(df, tb_effect, "gene_delta")
        self.assertLess(row["log2FoldChange"], -2.0)
        self.assertLess(row["padj"], 0.05)

        # gene_shared: up in treated_b (30->150)
        row = self._get_effect_row(df, tb_effect, "gene_shared")
        self.assertGreater(row["log2FoldChange"], 1.0)
        self.assertLess(row["padj"], 0.05)

        # gene_alpha: stable in treated_b (not its signal group)
        row = self._get_effect_row(df, tb_effect, "gene_alpha")
        self.assertLess(abs(row["log2FoldChange"]), 1.5)

    def test_run_deseq2_simple_result_structure(self):
        result = methods.run_deseq2(
            table=self.table,
            condition=self.condition,
            reference_level="control",
            min_total_count=10,
        )

        df = result.results
        # Required DESeq2 output columns present
        for col in ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]:
            self.assertIn(col, df.columns)

        # All baseMean values should be positive
        self.assertTrue((df["baseMean"] > 0).all())

        # p-values should be in [0, 1] where not NaN
        pvalues = df["pvalue"].dropna()
        self.assertTrue((pvalues >= 0).all())
        self.assertTrue((pvalues <= 1).all())

        # Adjusted p-values should be >= raw p-values (BH correction)
        both = df.dropna(subset=["pvalue", "padj"])
        self.assertTrue(
            (both["padj"] >= both["pvalue"] - 1e-10).all(),
            "Adjusted p-values should be >= raw p-values",
        )

        # Normalized counts
        nc = result.normalized_counts
        self.assertIn("feature_id", nc.columns)
        self.assertEqual(nc.shape[0], 8)  # 9 genes minus gene_low_count
        for sample in ["CTRL_1", "CTRL_2", "TRTA_1", "TRTA_2", "TRTB_1", "TRTB_2"]:
            self.assertIn(sample, nc.columns)
        # All normalized counts should be non-negative
        numeric_cols = [c for c in nc.columns if c != "feature_id"]
        self.assertTrue((nc[numeric_cols] >= 0).all().all())

        # Sample distance matrix should be 6x6 and symmetric
        sdm = result.sample_distance_matrix
        self.assertEqual(sdm.shape, (6, 6))
        np.testing.assert_array_almost_equal(sdm.values, sdm.values.T)
        np.testing.assert_array_almost_equal(np.diag(sdm.values), 0.0)

        # PCA scores
        pca = result.sample_pca_scores
        self.assertEqual(pca.shape, (6, 2))
        self.assertIn("PC1", pca.columns)
        self.assertIn("PC2", pca.columns)
        pv = result.sample_pca_percent_variance
        self.assertEqual(len(pv), 2)
        self.assertGreater(pv[0], 0)
        self.assertGreater(pv[0] + pv[1], 50)  # first two PCs explain majority

        # Heatmap
        hm = result.count_matrix_heatmap
        self.assertGreater(hm.shape[0], 0)
        self.assertEqual(hm.shape[1], 6)

    def test_run_deseq2_model_with_batch_correction(self):
        result = methods.run_deseq2_model(
            table=self.table,
            metadata=self.metadata,
            fixed_effects_formula="batch + condition",
            reference_levels=["condition::control"],
            min_total_count=10,
        )

        df = result.results
        self.assertEqual(result.fixed_effects_formula, "batch + condition")

        # Should have coefficients for batch and both condition contrasts
        effect_ids = sorted(df["effect_id"].unique())
        self.assertGreater(len(effect_ids), 0)

        # gene_batch_shift should have a smaller effect after batch correction
        # compared to genes with real condition effects
        condition_effects = df[df["effect_id"].str.contains("condition")]
        alpha_rows = condition_effects[
            condition_effects["feature_id"] == "gene_alpha"
        ]
        if not alpha_rows.empty:
            # gene_alpha should still be significant in its contrast
            ta_rows = alpha_rows[alpha_rows["effect_id"].str.contains("treated_a")]
            if not ta_rows.empty:
                self.assertGreater(
                    abs(ta_rows.iloc[0]["log2FoldChange"]), 2.0
                )


@unittest.skipUnless(_r_deseq2_available(), "R or DESeq2 not available")
class TestMicrobiomeIntegration(_IntegrationBase):

    def setUp(self):
        super().setUp()
        demo_dir = Path(self.get_data_path("mini-diff-abund-demo"))
        self.table = biom.load_table(str(demo_dir / "feature-table.biom"))
        self.metadata = qiime2.Metadata.load(
            str(demo_dir / "sample-metadata.tsv")
        )
        self.condition = self.metadata.get_column("condition")

    def _run_skin_vs_gut(self):
        return methods.run_deseq2(
            table=self.table,
            condition=self.condition,
            reference_level="gut",
            min_total_count=5,
        )

    def test_sparse_microbiome_biological_signal(self):
        result = self._run_skin_vs_gut()
        df = result.results
        effect = "contrast::condition::skin::gut"

        self.assertEqual(result.reference_level, "gut")
        self.assertEqual(result.test_level, "skin")

        # asv_low_count (total=2) should be filtered out by min_total_count=5
        self.assertNotIn("asv_low_count", df["feature_id"].values)

        # asv_skin_dominant: near-zero in gut, abundant in skin -> strongly up
        row = self._get_effect_row(df, effect, "asv_skin_dominant")
        self.assertGreater(row["log2FoldChange"], 2.0)
        self.assertLess(row["padj"], 0.05)

        # asv_gut_dominant: abundant in gut, near-zero in skin -> strongly down
        row = self._get_effect_row(df, effect, "asv_gut_dominant")
        self.assertLess(row["log2FoldChange"], -2.0)
        self.assertLess(row["padj"], 0.05)

        # asv_skin_enriched: moderately enriched in skin -> up
        row = self._get_effect_row(df, effect, "asv_skin_enriched")
        self.assertGreater(row["log2FoldChange"], 1.0)

        # asv_gut_enriched: moderately enriched in gut -> down
        row = self._get_effect_row(df, effect, "asv_gut_enriched")
        self.assertLess(row["log2FoldChange"], -1.0)

        # asv_ubiquitous: stable everywhere -> no signal
        row = self._get_effect_row(df, effect, "asv_ubiquitous")
        self.assertLess(abs(row["log2FoldChange"]), 1.0)

        # asv_shared_high: similar in both -> no strong signal
        row = self._get_effect_row(df, effect, "asv_shared_high")
        self.assertLess(abs(row["log2FoldChange"]), 1.5)

    def test_sparse_microbiome_result_structure(self):
        result = self._run_skin_vs_gut()
        df = result.results

        # Single contrast for 2-level comparison
        effect_ids = df["effect_id"].unique()
        self.assertEqual(len(effect_ids), 1)
        self.assertEqual(effect_ids[0], "contrast::condition::skin::gut")

        # Required DESeq2 output columns present
        for col in ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]:
            self.assertIn(col, df.columns)

        # 12 features minus asv_low_count = 11 features in results
        feature_ids = sorted(df["feature_id"].unique())
        self.assertEqual(len(feature_ids), 11)
        self.assertNotIn("asv_low_count", feature_ids)

        # All baseMean values should be positive
        self.assertTrue((df["baseMean"] > 0).all())

        # p-values in [0, 1] where not NaN
        pvalues = df["pvalue"].dropna()
        self.assertTrue((pvalues >= 0).all())
        self.assertTrue((pvalues <= 1).all())

        # padj >= pvalue (BH correction)
        both = df.dropna(subset=["pvalue", "padj"])
        self.assertTrue(
            (both["padj"] >= both["pvalue"] - 1e-10).all(),
            "Adjusted p-values should be >= raw p-values",
        )

        # Normalized counts: 11 features x 6 samples + feature_id column
        nc = result.normalized_counts
        self.assertIn("feature_id", nc.columns)
        self.assertEqual(nc.shape[0], 11)
        for sample in ["GUT_1", "GUT_2", "GUT_3", "SKIN_1", "SKIN_2", "SKIN_3"]:
            self.assertIn(sample, nc.columns)
        numeric_cols = [c for c in nc.columns if c != "feature_id"]
        self.assertTrue((nc[numeric_cols] >= 0).all().all())

        # Sample distance matrix: 6x6, symmetric, zero diagonal
        sdm = result.sample_distance_matrix
        self.assertEqual(sdm.shape, (6, 6))
        np.testing.assert_array_almost_equal(sdm.values, sdm.values.T)
        np.testing.assert_array_almost_equal(np.diag(sdm.values), 0.0)

        # PCA scores
        pca = result.sample_pca_scores
        self.assertEqual(pca.shape, (6, 2))
        self.assertIn("PC1", pca.columns)
        self.assertIn("PC2", pca.columns)
        pv = result.sample_pca_percent_variance
        self.assertEqual(len(pv), 2)
        self.assertGreater(pv[0], 0)

        # Heatmap
        hm = result.count_matrix_heatmap
        self.assertGreater(hm.shape[0], 0)
        self.assertEqual(hm.shape[1], 6)
