# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Lesser General Public License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from pathlib import Path

import biom
import numpy as np
import pandas as pd
import qiime2
from q2_deseq2.types import DESeq2RunDirectoryFormat
from q2_deseq2.utils.run_data import _parse_run_results
from q2_types.feature_table import FeatureTable, Frequency
from qiime2.plugin.testing import TestPluginBase


def _load_table_artifact(tsv_path: str) -> qiime2.Artifact:
    df = pd.read_csv(tsv_path, sep="\t", index_col=0)
    table = biom.Table(
        df.values,
        observation_ids=df.index.tolist(),
        sample_ids=df.columns.tolist(),
    )
    return qiime2.Artifact.import_data(FeatureTable[Frequency], table)


class _IntegrationBase(TestPluginBase):
    package = "q2_deseq2.tests"

    def _get_effect_row(self, results, effect_id, feature_id):
        mask = (results["effect_id"] == effect_id) & (
            results["feature_id"] == feature_id
        )
        rows = results.loc[mask]
        self.assertEqual(len(rows), 1, f"Expected 1 row for {effect_id}/{feature_id}")
        return rows.iloc[0]

    def _parse(self, action_result):
        run_data = action_result[1].view(DESeq2RunDirectoryFormat)
        return _parse_run_results(run_data)


class TestRNASeqIntegration(_IntegrationBase):

    def setUp(self):
        super().setUp()
        demo_dir = Path(self.get_data_path("mini-diff-expr-demo"))
        self.table = _load_table_artifact(str(demo_dir / "feature-table.tsv"))
        self.metadata = qiime2.Metadata.load(str(demo_dir / "sample-metadata.tsv"))
        self.condition = self.metadata.get_column("condition")

    def test_run_deseq2_simple_two_contrasts(self):
        run_result, _ = self._parse(
            self.plugin.actions["_estimate"](
                table=self.table,
                condition=self.condition,
                reference_level="control",
                min_total_count=10,
            )
        )
        df = run_result.results

        self.assertEqual(run_result.reference_level, "control")
        effect_ids = sorted(df["effect_id"].unique())
        self.assertEqual(len(effect_ids), 2)
        self.assertIn("contrast::condition::treated_a::control", effect_ids)
        self.assertIn("contrast::condition::treated_b::control", effect_ids)

        self.assertNotIn("gene_low_count", df["feature_id"].values)

        ta_effect = "contrast::condition::treated_a::control"

        row = self._get_effect_row(df, ta_effect, "gene_alpha")
        self.assertGreater(row["log2FoldChange"], 2.0)
        self.assertLess(row["padj"], 0.05)

        row = self._get_effect_row(df, ta_effect, "gene_beta")
        self.assertLess(row["log2FoldChange"], -2.0)
        self.assertLess(row["padj"], 0.05)

        row = self._get_effect_row(df, ta_effect, "gene_shared")
        self.assertGreater(row["log2FoldChange"], 1.0)
        self.assertLess(row["padj"], 0.05)

        row = self._get_effect_row(df, ta_effect, "gene_housekeeping")
        self.assertLess(abs(row["log2FoldChange"]), 1.0)

        row = self._get_effect_row(df, ta_effect, "gene_gamma")
        self.assertLess(abs(row["log2FoldChange"]), 1.5)

        tb_effect = "contrast::condition::treated_b::control"

        row = self._get_effect_row(df, tb_effect, "gene_gamma")
        self.assertGreater(row["log2FoldChange"], 2.0)
        self.assertLess(row["padj"], 0.05)

        row = self._get_effect_row(df, tb_effect, "gene_delta")
        self.assertLess(row["log2FoldChange"], -2.0)
        self.assertLess(row["padj"], 0.05)

        row = self._get_effect_row(df, tb_effect, "gene_shared")
        self.assertGreater(row["log2FoldChange"], 1.0)
        self.assertLess(row["padj"], 0.05)

        row = self._get_effect_row(df, tb_effect, "gene_alpha")
        self.assertLess(abs(row["log2FoldChange"]), 1.5)

    def test_run_deseq2_simple_result_structure(self):
        run_result, _ = self._parse(
            self.plugin.actions["_estimate"](
                table=self.table,
                condition=self.condition,
                reference_level="control",
                min_total_count=10,
            )
        )
        df = run_result.results

        for col in ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]:
            self.assertIn(col, df.columns)

        self.assertTrue((df["baseMean"] > 0).all())

        pvalues = df["pvalue"].dropna()
        self.assertTrue((pvalues >= 0).all())
        self.assertTrue((pvalues <= 1).all())

        both = df.dropna(subset=["pvalue", "padj"])
        self.assertTrue((both["padj"] >= both["pvalue"] - 1e-10).all())

        nc = run_result.normalized_counts
        self.assertIn("feature_id", nc.columns)
        self.assertEqual(nc.shape[0], 8)
        for sample in ["CTRL_1", "CTRL_2", "TRTA_1", "TRTA_2", "TRTB_1", "TRTB_2"]:
            self.assertIn(sample, nc.columns)
        numeric_cols = [c for c in nc.columns if c != "feature_id"]
        self.assertTrue((nc[numeric_cols] >= 0).all().all())

        sdm = run_result.sample_distance_matrix
        self.assertEqual(sdm.shape, (6, 6))
        np.testing.assert_array_almost_equal(sdm.values, sdm.values.T)
        np.testing.assert_array_almost_equal(np.diag(sdm.values), 0.0)

        pca = run_result.sample_pca_scores
        self.assertEqual(pca.shape, (6, 2))
        self.assertIn("PC1", pca.columns)
        self.assertIn("PC2", pca.columns)
        pv = run_result.sample_pca_percent_variance
        self.assertEqual(len(pv), 2)
        self.assertGreater(pv[0], 0)
        self.assertGreater(pv[0] + pv[1], 50)

        hm = run_result.count_matrix_heatmap
        self.assertGreater(hm.shape[0], 0)
        self.assertEqual(hm.shape[1], 6)

    def test_run_deseq2_model_with_batch_correction(self):
        run_result, _ = self._parse(
            self.plugin.actions["_estimate_model"](
                table=self.table,
                metadata=self.metadata,
                fixed_effects_formula="batch + condition",
                reference_levels=["condition::control"],
                min_total_count=10,
            )
        )
        df = run_result.results

        self.assertEqual(run_result.fixed_effects_formula, "batch + condition")
        self.assertGreater(len(df["effect_id"].unique()), 0)

        condition_effects = df[df["effect_id"].str.contains("condition")]
        alpha_rows = condition_effects[condition_effects["feature_id"] == "gene_alpha"]
        if not alpha_rows.empty:
            ta_rows = alpha_rows[alpha_rows["effect_id"].str.contains("treated_a")]
            if not ta_rows.empty:
                self.assertGreater(abs(ta_rows.iloc[0]["log2FoldChange"]), 2.0)


class TestMicrobiomeIntegration(_IntegrationBase):

    def setUp(self):
        super().setUp()
        demo_dir = Path(self.get_data_path("mini-diff-abund-demo"))
        self.table = _load_table_artifact(str(demo_dir / "feature-table.tsv"))
        self.metadata = qiime2.Metadata.load(str(demo_dir / "sample-metadata.tsv"))
        self.condition = self.metadata.get_column("condition")

    def _run_skin_vs_gut(self):
        run_result, _ = self._parse(
            self.plugin.actions["_estimate"](
                table=self.table,
                condition=self.condition,
                reference_level="gut",
                min_total_count=5,
            )
        )
        return run_result

    def test_sparse_microbiome_biological_signal(self):
        run_result = self._run_skin_vs_gut()
        df = run_result.results
        effect = "contrast::condition::skin::gut"

        self.assertEqual(run_result.reference_level, "gut")
        self.assertEqual(run_result.test_level, "skin")

        self.assertNotIn("asv_low_count", df["feature_id"].values)

        row = self._get_effect_row(df, effect, "asv_skin_dominant")
        self.assertGreater(row["log2FoldChange"], 2.0)
        self.assertLess(row["padj"], 0.05)

        row = self._get_effect_row(df, effect, "asv_gut_dominant")
        self.assertLess(row["log2FoldChange"], -2.0)
        self.assertLess(row["padj"], 0.05)

        row = self._get_effect_row(df, effect, "asv_skin_enriched")
        self.assertGreater(row["log2FoldChange"], 1.0)

        row = self._get_effect_row(df, effect, "asv_gut_enriched")
        self.assertLess(row["log2FoldChange"], -1.0)

        row = self._get_effect_row(df, effect, "asv_ubiquitous")
        self.assertLess(abs(row["log2FoldChange"]), 1.0)

        row = self._get_effect_row(df, effect, "asv_shared_high")
        self.assertLess(abs(row["log2FoldChange"]), 1.5)

    def test_sparse_microbiome_result_structure(self):
        run_result = self._run_skin_vs_gut()
        df = run_result.results

        effect_ids = df["effect_id"].unique()
        self.assertEqual(len(effect_ids), 1)
        self.assertEqual(effect_ids[0], "contrast::condition::skin::gut")

        for col in ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]:
            self.assertIn(col, df.columns)

        feature_ids = sorted(df["feature_id"].unique())
        self.assertEqual(len(feature_ids), 10)
        self.assertNotIn("asv_low_count", feature_ids)
        self.assertNotIn("asv_zero_heavy", feature_ids)

        self.assertTrue((df["baseMean"] > 0).all())

        pvalues = df["pvalue"].dropna()
        self.assertTrue((pvalues >= 0).all())
        self.assertTrue((pvalues <= 1).all())

        both = df.dropna(subset=["pvalue", "padj"])
        self.assertTrue((both["padj"] >= both["pvalue"] - 1e-10).all())

        nc = run_result.normalized_counts
        self.assertIn("feature_id", nc.columns)
        self.assertEqual(nc.shape[0], 10)
        for sample in ["GUT_1", "GUT_2", "GUT_3", "SKIN_1", "SKIN_2", "SKIN_3"]:
            self.assertIn(sample, nc.columns)
        numeric_cols = [c for c in nc.columns if c != "feature_id"]
        self.assertTrue((nc[numeric_cols] >= 0).all().all())

        sdm = run_result.sample_distance_matrix
        self.assertEqual(sdm.shape, (6, 6))
        np.testing.assert_array_almost_equal(sdm.values, sdm.values.T)
        np.testing.assert_array_almost_equal(np.diag(sdm.values), 0.0)

        pca = run_result.sample_pca_scores
        self.assertEqual(pca.shape, (6, 2))
        self.assertIn("PC1", pca.columns)
        self.assertIn("PC2", pca.columns)
        pv = run_result.sample_pca_percent_variance
        self.assertEqual(len(pv), 2)
        self.assertGreater(pv[0], 0)

        hm = run_result.count_matrix_heatmap
        self.assertGreater(hm.shape[0], 0)
        self.assertEqual(hm.shape[1], 6)
