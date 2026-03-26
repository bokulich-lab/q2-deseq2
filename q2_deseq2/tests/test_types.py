import shutil
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd
from pandas.testing import assert_frame_equal
from q2_types.feature_data import FeatureData
from qiime2.plugin import ValidationError
from qiime2.plugin.testing import TestPluginBase

from q2_deseq2.types import (
    DESeq2Run,
    DESeq2RunDirectoryFormat,
    DESeq2RunMetadataFormat,
    DESeq2Stats,
    DESeq2StatsDirectoryFormat,
    DESeq2StatsFormat,
)


class TestFormats(TestPluginBase):
    package = "q2_deseq2.tests"

    def test_deseq2_stats_format_validate(self):
        filepath = self.get_data_path("deseq2-stats-valid.tsv")
        format = DESeq2StatsFormat(filepath, mode="r")

        format.validate("min")
        format.validate("max")

    def test_deseq2_stats_format_rejects_missing_feature_id(self):
        filepath = self.get_data_path("deseq2-stats-missing-feature-id.tsv")
        format = DESeq2StatsFormat(filepath, mode="r")

        with self.assertRaisesRegex(ValidationError, 'include a "feature_id"'):
            format.validate()

    def test_deseq2_stats_format_rejects_mismatched_fields(self):
        filepath = self.get_data_path("deseq2-stats-row-mismatch.tsv")
        format = DESeq2StatsFormat(filepath, mode="r")

        with self.assertRaisesRegex(
            ValidationError, r"Line 3 has 2 fields, expected 3\."
        ):
            format.validate()

    def test_deseq2_stats_format_min_validation_skips_rows_past_limit(self):
        filepath = self.get_data_path("deseq2-stats-row-mismatch-after-min.tsv")
        format = DESeq2StatsFormat(filepath, mode="r")

        format.validate("min")
        with self.assertRaisesRegex(
            ValidationError, r"Line 7 has 2 fields, expected 3\."
        ):
            format.validate("max")

    def test_deseq2_run_metadata_format_validate_legacy_effect_fields(self):
        filepath = self.get_data_path("deseq2-run-metadata-legacy.json")
        format = DESeq2RunMetadataFormat(filepath, mode="r")

        format.validate()

    def test_deseq2_run_metadata_format_validate_default_effect_id(self):
        filepath = self.get_data_path("deseq2-run-metadata-default-effect.json")
        format = DESeq2RunMetadataFormat(filepath, mode="r")

        format.validate()

    def test_deseq2_run_metadata_format_rejects_invalid_json(self):
        filepath = self.get_data_path("deseq2-run-metadata-invalid.json")
        format = DESeq2RunMetadataFormat(filepath, mode="r")

        with self.assertRaisesRegex(ValidationError, "is not valid JSON"):
            format.validate()

    def test_deseq2_run_metadata_format_rejects_missing_alpha(self):
        filepath = self.get_data_path("deseq2-run-metadata-missing-alpha.json")
        format = DESeq2RunMetadataFormat(filepath, mode="r")

        with self.assertRaisesRegex(ValidationError, "missing required fields: alpha"):
            format.validate()

    def test_deseq2_run_metadata_format_rejects_missing_effect_fields(self):
        filepath = self.get_data_path("deseq2-run-metadata-missing-effect-fields.json")
        format = DESeq2RunMetadataFormat(filepath, mode="r")

        with self.assertRaisesRegex(
            ValidationError, 'must include either legacy comparison fields'
        ):
            format.validate()

    def test_deseq2_run_metadata_format_rejects_non_numeric_alpha(self):
        filepath = self.get_data_path("deseq2-run-metadata-nonnumeric-alpha.json")
        format = DESeq2RunMetadataFormat(filepath, mode="r")

        with self.assertRaisesRegex(ValidationError, "alpha must be numeric"):
            format.validate()

    def test_deseq2_run_metadata_format_rejects_out_of_range_alpha(self):
        filepath = self.get_data_path("deseq2-run-metadata-alpha-out-of-range.json")
        format = DESeq2RunMetadataFormat(filepath, mode="r")

        with self.assertRaisesRegex(
            ValidationError, "alpha must be between 0 and 1"
        ):
            format.validate()

    def test_deseq2_run_metadata_format_rejects_unknown_test(self):
        filepath = self.get_data_path("deseq2-run-metadata-invalid-test.json")
        format = DESeq2RunMetadataFormat(filepath, mode="r")

        with self.assertRaisesRegex(
            ValidationError, r'test must be "wald" or "lrt"'
        ):
            format.validate()

    def test_deseq2_run_directory_format_validate(self):
        filepath = self.get_data_path("deseq2-run-valid")
        format = DESeq2RunDirectoryFormat(filepath, mode="r")

        format.validate("min")
        format.validate("max")

    def test_deseq2_run_directory_format_accepts_extended_run_outputs(self):
        source = Path(self.get_data_path("deseq2-run-valid"))

        with TemporaryDirectory() as temp_dir:
            workdir = Path(temp_dir) / "deseq2-run-valid"
            shutil.copytree(source, workdir)
            (workdir / "sample_metadata.tsv").write_text(
                "sample_id\tcondition\tbatch\n"
                "Sample1\tcontrol\tA\n"
                "Sample2\ttreated\tB\n",
                encoding="utf-8",
            )
            (workdir / "sample_distances.tsv").write_text(
                "sample_id\tSample1\tSample2\n"
                "Sample1\t0.0\t1.5\n"
                "Sample2\t1.5\t0.0\n",
                encoding="utf-8",
            )
            (workdir / "sample_distance_order.txt").write_text(
                "Sample2\nSample1\n", encoding="utf-8"
            )
            (workdir / "sample_pca.tsv").write_text(
                "sample_id\tPC1\tPC2\n"
                "Sample1\t-2.1\t0.4\n"
                "Sample2\t2.1\t-0.4\n",
                encoding="utf-8",
            )
            (workdir / "count_matrix_heatmap.tsv").write_text(
                "feature_id\tSample2\tSample1\n"
                "GG_OTU_2\t7.5\t6.8\n"
                "GG_OTU_1\t5.1\t4.6\n",
                encoding="utf-8",
            )

            format = DESeq2RunDirectoryFormat(str(workdir), mode="r")
            format.validate("min")
            format.validate("max")

    def test_deseq2_run_directory_format_rejects_removed_png_outputs(self):
        source = Path(self.get_data_path("deseq2-run-valid"))

        with TemporaryDirectory() as temp_dir:
            workdir = Path(temp_dir) / "deseq2-run-valid"
            shutil.copytree(source, workdir)
            (workdir / "ma_plot.png").write_bytes(b"\x89PNG\r\n\x1a\n")

            format = DESeq2RunDirectoryFormat(str(workdir), mode="r")

            with self.assertRaisesRegex(ValidationError, "Unrecognized file"):
                format.validate()

    def test_deseq2_run_directory_format_rejects_invalid_sample_metadata(self):
        source = Path(self.get_data_path("deseq2-run-valid"))

        with TemporaryDirectory() as temp_dir:
            workdir = Path(temp_dir) / "deseq2-run-valid"
            shutil.copytree(source, workdir)
            (workdir / "sample_metadata.tsv").write_text(
                "condition\tbatch\ncontrol\tA\n",
                encoding="utf-8",
            )

            format = DESeq2RunDirectoryFormat(str(workdir), mode="r")

            with self.assertRaisesRegex(ValidationError, 'start with a "sample_id"'):
                format.validate()

    def test_deseq2_run_directory_format_rejects_invalid_metadata(self):
        filepath = self.get_data_path("deseq2-run-invalid-metadata")
        format = DESeq2RunDirectoryFormat(filepath, mode="r")

        with self.assertRaisesRegex(
            ValidationError, r'test must be "wald" or "lrt"'
        ):
            format.validate()


class TestTypes(TestPluginBase):
    package = "q2_deseq2.tests"

    def test_deseq2_stats_semantic_type_registration(self):
        self.assertRegisteredSemanticType(DESeq2Stats)

    def test_deseq2_run_semantic_type_registration(self):
        self.assertRegisteredSemanticType(DESeq2Run)

    def test_deseq2_stats_to_format_registration(self):
        self.assertSemanticTypeRegisteredToFormat(
            FeatureData[DESeq2Stats], DESeq2StatsDirectoryFormat
        )

    def test_deseq2_run_to_format_registration(self):
        self.assertSemanticTypeRegisteredToFormat(
            DESeq2Run, DESeq2RunDirectoryFormat
        )


class TestTransformers(TestPluginBase):
    package = "q2_deseq2.tests"

    def setUp(self):
        super().setUp()
        filepath = self.get_data_path("deseq2-stats-valid.tsv")
        self.stats_df = pd.read_csv(filepath, sep="\t")

    def test_dataframe_to_deseq2_stats(self):
        transformer = self.get_transformer(pd.DataFrame, DESeq2StatsFormat)
        observed = transformer(self.stats_df)

        self.assertIsInstance(observed, DESeq2StatsFormat)
        observed.validate("max")
        observed = pd.read_csv(str(observed), sep="\t")
        assert_frame_equal(observed, self.stats_df)

    def test_deseq2_stats_to_dataframe(self):
        _, observed = self.transform_format(
            DESeq2StatsFormat, pd.DataFrame, "deseq2-stats-valid.tsv"
        )

        self.assertIsInstance(observed, pd.DataFrame)
        assert_frame_equal(observed, self.stats_df)
