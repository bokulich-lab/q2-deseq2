# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Lesser General Public License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import csv
import json

from qiime2.plugin import model, ValidationError


class DESeq2StatsFormat(model.TextFileFormat):
    def _validate_(self, level):
        record_limit = None if level == "max" else 5

        with self.open() as fh:
            reader = csv.reader(fh, delimiter="\t")
            try:
                header = next(reader)
            except StopIteration as exc:
                raise ValidationError("DESeq2 stats file is empty.") from exc

            if "feature_id" not in header:
                raise ValidationError(
                    'DESeq2 stats table must include a "feature_id" header column.'
                )

            for line_no, row in enumerate(reader, start=2):
                if record_limit is not None and line_no > (record_limit + 1):
                    break

                if not row:
                    continue

                if len(row) != len(header):
                    raise ValidationError(
                        f"Line {line_no} has {len(row)} fields, expected {len(header)}."
                    )


DESeq2StatsDirectoryFormat = model.SingleFileDirectoryFormat(
    "DESeq2StatsDirectoryFormat", "deseq2_stats.tsv", DESeq2StatsFormat
)


class DESeq2RunMetadataFormat(model.TextFileFormat):
    def _validate_(self, level):
        with self.open() as fh:
            try:
                metadata = json.load(fh)
            except json.JSONDecodeError as exc:
                raise ValidationError("DESeq2 run metadata is not valid JSON.") from exc

        missing = {"alpha"}.difference(metadata)
        if missing:
            missing_list = ", ".join(sorted(missing))
            raise ValidationError(
                f"DESeq2 run metadata is missing required fields: {missing_list}."
            )

        has_legacy_effect_fields = {"test_level", "reference_level"} <= set(metadata)
        has_model_effect_field = "default_effect_id" in metadata
        if not has_legacy_effect_fields and not has_model_effect_field:
            raise ValidationError(
                "DESeq2 run metadata must include either legacy comparison fields "
                '("test_level" and "reference_level") or "default_effect_id".'
            )

        try:
            alpha = float(metadata["alpha"])
        except (TypeError, ValueError) as exc:
            raise ValidationError("DESeq2 run metadata alpha must be numeric.") from exc

        if not 0 < alpha < 1:
            raise ValidationError("DESeq2 run metadata alpha must be between 0 and 1.")

        if "test" in metadata and metadata["test"] not in {"wald", "lrt", ""}:
            raise ValidationError('DESeq2 run metadata test must be "wald" or "lrt".')


class DESeq2DistanceMatrixFormat(model.TextFileFormat):
    def _validate_(self, level):
        record_limit = None if level == "max" else 5

        with self.open() as fh:
            reader = csv.reader(fh, delimiter="\t")
            try:
                header = next(reader)
            except StopIteration as exc:
                raise ValidationError("DESeq2 sample distance file is empty.") from exc

            if not header or header[0] != "sample_id":
                raise ValidationError(
                    'DESeq2 sample distance table must start with a '
                    '"sample_id" header column.'
                )

            for line_no, row in enumerate(reader, start=2):
                if record_limit is not None and line_no > (record_limit + 1):
                    break

                if not row:
                    continue

                if len(row) != len(header):
                    raise ValidationError(
                        f"Line {line_no} has {len(row)} fields, expected {len(header)}."
                    )


class DESeq2DistanceOrderFormat(model.TextFileFormat):
    def _validate_(self, level):
        with self.open() as fh:
            lines = [line.strip() for line in fh]

        if lines and not any(lines):
            raise ValidationError(
                "DESeq2 sample distance order file must contain at least one sample ID."
            )


class DESeq2SampleMetadataFormat(model.TextFileFormat):
    def _validate_(self, level):
        record_limit = None if level == "max" else 5

        with self.open() as fh:
            reader = csv.reader(fh, delimiter="\t")
            try:
                header = next(reader)
            except StopIteration as exc:
                raise ValidationError("DESeq2 sample metadata file is empty.") from exc

            if not header or header[0] != "sample_id":
                raise ValidationError(
                    'DESeq2 sample metadata table must start with a '
                    '"sample_id" header column.'
                )

            for line_no, row in enumerate(reader, start=2):
                if record_limit is not None and line_no > (record_limit + 1):
                    break

                if not row:
                    continue

                if len(row) != len(header):
                    raise ValidationError(
                        f"Line {line_no} has {len(row)} fields, expected {len(header)}."
                    )


class DESeq2SamplePCAFormat(model.TextFileFormat):
    def _validate_(self, level):
        record_limit = None if level == "max" else 5

        with self.open() as fh:
            reader = csv.reader(fh, delimiter="\t")
            try:
                header = next(reader)
            except StopIteration as exc:
                raise ValidationError("DESeq2 sample PCA file is empty.") from exc

            required = {"sample_id", "PC1", "PC2"}
            missing = required.difference(header)
            if missing:
                missing_list = ", ".join(sorted(missing))
                raise ValidationError(
                    "DESeq2 sample PCA table is missing required columns: "
                    f"{missing_list}."
                )

            for line_no, row in enumerate(reader, start=2):
                if record_limit is not None and line_no > (record_limit + 1):
                    break

                if not row:
                    continue

                if len(row) != len(header):
                    raise ValidationError(
                        f"Line {line_no} has {len(row)} fields, expected {len(header)}."
                    )


class DESeq2CountMatrixHeatmapFormat(model.TextFileFormat):
    def _validate_(self, level):
        record_limit = None if level == "max" else 5

        with self.open() as fh:
            reader = csv.reader(fh, delimiter="\t")
            try:
                header = next(reader)
            except StopIteration as exc:
                raise ValidationError(
                    "DESeq2 count-matrix heatmap file is empty."
                ) from exc

            if not header or header[0] != "feature_id":
                raise ValidationError(
                    'DESeq2 count-matrix heatmap table must start with a '
                    '"feature_id" header column.'
                )

            for line_no, row in enumerate(reader, start=2):
                if record_limit is not None and line_no > (record_limit + 1):
                    break

                if not row:
                    continue

                if len(row) != len(header):
                    raise ValidationError(
                        f"Line {line_no} has {len(row)} fields, expected {len(header)}."
                    )


class DESeq2RunDirectoryFormat(model.DirectoryFormat):
    results = model.File("deseq2_results.tsv", format=DESeq2StatsFormat)
    normalized_counts = model.File("normalized_counts.tsv", format=DESeq2StatsFormat)
    metadata = model.File("metadata.json", format=DESeq2RunMetadataFormat)
    sample_metadata = model.File(
        "sample_metadata.tsv", format=DESeq2SampleMetadataFormat, optional=True
    )
    sample_distances = model.File(
        "sample_distances.tsv", format=DESeq2DistanceMatrixFormat, optional=True
    )
    sample_distance_order = model.File(
        "sample_distance_order.txt", format=DESeq2DistanceOrderFormat, optional=True
    )
    sample_pca = model.File(
        "sample_pca.tsv", format=DESeq2SamplePCAFormat, optional=True
    )
    count_matrix_heatmap = model.File(
        "count_matrix_heatmap.tsv",
        format=DESeq2CountMatrixHeatmapFormat,
        optional=True,
    )
