# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
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

        required_fields = {"test_level", "reference_level", "alpha"}
        missing = required_fields.difference(metadata)
        if missing:
            missing_list = ", ".join(sorted(missing))
            raise ValidationError(
                f"DESeq2 run metadata is missing required fields: {missing_list}."
            )

        try:
            alpha = float(metadata["alpha"])
        except (TypeError, ValueError) as exc:
            raise ValidationError("DESeq2 run metadata alpha must be numeric.") from exc

        if not 0 < alpha < 1:
            raise ValidationError("DESeq2 run metadata alpha must be between 0 and 1.")


class DESeq2ImageFormat(model.BinaryFileFormat):
    PNG_HEADER = b"\x89PNG\r\n\x1a\n"

    def _validate_(self, level): ...


class DESeq2RunDirectoryFormat(model.DirectoryFormat):
    results = model.File("deseq2_results.tsv", format=DESeq2StatsFormat)
    normalized_counts = model.File("normalized_counts.tsv", format=DESeq2StatsFormat)
    ma_plot = model.File("ma_plot.png", format=DESeq2ImageFormat)
    volcano_plot = model.File("volcano_plot.png", format=DESeq2ImageFormat)
    metadata = model.File("metadata.json", format=DESeq2RunMetadataFormat)
