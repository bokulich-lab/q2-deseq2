# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import csv

from qiime2.plugin import model, ValidationError


class DESeq2StatsFormat(model.TextFileFormat):
    def _validate_(self, level):
        record_limit = None if level == 'max' else 5

        with self.open() as fh:
            reader = csv.reader(fh, delimiter='\t')
            try:
                header = next(reader)
            except StopIteration as exc:
                raise ValidationError('DESeq2 stats file is empty.') from exc

            if 'feature_id' not in header:
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
                        f'Line {line_no} has {len(row)} fields, expected {len(header)}.'
                    )


DESeq2StatsDirectoryFormat = model.SingleFileDirectoryFormat(
    'DESeq2StatsDirectoryFormat', 'deseq2_stats.tsv', DESeq2StatsFormat
)
