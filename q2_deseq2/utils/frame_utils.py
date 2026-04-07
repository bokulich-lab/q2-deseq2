# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd


def _first_non_empty_string(value) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def _first_value_from_column(frame: pd.DataFrame, column_name: str) -> str:
    if column_name not in frame.columns:
        return ""
    for value in frame[column_name]:
        normalized = _first_non_empty_string(value)
        if normalized:
            return normalized
    return ""


def _unique_non_empty_values(values) -> tuple[str, ...]:
    ordered_values = []
    seen = set()
    for value in values:
        normalized = _first_non_empty_string(value)
        if normalized and normalized not in seen:
            ordered_values.append(normalized)
            seen.add(normalized)
    return tuple(ordered_values)
