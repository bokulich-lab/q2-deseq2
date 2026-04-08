# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Lesser General Public License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd


def _first_non_empty_string(value) -> str:
    """Return a stripped string for *value*, or ``""`` if it is None/NaN.

    Example::

        >>> _first_non_empty_string("  hello  ")
        'hello'
        >>> _first_non_empty_string(None)
        ''
    """
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def _first_value_from_column(frame: pd.DataFrame, column_name: str) -> str:
    """Return the first non-empty string value from *column_name* in *frame*.

    Returns ``""`` if the column does not exist or every value is empty/NaN.

    Example::

        >>> import pandas as pd
        >>> df = pd.DataFrame({"col": [None, "A", "B"]})
        >>> _first_value_from_column(df, "col")
        'A'
    """
    if column_name not in frame.columns:
        return ""
    for value in frame[column_name]:
        normalized = _first_non_empty_string(value)
        if normalized:
            return normalized
    return ""


def _unique_non_empty_values(values) -> tuple[str, ...]:
    """Return unique non-empty string values from *values*, preserving insertion order.

    Example::

        >>> _unique_non_empty_values(["b", "a", "b", None, ""])
        ('b', 'a')
    """
    ordered_values = []
    seen = set()
    for value in values:
        normalized = _first_non_empty_string(value)
        if normalized and normalized not in seen:
            ordered_values.append(normalized)
            seen.add(normalized)
    return tuple(ordered_values)
