# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Lesser General Public License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re

import biom
import pandas as pd
from pandas.api.types import is_numeric_dtype
from rachis import Metadata

_FORMULA_ALLOWED_CHARS = re.compile(r"^[A-Za-z0-9_:+*() \t]+$")
_FORMULA_TOKEN_RE = re.compile(r"\b[A-Za-z_][A-Za-z0-9_]*\b")
_COEF_SPEC_RE = re.compile(r"^coef::(.+)$")
_CONTRAST_SPEC_RE = re.compile(r"^contrast::([^:]+)::([^:]+)::([^:]+)$")
_SIMPLE_SPEC_RE = re.compile(
    r"^simple::([^:]+)::([^:]+)::([^:]+)\|within::([^:]+)::([^:]+)$"
)


def _filter_counts(counts: pd.DataFrame, min_total_count: int) -> pd.DataFrame:
    """Remove features whose total count across all samples is below *min_total_count*.

    Raises ``ValueError`` if no features remain after filtering.

    Example::

        >>> import pandas as pd
        >>> df = pd.DataFrame({"s1": [1, 10], "s2": [1, 5]}, index=["g1", "g2"])
        >>> _filter_counts(df, min_total_count=5)
            s1  s2
        g2  10   5
    """
    filtered = counts.loc[counts.sum(axis=1) >= min_total_count]
    if filtered.empty:
        raise ValueError(
            "No genes remain after filtering. Lower min_total_count or "
            "provide a denser feature table."
        )
    return filtered


def _collect_matching_samples(
    sample_ids: list[str], metadata_index: pd.Index, context: str
) -> list[str]:
    """Return the subset of *sample_ids* that are present in *metadata_index*.

    Raises ``ValueError`` if fewer than two samples match, including the
    *context* string in the error message to identify which metadata source
    failed to overlap.
    """
    matched_samples = [
        sample_id for sample_id in sample_ids if sample_id in metadata_index
    ]
    if len(matched_samples) < 2:
        raise ValueError(
            f"At least two samples must overlap between the feature table "
            f"and {context}."
        )
    return matched_samples


def _normalize_formula(formula: str, parameter_name: str) -> str:
    """Strip and validate an R-style model formula string.

    Allowed characters: letters, digits, ``_``, ``:``, ``+``, ``*``,
    parentheses, and whitespace.  Raises ``ValueError`` for empty or
    disallowed-character input, naming *parameter_name* in the message.

    Example::

        >>> _normalize_formula("  condition + batch  ", "fixed_effects_formula")
        'condition + batch'
    """
    normalized = str(formula).strip()
    if not normalized:
        raise ValueError(f"{parameter_name} is required.")
    if not _FORMULA_ALLOWED_CHARS.fullmatch(normalized):
        raise ValueError(
            f"{parameter_name} may only contain metadata column names, "
            f"spaces, '+', ':', '*', and parentheses."
        )
    return normalized


def _extract_formula_columns(formula: str, parameter_name: str) -> list[str]:
    """Return the unique metadata column names referenced in *formula*.

    Column names are identified by the token regex ``\\b[A-Za-z_][A-Za-z0-9_]*\\b``
    applied to the normalised formula.  Raises ``ValueError`` if no column
    names are found.

    Example::

        >>> _extract_formula_columns("group + batch", "fixed_effects_formula")
        ['group', 'batch']
    """
    normalized = _normalize_formula(formula, parameter_name)
    columns = list(dict.fromkeys(_FORMULA_TOKEN_RE.findall(normalized)))
    if not columns:
        raise ValueError(
            f"{parameter_name} must reference at least one metadata column."
        )
    return columns


def _parse_reference_level_spec(reference_level_spec: str) -> tuple[str, str]:
    """Parse a ``"column::level"`` reference-level specification.

    Returns ``(column, level)``.  Raises ``ValueError`` when the separator
    ``::`` is absent or either component is empty.

    Example::

        >>> _parse_reference_level_spec("treatment::control")
        ('treatment', 'control')
    """
    column, separator, level = str(reference_level_spec).strip().partition("::")
    if separator != "::" or not column or not level:
        raise ValueError("reference_levels entries must use the form 'column::level'.")
    return column, level


def _normalize_reference_levels(reference_levels: list[str] | None) -> list[str]:
    """Validate and deduplicate a list of ``"column::level"`` reference-level specs.

    Skips blank entries.  Raises ``ValueError`` if the same column appears
    more than once.

    Example::

        >>> _normalize_reference_levels(["treatment::control", "batch::A"])
        ['treatment::control', 'batch::A']
    """
    normalized = []
    seen_columns = set()
    for raw_spec in reference_levels or []:
        spec = str(raw_spec).strip()
        if not spec:
            continue
        column, level = _parse_reference_level_spec(spec)
        if column in seen_columns:
            raise ValueError(
                f"Duplicate reference_levels entry provided for metadata column "
                f'"{column}".'
            )
        normalized.append(f"{column}::{level}")
        seen_columns.add(column)
    return normalized


def _validate_effect_specs(effect_specs: list[str] | None, test: str) -> list[str]:
    """Validate a list of effect specification strings.

    Each non-empty entry must match one of:

    * ``coef::<resultsName>``
    * ``contrast::<factor>::<numerator>::<denominator>``
    * ``simple::<factor>::<numerator>::<denominator>|within::<factor>::<level>``

    Raises ``ValueError`` for malformed entries, or if any specs are provided
    when *test* is ``"lrt"``.
    """
    normalized = []
    for raw_spec in effect_specs or []:
        spec = str(raw_spec).strip()
        if not spec:
            continue
        if (
            _COEF_SPEC_RE.match(spec) is None
            and _CONTRAST_SPEC_RE.match(spec) is None
            and _SIMPLE_SPEC_RE.match(spec) is None
        ):
            raise ValueError(
                "effect_specs entries must use one of: "
                "'coef::<resultsName>', "
                "'contrast::<factor>::<numerator>::<denominator>', or "
                "'simple::<factor>::<numerator>::<denominator>|"
                "within::<factor>::<level>'."
            )
        normalized.append(spec)

    if test == "lrt" and normalized:
        raise ValueError("effect_specs are not supported when test='lrt'.")

    return normalized


def _coerce_metadata_column(column: pd.Series) -> pd.Series:
    """Return *column* cast to numeric (if already numeric dtype) or to ``str``.

    Numeric columns are coerced strictly via ``pd.to_numeric(..., errors='raise')``;
    all other columns are cast to ``str``.
    """
    if is_numeric_dtype(column):
        return pd.to_numeric(column, errors="raise")
    return column.astype(str)


def _coerce_metadata_frame(metadata: Metadata) -> pd.DataFrame:
    """Convert a QIIME 2 metadata object to a plain DataFrame with string index/columns.

    Calls ``.to_dataframe()`` then ensures both the index and column labels are
    Python strings.
    """
    metadata_df = metadata.to_dataframe()
    metadata_df.index = metadata_df.index.map(str)
    metadata_df.columns = metadata_df.columns.map(str)
    return metadata_df


def _prepare_inputs(
    table: biom.Table,
    condition,
    min_total_count: int,
    reference_level: str,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], str]:
    """Prepare count table and coldata for a simple two-condition DESeq2 run.

    Intersects the feature table with the condition metadata, filters
    low-count features, infers the reference level when not supplied, and
    returns the matched data ready to be passed to the R script.

    Returns:
        ``(counts, coldata, comparison_levels, reference_level)``

        * *counts*: genes × matched-samples count DataFrame.
        * *coldata*: matched-samples × 1 DataFrame with a ``"condition"`` column.
        * *comparison_levels*: non-reference condition levels, sorted.
        * *reference_level*: the baseline level used for all contrasts.
    """
    counts = table.to_dataframe(dense=True).round().astype(int)
    sample_ids = list(table.ids(axis="sample"))

    metadata = condition.to_series().dropna().astype(str)
    matched_samples = _collect_matching_samples(
        sample_ids, metadata.index, "condition metadata"
    )

    counts = counts.loc[:, matched_samples]
    metadata = metadata.loc[matched_samples]

    if metadata.nunique() < 2:
        raise ValueError("Condition metadata must contain at least two unique levels.")

    counts = _filter_counts(counts, min_total_count)

    levels = sorted(metadata.unique().tolist())
    if not reference_level:
        if len(levels) != 2:
            raise ValueError(
                "Condition metadata has more than two levels. "
                "Set reference_level to define the baseline for all pairwise contrasts."
            )
        reference_level = levels[0]
    elif reference_level not in levels:
        raise ValueError(
            f'reference_level "{reference_level}" is not present in the '
            f"condition metadata."
        )

    comparison_levels = [level for level in levels if level != reference_level]
    if not comparison_levels:
        raise ValueError(
            "Condition metadata must contain at least one non-reference level."
        )

    coldata = pd.DataFrame({"condition": metadata})
    return counts, coldata, comparison_levels, reference_level


def _prepare_model_inputs(
    table: biom.Table,
    metadata: Metadata,
    fixed_effects_formula: str,
    min_total_count: int,
    reference_levels: list[str] | None = None,
    reduced_formula: str = "",
) -> tuple[pd.DataFrame, pd.DataFrame, str, list[str], str]:
    """Prepare inputs for a general fixed-effects DESeq2 model.

    Validates the formula(s), checks that all referenced metadata columns
    exist, intersects the feature table with the metadata, drops samples with
    missing covariate values, and filters low-count features.

    Returns:
        ``(counts, coldata, normalized_formula, normalized_reference_levels,
        normalized_reduced_formula)``

        * *counts*: genes × matched-samples count DataFrame.
        * *coldata*: matched-samples × model-columns metadata DataFrame.
        * *normalized_formula*: validated ``fixed_effects_formula`` string.
        * *normalized_reference_levels*: validated ``"column::level"`` list.
        * *normalized_reduced_formula*: validated reduced formula (LRT), or ``""``.
    """
    normalized_formula = _normalize_formula(
        fixed_effects_formula, parameter_name="fixed_effects_formula"
    )
    normalized_reduced_formula = (
        _normalize_formula(reduced_formula, parameter_name="reduced_formula")
        if str(reduced_formula).strip()
        else ""
    )

    metadata_df = _coerce_metadata_frame(metadata)
    sample_ids = list(table.ids(axis="sample"))
    counts = table.to_dataframe(dense=True).round().astype(int)
    referenced_columns = _extract_formula_columns(
        normalized_formula, parameter_name="fixed_effects_formula"
    )
    if normalized_reduced_formula:
        referenced_columns.extend(
            _extract_formula_columns(
                normalized_reduced_formula, parameter_name="reduced_formula"
            )
        )
    referenced_columns = list(dict.fromkeys(referenced_columns))

    missing_columns = [
        column for column in referenced_columns if column not in metadata_df.columns
    ]
    if missing_columns:
        missing_display = ", ".join(sorted(missing_columns))
        raise ValueError(
            "The metadata is missing columns required by the model formula: "
            f"{missing_display}."
        )

    normalized_reference_levels = _normalize_reference_levels(reference_levels)
    for reference_level_spec in normalized_reference_levels:
        column, _ = _parse_reference_level_spec(reference_level_spec)
        if column not in referenced_columns:
            raise ValueError(
                f'reference_levels entry "{reference_level_spec}" refers to '
                f'"{column}", but that column is not used in the fitted model.'
            )

    matched_samples = _collect_matching_samples(
        sample_ids, metadata_df.index, "metadata"
    )
    coldata = metadata_df.loc[matched_samples, referenced_columns].copy()
    counts = counts.loc[:, matched_samples]

    coldata = coldata.loc[~coldata.isna().any(axis=1)].copy()
    if coldata.shape[0] < 2:
        raise ValueError(
            "At least two samples with complete metadata are required "
            "after dropping missing values."
        )

    counts = counts.loc[:, coldata.index]
    counts = _filter_counts(counts, min_total_count)
    for column in coldata.columns:
        coldata[column] = _coerce_metadata_column(coldata[column])

    for reference_level_spec in normalized_reference_levels:
        column, level = _parse_reference_level_spec(reference_level_spec)
        if is_numeric_dtype(coldata[column]):
            raise ValueError(
                f'reference_levels entry "{reference_level_spec}" refers to numeric '
                f'metadata column "{column}".'
            )
        observed_levels = set(coldata[column].astype(str))
        if level not in observed_levels:
            raise ValueError(
                f'reference_levels entry "{reference_level_spec}" refers to level '
                f'"{level}", which is not present in metadata column "{column}".'
            )

    return (
        counts,
        coldata,
        normalized_formula,
        normalized_reference_levels,
        normalized_reduced_formula,
    )
