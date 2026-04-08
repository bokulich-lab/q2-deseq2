# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Lesser General Public License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import euclidean_distances


def _compute_normalized_counts(
    counts_df: pd.DataFrame, size_factors: pd.Series
) -> pd.DataFrame:
    """Divide raw counts by per-sample DESeq2 size factors.

    Returns a DataFrame with a leading ``feature_id`` column followed by one
    column per sample.  Raises ``RuntimeError`` if any sample is missing a
    size factor.
    """
    aligned_size_factors = size_factors.reindex(counts_df.columns.map(str))
    missing_samples = aligned_size_factors[aligned_size_factors.isna()].index.tolist()
    if missing_samples:
        missing_display = ", ".join(missing_samples)
        raise RuntimeError(
            "DESeq2 completed but did not return size factors for sample(s): "
            f"{missing_display}."
        )

    normalized = counts_df.div(aligned_size_factors.astype(float), axis=1)
    normalized.index = normalized.index.map(str)
    normalized.columns = normalized.columns.map(str)
    normalized.insert(0, "feature_id", normalized.index)
    normalized = normalized.reset_index(drop=True)
    normalized.columns.name = None
    return normalized


def _compute_sample_distance_matrix(vst_counts: pd.DataFrame) -> pd.DataFrame:
    """Compute the pairwise Euclidean distance matrix across samples.

    Uses VST-normalised counts (genes × samples) via
    ``sklearn.metrics.pairwise.euclidean_distances``.  The returned DataFrame
    is labelled by sample ID on both axes.
    """
    sample_ids = vst_counts.columns.map(str)
    distances = euclidean_distances(vst_counts.to_numpy(dtype=float).T)
    return pd.DataFrame(distances, index=sample_ids, columns=sample_ids)


def _compute_sample_distance_order(
    sample_distance_matrix: pd.DataFrame,
) -> tuple[str, ...]:
    """Return a sample ordering produced by complete-linkage hierarchical clustering.

    Uses ``scipy.cluster.hierarchy.linkage`` (``method="complete"``) on the
    condensed form of the pre-computed distance matrix, then extracts the
    dendrogram leaf order via ``leaves_list``.

    Returns the original sample tuple unchanged when fewer than two samples
    are present.
    """
    sample_ids = tuple(sample_distance_matrix.index.map(str))
    if len(sample_ids) < 2:
        return sample_ids
    condensed = squareform(sample_distance_matrix.to_numpy(dtype=float))
    order = leaves_list(linkage(condensed, method="complete"))
    return tuple(sample_ids[i] for i in order)


def _compute_sample_pca(
    vst_counts: pd.DataFrame,
) -> tuple[pd.DataFrame, tuple[float, float]]:
    """Project samples onto the first two principal components of VST counts.

    Up to the 500 highest-variance genes are used (matching the DESeq2 R
    convention). PCA is performed by ``sklearn.decomposition.PCA``, which
    handles mean-centring internally.

    Returns:
        scores: DataFrame (samples × 2) with columns ``PC1`` and ``PC2``.
        percent_variance: ``(pct_PC1, pct_PC2)`` as floats in ``[0, 100]``.
            Empty tuple when no variance can be computed.
    """
    sample_ids = vst_counts.columns.map(str)
    sample_count = len(sample_ids)
    if sample_count == 0:
        return pd.DataFrame(columns=["PC1", "PC2"]), ()

    ntop = min(500, vst_counts.shape[0])
    if ntop == 0:
        empty_scores = pd.DataFrame(
            {"PC1": np.zeros(sample_count), "PC2": np.zeros(sample_count)},
            index=sample_ids,
        )
        empty_scores.index.name = None
        return empty_scores, ()

    matrix = vst_counts.to_numpy(dtype=float)
    ddof = 1 if sample_count > 1 else 0
    row_variances = np.var(matrix, axis=1, ddof=ddof)
    selected_rows = np.argsort(-row_variances, kind="mergesort")[:ntop]
    selected_matrix = matrix[selected_rows, :].T  # shape: (n_samples, ntop)

    n_components = min(2, sample_count, ntop)
    pca = PCA(n_components=n_components)
    scores = pca.fit_transform(selected_matrix)

    pc1_values = scores[:, 0] if n_components >= 1 else np.zeros(sample_count)
    pc2_values = scores[:, 1] if n_components >= 2 else np.zeros(sample_count)
    percent_variance = pca.explained_variance_ratio_
    percent_pc1 = (
        float(percent_variance[0] * 100) if len(percent_variance) >= 1 else 0.0
    )
    percent_pc2 = (
        float(percent_variance[1] * 100) if len(percent_variance) >= 2 else 0.0
    )

    sample_pca_scores = pd.DataFrame(
        {"PC1": pc1_values, "PC2": pc2_values},
        index=sample_ids,
    )
    sample_pca_scores.index.name = None
    return sample_pca_scores, (percent_pc1, percent_pc2)


def _compute_count_matrix_heatmap(
    vst_counts: pd.DataFrame,
    sample_distance_order: tuple[str, ...],
) -> pd.DataFrame:
    """Select the top 100 highest-variance features and reorder columns by clustering.

    Features are ranked by per-row variance of VST counts.  Columns are
    arranged according to *sample_distance_order* (as returned by
    :func:`_compute_sample_distance_order`), falling back to the original
    column order for any sample absent from the ordering.

    Returns a (features × samples) DataFrame ready to be rendered as a heatmap.
    """
    ordered_sample_ids = [
        sample_id
        for sample_id in sample_distance_order
        if sample_id in vst_counts.columns
    ]
    if not ordered_sample_ids:
        ordered_sample_ids = list(vst_counts.columns.map(str))

    sample_count = vst_counts.shape[1]
    ddof = 1 if sample_count > 1 else 0
    feature_variance = vst_counts.var(axis=1, ddof=ddof)
    top_heatmap_n = min(100, len(feature_variance))
    selected_feature_ids = feature_variance.sort_values(
        ascending=False, kind="mergesort"
    ).index[:top_heatmap_n]
    heatmap = vst_counts.loc[selected_feature_ids, ordered_sample_ids].copy()
    heatmap.index = heatmap.index.map(str)
    heatmap.columns = heatmap.columns.map(str)
    heatmap.index.name = None
    heatmap.columns.name = None
    return heatmap


def _compute_run_analytics(
    counts_df: pd.DataFrame,
    size_factors: pd.Series,
    vst_counts: pd.DataFrame,
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    tuple[str, ...],
    pd.DataFrame,
    tuple[float, float],
    pd.DataFrame,
]:
    """Orchestrate all post-DESeq2 analytics in a single call.

    Validates that every sample and feature in *counts_df* is present in
    *vst_counts*, then delegates to the individual compute helpers.

    Returns a 6-tuple:
        ``(normalized_counts, sample_distance_matrix, sample_distance_order,
        sample_pca_scores, sample_pca_percent_variance, count_matrix_heatmap)``

    Raises ``RuntimeError`` if *vst_counts* is missing any expected sample or
    feature.
    """
    expected_sample_ids = list(counts_df.columns.map(str))
    missing_vst_samples = [
        sample_id
        for sample_id in expected_sample_ids
        if sample_id not in vst_counts.columns
    ]
    if missing_vst_samples:
        missing_display = ", ".join(missing_vst_samples)
        raise RuntimeError(
            "DESeq2 completed but the VST matrix is missing sample(s): "
            f"{missing_display}."
        )

    expected_feature_ids = list(counts_df.index.map(str))
    missing_vst_features = [
        feature_id
        for feature_id in expected_feature_ids
        if feature_id not in vst_counts.index
    ]
    if missing_vst_features:
        missing_display = ", ".join(missing_vst_features[:10])
        if len(missing_vst_features) > 10:
            missing_display += ", ..."
        raise RuntimeError(
            "DESeq2 completed but the VST matrix is missing feature(s): "
            f"{missing_display}."
        )

    vst_counts = vst_counts.loc[expected_feature_ids, expected_sample_ids]
    normalized_counts = _compute_normalized_counts(counts_df, size_factors)
    sample_distance_matrix = _compute_sample_distance_matrix(vst_counts)
    sample_distance_matrix.index.name = None
    sample_distance_matrix.columns.name = None
    sample_distance_order = _compute_sample_distance_order(sample_distance_matrix)
    sample_pca_scores, sample_pca_percent_variance = _compute_sample_pca(vst_counts)
    count_matrix_heatmap = _compute_count_matrix_heatmap(
        vst_counts,
        sample_distance_order,
    )
    return (
        normalized_counts,
        sample_distance_matrix,
        sample_distance_order,
        sample_pca_scores,
        sample_pca_percent_variance,
        count_matrix_heatmap,
    )
