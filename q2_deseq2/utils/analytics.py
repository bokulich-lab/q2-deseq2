# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Lesser General Public License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd


def _compute_normalized_counts(
    counts_df: pd.DataFrame, size_factors: pd.Series
) -> pd.DataFrame:
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
    sample_ids = vst_counts.columns.map(str)
    sample_matrix = vst_counts.to_numpy(dtype=float).T
    squared_norms = np.sum(sample_matrix * sample_matrix, axis=1, keepdims=True)
    squared_distances = (
        squared_norms + squared_norms.T - 2 * (sample_matrix @ sample_matrix.T)
    )
    squared_distances = np.maximum(squared_distances, 0.0)
    distances = np.sqrt(squared_distances)
    return pd.DataFrame(distances, index=sample_ids, columns=sample_ids)


def _compute_sample_distance_order(
    sample_distance_matrix: pd.DataFrame,
) -> tuple[str, ...]:
    sample_ids = tuple(sample_distance_matrix.index.map(str))
    sample_count = len(sample_ids)
    if sample_count < 2:
        return sample_ids

    distances = sample_distance_matrix.to_numpy(dtype=float)
    clusters = {
        index: {"order": (index,), "height": 0.0} for index in range(sample_count)
    }
    cluster_distances = {
        (left, right): float(distances[left, right])
        for left in range(sample_count)
        for right in range(left + 1, sample_count)
    }

    active_clusters = list(range(sample_count))
    next_cluster_id = sample_count

    while len(active_clusters) > 1:
        best_pair = None
        best_key = None
        active_clusters_sorted = sorted(
            active_clusters, key=lambda cluster_id: clusters[cluster_id]["order"]
        )
        for left_index, left_cluster in enumerate(active_clusters_sorted[:-1]):
            for right_cluster in active_clusters_sorted[left_index + 1 :]:
                pair = (
                    (left_cluster, right_cluster)
                    if left_cluster < right_cluster
                    else (right_cluster, left_cluster)
                )
                candidate_key = (
                    cluster_distances[pair],
                    clusters[left_cluster]["order"],
                    clusters[right_cluster]["order"],
                )
                if best_key is None or candidate_key < best_key:
                    best_key = candidate_key
                    best_pair = pair

        if best_pair is None or best_key is None:
            break

        left_cluster, right_cluster = sorted(
            best_pair,
            key=lambda cluster_id: (
                clusters[cluster_id]["height"],
                clusters[cluster_id]["order"],
            ),
        )
        merged_order = (
            clusters[left_cluster]["order"] + clusters[right_cluster]["order"]
        )
        merged_height = float(best_key[0])
        clusters[next_cluster_id] = {
            "order": merged_order,
            "height": merged_height,
        }

        remaining_clusters = [
            cluster_id for cluster_id in active_clusters if cluster_id not in best_pair
        ]
        for other_cluster in remaining_clusters:
            left_distance = cluster_distances[
                (
                    (left_cluster, other_cluster)
                    if left_cluster < other_cluster
                    else (other_cluster, left_cluster)
                )
            ]
            right_distance = cluster_distances[
                (
                    (right_cluster, other_cluster)
                    if right_cluster < other_cluster
                    else (other_cluster, right_cluster)
                )
            ]
            pair = (
                (other_cluster, next_cluster_id)
                if other_cluster < next_cluster_id
                else (next_cluster_id, other_cluster)
            )
            cluster_distances[pair] = max(left_distance, right_distance)

        active_clusters = remaining_clusters + [next_cluster_id]
        next_cluster_id += 1

    final_order = clusters[active_clusters[0]]["order"]
    return tuple(sample_ids[index] for index in final_order)


def _compute_sample_pca(
    vst_counts: pd.DataFrame,
) -> tuple[pd.DataFrame, tuple[float, float]]:
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
    selected_matrix = matrix[selected_rows, :].T
    centered_matrix = selected_matrix - selected_matrix.mean(axis=0, keepdims=True)

    if centered_matrix.size == 0:
        singular_values = np.array([], dtype=float)
        sample_scores = np.zeros((sample_count, 0), dtype=float)
    else:
        left_singular_vectors, singular_values, _ = np.linalg.svd(
            centered_matrix, full_matrices=False
        )
        sample_scores = left_singular_vectors * singular_values

    explained_variance = singular_values * singular_values
    total_variance = float(explained_variance.sum())
    if total_variance > 0:
        percent_variance = explained_variance / total_variance
    else:
        percent_variance = np.zeros(max(1, len(singular_values)), dtype=float)
        percent_variance[0] = 1.0

    pc1_values = (
        sample_scores[:, 0]
        if sample_scores.shape[1] >= 1
        else np.zeros(sample_count, dtype=float)
    )
    pc2_values = (
        sample_scores[:, 1]
        if sample_scores.shape[1] >= 2
        else np.zeros(sample_count, dtype=float)
    )
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
