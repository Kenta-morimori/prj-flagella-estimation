"""Deterministic feature baseline for Phase 4 clip classification."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


FEATURE_NAMES = (
    "intensity_mean",
    "intensity_std",
    "intensity_q90",
    "foreground_fraction_0p1",
    "foreground_fraction_0p5",
    "frame_mean_std",
    "temporal_diff_mean",
    "temporal_diff_std",
    "temporal_diff_max",
    "centroid_step_mean",
    "centroid_step_max",
    "radial_spread_mean",
)


@dataclass(frozen=True)
class NearestCentroidModel:
    classes: np.ndarray
    feature_mean: np.ndarray
    feature_scale: np.ndarray
    centroids: np.ndarray


def extract_clip_features(clip: np.ndarray) -> np.ndarray:
    """Extract compact spatial and temporal summary features from a uint8 clip."""

    if clip.dtype != np.uint8 or clip.ndim != 3 or clip.shape[0] == 0:
        raise ValueError("clip must be a non-empty uint8 array with shape (T, H, W)")

    frames = clip.astype(np.float64) / 255.0
    frame_means = frames.mean(axis=(1, 2))
    temporal_diff = (
        np.abs(np.diff(frames, axis=0))
        if frames.shape[0] > 1
        else np.zeros((1, *frames.shape[1:]), dtype=np.float64)
    )
    centroids, radial_spreads = _weighted_geometry(frames)
    centroid_steps = (
        np.linalg.norm(np.diff(centroids, axis=0), axis=1)
        if len(centroids) > 1
        else np.zeros(1, dtype=np.float64)
    )

    return np.asarray(
        [
            frames.mean(),
            frames.std(),
            np.quantile(frames, 0.9),
            np.mean(frames > 0.1),
            np.mean(frames > 0.5),
            frame_means.std(),
            temporal_diff.mean(),
            temporal_diff.std(),
            temporal_diff.max(),
            centroid_steps.mean(),
            centroid_steps.max(),
            radial_spreads.mean(),
        ],
        dtype=np.float64,
    )


def fit_nearest_centroid(
    features: np.ndarray, labels: np.ndarray
) -> NearestCentroidModel:
    """Fit a standardized nearest-centroid classifier."""

    features = np.asarray(features, dtype=np.float64)
    labels = np.asarray(labels, dtype=np.int64)
    if features.ndim != 2 or len(features) != len(labels) or len(features) == 0:
        raise ValueError("features must be a non-empty 2D array aligned with labels")
    if not np.isfinite(features).all():
        raise ValueError("features contain non-finite values")

    classes = np.unique(labels)
    if len(classes) < 2:
        raise ValueError("training data must contain at least two classes")
    feature_mean = features.mean(axis=0)
    feature_scale = features.std(axis=0)
    feature_scale = np.where(feature_scale < 1.0e-12, 1.0, feature_scale)
    standardized = (features - feature_mean) / feature_scale
    centroids = np.stack(
        [standardized[labels == class_id].mean(axis=0) for class_id in classes]
    )
    return NearestCentroidModel(
        classes=classes,
        feature_mean=feature_mean,
        feature_scale=feature_scale,
        centroids=centroids,
    )


def predict_nearest_centroid(
    model: NearestCentroidModel, features: np.ndarray
) -> np.ndarray:
    """Predict labels by squared Euclidean distance to standardized centroids."""

    features = np.asarray(features, dtype=np.float64)
    if features.ndim != 2 or features.shape[1] != model.feature_mean.shape[0]:
        raise ValueError("feature shape does not match the fitted model")
    standardized = (features - model.feature_mean) / model.feature_scale
    distances = np.sum(
        (standardized[:, np.newaxis, :] - model.centroids[np.newaxis, :, :]) ** 2,
        axis=2,
    )
    return model.classes[np.argmin(distances, axis=1)]


def confusion_matrix(
    y_true: np.ndarray, y_pred: np.ndarray, classes: np.ndarray
) -> np.ndarray:
    """Build a row=true, column=predicted confusion matrix."""

    y_true = np.asarray(y_true, dtype=np.int64)
    y_pred = np.asarray(y_pred, dtype=np.int64)
    classes = np.asarray(classes, dtype=np.int64)
    class_to_index = {int(value): index for index, value in enumerate(classes)}
    matrix = np.zeros((len(classes), len(classes)), dtype=np.int64)
    for actual, predicted in zip(y_true, y_pred, strict=True):
        if int(actual) not in class_to_index or int(predicted) not in class_to_index:
            raise ValueError("labels must be present in classes")
        matrix[class_to_index[int(actual)], class_to_index[int(predicted)]] += 1
    return matrix


def classification_metrics(
    y_true: np.ndarray, y_pred: np.ndarray, classes: np.ndarray
) -> dict[str, float | int]:
    """Return baseline classification metrics without external ML dependencies."""

    matrix = confusion_matrix(y_true, y_pred, classes)
    total = int(matrix.sum())
    if total == 0:
        raise ValueError("cannot evaluate an empty split")
    recalls: list[float] = []
    f1_scores: list[float] = []
    for index in range(len(classes)):
        true_positive = int(matrix[index, index])
        actual_count = int(matrix[index, :].sum())
        predicted_count = int(matrix[:, index].sum())
        if actual_count > 0:
            recalls.append(true_positive / actual_count)
        precision = true_positive / predicted_count if predicted_count else 0.0
        recall = true_positive / actual_count if actual_count else 0.0
        f1_scores.append(
            2.0 * precision * recall / (precision + recall)
            if precision + recall
            else 0.0
        )
    return {
        "sample_count": total,
        "accuracy": float(np.trace(matrix) / total),
        "balanced_accuracy": float(np.mean(recalls)),
        "macro_f1": float(np.mean(f1_scores)),
    }


def _weighted_geometry(frames: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    height, width = frames.shape[1:]
    yy, xx = np.mgrid[:height, :width]
    centroids: list[tuple[float, float]] = []
    radial_spreads: list[float] = []
    for frame in frames:
        weight = frame.sum()
        if weight <= 0.0:
            centroid_x = (width - 1) / 2.0
            centroid_y = (height - 1) / 2.0
            radial_spread = 0.0
        else:
            centroid_x = float(np.sum(frame * xx) / weight)
            centroid_y = float(np.sum(frame * yy) / weight)
            squared_radius = (xx - centroid_x) ** 2 + (yy - centroid_y) ** 2
            radial_spread = float(np.sqrt(np.sum(frame * squared_radius) / weight))
        centroids.append((centroid_x / width, centroid_y / height))
        radial_spreads.append(radial_spread / max(height, width))
    return np.asarray(centroids), np.asarray(radial_spreads)
