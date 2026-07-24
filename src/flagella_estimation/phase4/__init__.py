"""Phase 4 dataset loading, training, and evaluation helpers."""

from flagella_estimation.phase4.freeze import (
    DatasetFreezeAudit,
    DatasetFreezePolicy,
    audit_phase4_dataset_freeze,
)
from flagella_estimation.phase4.training import (
    Phase4BaselineConfig,
    load_baseline_config,
    train_baseline_classifier,
)
from flagella_estimation.phase4.learning_curve import (
    Phase4LearningCurveConfig,
    evaluate_grouped_learning_curve,
    load_learning_curve_config,
)

__all__ = [
    "DatasetFreezeAudit",
    "DatasetFreezePolicy",
    "Phase4BaselineConfig",
    "Phase4LearningCurveConfig",
    "audit_phase4_dataset_freeze",
    "evaluate_grouped_learning_curve",
    "load_baseline_config",
    "load_learning_curve_config",
    "train_baseline_classifier",
]
