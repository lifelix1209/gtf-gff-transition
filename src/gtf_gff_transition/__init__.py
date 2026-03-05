"""Utilities for GTF/GFF3 conversion and Cell Ranger-oriented workflows."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from .convert import convert_gff3_to_gtf, convert_gtf_to_gff3, convert_text

if TYPE_CHECKING:  # pragma: no cover
    from .cellranger import BuildStats
    from .validate_cellranger import ValidationReport

__all__ = [
    "BuildStats",
    "ValidationReport",
    "build_cellranger_gtf",
    "convert_gff3_to_gtf",
    "convert_gtf_to_gff3",
    "convert_text",
    "validate_cellranger_gtf",
]

__version__ = "0.1.0"


def __getattr__(name: str) -> Any:
    if name in {"BuildStats", "build_cellranger_gtf"}:
        from .cellranger import BuildStats, build_cellranger_gtf

        return {"BuildStats": BuildStats, "build_cellranger_gtf": build_cellranger_gtf}[name]
    if name in {"ValidationReport", "validate_cellranger_gtf"}:
        from .validate_cellranger import ValidationReport, validate_cellranger_gtf

        return {"ValidationReport": ValidationReport, "validate_cellranger_gtf": validate_cellranger_gtf}[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
