from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass


@dataclass
class FeatureRecord:
    """A parsed 9-column annotation record."""

    seqid: str
    source: str
    feature_type: str
    start: str
    end: str
    score: str
    strand: str
    phase: str
    attributes: OrderedDict[str, str]


def is_transcript_feature(feature_type: str) -> bool:
    return feature_type.lower() in {"transcript", "mrna", "rna", "lnc_rna", "ncrna"}


def is_gene_feature(feature_type: str) -> bool:
    return feature_type.lower() == "gene"
