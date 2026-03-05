from __future__ import annotations

from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .model import FeatureRecord, is_gene_feature, is_transcript_feature
from .parsing import (
    format_feature_line,
    format_gff3_attributes,
    format_gtf_attributes,
    parse_feature_line,
    parse_gff3_attributes,
    parse_gtf_attributes,
)


@dataclass
class _LineEntry:
    kind: str
    value: str | FeatureRecord


class _IdFactory:
    def __init__(self) -> None:
        self._counters: defaultdict[str, int] = defaultdict(int)

    def next(self, prefix: str) -> str:
        self._counters[prefix] += 1
        return f"{prefix}_{self._counters[prefix]}"


def _normalize_format(name: str) -> str:
    value = name.strip().lower()
    if value == "gff":
        return "gff3"
    if value not in {"gtf", "gff3"}:
        raise ValueError(f"Unsupported format: {name}")
    return value


def _first_parent(parent_value: str | None) -> str | None:
    if not parent_value:
        return None
    for part in parent_value.split(","):
        part = part.strip()
        if part:
            return part
    return None


def _walk_ancestors(start_id: str, id_to_record: dict[str, FeatureRecord]):
    current = start_id
    seen: set[str] = set()

    while current and current not in seen:
        seen.add(current)
        record = id_to_record.get(current)
        if record is None:
            return
        yield record
        current = _first_parent(record.attributes.get("Parent"))


def _to_gff3_feature_type(feature_type: str) -> str:
    return "mRNA" if feature_type.lower() == "transcript" else feature_type


def _to_gtf_feature_type(feature_type: str) -> str:
    return "transcript" if is_transcript_feature(feature_type) else feature_type


def _convert_gtf_record_to_gff3(record: FeatureRecord, id_factory: _IdFactory) -> FeatureRecord:
    attrs_in = record.attributes
    attrs_out: OrderedDict[str, str] = OrderedDict()

    gene_id = attrs_in.get("gene_id")
    transcript_id = attrs_in.get("transcript_id")
    record_id = attrs_in.get("ID")
    parent_id = attrs_in.get("Parent")

    if is_gene_feature(record.feature_type):
        attrs_out["ID"] = record_id or gene_id or id_factory.next("gene")
    elif is_transcript_feature(record.feature_type):
        attrs_out["ID"] = record_id or transcript_id or id_factory.next("transcript")
        parent = parent_id or gene_id
        if parent:
            attrs_out["Parent"] = parent
    else:
        parent = parent_id or transcript_id or gene_id
        if parent:
            attrs_out["Parent"] = parent
        if record_id:
            attrs_out["ID"] = record_id

    if gene_id:
        attrs_out.setdefault("gene_id", gene_id)
    if transcript_id:
        attrs_out.setdefault("transcript_id", transcript_id)

    for key, value in attrs_in.items():
        if key in {"ID", "Parent", "gene_id", "transcript_id"}:
            continue
        attrs_out[key] = value

    return FeatureRecord(
        seqid=record.seqid,
        source=record.source,
        feature_type=_to_gff3_feature_type(record.feature_type),
        start=record.start,
        end=record.end,
        score=record.score,
        strand=record.strand,
        phase=record.phase,
        attributes=attrs_out,
    )


def _infer_gene_transcript_ids(
    record: FeatureRecord,
    id_to_record: dict[str, FeatureRecord],
    id_factory: _IdFactory,
) -> tuple[str | None, str | None]:
    attrs = record.attributes

    gene_id = attrs.get("gene_id")
    transcript_id = attrs.get("transcript_id")

    this_id = attrs.get("ID")
    parent_id = _first_parent(attrs.get("Parent"))

    if is_gene_feature(record.feature_type):
        gene_id = gene_id or this_id or id_factory.next("gene")
        return gene_id, None

    if is_transcript_feature(record.feature_type):
        transcript_id = transcript_id or this_id or id_factory.next("transcript")
        if not gene_id and parent_id:
            gene_id = parent_id
            parent_record = id_to_record.get(parent_id)
            if parent_record and is_gene_feature(parent_record.feature_type):
                gene_id = parent_record.attributes.get("gene_id") or parent_record.attributes.get("ID") or parent_id
        if not gene_id:
            gene_id = id_factory.next("gene")
        return gene_id, transcript_id

    if not transcript_id and parent_id:
        parent_record = id_to_record.get(parent_id)
        if parent_record is None:
            transcript_id = parent_id
        elif is_transcript_feature(parent_record.feature_type):
            transcript_id = parent_record.attributes.get("transcript_id") or parent_record.attributes.get("ID") or parent_id
            if not gene_id:
                parent_parent = _first_parent(parent_record.attributes.get("Parent"))
                gene_id = parent_record.attributes.get("gene_id") or parent_parent
        elif is_gene_feature(parent_record.feature_type):
            gene_id = parent_record.attributes.get("gene_id") or parent_record.attributes.get("ID") or parent_id

    if not transcript_id and parent_id:
        for ancestor in _walk_ancestors(parent_id, id_to_record):
            if is_transcript_feature(ancestor.feature_type):
                transcript_id = ancestor.attributes.get("transcript_id") or ancestor.attributes.get("ID")
                if transcript_id:
                    break

    if not gene_id and parent_id:
        for ancestor in _walk_ancestors(parent_id, id_to_record):
            if is_gene_feature(ancestor.feature_type):
                gene_id = ancestor.attributes.get("gene_id") or ancestor.attributes.get("ID")
                if gene_id:
                    break

    if not transcript_id:
        transcript_id = this_id or id_factory.next("transcript")
    if not gene_id:
        gene_id = id_factory.next("gene")

    return gene_id, transcript_id


def _convert_gff3_record_to_gtf(
    record: FeatureRecord,
    id_to_record: dict[str, FeatureRecord],
    id_factory: _IdFactory,
) -> FeatureRecord:
    gene_id, transcript_id = _infer_gene_transcript_ids(record, id_to_record, id_factory)

    attrs_out: OrderedDict[str, str] = OrderedDict()
    if gene_id:
        attrs_out["gene_id"] = gene_id
    if not is_gene_feature(record.feature_type) and transcript_id:
        attrs_out["transcript_id"] = transcript_id

    for key, value in record.attributes.items():
        if key in {"ID", "Parent", "gene_id", "transcript_id"}:
            continue
        attrs_out[key] = value

    return FeatureRecord(
        seqid=record.seqid,
        source=record.source,
        feature_type=_to_gtf_feature_type(record.feature_type),
        start=record.start,
        end=record.end,
        score=record.score,
        strand=record.strand,
        phase=record.phase,
        attributes=attrs_out,
    )


def _classify_lines(lines: Iterable[str]) -> list[_LineEntry]:
    entries: list[_LineEntry] = []
    for raw_line in lines:
        line = raw_line.rstrip("\n")
        if not line:
            entries.append(_LineEntry("blank", ""))
        elif line.startswith("#"):
            entries.append(_LineEntry("comment", line))
        else:
            entries.append(_LineEntry("data", line))
    return entries


def convert_gtf_to_gff3(lines: Iterable[str]) -> list[str]:
    entries = _classify_lines(lines)
    has_header = any(
        entry.kind == "comment" and str(entry.value).strip().lower() == "##gff-version 3"
        for entry in entries
    )

    id_factory = _IdFactory()
    out: list[str] = []
    first_data_written = False

    for entry in entries:
        if entry.kind in {"blank", "comment"}:
            out.append(str(entry.value))
            continue

        if not first_data_written:
            if not has_header:
                out.append("##gff-version 3")
            first_data_written = True

        in_record = parse_feature_line(str(entry.value), parse_gtf_attributes)
        out_record = _convert_gtf_record_to_gff3(in_record, id_factory)
        out.append(format_feature_line(out_record, format_gff3_attributes))

    return out


def convert_gff3_to_gtf(lines: Iterable[str]) -> list[str]:
    entries = _classify_lines(lines)

    parsed_entries: list[_LineEntry] = []
    id_to_record: dict[str, FeatureRecord] = {}

    for entry in entries:
        if entry.kind != "data":
            parsed_entries.append(entry)
            continue

        record = parse_feature_line(str(entry.value), parse_gff3_attributes)
        parsed_entries.append(_LineEntry("data", record))

        record_id = record.attributes.get("ID")
        if record_id and record_id not in id_to_record:
            id_to_record[record_id] = record

    id_factory = _IdFactory()
    out: list[str] = []

    for entry in parsed_entries:
        if entry.kind in {"blank", "comment"}:
            out.append(str(entry.value))
            continue

        in_record = entry.value
        if not isinstance(in_record, FeatureRecord):
            raise TypeError("Expected parsed FeatureRecord")

        out_record = _convert_gff3_record_to_gtf(in_record, id_to_record, id_factory)
        out.append(format_feature_line(out_record, format_gtf_attributes))

    return out


def convert_text(text: str, source_format: str, target_format: str) -> str:
    src = _normalize_format(source_format)
    dst = _normalize_format(target_format)

    if src == dst:
        raise ValueError("Source and target formats are identical")

    lines = text.splitlines()
    if src == "gtf" and dst == "gff3":
        converted = convert_gtf_to_gff3(lines)
    elif src == "gff3" and dst == "gtf":
        converted = convert_gff3_to_gtf(lines)
    else:
        raise ValueError(f"Unsupported conversion: {source_format} -> {target_format}")

    return "\n".join(converted) + ("\n" if text.endswith("\n") else "")


def convert_file(
    input_path: str | Path,
    output_path: str | Path,
    source_format: str,
    target_format: str,
) -> None:
    input_path = Path(input_path)
    output_path = Path(output_path)

    text = input_path.read_text(encoding="utf-8")
    converted = convert_text(text, source_format=source_format, target_format=target_format)
    output_path.write_text(converted, encoding="utf-8")
