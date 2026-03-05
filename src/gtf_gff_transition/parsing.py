from __future__ import annotations

from collections import OrderedDict
from typing import Callable
from urllib.parse import quote, unquote

from .model import FeatureRecord


def parse_feature_line(
    line: str,
    attribute_parser: Callable[[str], OrderedDict[str, str]],
) -> FeatureRecord:
    cols = line.rstrip("\n").split("\t")
    if len(cols) != 9:
        raise ValueError(f"Expected 9 columns, got {len(cols)}: {line.rstrip()}")

    return FeatureRecord(
        seqid=cols[0],
        source=cols[1],
        feature_type=cols[2],
        start=cols[3],
        end=cols[4],
        score=cols[5],
        strand=cols[6],
        phase=cols[7],
        attributes=attribute_parser(cols[8]),
    )


def format_feature_line(
    record: FeatureRecord,
    attribute_formatter: Callable[[OrderedDict[str, str]], str],
) -> str:
    return "\t".join(
        [
            record.seqid,
            record.source,
            record.feature_type,
            record.start,
            record.end,
            record.score,
            record.strand,
            record.phase,
            attribute_formatter(record.attributes),
        ]
    )


def parse_gtf_attributes(text: str) -> OrderedDict[str, str]:
    attrs: OrderedDict[str, str] = OrderedDict()
    text = text.strip()
    if not text or text == ".":
        return attrs

    # GTF convention: key "value"; key2 "value2";
    for part in text.split(";"):
        part = part.strip()
        if not part:
            continue
        if " " not in part:
            attrs[part] = ""
            continue

        key, raw_value = part.split(" ", 1)
        value = raw_value.strip()
        if value.startswith('"') and value.endswith('"') and len(value) >= 2:
            value = value[1:-1]
        value = value.replace(r"\"", '"').replace(r"\\", "\\")
        attrs[key] = value

    return attrs


def format_gtf_attributes(attrs: OrderedDict[str, str]) -> str:
    if not attrs:
        return "."

    rendered = []
    for key, value in attrs.items():
        escaped = str(value).replace("\\", r"\\").replace('"', r'\"')
        rendered.append(f'{key} "{escaped}";')
    return " ".join(rendered)


def parse_gff3_attributes(text: str) -> OrderedDict[str, str]:
    attrs: OrderedDict[str, str] = OrderedDict()
    text = text.strip()
    if not text or text == ".":
        return attrs

    for part in text.split(";"):
        part = part.strip()
        if not part:
            continue

        if "=" in part:
            key, raw_value = part.split("=", 1)
            attrs[unquote(key)] = unquote(raw_value)
        else:
            attrs[unquote(part)] = ""

    return attrs


def format_gff3_attributes(attrs: OrderedDict[str, str]) -> str:
    if not attrs:
        return "."

    rendered = []
    for key, value in attrs.items():
        enc_key = quote(str(key), safe="._:-")
        enc_value = quote(str(value), safe="._:-,|")
        rendered.append(f"{enc_key}={enc_value}")
    return ";".join(rendered)
