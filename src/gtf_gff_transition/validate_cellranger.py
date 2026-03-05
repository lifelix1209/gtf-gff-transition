from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from .parsing import parse_gtf_attributes

PASS_ASCII = r"""
 ____    _    ____ ____  
|  _ \  / \  / ___/ ___| 
| |_) |/ _ \ \___ \___ \ 
|  __// ___ \ ___) |__) |
|_|  /_/   \_\____/____/ 
"""


@dataclass
class ValidationReport:
    is_pass: bool
    errors: list[str]
    warnings: list[str]
    stats: dict[str, int]


def _parse_gtf_line(line: str) -> tuple[list[str], dict[str, str]]:
    cols = line.rstrip("\n").split("\t")
    if len(cols) != 9:
        raise ValueError(f"Expected 9 columns, got {len(cols)}")
    attrs = parse_gtf_attributes(cols[8])
    return cols, attrs


def validate_cellranger_gtf(
    gtf_path: str | Path,
    *,
    require_hierarchy: bool = False,
    strict_recommended: bool = False,
) -> ValidationReport:
    path = Path(gtf_path)

    feature_counts: Counter[str] = Counter()
    errors: list[str] = []
    warnings: list[str] = []

    gene_ids: set[str] = set()
    transcript_ids: set[str] = set()

    total_data_lines = 0
    malformed_lines = 0
    invalid_coords = 0
    invalid_strand = 0

    missing_gene_id_by_feature: Counter[str] = Counter()
    missing_transcript_id_by_feature: Counter[str] = Counter()
    missing_gene_name_gene = 0
    missing_gene_biotype_gene = 0

    # Keep links observed in exon/transcript records for hierarchy checks.
    transcript_ids_seen_in_exon: set[str] = set()
    gene_ids_seen_in_exon: set[str] = set()
    gene_ids_seen_in_transcript: set[str] = set()

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            total_data_lines += 1

            try:
                cols, attrs = _parse_gtf_line(line)
            except Exception:
                malformed_lines += 1
                continue

            seqname, _source, feature, start, end, _score, strand, _frame, _attr_text = cols
            feature_lower = feature.lower()
            feature_counts[feature_lower] += 1

            if not seqname:
                malformed_lines += 1

            try:
                start_i = int(start)
                end_i = int(end)
                if start_i <= 0 or end_i <= 0 or start_i > end_i:
                    invalid_coords += 1
            except ValueError:
                invalid_coords += 1

            if strand not in {"+", "-", "."}:
                invalid_strand += 1

            gene_id = attrs.get("gene_id")
            transcript_id = attrs.get("transcript_id")

            if feature_lower == "gene":
                if gene_id:
                    gene_ids.add(gene_id)
                else:
                    missing_gene_id_by_feature["gene"] += 1

                if "gene_name" not in attrs:
                    missing_gene_name_gene += 1
                if "gene_biotype" not in attrs and "gene_type" not in attrs:
                    missing_gene_biotype_gene += 1

            elif feature_lower == "transcript":
                if gene_id:
                    gene_ids_seen_in_transcript.add(gene_id)
                else:
                    missing_gene_id_by_feature["transcript"] += 1

                if transcript_id:
                    transcript_ids.add(transcript_id)
                else:
                    missing_transcript_id_by_feature["transcript"] += 1

            elif feature_lower == "exon":
                if gene_id:
                    gene_ids_seen_in_exon.add(gene_id)
                else:
                    missing_gene_id_by_feature["exon"] += 1

                if transcript_id:
                    transcript_ids_seen_in_exon.add(transcript_id)
                else:
                    missing_transcript_id_by_feature["exon"] += 1

            else:
                # For other features we do not enforce hard requirements.
                pass

    if total_data_lines == 0:
        errors.append("GTF 中没有可用数据行。")

    if feature_counts.get("exon", 0) == 0:
        errors.append("缺少 exon 特征；Cell Ranger 需要 exon 注释。")

    if malformed_lines > 0:
        errors.append(f"存在格式错误的行：{malformed_lines}（非 9 列或关键列非法）。")

    if invalid_coords > 0:
        errors.append(f"存在非法坐标行：{invalid_coords}（start/end 非正整数或 start > end）。")

    if invalid_strand > 0:
        errors.append(f"存在非法 strand 行：{invalid_strand}（strand 必须是 + / - / .）。")

    for feature in sorted(missing_gene_id_by_feature):
        count = missing_gene_id_by_feature[feature]
        if count > 0:
            errors.append(f"{feature} 缺少 gene_id：{count} 行。")

    for feature in sorted(missing_transcript_id_by_feature):
        count = missing_transcript_id_by_feature[feature]
        if count > 0:
            errors.append(f"{feature} 缺少 transcript_id：{count} 行。")

    # Hierarchy checks:
    # - If transcript lines exist, exon transcript_id should map to transcript lines.
    # - If gene lines exist, transcript/exon gene_id should map to gene lines.
    if feature_counts.get("transcript", 0) > 0:
        orphan_exon_transcript = transcript_ids_seen_in_exon - transcript_ids
        if orphan_exon_transcript:
            errors.append(
                f"存在 exon 的 transcript_id 在 transcript 行中找不到：{len(orphan_exon_transcript)} 个。"
            )
    elif require_hierarchy:
        errors.append("缺少 transcript 特征行（启用了 --require-hierarchy）。")
    else:
        warnings.append("未检测到 transcript 特征行；若使用 ARC 等流程建议补齐 gene/transcript 层级。")

    if feature_counts.get("gene", 0) > 0:
        orphan_tx_gene = gene_ids_seen_in_transcript - gene_ids
        orphan_exon_gene = gene_ids_seen_in_exon - gene_ids
        if orphan_tx_gene:
            errors.append(f"存在 transcript 的 gene_id 在 gene 行中找不到：{len(orphan_tx_gene)} 个。")
        if orphan_exon_gene:
            errors.append(f"存在 exon 的 gene_id 在 gene 行中找不到：{len(orphan_exon_gene)} 个。")
    elif require_hierarchy:
        errors.append("缺少 gene 特征行（启用了 --require-hierarchy）。")
    else:
        warnings.append("未检测到 gene 特征行；建议补齐以便下游层级一致性检查。")

    if missing_gene_name_gene > 0:
        warnings.append(f"gene 缺少 gene_name：{missing_gene_name_gene} 行。")
    if missing_gene_biotype_gene > 0:
        warnings.append(f"gene 缺少 gene_biotype/gene_type：{missing_gene_biotype_gene} 行。")

    if strict_recommended and warnings:
        errors.extend(f"[strict] {item}" for item in warnings)
        warnings = []

    report = ValidationReport(
        is_pass=(len(errors) == 0),
        errors=errors,
        warnings=warnings,
        stats={
            "data_lines": total_data_lines,
            "genes": feature_counts.get("gene", 0),
            "transcripts": feature_counts.get("transcript", 0),
            "exons": feature_counts.get("exon", 0),
        },
    )
    return report


def _format_report(report: ValidationReport, gtf_path: str | Path) -> str:
    lines: list[str] = []
    lines.append(f"GTF: {Path(gtf_path)}")
    lines.append(
        "stats: "
        f"data_lines={report.stats['data_lines']} "
        f"genes={report.stats['genes']} "
        f"transcripts={report.stats['transcripts']} "
        f"exons={report.stats['exons']}"
    )

    if report.is_pass:
        lines.append(PASS_ASCII.rstrip("\n"))
        lines.append("Cell Ranger GTF check: PASS")
    else:
        lines.append("Cell Ranger GTF check: FAIL")
        lines.append("失败原因：")
        for item in report.errors:
            lines.append(f"- {item}")

    if report.warnings:
        lines.append("警告：")
        for item in report.warnings:
            lines.append(f"- {item}")

    return "\n".join(lines) + "\n"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="gtf-gff-validate-cellranger",
        description="Validate whether a GTF can be used as Cell Ranger input.",
    )
    parser.add_argument("-i", "--input", required=True, help="Input GTF path")
    parser.add_argument(
        "--require-hierarchy",
        action="store_true",
        help="Require gene/transcript feature lines and enforce complete hierarchy mapping.",
    )
    parser.add_argument(
        "--strict-recommended",
        action="store_true",
        help="Treat recommended-field warnings as failures.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    report = validate_cellranger_gtf(
        args.input,
        require_hierarchy=args.require_hierarchy,
        strict_recommended=args.strict_recommended,
    )
    print(_format_report(report, args.input), end="")
    return 0 if report.is_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
