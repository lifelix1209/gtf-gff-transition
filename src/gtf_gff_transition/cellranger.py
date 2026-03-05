from __future__ import annotations

import argparse
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .parsing import format_gtf_attributes, parse_gff3_attributes, parse_gtf_attributes

GENE_FEATURE_TYPES = {
    "gene",
    "ncrna_gene",
    "pseudogene",
    "j_gene_segment",
    "v_gene_segment",
}

TRANSCRIPT_FEATURE_TYPES = {
    "transcript",
    "mrna",
    "rna",
    "lnc_rna",
    "ncrna",
    "pseudogenic_transcript",
    "mirna",
    "rrna",
    "scrna",
    "snorna",
    "snrna",
    "trna",
    "y_rna",
}

CELLRANGER_DEFAULT_BIOTYPES = {
    "protein_coding",
    "lncrna",
    "antisense",
    "ig_lv_gene",
    "ig_v_gene",
    "ig_v_pseudogene",
    "ig_d_gene",
    "ig_j_gene",
    "ig_j_pseudogene",
    "ig_c_gene",
    "ig_c_pseudogene",
    "tr_v_gene",
    "tr_v_pseudogene",
    "tr_d_gene",
    "tr_j_gene",
    "tr_j_pseudogene",
    "tr_c_gene",
}


@dataclass
class GeneMeta:
    gene_id: str
    gene_name: str
    gene_biotype: str


@dataclass
class TranscriptMeta:
    transcript_id: str
    gene_key: str
    transcript_name: str
    transcript_biotype: str


@dataclass
class BuildStats:
    genes_written: int
    transcripts_written: int
    exons_written: int
    input_lines: int
    skipped_exons_missing_transcript: int


@dataclass
class _IdHints:
    gene_style_by_key: dict[str, str]
    transcript_style_by_key: dict[str, str]
    transcript_to_gene_key: dict[str, str]


def _normalize_format_id(value: str | None) -> str | None:
    if value is None:
        return None
    value = value.strip()
    if not value:
        return None
    if value.startswith("gene:"):
        tail = value.split(":", 1)[1]
        return tail or value
    if value.startswith("transcript:"):
        tail = value.split(":", 1)[1]
        return tail or value
    return value


def _normalize_biotype(value: str | None) -> str:
    if not value:
        return "unknown"
    return value.strip() or "unknown"


def _first_parent(parent_value: str | None) -> str | None:
    if not parent_value:
        return None
    for part in parent_value.split(","):
        part = part.strip()
        if part:
            return part
    return None


def _is_gene_feature(feature_type: str, attrs: OrderedDict[str, str]) -> bool:
    normalized = feature_type.lower()
    if normalized in GENE_FEATURE_TYPES:
        return True

    feature_id = attrs.get("ID", "")
    if feature_id.startswith("gene:"):
        return True

    return bool(attrs.get("gene_id") and not attrs.get("Parent"))


def _is_transcript_feature(feature_type: str, attrs: OrderedDict[str, str]) -> bool:
    normalized = feature_type.lower()
    if normalized in TRANSCRIPT_FEATURE_TYPES:
        return True

    if normalized in {"exon", "cds", "five_prime_utr", "three_prime_utr", "start_codon", "stop_codon"}:
        return False

    feature_id = attrs.get("ID", "")
    if feature_id.startswith("transcript:"):
        return True

    return "transcript_id" in attrs


def _parse_data_line(line: str) -> tuple[list[str], OrderedDict[str, str]]:
    cols = line.rstrip("\n").split("\t")
    if len(cols) != 9:
        raise ValueError(f"Expected 9 columns in annotation line, got {len(cols)}")
    attrs = parse_gff3_attributes(cols[8])
    return cols, attrs


def _parse_gtf_data_line(line: str) -> OrderedDict[str, str]:
    cols = line.rstrip("\n").split("\t")
    if len(cols) != 9:
        raise ValueError(f"Expected 9 columns in GTF line, got {len(cols)}")
    return parse_gtf_attributes(cols[8])


def _pick_gene_key(attrs: OrderedDict[str, str]) -> str | None:
    return _normalize_format_id(attrs.get("gene_id") or attrs.get("ID"))


def _pick_transcript_key(attrs: OrderedDict[str, str]) -> str | None:
    return _normalize_format_id(attrs.get("transcript_id") or attrs.get("ID"))


def _pick_parent_gene_key(attrs: OrderedDict[str, str]) -> str | None:
    if attrs.get("gene_id"):
        return _normalize_format_id(attrs.get("gene_id"))
    return _normalize_format_id(_first_parent(attrs.get("Parent")))


def _pick_parent_transcript_key(attrs: OrderedDict[str, str]) -> str | None:
    if attrs.get("transcript_id"):
        return _normalize_format_id(attrs.get("transcript_id"))
    return _normalize_format_id(_first_parent(attrs.get("Parent")))


def _pick_gene_biotype(attrs: OrderedDict[str, str], feature_type: str) -> str:
    value = attrs.get("gene_biotype") or attrs.get("gene_type") or attrs.get("biotype")
    if value:
        return _normalize_biotype(value)
    normalized = feature_type.lower()
    if normalized in {"ncrna_gene", "lnc_rna"}:
        return "lncRNA"
    return "unknown"


def _pick_transcript_biotype(attrs: OrderedDict[str, str], gene_biotype: str) -> str:
    value = attrs.get("transcript_biotype") or attrs.get("transcript_type") or attrs.get("biotype")
    if value:
        return _normalize_biotype(value)
    return gene_biotype


def _build_gtf_attrs_for_gene(meta: GeneMeta) -> OrderedDict[str, str]:
    attrs: OrderedDict[str, str] = OrderedDict()
    attrs["gene_id"] = meta.gene_id
    attrs["gene_name"] = meta.gene_name
    attrs["gene_biotype"] = meta.gene_biotype
    attrs["gene_type"] = meta.gene_biotype
    return attrs


def _build_gtf_attrs_for_transcript(meta: TranscriptMeta, gene: GeneMeta) -> OrderedDict[str, str]:
    attrs: OrderedDict[str, str] = OrderedDict()
    attrs["gene_id"] = gene.gene_id
    attrs["transcript_id"] = meta.transcript_id
    attrs["gene_name"] = gene.gene_name
    attrs["gene_biotype"] = gene.gene_biotype
    attrs["gene_type"] = gene.gene_biotype
    attrs["transcript_biotype"] = meta.transcript_biotype
    attrs["transcript_type"] = meta.transcript_biotype
    attrs["transcript_name"] = meta.transcript_name
    return attrs


def _build_gtf_attrs_for_exon(
    in_attrs: OrderedDict[str, str],
    transcript: TranscriptMeta,
    gene: GeneMeta,
) -> OrderedDict[str, str]:
    attrs: OrderedDict[str, str] = OrderedDict()
    attrs["gene_id"] = gene.gene_id
    attrs["transcript_id"] = transcript.transcript_id
    attrs["gene_name"] = gene.gene_name
    attrs["gene_biotype"] = gene.gene_biotype
    attrs["gene_type"] = gene.gene_biotype
    attrs["transcript_biotype"] = transcript.transcript_biotype
    attrs["transcript_type"] = transcript.transcript_biotype

    exon_number = in_attrs.get("exon_number") or in_attrs.get("rank")
    if exon_number:
        attrs["exon_number"] = exon_number
    exon_id = in_attrs.get("exon_id")
    if exon_id:
        attrs["exon_id"] = exon_id
    return attrs


def _render_gtf_line(cols: list[str], feature_type: str, attrs: OrderedDict[str, str]) -> str:
    fields = [
        cols[0],
        cols[1],
        feature_type,
        cols[3],
        cols[4],
        cols[5],
        cols[6],
        cols[7],
        format_gtf_attributes(attrs),
    ]
    return "\t".join(fields)


def _read_id_hints_from_gtf(gtf_path: Path) -> _IdHints:
    gene_style_by_key: dict[str, str] = {}
    transcript_style_by_key: dict[str, str] = {}
    transcript_to_gene_key: dict[str, str] = {}

    with gtf_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            attrs = _parse_gtf_data_line(line)
            gene_id = attrs.get("gene_id")
            transcript_id = attrs.get("transcript_id")
            gene_key = _normalize_format_id(gene_id)
            transcript_key = _normalize_format_id(transcript_id)

            if gene_key and gene_id and gene_key not in gene_style_by_key:
                gene_style_by_key[gene_key] = gene_id
            if transcript_key and transcript_id and transcript_key not in transcript_style_by_key:
                transcript_style_by_key[transcript_key] = transcript_id
            if gene_key and transcript_key and transcript_key not in transcript_to_gene_key:
                transcript_to_gene_key[transcript_key] = gene_key

    return _IdHints(
        gene_style_by_key=gene_style_by_key,
        transcript_style_by_key=transcript_style_by_key,
        transcript_to_gene_key=transcript_to_gene_key,
    )


def _merge_gene_meta(existing: GeneMeta | None, incoming: GeneMeta) -> GeneMeta:
    if existing is None:
        return incoming
    if existing.gene_name == existing.gene_id and incoming.gene_name != incoming.gene_id:
        existing.gene_name = incoming.gene_name
    if existing.gene_biotype == "unknown" and incoming.gene_biotype != "unknown":
        existing.gene_biotype = incoming.gene_biotype
    return existing


def _merge_transcript_meta(existing: TranscriptMeta | None, incoming: TranscriptMeta) -> TranscriptMeta:
    if existing is None:
        return incoming
    if not existing.gene_key and incoming.gene_key:
        existing.gene_key = incoming.gene_key
    if existing.transcript_name == existing.transcript_id and incoming.transcript_name != incoming.transcript_id:
        existing.transcript_name = incoming.transcript_name
    if existing.transcript_biotype == "unknown" and incoming.transcript_biotype != "unknown":
        existing.transcript_biotype = incoming.transcript_biotype
    return existing


def _collect_metadata(
    gff3_path: Path,
    hints: _IdHints,
) -> tuple[dict[str, GeneMeta], dict[str, TranscriptMeta], set[str], set[str], int]:
    genes: dict[str, GeneMeta] = {}
    transcripts: dict[str, TranscriptMeta] = {}
    transcript_with_exon: set[str] = set()
    transcript_with_feature: set[str] = set()
    input_lines = 0

    with gff3_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            input_lines += 1

            cols, attrs = _parse_data_line(line)
            feature_type = cols[2]

            if _is_gene_feature(feature_type, attrs):
                gene_key = _pick_gene_key(attrs)
                if not gene_key:
                    continue

                gene_id = hints.gene_style_by_key.get(gene_key) or attrs.get("gene_id") or gene_key
                gene_name = attrs.get("gene_name") or attrs.get("Name") or gene_id
                gene_biotype = _pick_gene_biotype(attrs, feature_type)

                incoming = GeneMeta(gene_id=gene_id, gene_name=gene_name, gene_biotype=gene_biotype)
                genes[gene_key] = _merge_gene_meta(genes.get(gene_key), incoming)
                continue

            if _is_transcript_feature(feature_type, attrs):
                transcript_key = _pick_transcript_key(attrs)
                if not transcript_key:
                    continue

                gene_key = _pick_parent_gene_key(attrs) or hints.transcript_to_gene_key.get(transcript_key)
                if not gene_key:
                    continue

                transcript_id = (
                    hints.transcript_style_by_key.get(transcript_key) or attrs.get("transcript_id") or transcript_key
                )

                gene = genes.get(gene_key)
                gene_biotype = gene.gene_biotype if gene else _pick_gene_biotype(attrs, feature_type)
                transcript_biotype = _pick_transcript_biotype(attrs, gene_biotype)
                transcript_name = attrs.get("transcript_name") or attrs.get("Name") or transcript_id

                incoming = TranscriptMeta(
                    transcript_id=transcript_id,
                    gene_key=gene_key,
                    transcript_name=transcript_name,
                    transcript_biotype=transcript_biotype,
                )
                transcripts[transcript_key] = _merge_transcript_meta(transcripts.get(transcript_key), incoming)
                transcript_with_feature.add(transcript_key)

                if gene_key not in genes:
                    fallback_gene_id = hints.gene_style_by_key.get(gene_key) or gene_key
                    genes[gene_key] = GeneMeta(
                        gene_id=fallback_gene_id,
                        gene_name=fallback_gene_id,
                        gene_biotype=_normalize_biotype(transcript_biotype),
                    )
                continue

            if cols[2].lower() == "exon":
                parent_tx_key = _pick_parent_transcript_key(attrs)
                if parent_tx_key:
                    transcript_with_exon.add(parent_tx_key)
                    if parent_tx_key not in transcripts and parent_tx_key in hints.transcript_to_gene_key:
                        gene_key = hints.transcript_to_gene_key[parent_tx_key]
                        transcript_id = hints.transcript_style_by_key.get(parent_tx_key) or parent_tx_key
                        transcripts[parent_tx_key] = TranscriptMeta(
                            transcript_id=transcript_id,
                            gene_key=gene_key,
                            transcript_name=transcript_id,
                            transcript_biotype="unknown",
                        )
                        if gene_key not in genes:
                            fallback_gene_id = hints.gene_style_by_key.get(gene_key) or gene_key
                            genes[gene_key] = GeneMeta(
                                gene_id=fallback_gene_id,
                                gene_name=fallback_gene_id,
                                gene_biotype="unknown",
                            )

    return genes, transcripts, transcript_with_exon, transcript_with_feature, input_lines


def _parse_biotype_filter(value: str | None) -> set[str] | None:
    if value is None:
        return None
    parts = [item.strip() for item in value.split(",")]
    normalized = {_normalize_biotype(part).lower() for part in parts if part}
    if not normalized:
        return None
    return normalized


def build_cellranger_gtf(
    gff3_path: str | Path,
    gtf_path: str | Path,
    output_path: str | Path,
    *,
    cellranger_biotypes_only: bool = False,
    allowed_biotypes: Iterable[str] | None = None,
) -> BuildStats:
    gff3_path = Path(gff3_path)
    gtf_path = Path(gtf_path)
    output_path = Path(output_path)

    hints = _read_id_hints_from_gtf(gtf_path)
    genes, transcripts, transcript_with_exon, transcript_with_feature, input_lines = _collect_metadata(
        gff3_path, hints
    )

    if allowed_biotypes is None:
        keep_biotypes = set(CELLRANGER_DEFAULT_BIOTYPES) if cellranger_biotypes_only else None
    else:
        keep_biotypes = {_normalize_biotype(value).lower() for value in allowed_biotypes}

    valid_transcripts: set[str] = set()
    for tx_key, meta in transcripts.items():
        if tx_key not in transcript_with_exon:
            continue
        if tx_key not in transcript_with_feature:
            continue
        if not meta.gene_key or meta.gene_key not in genes:
            continue

        gene_biotype = _normalize_biotype(genes[meta.gene_key].gene_biotype).lower()
        if keep_biotypes is not None and gene_biotype not in keep_biotypes:
            continue
        valid_transcripts.add(tx_key)

    valid_genes = {transcripts[tx_key].gene_key for tx_key in valid_transcripts}

    genes_written = 0
    transcripts_written = 0
    exons_written = 0
    skipped_exons_missing_transcript = 0

    emitted_genes: set[str] = set()
    emitted_transcripts: set[str] = set()

    with gff3_path.open("r", encoding="utf-8") as in_handle, output_path.open("w", encoding="utf-8") as out_handle:
        out_handle.write(
            "# generated_by gtf-gff-transition cellranger builder\n"
            f"# source_gff3 {gff3_path.name}\n"
            f"# source_gtf {gtf_path.name}\n"
        )

        for raw_line in in_handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            cols, attrs = _parse_data_line(line)
            feature_type = cols[2]

            if _is_gene_feature(feature_type, attrs):
                gene_key = _pick_gene_key(attrs)
                if not gene_key or gene_key not in valid_genes or gene_key in emitted_genes:
                    continue
                out_attrs = _build_gtf_attrs_for_gene(genes[gene_key])
                out_handle.write(_render_gtf_line(cols, "gene", out_attrs) + "\n")
                emitted_genes.add(gene_key)
                genes_written += 1
                continue

            if _is_transcript_feature(feature_type, attrs):
                tx_key = _pick_transcript_key(attrs)
                if not tx_key or tx_key not in valid_transcripts or tx_key in emitted_transcripts:
                    continue
                tx_meta = transcripts[tx_key]
                gene_meta = genes.get(tx_meta.gene_key)
                if gene_meta is None:
                    continue
                out_attrs = _build_gtf_attrs_for_transcript(tx_meta, gene_meta)
                out_handle.write(_render_gtf_line(cols, "transcript", out_attrs) + "\n")
                emitted_transcripts.add(tx_key)
                transcripts_written += 1
                continue

            if feature_type.lower() == "exon":
                tx_key = _pick_parent_transcript_key(attrs)
                if not tx_key or tx_key not in valid_transcripts:
                    skipped_exons_missing_transcript += 1
                    continue
                tx_meta = transcripts.get(tx_key)
                if tx_meta is None:
                    skipped_exons_missing_transcript += 1
                    continue
                gene_meta = genes.get(tx_meta.gene_key)
                if gene_meta is None:
                    skipped_exons_missing_transcript += 1
                    continue

                out_attrs = _build_gtf_attrs_for_exon(attrs, tx_meta, gene_meta)
                out_handle.write(_render_gtf_line(cols, "exon", out_attrs) + "\n")
                exons_written += 1

    return BuildStats(
        genes_written=genes_written,
        transcripts_written=transcripts_written,
        exons_written=exons_written,
        input_lines=input_lines,
        skipped_exons_missing_transcript=skipped_exons_missing_transcript,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="gtf-gff-cellranger",
        description="Merge GFF3 and GTF information into a Cell Ranger-compatible GTF.",
    )
    parser.add_argument("--gff3", required=True, help="Input GFF3 file path")
    parser.add_argument("--gtf", required=True, help="Input GTF file path")
    parser.add_argument("-o", "--output", required=True, help="Output GTF path")
    parser.add_argument(
        "--cellranger-biotypes-only",
        action="store_true",
        help="Keep only the recommended 10x biotypes for Cell Ranger references.",
    )
    parser.add_argument(
        "--allowed-biotypes",
        default=None,
        help="Optional comma-separated biotypes whitelist; overrides --cellranger-biotypes-only.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    allowed = _parse_biotype_filter(args.allowed_biotypes)
    stats = build_cellranger_gtf(
        gff3_path=args.gff3,
        gtf_path=args.gtf,
        output_path=args.output,
        cellranger_biotypes_only=args.cellranger_biotypes_only,
        allowed_biotypes=allowed,
    )

    print(
        "done",
        f"genes={stats.genes_written}",
        f"transcripts={stats.transcripts_written}",
        f"exons={stats.exons_written}",
        f"input_lines={stats.input_lines}",
        f"skipped_exons={stats.skipped_exons_missing_transcript}",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
