from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .convert import convert_text


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="gtf-gff-transition",
        description="Convert between GTF and GFF3 annotation formats.",
    )
    parser.add_argument("-i", "--input", required=True, help="Input file path, or '-' for stdin")
    parser.add_argument("-o", "--output", required=True, help="Output file path, or '-' for stdout")
    parser.add_argument(
        "--from",
        dest="source_format",
        required=True,
        choices=["gtf", "gff", "gff3"],
        help="Source format",
    )
    parser.add_argument(
        "--to",
        dest="target_format",
        required=True,
        choices=["gtf", "gff", "gff3"],
        help="Target format",
    )
    return parser


def _read_input(path: str) -> str:
    if path == "-":
        return sys.stdin.read()
    return Path(path).read_text(encoding="utf-8")


def _write_output(path: str, content: str) -> None:
    if path == "-":
        sys.stdout.write(content)
        return
    Path(path).write_text(content, encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        text = _read_input(args.input)
        converted = convert_text(
            text,
            source_format=args.source_format,
            target_format=args.target_format,
        )
        _write_output(args.output, converted)
    except Exception as exc:  # pragma: no cover - handled by CLI behavior
        parser.exit(status=1, message=f"error: {exc}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
