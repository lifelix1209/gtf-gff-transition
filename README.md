# gtf-gff-transition

`gtf-gff-transition` is a Python package for annotation format conversion and Cell Ranger-oriented preprocessing.

It provides:

- Bidirectional conversion between `GTF` and `GFF3`
- A merge workflow that combines `GFF3 + GTF` into a Cell Ranger-compatible `GTF`
- A validator that checks whether a `GTF` is suitable as Cell Ranger input (`PASS/FAIL` with detailed reasons)

## Installation

Install from package index:

```bash
pip install gtf-gff-transition
```

Install for local development:

```bash
pip install -e .[dev]
```

## CLI Commands

The package installs three commands:

1. `gtf-gff-transition`
2. `gtf-gff-cellranger`
3. `gtf-gff-validate-cellranger`

### 1) GTF/GFF3 Conversion

```bash
gtf-gff-transition -i input.gtf -o output.gff3 --from gtf --to gff3
gtf-gff-transition -i input.gff3 -o output.gtf --from gff3 --to gtf
```

Standard input/output is supported via `-i -` and `-o -`.

### 2) Build Cell Ranger-Compatible GTF

```bash
gtf-gff-cellranger \
  --gff3 input.gff3 \
  --gtf input.gtf \
  -o output.cellranger.gtf
```

Filter to the default 10x-recommended biotypes:

```bash
gtf-gff-cellranger \
  --gff3 input.gff3 \
  --gtf input.gtf \
  --cellranger-biotypes-only \
  -o output.cellranger.recommended.gtf
```

### 3) Validate GTF for Cell Ranger

```bash
gtf-gff-validate-cellranger -i output.cellranger.gtf
```

Behavior:

- PASS: prints an ASCII `PASS` banner and exits with code `0`
- FAIL: prints failure reasons (missing fields, hierarchy breaks, malformed lines, etc.) and exits with code `1`

Require strict `gene/transcript/exon` hierarchy:

```bash
gtf-gff-validate-cellranger -i output.cellranger.gtf --require-hierarchy
```

Treat recommended-field warnings as failures:

```bash
gtf-gff-validate-cellranger -i output.cellranger.gtf --strict-recommended
```

## Library API

```python
from gtf_gff_transition import (
    build_cellranger_gtf,
    convert_text,
    validate_cellranger_gtf,
)

converted = convert_text(
    "chr1\tsrc\texon\t1\t10\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n",
    source_format="gtf",
    target_format="gff3",
)

stats = build_cellranger_gtf("input.gff3", "input.gtf", "output.cellranger.gtf")
report = validate_cellranger_gtf("output.cellranger.gtf", require_hierarchy=True)
print(stats, report.is_pass)
```

## Implementation Logic

### GFF3 + GTF Merge Strategy (Cell Ranger Builder)

Core implementation:

- `src/gtf_gff_transition/cellranger.py`
- `tests/test_cellranger.py`

Pipeline:

1. Read ID-style hints from the input GTF (`gene_id`, `transcript_id`) to preserve ID style.
2. Scan GFF3 to collect gene/transcript metadata.
3. Normalize feature types:
   - gene-like features (`gene`, `ncRNA_gene`, `pseudogene`, etc.) -> `gene`
   - transcript-like features (`mRNA`, `lnc_RNA`, `miRNA`, `transcript`, etc.) -> `transcript`
   - structural output keeps `exon` and `CDS`
4. Merge attributes:
   - `gene`: `gene_id`, `gene_name`, `gene_biotype`, `gene_type`
   - `transcript`: `gene_id`, `transcript_id`, `gene_name`, `gene_biotype`, `gene_type`, `transcript_biotype`, `transcript_type`, `transcript_name`
   - `exon`: `gene_id`, `transcript_id`, and optional `exon_number`/`exon_id`
   - `CDS`: `gene_id`, `transcript_id`, and optional `protein_id`
5. Enforce hierarchy closure:
   - keep transcripts only if transcript features exist in GFF3 and have exons
   - keep exons/CDS only when their transcript is retained
6. Optionally filter by biotype:
   - `--cellranger-biotypes-only` uses a default 10x-oriented whitelist
   - `--allowed-biotypes` applies a custom whitelist

### Cell Ranger Input Validator Strategy

Core implementation:

- `src/gtf_gff_transition/validate_cellranger.py`
- `tests/test_validate_cellranger.py`

Validation checks:

- 9-column GTF structure
- coordinate sanity (`start/end`)
- strand sanity (`+`, `-`, `.`)
- required IDs:
  - `gene_id` on gene/transcript/exon where applicable
  - `transcript_id` on transcript/exon
- hierarchy consistency:
  - exon transcript IDs map to transcript records
  - transcript/exon gene IDs map to gene records (when gene records exist)
- recommended metadata:
  - `gene_name`
  - `gene_biotype`/`gene_type`

## Reproducible Example with Included Data

Input files:

- `data/astatotilapia_calliptera_gca964374335v1.gff3`
- `data/astatotilapia_calliptera_gca964374335v1.gtf`

Build recommended output:

```bash
PYTHONPATH=src python -m gtf_gff_transition.cellranger \
  --gff3 data/astatotilapia_calliptera_gca964374335v1.gff3 \
  --gtf data/astatotilapia_calliptera_gca964374335v1.gtf \
  --cellranger-biotypes-only \
  --output data/astatotilapia_calliptera_cellranger_recommended.gtf
```

Build full output:

```bash
PYTHONPATH=src python -m gtf_gff_transition.cellranger \
  --gff3 data/astatotilapia_calliptera_gca964374335v1.gff3 \
  --gtf data/astatotilapia_calliptera_gca964374335v1.gtf \
  --output data/astatotilapia_calliptera_cellranger_full.gtf
```

Validate output:

```bash
PYTHONPATH=src python -m gtf_gff_transition.validate_cellranger \
  -i data/astatotilapia_calliptera_cellranger_recommended.gtf \
  --require-hierarchy
```

Observed dataset-level output stats:

- `astatotilapia_calliptera_cellranger_full.gtf`: total `2472066`, gene `36129`, transcript `100860`, exon `1210270`, CDS `1124804`
- `astatotilapia_calliptera_cellranger_recommended.gtf`: total `2462622`, gene `33376`, transcript `98107`, exon `1206332`, CDS `1124804`

## References

References used for format and Cell Ranger compatibility guidance (retrieved on 2026-03-05):

1. Ensembl format documentation (GFF/GTF):  
   https://www.ensembl.org/info/website/upload/gff.html
2. GENCODE GTF data format notes:  
   https://www.gencodegenes.org/pages/data_format.html
3. 10x Cell Ranger tutorial (reference build workflow):  
   https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-mr
4. 10x Cell Ranger reference / `mkgtf` guidance:  
   https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/cr-3p-references
5. 10x Cell Ranger ARC custom reference requirements:  
   https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/tutorials/create-a-custom-reference-genome-for-cell-ranger-arc
