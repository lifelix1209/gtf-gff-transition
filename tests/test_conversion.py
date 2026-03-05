from gtf_gff_transition.convert import convert_text
from gtf_gff_transition.parsing import parse_feature_line, parse_gtf_attributes


def _data_records(gtf_text: str):
    records = []
    for line in gtf_text.splitlines():
        if not line or line.startswith("#"):
            continue
        records.append(parse_feature_line(line, parse_gtf_attributes))
    return records


def test_gtf_to_gff3_basic_mapping():
    gtf_text = """chr1\tsrc\tgene\t1\t1000\t.\t+\t.\tgene_id \"geneA\"; gene_name \"GA\";
chr1\tsrc\ttranscript\t1\t1000\t.\t+\t.\tgene_id \"geneA\"; transcript_id \"tx1\";
chr1\tsrc\texon\t1\t100\t.\t+\t.\tgene_id \"geneA\"; transcript_id \"tx1\"; exon_number \"1\";
"""

    out = convert_text(gtf_text, source_format="gtf", target_format="gff3")
    lines = [line for line in out.splitlines() if line]

    assert lines[0] == "##gff-version 3"
    assert "\tgene\t" in lines[1]
    assert "ID=geneA" in lines[1]
    assert "gene_name=GA" in lines[1]

    assert "\tmRNA\t" in lines[2]
    assert "ID=tx1" in lines[2]
    assert "Parent=geneA" in lines[2]

    assert "\texon\t" in lines[3]
    assert "Parent=tx1" in lines[3]


def test_gff3_to_gtf_parent_inference():
    gff3_text = """##gff-version 3
chr1\tsrc\tgene\t1\t1000\t.\t+\t.\tID=geneA
chr1\tsrc\tmRNA\t1\t1000\t.\t+\t.\tID=tx1;Parent=geneA
chr1\tsrc\texon\t1\t100\t.\t+\t.\tID=exon1;Parent=tx1;exon_number=1
"""

    out = convert_text(gff3_text, source_format="gff3", target_format="gtf")
    records = _data_records(out)

    assert records[0].feature_type == "gene"
    assert records[0].attributes["gene_id"] == "geneA"

    assert records[1].feature_type == "transcript"
    assert records[1].attributes["gene_id"] == "geneA"
    assert records[1].attributes["transcript_id"] == "tx1"

    assert records[2].feature_type == "exon"
    assert records[2].attributes["gene_id"] == "geneA"
    assert records[2].attributes["transcript_id"] == "tx1"
    assert records[2].attributes["exon_number"] == "1"


def test_round_trip_single_exon_keeps_ids():
    gtf_text = (
        "chr1\tsrc\texon\t1\t100\t.\t+\t.\t"
        'gene_id "geneZ"; transcript_id "txZ"; exon_number "1";\n'
    )

    gff3_text = convert_text(gtf_text, source_format="gtf", target_format="gff3")
    back_to_gtf = convert_text(gff3_text, source_format="gff3", target_format="gtf")
    record = _data_records(back_to_gtf)[0]

    assert record.attributes["gene_id"] == "geneZ"
    assert record.attributes["transcript_id"] == "txZ"
    assert record.attributes["exon_number"] == "1"
