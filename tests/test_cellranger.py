from gtf_gff_transition.cellranger import build_cellranger_gtf
from gtf_gff_transition.parsing import parse_feature_line, parse_gtf_attributes


def _records(text: str):
    return [
        parse_feature_line(line, parse_gtf_attributes)
        for line in text.splitlines()
        if line and not line.startswith("#")
    ]


def test_build_cellranger_gtf_merges_gff3_gtf_and_normalizes_features(tmp_path):
    gff3 = tmp_path / "input.gff3"
    gtf = tmp_path / "input.gtf"
    out = tmp_path / "output.gtf"

    gff3.write_text(
        "\n".join(
            [
                "##gff-version 3",
                "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=gene:GENE1;gene_id=GENE1;biotype=protein_coding",
                "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=transcript:TX1;Parent=gene:GENE1;transcript_id=TX1;biotype=protein_coding",
                "chr1\tsrc\texon\t1\t50\t.\t+\t.\tParent=transcript:TX1;rank=1;exon_id=EX1",
                "chr1\tsrc\tCDS\t1\t50\t.\t+\t0\tID=CDS:PX1;Parent=transcript:TX1;protein_id=PX1",
                "chr1\tsrc\texon\t80\t100\t.\t+\t.\tParent=transcript:TX1;rank=2;exon_id=EX2",
                "chr1\tsrc\tCDS\t80\t100\t.\t+\t2\tID=CDS:PX1;Parent=transcript:TX1;protein_id=PX1",
                "chr1\tsrc\tncRNA_gene\t200\t300\t.\t+\t.\tID=gene:GENE2;gene_id=GENE2;biotype=miRNA",
                "chr1\tsrc\tmiRNA\t200\t300\t.\t+\t.\tID=transcript:TX2;Parent=gene:GENE2;transcript_id=TX2;biotype=miRNA",
                "chr1\tsrc\texon\t200\t300\t.\t+\t.\tParent=transcript:TX2;rank=1;exon_id=EX3",
                "chr1\tsrc\tCDS\t220\t280\t.\t+\t1\tID=CDS:PX2;Parent=transcript:TX2;protein_id=PX2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    gtf.write_text(
        "\n".join(
            [
                'chr1\tsrc\ttranscript\t1\t100\t.\t+\t.\tgene_id "gene:GENE1"; transcript_id "transcript:TX1";',
                'chr1\tsrc\texon\t1\t50\t.\t+\t.\tgene_id "gene:GENE1"; transcript_id "transcript:TX1";',
                'chr1\tsrc\ttranscript\t200\t300\t.\t+\t.\tgene_id "gene:GENE2"; transcript_id "transcript:TX2";',
                'chr1\tsrc\texon\t200\t300\t.\t+\t.\tgene_id "gene:GENE2"; transcript_id "transcript:TX2";',
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    stats = build_cellranger_gtf(gff3, gtf, out)
    assert stats.genes_written == 2
    assert stats.transcripts_written == 2
    assert stats.exons_written == 3
    assert stats.cds_written == 3

    records = _records(out.read_text(encoding="utf-8"))
    assert records[0].feature_type == "gene"
    assert records[0].attributes["gene_id"] == "gene:GENE1"
    assert records[0].attributes["gene_biotype"] == "protein_coding"

    assert records[1].feature_type == "transcript"
    assert records[1].attributes["transcript_id"] == "transcript:TX1"
    assert records[1].attributes["gene_id"] == "gene:GENE1"

    exon_record = [record for record in records if record.feature_type == "exon"][0]
    assert exon_record.attributes["exon_number"] == "1"
    assert exon_record.attributes["exon_id"] == "EX1"

    cds_record = [record for record in records if record.feature_type == "CDS"][0]
    assert cds_record.attributes["transcript_id"] == "transcript:TX1"
    assert cds_record.attributes["protein_id"] == "PX1"


def test_build_cellranger_gtf_biotype_filter_keeps_recommended_types(tmp_path):
    gff3 = tmp_path / "input.gff3"
    gtf = tmp_path / "input.gtf"
    out = tmp_path / "output.gtf"

    gff3.write_text(
        "\n".join(
            [
                "##gff-version 3",
                "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=gene:GENE1;gene_id=GENE1;biotype=protein_coding",
                "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=transcript:TX1;Parent=gene:GENE1;transcript_id=TX1;biotype=protein_coding",
                "chr1\tsrc\texon\t1\t100\t.\t+\t.\tParent=transcript:TX1;rank=1",
                "chr1\tsrc\tCDS\t1\t90\t.\t+\t0\tID=CDS:PX1;Parent=transcript:TX1;protein_id=PX1",
                "chr1\tsrc\tncRNA_gene\t200\t300\t.\t+\t.\tID=gene:GENE2;gene_id=GENE2;biotype=miRNA",
                "chr1\tsrc\tmiRNA\t200\t300\t.\t+\t.\tID=transcript:TX2;Parent=gene:GENE2;transcript_id=TX2;biotype=miRNA",
                "chr1\tsrc\texon\t200\t300\t.\t+\t.\tParent=transcript:TX2;rank=1",
                "chr1\tsrc\tCDS\t210\t290\t.\t+\t1\tID=CDS:PX2;Parent=transcript:TX2;protein_id=PX2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    gtf.write_text(
        "\n".join(
            [
                'chr1\tsrc\ttranscript\t1\t100\t.\t+\t.\tgene_id "gene:GENE1"; transcript_id "transcript:TX1";',
                'chr1\tsrc\ttranscript\t200\t300\t.\t+\t.\tgene_id "gene:GENE2"; transcript_id "transcript:TX2";',
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    stats = build_cellranger_gtf(gff3, gtf, out, cellranger_biotypes_only=True)
    assert stats.genes_written == 1
    assert stats.transcripts_written == 1
    assert stats.exons_written == 1
    assert stats.cds_written == 1

    records = _records(out.read_text(encoding="utf-8"))
    assert len(records) == 4
    assert all(record.attributes["gene_id"] == "gene:GENE1" for record in records)
