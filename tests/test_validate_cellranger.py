from gtf_gff_transition.validate_cellranger import _format_report, validate_cellranger_gtf


def test_validate_cellranger_pass_with_hierarchy(tmp_path):
    gtf = tmp_path / "ok.gtf"
    gtf.write_text(
        "\n".join(
            [
                'chr1\tsrc\tgene\t1\t100\t.\t+\t.\tgene_id "g1"; gene_name "G1"; gene_biotype "protein_coding";',
                'chr1\tsrc\ttranscript\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; transcript_biotype "protein_coding";',
                'chr1\tsrc\texon\t1\t50\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; exon_number "1";',
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    report = validate_cellranger_gtf(gtf, require_hierarchy=True)
    assert report.is_pass
    text = _format_report(report, gtf)
    assert "Cell Ranger GTF check: PASS" in text
    assert "____" in text


def test_validate_cellranger_fail_missing_transcript_id_in_exon(tmp_path):
    gtf = tmp_path / "bad.gtf"
    gtf.write_text(
        "\n".join(
            [
                'chr1\tsrc\tgene\t1\t100\t.\t+\t.\tgene_id "g1";',
                'chr1\tsrc\ttranscript\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
                'chr1\tsrc\texon\t1\t50\t.\t+\t.\tgene_id "g1";',
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    report = validate_cellranger_gtf(gtf, require_hierarchy=True)
    assert not report.is_pass
    assert any("exon 缺少 transcript_id" in item for item in report.errors)


def test_validate_cellranger_fail_when_hierarchy_required_but_missing_levels(tmp_path):
    gtf = tmp_path / "flat.gtf"
    gtf.write_text(
        "\n".join(
            [
                'chr1\tsrc\texon\t1\t50\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
                'chr1\tsrc\texon\t70\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    report = validate_cellranger_gtf(gtf, require_hierarchy=True)
    assert not report.is_pass
    assert any("缺少 transcript 特征行" in item for item in report.errors)
    assert any("缺少 gene 特征行" in item for item in report.errors)


def test_validate_cellranger_warn_without_gene_name_by_default(tmp_path):
    gtf = tmp_path / "warn.gtf"
    gtf.write_text(
        "\n".join(
            [
                'chr1\tsrc\tgene\t1\t100\t.\t+\t.\tgene_id "g1"; gene_biotype "protein_coding";',
                'chr1\tsrc\ttranscript\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
                'chr1\tsrc\texon\t1\t50\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    report = validate_cellranger_gtf(gtf)
    assert report.is_pass
    assert any("gene 缺少 gene_name" in item for item in report.warnings)
