from gtf_gff_transition.cli import main


def test_cli_convert_file_to_file(tmp_path):
    input_file = tmp_path / "input.gtf"
    output_file = tmp_path / "output.gff3"

    input_file.write_text(
        "chr1\tsrc\texon\t1\t100\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n",
        encoding="utf-8",
    )

    rc = main(
        [
            "-i",
            str(input_file),
            "-o",
            str(output_file),
            "--from",
            "gtf",
            "--to",
            "gff3",
        ]
    )

    assert rc == 0
    content = output_file.read_text(encoding="utf-8")
    assert content.startswith("##gff-version 3\n")
    assert "Parent=t1" in content
