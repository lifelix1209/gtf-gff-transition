import gtf_gff_transition as pkg


def test_public_api_exports_expected_symbols():
    assert hasattr(pkg, "__version__")
    assert callable(pkg.convert_text)
    assert callable(pkg.convert_gtf_to_gff3)
    assert callable(pkg.convert_gff3_to_gtf)
    assert callable(pkg.build_cellranger_gtf)
    assert callable(pkg.validate_cellranger_gtf)
