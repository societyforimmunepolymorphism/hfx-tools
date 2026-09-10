import json

from hfx_tools.cli import main


def test_generic_inspect_subcommand_accepts_an_hfx_document(tmp_path, capsys):
    metadata_path = tmp_path / "metadata.json"
    metadata_path.write_text(
        json.dumps({"version": "0.1.1", "metadata": {"frequencyLocation": "inline"}}),
        encoding="utf-8",
    )

    main(["inspect", str(metadata_path)])

    assert "frequencyLocation: inline" in capsys.readouterr().out


def test_reduce_loci_subcommand_writes_reduced_bundle(tmp_path, capsys):
    from hfx_tools.io import read_hfx_bundle
    from hfx_tools.pack import pack_hfx

    folder = tmp_path / "src"
    folder.mkdir()
    (folder / "freqs.csv").write_text(
        "haplotype,frequency\n"
        "A*01:01~B*08:01~DRB1*03:01,0.6\n"
        "A*01:01~B*08:01~DRB1*04:01,0.4\n",
        encoding="utf-8",
    )
    meta = {
        "version": "0.1.1",
        "metadata": {
            "outputResolution": [
                {"locus": "A", "resolution": "g"},
                {"locus": "B", "resolution": "g"},
                {"locus": "DRB1", "resolution": "g"},
            ],
            "hfeMethod": {"method": "EM", "parameters": []},
            "cohortDescription": {
                "species": "Homo sapiens",
                "population": [{"name": "TEST", "geoLocation": {"ISO3166": "US"}}],
                "cohortSize": 10,
            },
            "nomenclatureUsed": {"database": "IPD-IMGT/HLA", "version": "3.57.0"},
            "frequencyLocation": "file://freqs.csv",
        },
    }
    meta_path = folder / "metadata.json"
    meta_path.write_text(json.dumps(meta), encoding="utf-8")
    src = tmp_path / "src.hfx"
    pack_hfx(meta_path, src, write_manifest=True, hash_alg="sha256")

    out = tmp_path / "reduced.hfx"
    main(["reduce-loci", str(src), "--loci", "A", "B", "-o", str(out)])

    assert out.exists()
    assert "Wrote:" in capsys.readouterr().out
    bundle = read_hfx_bundle(out)
    # both source rows share A*01:01~B*08:01 -> single reduced haplotype summing to 1.0
    assert bundle.rows == [("A*01:01~B*08:01", 1.0)]
