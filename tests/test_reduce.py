"""Tests for locus reduction (marginalization / projection).

Split into pure-math unit tests (no I/O) and HFX-in/HFX-out integration tests
covering CSV, Parquet, and inline storage styles.
"""

import json
import zipfile

import pytest

from hfx_tools.io import read_hfx_bundle
from hfx_tools.pack import pack_hfx
from hfx_tools.reduce import (
    frequency_conservation_error,
    haplotype_loci,
    parse_haplotype,
    project_haplotype,
    reduce_frequency_rows,
    reduce_loci_hfx,
    validate_selected_loci,
)

NINE_LOCUS = (
    "A*01:01~C*07:01~B*08:01~DRB3*01:01~DRB1*03:01~DQA1*05:01~DQB1*02:01~DPA1*01:03~DPB1*04:01"
)
NINE_LOCI = ["A", "C", "B", "DRB3", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]

# Variant that differs from NINE_LOCUS only at omitted loci (C, DRB3, DQA1, DPA1,
# DPB1) -> projects to the same A/B/DRB1/DQB1 haplotype and therefore merges.
NINE_LOCUS_MERGE = (
    "A*01:01~C*99:99~B*08:01~DRB3*02:02~DRB1*03:01~DQA1*01:01~DQB1*02:01~DPA1*02:02~DPB1*01:01"
)
# Variant that differs at a kept locus (B) -> stays distinct after projection.
NINE_LOCUS_ALT_B = (
    "A*01:01~C*07:01~B*07:02~DRB3*01:01~DRB1*03:01~DQA1*05:01~DQB1*02:01~DPA1*01:03~DPB1*04:01"
)
# A distinct second-population haplotype.
NINE_LOCUS_POP2 = (
    "A*02:01~C*05:01~B*44:02~DRB3*02:02~DRB1*04:01~DQA1*03:01~DQB1*03:01~DPA1*01:03~DPB1*04:01"
)


# ---------------------------------------------------------------------------
# Pure projection unit tests (no archive I/O)
# ---------------------------------------------------------------------------


def test_parse_haplotype_returns_ordered_locus_allele_pairs():
    assert parse_haplotype("A*01:01~B*08:01~DRB1*03:01") == [
        ("A", "A*01:01"),
        ("B", "B*08:01"),
        ("DRB1", "DRB1*03:01"),
    ]


def test_parse_haplotype_rejects_malformed_and_duplicate():
    with pytest.raises(ValueError):
        parse_haplotype("A0101~B*08:01")  # missing '*'
    with pytest.raises(ValueError):
        parse_haplotype("A*01:01~A*02:01")  # duplicate locus
    with pytest.raises(ValueError):
        parse_haplotype("")  # empty


def test_haplotype_loci():
    assert haplotype_loci(NINE_LOCUS) == NINE_LOCI


def test_project_preserves_requested_order_and_allele_values():
    # request an order different from the source order
    assert (
        project_haplotype(NINE_LOCUS, ["DQB1", "A", "DRB1"])
        == "DQB1*02:01~A*01:01~DRB1*03:01"
    )


def test_project_missing_locus_raises():
    with pytest.raises(ValueError):
        project_haplotype("A*01:01~B*08:01", ["A", "DRB1"])


def test_reduce_aggregates_identical_projections():
    rows = [
        ("A*01:01~C*07:01~B*08:01~DRB1*03:01~DQB1*02:01", 0.040),
        ("A*01:01~C*07:02~B*08:01~DRB1*03:01~DQB1*02:01", 0.025),
    ]
    out = reduce_frequency_rows(rows, ["A", "B", "DRB1", "DQB1"])
    assert out == [("A*01:01~B*08:01~DRB1*03:01~DQB1*02:01", 0.065)]


def test_reduce_keeps_distinct_haplotypes_separate():
    rows = [
        ("A*01:01~B*08:01~DRB1*03:01", 0.04),
        ("A*02:01~B*44:02~DRB1*04:01", 0.03),
    ]
    out = reduce_frequency_rows(rows, ["A", "B", "DRB1"])
    assert dict(out) == {
        "A*01:01~B*08:01~DRB1*03:01": 0.04,
        "A*02:01~B*44:02~DRB1*04:01": 0.03,
    }


def test_reduce_9_to_4():
    rows = [
        (NINE_LOCUS, 0.5),
        # differs only at omitted loci (C, DRB3, DQA1, DPA1, DPB1) -> merges
        (NINE_LOCUS_MERGE, 0.3),
        # differs at a kept locus (B) -> stays separate
        (NINE_LOCUS_ALT_B, 0.2),
    ]
    out = reduce_frequency_rows(rows, ["A", "B", "DRB1", "DQB1"])
    assert dict(out) == {
        "A*01:01~B*08:01~DRB1*03:01~DQB1*02:01": pytest.approx(0.8),
        "A*01:01~B*07:02~DRB1*03:01~DQB1*02:01": pytest.approx(0.2),
    }


def test_reduce_another_subset_9_to_6():
    out = reduce_frequency_rows([(NINE_LOCUS, 1.0)], ["A", "C", "B", "DRB1", "DQB1", "DPB1"])
    assert out == [("A*01:01~C*07:01~B*08:01~DRB1*03:01~DQB1*02:01~DPB1*04:01", 1.0)]


def test_reduce_all_loci_is_scientific_noop():
    rows = [(NINE_LOCUS, 0.6), (NINE_LOCUS_POP2, 0.4)]
    out = reduce_frequency_rows(rows, NINE_LOCI)
    assert out == rows  # order and values unchanged


def test_validate_rejects_duplicate_loci():
    with pytest.raises(ValueError, match="Duplicate"):
        validate_selected_loci(["A", "B", "A"], NINE_LOCI)


def test_validate_rejects_unknown_locus():
    with pytest.raises(ValueError, match="not present"):
        validate_selected_loci(["A", "DRB9"], NINE_LOCI)


def test_validate_rejects_empty():
    with pytest.raises(ValueError):
        validate_selected_loci([], NINE_LOCI)


def test_frequency_conservation_error():
    rows_in = [(NINE_LOCUS, 0.3), (NINE_LOCUS_POP2, 0.7)]
    rows_out = reduce_frequency_rows(rows_in, ["A", "DRB1"])
    assert frequency_conservation_error(rows_in, rows_out) == pytest.approx(0.0, abs=1e-12)


# ---------------------------------------------------------------------------
# Integration fixtures & tests (HFX in -> HFX out)
# ---------------------------------------------------------------------------


def _base_metadata(freq_location, header_map=None):
    md = {
        "outputResolution": [{"locus": locus, "resolution": "g"} for locus in NINE_LOCI],
        "hfeMethod": {"method": "EM", "parameters": [{"parameter": "iterations", "value": "100"}]},
        "cohortDescription": {
            "species": "Homo sapiens",
            "population": [{"name": "TEST", "geoLocation": {"ISO3166": "US"}}],
            "cohortSize": 1000,
        },
        "nomenclatureUsed": {"database": "IPD-IMGT/HLA", "version": "3.57.0"},
        "frequencyLocation": freq_location,
    }
    if header_map is not None:
        md["frequencyFileHeader"] = header_map
    return {"version": "0.1.1", "metadata": md}


SAMPLE_ROWS = [
    (NINE_LOCUS, 0.5),
    (NINE_LOCUS_MERGE, 0.3),
    (NINE_LOCUS_POP2, 0.2),
]


def _build_csv_hfx(tmp_path, header_map=None):
    """Build a CSV-backed .hfx, optionally with non-standard headers."""
    folder = tmp_path / "src_csv"
    folder.mkdir()
    haplo_col, freq_col = "haplotype", "frequency"
    if header_map:
        # header_map maps original->canonical; invert to get output column names
        for orig, canon in header_map.items():
            if canon == "haplotype":
                haplo_col = orig
            elif canon == "frequency":
                freq_col = orig
    csv_path = folder / "freqs.csv"
    lines = [f"{haplo_col},{freq_col}"]
    lines += [f"{h},{f}" for h, f in SAMPLE_ROWS]
    csv_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    meta = _base_metadata("file://freqs.csv", header_map=header_map)
    meta_path = folder / "metadata.json"
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")

    out = tmp_path / "src_csv.hfx"
    pack_hfx(meta_path, out, write_manifest=True, hash_alg="sha256")
    return out


def _build_parquet_hfx(tmp_path):
    pd = pytest.importorskip("pandas")
    pytest.importorskip("pyarrow")
    folder = tmp_path / "src_parquet"
    folder.mkdir()
    pq_path = folder / "freqs.parquet"
    pd.DataFrame(
        {"haplotype": [h for h, _ in SAMPLE_ROWS], "frequency": [f for _, f in SAMPLE_ROWS]}
    ).to_parquet(pq_path, index=False)

    meta = _base_metadata("file://freqs.parquet")
    meta_path = folder / "metadata.json"
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")

    out = tmp_path / "src_parquet.hfx"
    pack_hfx(meta_path, out, write_manifest=True, hash_alg="sha256")
    return out


def _build_inline_hfx(tmp_path):
    folder = tmp_path / "src_inline"
    folder.mkdir()
    meta = _base_metadata("inline")
    meta["frequencyData"] = [{"haplotype": h, "frequency": f} for h, f in SAMPLE_ROWS]
    meta_path = folder / "metadata.json"
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")

    out = tmp_path / "src_inline.hfx"
    pack_hfx(meta_path, out, write_manifest=True, hash_alg="sha256")
    return out


def test_reduce_csv_hfx_roundtrip(tmp_path):
    src = _build_csv_hfx(tmp_path)
    out = tmp_path / "reduced.hfx"
    result = reduce_loci_hfx(src, ["A", "B", "DRB1", "DQB1"], out)

    assert result["storage"] == "csv"
    assert result["frequency_conservation_error"] == pytest.approx(0.0, abs=1e-9)
    assert result["n_output_haplotypes"] == 2  # first two rows merge

    # readback validation
    bundle = read_hfx_bundle(out)
    assert bundle.storage == "csv"
    assert sum(f for _, f in bundle.rows) == pytest.approx(1.0)
    loci = [e["locus"] for e in bundle.metadata["outputResolution"]]
    assert loci == ["A", "B", "DRB1", "DQB1"]

    # provenance recorded in hfeMethod.parameters (schema-legal)
    params = {p["parameter"]: p["value"] for p in bundle.metadata["hfeMethod"]["parameters"]}
    assert params["derivation.selectedLoci"] == "A,B,DRB1,DQB1"
    assert "iterations" in params  # original method params preserved

    # newly generated manifest/checksums present
    with zipfile.ZipFile(out) as z:
        assert "MANIFEST.json" in z.namelist()
        assert "SHA256SUMS" in z.namelist()


def test_reduce_preserves_nonstandard_csv_headers(tmp_path):
    header_map = {"Haplo": "haplotype", "Freq": "frequency"}
    src = _build_csv_hfx(tmp_path, header_map=header_map)
    out = tmp_path / "reduced_hdr.hfx"
    reduce_loci_hfx(src, ["A", "B", "DRB1"], out)

    # The output CSV must use the original non-standard column names.
    with zipfile.ZipFile(out) as z:
        data_name = [n for n in z.namelist() if n.endswith(".csv")][0]
        header_line = z.read(data_name).decode("utf-8").splitlines()[0]
    assert header_line == "Haplo,Freq"

    # And it must still read back correctly via the header mapping.
    bundle = read_hfx_bundle(out)
    assert bundle.metadata["frequencyFileHeader"] == header_map
    assert sum(f for _, f in bundle.rows) == pytest.approx(1.0)


def test_reduce_parquet_hfx_roundtrip(tmp_path):
    src = _build_parquet_hfx(tmp_path)
    out = tmp_path / "reduced.hfx"
    result = reduce_loci_hfx(src, ["A", "B", "DRB1", "DQB1"], out)

    assert result["storage"] == "parquet"
    with zipfile.ZipFile(out) as z:
        assert any(n.endswith(".parquet") for n in z.namelist())
        assert not any(n.endswith(".csv") for n in z.namelist())  # no silent conversion

    bundle = read_hfx_bundle(out)
    assert bundle.storage == "parquet"
    assert sum(f for _, f in bundle.rows) == pytest.approx(1.0)


def test_reduce_inline_hfx_roundtrip(tmp_path):
    src = _build_inline_hfx(tmp_path)
    out = tmp_path / "reduced.hfx"
    result = reduce_loci_hfx(src, ["A", "B", "DRB1", "DQB1"], out)

    assert result["storage"] == "inline"
    with zipfile.ZipFile(out) as z:
        meta = json.loads(z.read("metadata.json"))
    assert meta["metadata"]["frequencyLocation"] == "inline"
    assert "frequencyData" in meta
    assert sum(r["frequency"] for r in meta["frequencyData"]) == pytest.approx(1.0)


def test_reduce_all_loci_noop_produces_valid_bundle(tmp_path):
    src = _build_csv_hfx(tmp_path)
    out = tmp_path / "noop.hfx"
    result = reduce_loci_hfx(src, NINE_LOCI, out)
    assert result["n_input_haplotypes"] == result["n_output_haplotypes"] == 3
    bundle = read_hfx_bundle(out)
    assert sum(f for _, f in bundle.rows) == pytest.approx(1.0)


def test_reduce_invalid_locus_raises(tmp_path):
    src = _build_csv_hfx(tmp_path)
    with pytest.raises(ValueError, match="not present"):
        reduce_loci_hfx(src, ["A", "DRB9"], tmp_path / "bad.hfx")


def test_reduce_duplicate_locus_raises(tmp_path):
    src = _build_csv_hfx(tmp_path)
    with pytest.raises(ValueError, match="Duplicate"):
        reduce_loci_hfx(src, ["A", "A"], tmp_path / "bad.hfx")


def test_reduce_total_frequency_is_conserved(tmp_path):
    src = _build_csv_hfx(tmp_path)
    out = tmp_path / "reduced.hfx"
    result = reduce_loci_hfx(src, ["A", "DRB1"], out)
    assert result["input_frequency_sum"] == pytest.approx(result["output_frequency_sum"])
    assert result["frequency_conservation_error"] <= 1e-9
