from __future__ import annotations

import csv
import json
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

from .util import safe_relpath


def read_hfx_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def write_hfx_json(path: Path, obj: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, sort_keys=False)
        f.write("\n")


def parse_frequency_location(freq_loc: str) -> tuple[str, str | None]:
    """
    Returns (kind, value)
      - ("inline", None)
      - ("file", relative_path) for file://...
      - ("http", uri) for http(s)://... (not bundled in MVP)
    """
    if freq_loc == "inline":
        return ("inline", None)

    u = urlparse(freq_loc)
    if u.scheme in ("http", "https"):
        return ("http", freq_loc)
    if u.scheme == "file":
        # file://data/f.csv -> path = "data/f.csv"
        # urlparse gives netloc + path; handle both
        raw = (u.netloc + u.path).lstrip("/")
        return ("file", safe_relpath(raw))

    # allow plain relative path as a convenience (not strictly per schema format=uri)
    return ("file", safe_relpath(freq_loc))


def _resolve_header_mapping(hfx_obj: dict[str, Any]) -> dict[str, str]:
    """Return {original_col: canonical_col} from metadata.frequencyFileHeader.

    The schema defines frequencyFileHeader as an object mapping CSV header
    names to expected field names (e.g. {"Haplo": "haplotype", "Freq": "frequency"}).
    """
    return hfx_obj.get("metadata", {}).get("frequencyFileHeader", {})


def load_frequency_rows(hfx_path: Path, hfx_obj: dict[str, Any]) -> list[tuple[str, float]]:
    md = hfx_obj.get("metadata", {})
    freq_loc = md.get("frequencyLocation")
    if not freq_loc:
        raise ValueError("metadata.frequencyLocation is required")

    kind, val = parse_frequency_location(freq_loc)
    header_map = _resolve_header_mapping(hfx_obj)

    if kind == "inline":
        rows = hfx_obj.get("frequencyData")
        if rows is None:
            raise ValueError("frequencyLocation is 'inline' but top-level frequencyData is missing")
        out = []
        for r in rows:
            out.append((r["haplotype"], float(r["frequency"])))
        return out

    if kind == "http":
        raise ValueError(
            "http(s) frequencyLocation not supported in MVP loader; "
            "please download locally or bundle with file://"
        )

    # file
    rel = val
    assert rel is not None
    # Resolve relative to parent of metadata/ if inside build folder structure
    base = hfx_path.parent
    if base.name == "metadata":
        base = base.parent
    freq_file = (base / rel).resolve()
    if not freq_file.exists():
        raise FileNotFoundError(f"Referenced frequency file not found: {freq_file}")

    if freq_file.suffix.lower() == ".csv":
        return load_csv(freq_file, header_map=header_map)
    if freq_file.suffix.lower() == ".parquet":
        return load_parquet(freq_file, header_map=header_map)

    raise ValueError(f"Unsupported frequency file type: {freq_file.suffix}")


def _apply_header_map(fieldnames: list, header_map: dict[str, str]) -> dict[str, str]:
    """Build a reverse lookup: {original_col: canonical_col} for the columns we need."""
    # header_map is {csv_col: canonical_col}, e.g. {"Haplo": "haplotype"}
    reverse: dict[str, str] = {}
    for orig, canon in header_map.items():
        if orig in fieldnames:
            reverse[orig] = canon
    return reverse


def load_csv(path: Path, header_map: dict[str, str] | None = None) -> list[tuple[str, float]]:
    out: list[tuple[str, float]] = []
    header_map = header_map or {}
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        mapping = _apply_header_map(reader.fieldnames or [], header_map)
        # Determine which column names to use for haplotype and frequency
        haplo_col = "haplotype"
        freq_col = "frequency"
        for orig, canon in mapping.items():
            if canon == "haplotype":
                haplo_col = orig
            elif canon == "frequency":
                freq_col = orig
        if haplo_col not in reader.fieldnames or freq_col not in reader.fieldnames:
            raise ValueError(
                f"CSV must have columns haplotype,frequency (or mapped via frequencyFileHeader); "
                f"found {reader.fieldnames}"
            )
        for row in reader:
            out.append((row[haplo_col], float(row[freq_col])))
    return out


def load_parquet(path: Path, header_map: dict[str, str] | None = None) -> list[tuple[str, float]]:
    try:
        import pandas as pd  # type: ignore
    except Exception as e:
        raise ImportError(
            "Parquet support requires pandas + pyarrow. "
            "Install with: pip install 'hfx-tools[parquet]'"
        ) from e

    header_map = header_map or {}
    df = pd.read_parquet(path)
    # Apply header mapping: rename columns to canonical names
    rename = {orig: canon for orig, canon in header_map.items() if orig in df.columns}
    if rename:
        df = df.rename(columns=rename)
    if "haplotype" not in df.columns or "frequency" not in df.columns:
        raise ValueError(
            f"Parquet must have columns haplotype,frequency (or mapped via frequencyFileHeader); "
            f"found {list(df.columns)}"
        )
    return [(str(h), float(f)) for h, f in zip(df["haplotype"], df["frequency"], strict=True)]


# ---------------------------------------------------------------------------
# Reusable HFX bundle access
#
# The helpers below centralise reading a bundled ``.hfx`` archive (or a bare
# metadata JSON) into a normalised in-memory representation, and writing
# frequency rows back out in a chosen storage style. They are intended to be
# shared infrastructure for HFX transformation tools (e.g. locus reduction),
# so that data access stays in one place rather than being reimplemented.
# ---------------------------------------------------------------------------

# Storage styles for frequency data.
STORAGE_INLINE = "inline"
STORAGE_CSV = "csv"
STORAGE_PARQUET = "parquet"


@dataclass
class HfxBundle:
    """Normalised, in-memory view of an HFX document and its frequency data.

    Attributes:
        hfx_obj: The parsed HFX JSON document (top-level, includes ``metadata``).
        rows: Frequency rows as ``(haplotype, frequency)`` tuples.
        storage: One of ``STORAGE_INLINE``, ``STORAGE_CSV``, ``STORAGE_PARQUET``.
        data_filename: For file-backed data, the frequency file name inside the
            bundle (e.g. ``"HEGL_AFA.csv"``); ``None`` for inline data.
        header_map: The ``metadata.frequencyFileHeader`` mapping (may be empty).
        source_path: Path the bundle/metadata was read from.
    """

    hfx_obj: dict[str, Any]
    rows: list[tuple[str, float]]
    storage: str
    data_filename: str | None
    header_map: dict[str, str] = field(default_factory=dict)
    source_path: Path | None = None

    @property
    def metadata(self) -> dict[str, Any]:
        return self.hfx_obj.setdefault("metadata", {})


def _storage_from_location(freq_loc: str, data_filename: str | None) -> str:
    """Classify the storage style from a frequencyLocation and file name."""
    kind, _ = parse_frequency_location(freq_loc)
    if kind == "inline":
        return STORAGE_INLINE
    if kind == "file" and data_filename:
        suffix = Path(data_filename).suffix.lower()
        if suffix == ".parquet":
            return STORAGE_PARQUET
        if suffix == ".csv":
            return STORAGE_CSV
    raise ValueError(
        f"Unsupported or unresolved frequency storage for location '{freq_loc}'"
        f" (file '{data_filename}')"
    )


def read_hfx_bundle(path: Path) -> HfxBundle:
    """Read a ``.hfx`` archive or a metadata JSON into an :class:`HfxBundle`.

    This is the generic entry point for HFX transformations that need access to
    both the metadata document and the underlying frequency rows, regardless of
    whether frequencies are stored inline, as CSV, or as Parquet.

    Args:
        path: Path to a ``.hfx`` bundle or a metadata JSON file.

    Returns:
        An :class:`HfxBundle` with the parsed document, frequency rows, and the
        detected storage style.
    """
    path = Path(path)

    if path.suffix.lower() == ".hfx":
        return _read_hfx_archive(path)
    return _read_hfx_metadata_file(path)


def _read_hfx_metadata_file(path: Path) -> HfxBundle:
    hfx_obj = read_hfx_json(path)
    md = hfx_obj.get("metadata", {})
    freq_loc = md.get("frequencyLocation")
    if not freq_loc:
        raise ValueError("metadata.frequencyLocation is required")
    header_map = _resolve_header_mapping(hfx_obj)
    rows = load_frequency_rows(path, hfx_obj)

    kind, rel = parse_frequency_location(freq_loc)
    data_filename = Path(rel).name if (kind == "file" and rel) else None
    storage = _storage_from_location(freq_loc, data_filename)
    return HfxBundle(
        hfx_obj=hfx_obj,
        rows=rows,
        storage=storage,
        data_filename=data_filename,
        header_map=header_map,
        source_path=path,
    )


def _read_hfx_archive(path: Path) -> HfxBundle:
    with zipfile.ZipFile(path, "r") as z:
        names = z.namelist()
        if "metadata.json" not in names:
            raise ValueError(f"No metadata.json found inside {path}")
        hfx_obj = json.loads(z.read("metadata.json").decode("utf-8"))
        md = hfx_obj.get("metadata", {})
        freq_loc = md.get("frequencyLocation")
        if not freq_loc:
            raise ValueError("metadata.frequencyLocation is required")
        header_map = _resolve_header_mapping(hfx_obj)
        kind, rel = parse_frequency_location(freq_loc)

        if kind == "inline":
            rows = hfx_obj.get("frequencyData")
            if rows is None:
                raise ValueError(
                    "frequencyLocation is 'inline' but top-level frequencyData is missing"
                )
            out_rows = [(r["haplotype"], float(r["frequency"])) for r in rows]
            return HfxBundle(
                hfx_obj=hfx_obj,
                rows=out_rows,
                storage=STORAGE_INLINE,
                data_filename=None,
                header_map=header_map,
                source_path=path,
            )

        if kind == "http":
            raise ValueError(
                "http(s) frequencyLocation is not supported for bundled transformation; "
                "download locally or bundle with file://"
            )

        # file-backed: the data file is stored at the top level of the archive
        data_filename = Path(rel).name if rel else None
        if not data_filename or data_filename not in names:
            raise FileNotFoundError(
                f"Referenced frequency file '{rel}' not found inside bundle {path}"
            )
        raw = z.read(data_filename)
        suffix = Path(data_filename).suffix.lower()
        if suffix == ".csv":
            rows = _rows_from_csv_bytes(raw, header_map)
            storage = STORAGE_CSV
        elif suffix == ".parquet":
            rows = _rows_from_parquet_bytes(raw, header_map)
            storage = STORAGE_PARQUET
        else:
            raise ValueError(f"Unsupported frequency file type in bundle: {suffix}")

        return HfxBundle(
            hfx_obj=hfx_obj,
            rows=rows,
            storage=storage,
            data_filename=data_filename,
            header_map=header_map,
            source_path=path,
        )


def _rows_from_csv_bytes(raw: bytes, header_map: dict[str, str]) -> list[tuple[str, float]]:
    text = raw.decode("utf-8")
    reader = csv.DictReader(text.splitlines())
    mapping = _apply_header_map(reader.fieldnames or [], header_map)
    haplo_col, freq_col = "haplotype", "frequency"
    for orig, canon in mapping.items():
        if canon == "haplotype":
            haplo_col = orig
        elif canon == "frequency":
            freq_col = orig
    if haplo_col not in (reader.fieldnames or []) or freq_col not in (reader.fieldnames or []):
        raise ValueError(
            f"CSV must have columns haplotype,frequency (or mapped via frequencyFileHeader); "
            f"found {reader.fieldnames}"
        )
    return [(row[haplo_col], float(row[freq_col])) for row in reader]


def _rows_from_parquet_bytes(raw: bytes, header_map: dict[str, str]) -> list[tuple[str, float]]:
    try:
        import pandas as pd  # type: ignore
    except Exception as e:  # pragma: no cover - exercised only without pandas
        raise ImportError(
            "Parquet support requires pandas + pyarrow. Install with: pip install -e '.[parquet]'"
        ) from e
    import io as _io

    df = pd.read_parquet(_io.BytesIO(raw))
    rename = {orig: canon for orig, canon in header_map.items() if orig in df.columns}
    if rename:
        df = df.rename(columns=rename)
    if "haplotype" not in df.columns or "frequency" not in df.columns:
        raise ValueError(
            f"Parquet must have columns haplotype,frequency (or mapped via frequencyFileHeader); "
            f"found {list(df.columns)}"
        )
    return [(str(h), float(f)) for h, f in zip(df["haplotype"], df["frequency"], strict=False)]


def _canonical_columns(header_map: dict[str, str]) -> tuple[str, str]:
    """Return the (haplotype_column, frequency_column) names to write.

    If a ``frequencyFileHeader`` mapping is present it maps original column
    names to canonical fields; we invert it so that output files preserve the
    original (non-standard) column headers exactly.
    """
    haplo_col, freq_col = "haplotype", "frequency"
    for orig, canon in header_map.items():
        if canon == "haplotype":
            haplo_col = orig
        elif canon == "frequency":
            freq_col = orig
    return haplo_col, freq_col


def write_frequency_file(
    rows: list[tuple[str, float]],
    dest: Path,
    storage: str,
    header_map: dict[str, str] | None = None,
) -> None:
    """Write frequency rows to ``dest`` in the requested storage style.

    Preserves the original column headers implied by ``header_map`` so that the
    ``metadata.frequencyFileHeader`` contract continues to hold on readback.
    """
    header_map = header_map or {}
    haplo_col, freq_col = _canonical_columns(header_map)
    dest.parent.mkdir(parents=True, exist_ok=True)

    if storage == STORAGE_CSV:
        with dest.open("w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([haplo_col, freq_col])
            for haplo, freq in rows:
                writer.writerow([haplo, freq])
    elif storage == STORAGE_PARQUET:
        try:
            import pandas as pd  # type: ignore
        except Exception as e:  # pragma: no cover - exercised only without pandas
            raise ImportError(
                "Parquet support requires pandas + pyarrow. "
                "Install with: pip install -e '.[parquet]'"
            ) from e
        df = pd.DataFrame(
            {haplo_col: [h for h, _ in rows], freq_col: [f for _, f in rows]}
        )
        df.to_parquet(dest, index=False)
    else:
        raise ValueError(f"write_frequency_file does not support storage style '{storage}'")


def set_inline_frequency_data(hfx_obj: dict[str, Any], rows: list[tuple[str, float]]) -> None:
    """Store frequency rows inline on the HFX document per the schema.

    Sets ``metadata.frequencyLocation = "inline"`` and populates top-level
    ``frequencyData`` with ``{haplotype, frequency}`` objects.
    """
    hfx_obj.setdefault("metadata", {})["frequencyLocation"] = "inline"
    hfx_obj["frequencyData"] = [
        {"haplotype": haplo, "frequency": freq} for haplo, freq in rows
    ]

