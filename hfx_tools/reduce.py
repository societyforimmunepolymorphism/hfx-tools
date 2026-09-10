"""Locus reduction (marginalization / projection) for HFX frequency data.

This module projects a multilocus haplotype-frequency dataset onto a selected
subset of loci by *marginalizing out* the omitted loci: the frequency of a
reduced haplotype is the sum of the frequencies of all source haplotypes whose
projection onto the selected loci is identical.

Mathematically, for a set of selected loci ``L``::

    f(short) = sum{ f(long) : project(long, L) == short }

This is locus projection, not allele-resolution reduction: allele values are
preserved exactly.

The module is split into two layers:

* Pure functions (:func:`parse_haplotype`, :func:`project_haplotype`,
  :func:`reduce_frequency_rows`) implement the projection math and are fully
  unit-testable without any archive I/O.
* :func:`reduce_loci_hfx` orchestrates an HFX-in/HFX-out transformation, reusing
  the shared :mod:`hfx_tools.io` bundle access and :mod:`hfx_tools.pack` packing
  infrastructure, and updating metadata (locus definition + provenance) so the
  result is a valid, self-describing HFX bundle.
"""

from __future__ import annotations

import tempfile
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .io import (
    STORAGE_CSV,
    STORAGE_INLINE,
    STORAGE_PARQUET,
    HfxBundle,
    read_hfx_bundle,
    set_inline_frequency_data,
    write_frequency_file,
    write_hfx_json,
)
from .pack import pack_hfx

# Delimiter between alleles in an HFX/GL-string haplotype.
ALLELE_DELIM = "~"
# Delimiter between a locus and its allele designation, e.g. ``A*01:01``.
LOCUS_DELIM = "*"

# Default tolerance for the frequency-conservation check.
DEFAULT_FREQ_TOLERANCE = 1e-6


# ---------------------------------------------------------------------------
# Pure projection logic (no I/O)
# ---------------------------------------------------------------------------


def parse_haplotype(haplotype: str) -> list[tuple[str, str]]:
    """Parse an HFX haplotype string into ordered ``(locus, allele_token)`` pairs.

    A haplotype is a ``~``-delimited list of allele designations, each of the
    form ``LOCUS*ALLELE`` (e.g. ``A*01:01~B*08:01``). The locus is the substring
    before the first ``*``; the full token (including the locus prefix) is
    retained as the allele value so that projection preserves the exact source
    representation.

    Args:
        haplotype: The haplotype string.

    Returns:
        A list of ``(locus, allele_token)`` pairs in source order.

    Raises:
        ValueError: If the haplotype is empty, an allele lacks a ``*`` locus
            delimiter, or the same locus appears more than once.
    """
    if not haplotype or not haplotype.strip():
        raise ValueError("Empty haplotype string")

    pairs: list[tuple[str, str]] = []
    seen: set[str] = set()
    for token in haplotype.split(ALLELE_DELIM):
        if LOCUS_DELIM not in token:
            raise ValueError(
                f"Malformed allele '{token}' in haplotype '{haplotype}': "
                f"expected LOCUS{LOCUS_DELIM}ALLELE"
            )
        locus = token.split(LOCUS_DELIM, 1)[0]
        if not locus:
            raise ValueError(
                f"Malformed allele '{token}' in haplotype '{haplotype}': empty locus"
            )
        if locus in seen:
            raise ValueError(
                f"Locus '{locus}' appears more than once in haplotype '{haplotype}'"
            )
        seen.add(locus)
        pairs.append((locus, token))
    return pairs


def haplotype_loci(haplotype: str) -> list[str]:
    """Return the ordered list of loci present in a haplotype string."""
    return [locus for locus, _ in parse_haplotype(haplotype)]


def project_haplotype(haplotype: str, selected_loci: Sequence[str]) -> str:
    """Project a haplotype onto ``selected_loci``, preserving that locus order.

    Args:
        haplotype: The source haplotype string (``LOCUS*ALLELE`` tokens joined
            by ``~``).
        selected_loci: The loci to retain, in the desired output order.

    Returns:
        The reduced haplotype string containing only the selected loci, in the
        requested order.

    Raises:
        ValueError: If any requested locus is absent from the haplotype.
    """
    locus_to_token = dict(parse_haplotype(haplotype))
    projected: list[str] = []
    for locus in selected_loci:
        if locus not in locus_to_token:
            raise ValueError(
                f"Requested locus '{locus}' not present in haplotype '{haplotype}' "
                f"(available: {sorted(locus_to_token)})"
            )
        projected.append(locus_to_token[locus])
    return ALLELE_DELIM.join(projected)


def validate_selected_loci(
    selected_loci: Sequence[str], available_loci: Sequence[str]
) -> list[str]:
    """Validate a requested locus selection against the available loci.

    Args:
        selected_loci: The user-requested loci (order preserved).
        available_loci: The loci available in the source dataset.

    Returns:
        The selected loci as a list (unchanged order).

    Raises:
        ValueError: If ``selected_loci`` is empty, contains duplicates, or
            requests a locus not present in ``available_loci``.
    """
    if not selected_loci:
        raise ValueError("At least one locus must be selected")

    seen: set[str] = set()
    duplicates: list[str] = []
    for locus in selected_loci:
        if locus in seen:
            duplicates.append(locus)
        seen.add(locus)
    if duplicates:
        raise ValueError(f"Duplicate loci requested: {sorted(set(duplicates))}")

    available = set(available_loci)
    missing = [locus for locus in selected_loci if locus not in available]
    if missing:
        raise ValueError(
            f"Requested loci not present in source dataset: {missing} "
            f"(available: {sorted(available)})"
        )
    return list(selected_loci)


def reduce_frequency_rows(
    rows: Sequence[tuple[str, float]], selected_loci: Sequence[str]
) -> list[tuple[str, float]]:
    """Marginalize frequency rows onto ``selected_loci``.

    Every source haplotype is projected onto the selected loci; rows that map to
    the same reduced haplotype have their frequencies summed. The output order
    is deterministic: reduced haplotypes appear in the order their first
    contributing source row was encountered.

    Args:
        rows: Source ``(haplotype, frequency)`` rows.
        selected_loci: The loci to retain, in output order.

    Returns:
        The reduced ``(haplotype, frequency)`` rows.

    Raises:
        ValueError: If a requested locus is absent from any source haplotype
            (surfaced from :func:`project_haplotype`).
    """
    aggregated: dict[str, float] = {}
    order: list[str] = []
    for haplotype, freq in rows:
        short = project_haplotype(haplotype, selected_loci)
        if short not in aggregated:
            aggregated[short] = 0.0
            order.append(short)
        aggregated[short] += float(freq)
    return [(short, aggregated[short]) for short in order]


def frequency_conservation_error(
    input_rows: Sequence[tuple[str, float]],
    output_rows: Sequence[tuple[str, float]],
) -> float:
    """Return the absolute difference between input and output frequency sums."""
    return abs(
        sum(float(f) for _, f in input_rows) - sum(float(f) for _, f in output_rows)
    )


# ---------------------------------------------------------------------------
# Metadata handling
# ---------------------------------------------------------------------------


def _update_output_resolution(
    metadata: dict[str, Any], selected_loci: Sequence[str]
) -> None:
    """Restrict ``metadata.outputResolution`` to the selected loci, in order.

    Existing per-locus resolution entries are preserved for retained loci. Loci
    without an existing entry are omitted from ``outputResolution`` (the entry
    is only meaningful when a resolution was declared for that locus).
    """
    existing = metadata.get("outputResolution", []) or []
    by_locus = {
        entry.get("locus"): entry
        for entry in existing
        if isinstance(entry, dict) and entry.get("locus")
    }
    new_entries = [by_locus[locus] for locus in selected_loci if locus in by_locus]
    if new_entries:
        metadata["outputResolution"] = new_entries


def _record_projection_provenance(
    metadata: dict[str, Any],
    selected_loci: Sequence[str],
    source_loci: Sequence[str],
    source_identifier: str | None,
) -> None:
    """Record the locus-projection derivation in ``metadata.hfeMethod``.

    The HFX schema forbids additional metadata properties and provides no
    dedicated provenance field, so the standard ``hfeMethod.parameters`` array
    (``{parameter, value}`` items) is used to record that this dataset was
    derived by locus projection. The original ``hfeMethod.method`` is preserved.
    """
    hfe = metadata.get("hfeMethod")
    if not isinstance(hfe, dict):
        hfe = {"method": "unspecified"}
    hfe.setdefault("method", "unspecified")
    params = list(hfe.get("parameters", []) or [])

    params.append(
        {
            "parameter": "derivation",
            "value": "locus-projection (marginalization over omitted loci)",
        }
    )
    params.append(
        {"parameter": "derivation.selectedLoci", "value": ",".join(selected_loci)}
    )
    params.append(
        {"parameter": "derivation.sourceLoci", "value": ",".join(source_loci)}
    )
    params.append(
        {
            "parameter": "derivation.timestamp",
            "value": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        }
    )
    if source_identifier:
        params.append(
            {"parameter": "derivation.source", "value": source_identifier}
        )

    hfe["parameters"] = params
    metadata["hfeMethod"] = hfe


# ---------------------------------------------------------------------------
# Orchestration: HFX in -> HFX out
# ---------------------------------------------------------------------------


class ReduceResult(dict):
    """Result mapping returned by :func:`reduce_loci_hfx`.

    Keys:
        output_path: Path to the written ``.hfx`` bundle.
        selected_loci: The loci retained (output order).
        source_loci: The loci detected in the source dataset.
        n_input_haplotypes / n_output_haplotypes: Row counts.
        input_frequency_sum / output_frequency_sum: Total frequencies.
        frequency_conservation_error: Absolute difference of the two sums.
        storage: The storage style preserved (inline/csv/parquet).
    """


def _detect_source_loci(bundle: HfxBundle) -> list[str]:
    """Determine the ordered locus set for the source dataset.

    Prefers ``metadata.outputResolution`` when it is present and consistent;
    otherwise falls back to the loci observed in the first frequency row.
    """
    declared = [
        entry.get("locus")
        for entry in (bundle.metadata.get("outputResolution", []) or [])
        if isinstance(entry, dict) and entry.get("locus")
    ]
    if declared:
        return declared
    if bundle.rows:
        return haplotype_loci(bundle.rows[0][0])
    raise ValueError("Cannot determine source loci: no outputResolution and no data")


def reduce_loci_hfx(
    input_path: Path,
    selected_loci: Sequence[str],
    output_path: Path,
    write_manifest: bool = True,
    hash_alg: str | None = "sha256",
    tolerance: float = DEFAULT_FREQ_TOLERANCE,
    submission_identity: Mapping[str, Any] | None = None,
) -> ReduceResult:
    """Reduce an HFX bundle to a subset of loci and write a new HFX bundle.

    The transformation reads ``input_path`` (a ``.hfx`` bundle or a metadata
    JSON), marginalizes frequencies onto ``selected_loci`` (see
    :func:`reduce_frequency_rows`), preserves the source storage style
    (inline/CSV/Parquet) and header conventions, updates the locus definition
    and records the projection in provenance, then packs a fresh, validated
    ``.hfx`` with newly generated manifests/checksums.

    Args:
        input_path: Source ``.hfx`` bundle or metadata JSON.
        selected_loci: Loci to retain, in the desired output order.
        output_path: Destination ``.hfx`` path.
        write_manifest: Whether to embed ``MANIFEST.json`` in the new bundle.
        hash_alg: Checksum algorithm for the new bundle (``"md5"``/``"sha256"``/``None``).
        tolerance: Maximum allowed absolute deviation between input and output
            total frequency (frequency-conservation check).
        submission_identity: Optional identity embedded into the output bundle
            (passed through to :func:`hfx_tools.pack.pack_hfx`).

    Returns:
        A :class:`ReduceResult` describing the transformation.

    Raises:
        ValueError: On invalid/duplicate/missing loci, or if frequency is not
            conserved within ``tolerance``.
    """
    input_path = Path(input_path)
    output_path = Path(output_path)

    bundle = read_hfx_bundle(input_path)
    source_loci = _detect_source_loci(bundle)
    selected = validate_selected_loci(selected_loci, source_loci)

    reduced_rows = reduce_frequency_rows(bundle.rows, selected)

    conservation_error = frequency_conservation_error(bundle.rows, reduced_rows)
    if conservation_error > tolerance:
        raise ValueError(
            f"Frequency not conserved: input sum "
            f"{sum(f for _, f in bundle.rows):.10f} vs output sum "
            f"{sum(f for _, f in reduced_rows):.10f} "
            f"(deviation {conservation_error:.3g} > tolerance {tolerance:g})"
        )

    # Update metadata: preserve everything valid, adjust the locus definition,
    # and record the projection provenance in a schema-legal way.
    metadata = bundle.metadata
    _update_output_resolution(metadata, selected)
    _record_projection_provenance(
        metadata,
        selected_loci=selected,
        source_loci=source_loci,
        source_identifier=input_path.name,
    )

    # Write the transformed bundle, preserving the source storage style.
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        work = Path(tmpdir)
        header_map = bundle.header_map

        if bundle.storage == STORAGE_INLINE:
            set_inline_frequency_data(bundle.hfx_obj, reduced_rows)
        elif bundle.storage in (STORAGE_CSV, STORAGE_PARQUET):
            # Keep the original data file name so downstream references are stable.
            data_name = bundle.data_filename or (
                "frequencies.csv" if bundle.storage == STORAGE_CSV else "frequencies.parquet"
            )
            data_path = work / data_name
            write_frequency_file(reduced_rows, data_path, bundle.storage, header_map)
            metadata["frequencyLocation"] = f"file://{data_name}"
        else:  # pragma: no cover - defensive
            raise ValueError(f"Unsupported storage style: {bundle.storage}")

        # pack_hfx reads a metadata JSON on disk and bundles any file:// data
        # located next to it; write the (mutated) metadata into the work dir.
        meta_path = work / "metadata.json"
        write_hfx_json(meta_path, bundle.hfx_obj)

        pack_hfx(
            metadata_json=meta_path,
            out_path=output_path,
            write_manifest=write_manifest,
            hash_alg=hash_alg,
            submission_identity=submission_identity,
        )

    result = ReduceResult()
    result.update(
        {
            "output_path": str(output_path),
            "selected_loci": list(selected),
            "source_loci": list(source_loci),
            "n_input_haplotypes": len(bundle.rows),
            "n_output_haplotypes": len(reduced_rows),
            "input_frequency_sum": float(sum(f for _, f in bundle.rows)),
            "output_frequency_sum": float(sum(f for _, f in reduced_rows)),
            "frequency_conservation_error": float(conservation_error),
            "storage": bundle.storage,
        }
    )
    return result
