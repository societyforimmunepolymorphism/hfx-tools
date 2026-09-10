# hfx-tools

Tools for working with HFX submissions (Haplotype Frequency Exchange).

This repo provides composable command line tools and a Streamlit app for building, packing,
inspecting, and validating HFX documents, implementing the [HFX specification](https://github.com/societyforimmunepolymorphism/hfx). Key features include:

- **`build`** - Build HFX bundles from a folder with automatic validation
- **`pack`** - Pack HFX archives from metadata.json with optional manifests and checksums
- **`qc`** - Compute quality control statistics
- **`inspect`** - Inspect metadata or bundled HFX files
- **`reduce-loci`** - Project an HFX onto a subset of loci by marginalizing out omitted loci
- **Validation framework** - Extensible validation with built-in validators
- **Streamlit UI** - Web-based interface for building HFX files

## Key schema facts

- `metadata.frequencyLocation` controls where frequencies are stored: either `"inline"`
  or a URI (e.g., `file://frequencies.csv`) (see [HFX specification](https://github.com/societyforimmunepolymorphism/hfx)).

- If inline, the JSON may include `frequencyData` (array of `{haplotype, frequency}`).

- `metadata.frequencyFileHeader` maps CSV column names to the expected `haplotype`/`frequency`
  field names when the data file uses non-standard headers.

## Install

### Basic installation

```bash
pip install hfx-tools
```

### With optional dependencies

```bash
# For Parquet support
pip install "hfx-tools[parquet]"

# For Streamlit web UI
pip install "hfx-tools[streamlit]"
```

## Quick Start

**5-minute walkthrough** for the most common workflow:

```bash
# 1. Create input folder with metadata and data at top level
mkdir my_submission
cp my_metadata.json my_submission/
cp my_frequencies.csv my_submission/

# 2. Build and validate
hfx-build my_submission -n my_hfx_file

# 3. Done! Check output
ls -la my_submission/my_hfx_file.hfx
cat my_submission/my_hfx_file.build.log
```

For a guided interactive experience, launch the Streamlit web UI:
```bash
streamlit run "$(python -c 'import hfx_tools.streamlit_app as app; print(app.__file__)')"
```

## Architecture

**hfx-tools** follows a layered architecture:

```
CLI / Streamlit UI (user-facing)
    ↓
build.py (orchestration)
    ↓
validators.py (validation rules) ← pack.py (packing logic)
    ↓
io.py (file I/O, JSON parsing)
```

- **CLI layer** (`cli.py`) - Parses command-line arguments and delegates to build/pack/inspect/qc
- **Build orchestration** (`build.py`) - High-level workflow: reads metadata → detects files → validates → packs
- **Validation framework** (`validators.py`) - Pluggable validators for extensibility
- **Packing logic** (`pack.py`) - Low-level archive creation (ZIP with metadata, data, optional manifest)
- **I/O utilities** (`io.py`) - JSON parsing, file reading, consistent error handling

This design allows hackathon participants to:
1. Use the CLI for quick workflows
2. Call `build()` directly from Python for programmatic use
3. Register custom validators without modifying core code
4. Extend with custom QC statistics

## Usage

### Frequency Location Types

The HFX standard supports four types of frequency data locations:

1. **Inline** - `"frequencyLocation": "inline"` with `frequencyData` array in same JSON
2. **Remote** - `"frequencyLocation": "https://zenodo.org/.../data.csv"` or S3 URL
3. **File (CSV)** - `"frequencyLocation": "file://frequencies.csv"` pointing to file within HFX bundle
4. **File (Parquet)** - Same as above but with `.parquet` extension

### CLI: Build from folder

The most common workflow:

```bash
hfx-build /path/to/input_folder -n output_name
```

Expected folder structure:
```
input_folder/
├── metadata.json      # Required: HFX metadata (optionally with inline frequencyData)
└── frequencies.csv    # Optional: if frequencyLocation = "file://frequencies.csv"
```

This will:
1. Read `metadata.json` from `input_folder/`
2. **Auto-detect** a frequency data file in the same folder
3. **Auto-update** `metadata.frequencyLocation` to `file://<filename>` (unless already set to remote or inline)
4. **Validate** all data with built-in validators
5. **Pack** into a single `output_name.hfx` file
6. **Log** all validation results to `output_name.build.log`

Example:
```bash
cp metadata.json example/
cp frequencies.csv example/
hfx-build example -n my_submission
# Output: example/my_submission.hfx
```

Options:
- `-n, --name NAME` - Output filename (required, without .hfx)
- `-o, --out DIR` - Output directory (defaults to input folder)
- `--no-manifest` - Skip MANIFEST.json in archive
- `--hash {md5,sha256,none}` - Hash algorithm (default: sha256)
- `--no-auto-update-location` - Don't auto-update `metadata.frequencyLocation` (advanced)

### CLI: Pack (low-level)

For direct packing when you already have a metadata.json:

```bash
hfx-pack metadata.json -o dist/example.hfx --manifest --hash sha256
```

### CLI: Inspect

```bash
hfx-inspect metadata.json       # Inspect a metadata.json file
hfx-inspect example.hfx         # Inspect a bundled .hfx archive
hfx-tools inspect example.hfx   # Equivalent generic command
```

### CLI: QC

```bash
hfx-qc metadata.json --write-metadata --topk 10 100 1000
```

### CLI: Reduce loci (marginalization)

Project a multilocus HFX onto a subset of loci. Omitted loci are
**marginalized out**: the frequency of each reduced haplotype is the **sum** of
the frequencies of all source haplotypes that map to the same reduced
haplotype.

```bash
hfx-tools reduce-loci nine_locus.hfx \
    --loci A B DRB1 DQB1 \
    -o four_locus.hfx

# standalone convenience command (equivalent):
hfx-reduce-loci nine_locus.hfx --loci A B DRB1 DQB1 -o four_locus.hfx
```

For example, two nine-locus haplotypes that differ only at omitted loci collapse
into one four-locus haplotype whose frequency is the sum:

```
A*01:01~C*07:01~B*08:01~DRB1*03:01~DQB1*02:01   0.040
A*01:01~C*07:02~B*08:01~DRB1*03:01~DQB1*02:01   0.025
```

reduced with `--loci A B DRB1 DQB1` becomes:

```
A*01:01~B*08:01~DRB1*03:01~DQB1*02:01           0.065
```

Key properties:

- **Arbitrary reductions** - any subset of the source loci is supported
  (9→8, 9→6, 9→4, 6→3, 4→2, ...). Requesting the full source locus set is a
  scientific no-op (a new bundle is still produced).
- **Locus order** - the output haplotypes follow the order given in `--loci`.
- **Frequency conservation** - total frequency is preserved and checked; the
  build fails if the input and output sums differ by more than `--tolerance`
  (default `1e-6`). Frequencies are **not** renormalized.
- **Storage style preserved** - CSV→CSV, Parquet→Parquet, inline→inline. There
  is no silent conversion between storage styles, and non-standard column
  headers declared via `metadata.frequencyFileHeader` are preserved.
- **Metadata & provenance** - cohort, nomenclature, methodology, producer,
  license, and header mappings are preserved. `metadata.outputResolution` is
  restricted to the selected loci, and the projection is recorded in
  `metadata.hfeMethod.parameters` (the schema provides no dedicated provenance
  field, and forbids additional properties).
- **Fresh manifests** - the new bundle gets newly generated `MANIFEST.json` and
  checksums describing the reduced files.

Options:
- `--loci LOCUS [LOCUS ...]` - loci to retain, in output order (required)
- `-o, --out PATH` - output `.hfx` path (required)
- `--no-manifest` - skip `MANIFEST.json`
- `--hash {md5,sha256,none}` - checksum algorithm (default: sha256)
- `--tolerance FLOAT` - max allowed deviation between input/output total frequency

Programmatic use:

```python
from hfx_tools.reduce import reduce_loci_hfx

result = reduce_loci_hfx("nine_locus.hfx", ["A", "B", "DRB1", "DQB1"], "four_locus.hfx")
print(result["n_input_haplotypes"], "->", result["n_output_haplotypes"])
print("freq conservation error:", result["frequency_conservation_error"])
```

The projection math is also available as pure, I/O-free functions
(`hfx_tools.reduce.project_haplotype`, `reduce_frequency_rows`) for direct reuse.

### Streamlit: Web UI

Launch the interactive web interface:

```bash
streamlit run "$(python -c 'import hfx_tools.streamlit_app as app; print(app.__file__)')"
```

The Streamlit app provides:
- **Folder browser** - Select local folders containing metadata.json and data files
- **File upload** - Upload metadata.json and data files directly
- **Auto-update mode** - Automatically sets `metadata.frequencyLocation` to point to uploaded data
- **Metadata preview** - View JSON structure and what will be auto-updated before building
- **Validation preview** - Run validators and see results
- **HFX download** - Download the built .hfx file
- **Build logs** - View detailed validation and packing logs
- **HFX inspector** - Browse an existing `.hfx` archive and view its bundled metadata

For folder-based builds, the **Output folder** defaults to `output`, keeping generated archives
and build logs separate from example or source input folders.

### GitHub sign-in (optional)

The Streamlit app always requires a contributor name, affiliation, and GitHub identity before
building. GitHub usernames can be entered manually. To also offer **Sign in with GitHub**:

1. Create a GitHub OAuth App.
2. Set its authorization callback URL to your public Streamlit app URL.
3. Store these values in your deployment secrets, or in an untracked
   `.streamlit/secrets.toml` for local use:

   ```toml
   [github_oauth]
   client_id = "..."
   client_secret = "..."
   redirect_uri = "https://your-app-url/"
   ```

OAuth is enabled only when all three values are present. It requests the `read:user` scope and
uses the authenticated GitHub login as the identity. Manual GitHub identity remains available
when OAuth is not configured. Never commit this secrets file.

## Validation Framework

The build process includes an extensible validation framework with built-in validators:

1. **Schema version** - Ensures top-level `version` matches the current HFX schema (`0.1.1`)
2. **Metadata required fields** - Checks `outputResolution`, `hfeMethod`, `cohortDescription`, `nomenclatureUsed`, `frequencyLocation`
3. **Frequency location** - Validates frequency location format (inline, file://, http://)
4. **Frequency data format** - Checks inline frequency data structure, types, and duplicates
5. **File references** - Verifies that referenced data files exist

Validation results are logged and returned with error/warning levels. The build fails if any
error-level validations fail.

### Custom validators

```python
from hfx_tools.validators import ValidationFramework, ValidationResult


def my_custom_validator(metadata_json, hfx_obj, data_folder):
    return ValidationResult(
        validator_name="my_validator",
        passed=True,
        message="My validation passed",
        level="info",  # or "warning", "error"
    )


validator_framework = ValidationFramework()
validator_framework.register_validator("my_validator", my_custom_validator)
```

## Common Use Cases

### Scenario 1: Batch submission from local folders

```bash
for dir in submissions/*/; do
  hfx-build "$dir" -n "$(basename $dir)" -o dist/
done
```

### Scenario 2: Remote frequency data

Point to frequencies hosted on Zenodo or S3 without bundling:

```json
{
  "frequencyLocation": "https://zenodo.org/record/12345/files/data.csv"
}
```

```bash
hfx-build my_submission -n my_file --no-auto-update-location
```

### Scenario 3: Inline small frequencies

```json
{
  "frequencyLocation": "inline",
  "frequencyData": [
    {"haplotype": "A*01:01", "frequency": 0.123},
    {"haplotype": "A*01:02", "frequency": 0.456}
  ]
}
```

### Scenario 4: Non-standard CSV headers

If your CSV uses column names other than `haplotype`/`frequency`, map them in metadata:

```json
{
  "frequencyFileHeader": {
    "Haplo": "haplotype",
    "Freq": "frequency"
  }
}
```

### Scenario 5: Programmatic use in Python

```python
from hfx_tools.build import build_hfx_from_folder

result = build_hfx_from_folder(
    input_folder="my_data/",
    output_name="my_submission",
    output_dir="dist/",
    hash_alg="sha256",
    write_manifest=True,
)
print(f"Build {'succeeded' if result['success'] else 'failed'}")
for v in result["validation_results"]:
    print(f"  {v.level}: {v.message}")
```

## Developer API

### Using the Validation Framework

```python
from hfx_tools.validators import ValidationFramework, ValidationResult

validator = ValidationFramework()


def check_cohort_size(metadata_json, hfx_obj, data_folder):
    size = hfx_obj.get("metadata", {}).get("cohortDescription", {}).get("cohortSize", 0)
    if size < 100:
        return ValidationResult(
            validator_name="cohort_size",
            passed=False,
            message=f"Cohort too small: {size} < 100",
            level="warning",
        )
    return ValidationResult(
        validator_name="cohort_size", passed=True, message=f"Cohort size OK: {size}", level="info"
    )


validator.register_validator("cohort_size", check_cohort_size)
results = validator.validate(metadata_path, hfx_obj, data_folder)
```

## Package contents

```
hfx_tools/
├── __init__.py
├── build.py           # Build orchestration
├── cli.py             # Command-line interface
├── inspect.py         # HFX inspection tools
├── io.py              # JSON and file I/O
├── pack.py            # Low-level packing
├── qc.py              # Quality control
├── reduce.py          # Locus reduction (marginalization / projection)
├── streamlit_app.py   # Web UI
├── util.py            # Utilities
└── validators.py      # Validation framework
```

## Development & Contributing

### Development setup

```bash
git clone https://github.com/societyforimmunepolymorphism/hfx-tools
cd hfx-tools
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --editable .
python -m pip install pytest ruff build
```

### Running tests and linting

```bash
make fmt       # Format code
make lint      # Check code style
make test      # Run test suite
make build     # Build distribution
```

### Submitting changes

1. Fork the repo
2. Create a feature branch (`git checkout -b feature/my-feature`)
3. Add tests for new functionality
4. Run `make lint test` to verify
5. Submit a pull request

## Troubleshooting

### Issue: "Missing required field: metadata.frequencyLocation"

Add to your metadata.json:
```json
{
  "frequencyLocation": "file://frequencies.csv"
}
```

Or for inline data:
```json
{
  "frequencyLocation": "inline",
  "frequencyData": [...]
}
```

### Issue: Validation errors but can't see why

```bash
hfx-build my_data -n output
cat my_data/output.build.log
```

### Issue: File not found in bundle

Ensure the filename in `frequencyLocation` matches the actual file in your folder:
```
my_data/
├── metadata.json    # frequencyLocation: "file://my_file.csv"
└── my_file.csv      # ← must match
```

### Issue: CSV columns not recognized

Add `frequencyFileHeader` to your metadata to map your column names:
```json
{
  "frequencyFileHeader": {
    "Haplo": "haplotype",
    "Freq": "frequency"
  }
}
```

### Issue: Permission denied when creating .venv

```bash
mkdir -p ~/.hfx-tools
make sync VENV=~/.hfx-tools/.venv
```

## Resources

- [HFX Specification](https://github.com/societyforimmunepolymorphism/hfx) - Authoritative format specification and schema
- [phycus](https://github.com/nmdp-bioinformatics/phycus) - Related NMDP bioinformatics tools
- [Issues & Discussions](https://github.com/societyforimmunepolymorphism/hfx-tools/issues) - Report bugs or suggest features
- [HFX Spec Issues](https://github.com/societyforimmunepolymorphism/hfx/issues) - Discuss spec-related questions
