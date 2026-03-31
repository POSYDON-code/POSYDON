Validation suite for POSYDON binary evolution. Evolves a fixed set of test binaries on a candidate branch and compares results against a stored baseline to catch regressions.

## Quick Start

```bash
# 1. Generate a baseline from the main branch (once)
./generate_baseline.sh main

# 2. Validate a candidate branch against that baseline
./validate_binaries.sh feature/my-branch
```

Results are written to `outputs/<branch>/`. After validation, check:

- `outputs/<branch>/comparison_summary.txt` for a pass/fail overview across all metallicities
- `outputs/<branch>/comparison_<Z>Zsun.txt` for detailed per-metallicity diff reports
- `logs/<branch>/evolve_<Z>Zsun.log` for the full evolution output of each metallicity

By default, all eight POSYDON metallicities are run. To validate only a subset, pass a quoted space-separated list as the third argument:

```bash
./generate_baseline.sh main "" "1 0.45"
./validate_binaries.sh feature/my-branch main "1 0.45"
```

To re-run comparison with different tolerances without re-evolving:

```bash
./validate_binaries.sh feature/my-branch main "1 0.45" --skip-evolve --loose
```

## Scripts

### `validate_binaries.sh`

Top-level entry point. Evolves test binaries on a candidate branch, then compares results against an existing baseline.

```bash
./validate_binaries.sh <candidate_branch> [baseline_branch] [metallicities] [--loose] [--rtol VALUE] [--atol VALUE] [--skip-evolve]
```

By default, comparison is exact (rtol=0, atol=0). Use `--loose` for relaxed floating-point tolerances (rtol=1e-12, atol=1e-15), or set `--rtol`/`--atol` explicitly as per [np.allclose](https://numpy.org/devdocs/reference/generated/numpy.allclose.html). Use `--skip-evolve` to skip the evolution step and compare existing candidate outputs against the baseline.

### `generate_baseline.sh`

Generates baseline HDF5 files from a designated branch or tag. Can also promote existing outputs to baseline with `--promote`.

```bash
./generate_baseline.sh <branch> [sha] [metallicities]
./generate_baseline.sh --promote <branch> [metallicities]
```

### `evolve_binaries.sh`

Clones a POSYDON branch, creates a conda environment, installs POSYDON, and runs the binary suite at all requested metallicities. Called by `validate_binaries.sh` and `generate_baseline.sh`; can also be run standalone. Records the resolved commit SHA and branch name in each HDF5 file's metadata for provenance tracking.

```bash
./evolve_binaries.sh <branch> [sha] [metallicities]
```

### `binaries_suite.py`

Defines and evolves the set of 44 test binaries at a given metallicity. Each binary targets a specific edge case or past bug fix (e.g., matching failures, oRLO2 looping, SN type errors, NaN spins). Results are saved to an HDF5 file with a `/metadata` table that records metallicity, binary counts, `PATH_TO_POSYDON_DATA`, and optionally branch name, commit SHA, and generation timestamp (via `--branch`/`--sha`).

```bash
python binaries_suite.py --output results.h5 --metallicity 1
python binaries_suite.py --output results.h5 --metallicity 1 --branch main --sha abc123f
```

### `compare_runs.py`

Compares two HDF5 files produced by `binaries_suite.py` and reports differences in three categories:

- **Structural**: missing or extra binaries, evolution step count changes, binaries that newly fail or newly pass, missing HDF5 tables.
- **Qualitative**: changes to categorical columns such as state, event, step name, SN type, interpolation class, and mass transfer history.
- **Quantitative**: changes to any numeric column. By default, comparison is exact (bitwise identical floats). Use `--loose` for slightly relaxed tolerances (rtol=1e-12, atol=1e-15), or set `--rtol`/`--atol` explicitly as per [np.allclose](https://numpy.org/devdocs/reference/generated/numpy.allclose.html).

The script also compares warning and error tables, reporting new, removed, or changed warnings per binary. The report header includes source metadata (branch, commit SHA, generation time, POSYDON data path) read from each file when available.

```bash
python compare_runs.py baseline.h5 candidate.h5 [--loose] [--rtol VALUE] [--atol VALUE] [--verbose]
```

### `binaries_params.ini`

Configuration file for `SimulationProperties`. Defines the POSYDON evolution steps, supernova prescriptions, common envelope parameters, and output column selections. Metallicity is overridden at runtime by `binaries_suite.py`.

## Running Scripts Manually

The shell scripts handle cloning, environment setup, and orchestration. If you already have POSYDON installed in your current environment, you can run the Python scripts directly.

### Evolving binaries

```bash
# Evolve all 44 test binaries at solar metallicity
python binaries_suite.py --output my_results.h5 --metallicity 1

# Evolve at a specific metallicity with verbose output
python binaries_suite.py --output my_results.h5 --metallicity 0.01 --verbose

# Use a custom ini file
python binaries_suite.py --output my_results.h5 --metallicity 1 --ini /path/to/custom.ini

# Record branch/SHA provenance in HDF5 metadata (done automatically by evolve_binaries.sh)
python binaries_suite.py --output my_results.h5 --metallicity 1 --branch main --sha abc123f
```

The output HDF5 contains three tables: `evolution` (per-step binary data), `errors` (binaries that failed), and `warnings` (warnings raised during evolution). The `/metadata` table records metallicity, binary counts, `PATH_TO_POSYDON_DATA`, and optionally branch, commit SHA, and generation timestamp.

### Comparing two result files

```bash
# Exact comparison
python compare_runs.py file_a.h5 file_b.h5

# Relaxed tolerances
python compare_runs.py file_a.h5 file_b.h5 --loose

# Custom tolerances with verbose diagnostics
python compare_runs.py file_a.h5 file_b.h5 --rtol 1e-8 --atol 1e-12 --verbose
```

The two files do not need to come from the shell pipeline; any pair of HDF5 files produced by `binaries_suite.py` can be compared.

## Directory Structure

```
dev-tools/
├── README.md
├── validate_binaries.sh      # full validation pipeline
├── generate_baseline.sh      # create or promote baselines
├── evolve_binaries.sh        # clone, install, and run suite
├── binaries_suite.py         # test binary definitions and evolution
├── binaries_params.ini       # SimulationProperties configuration
├── compare_runs.py           # diff two HDF5 result files
├── baselines/                # stored baseline HDF5 files (per branch)
├── outputs/                  # candidate evolution results (per branch)
├── logs/                     # per-metallicity evolution logs (per branch)
└── workdirs/                 # cloned repos and conda environments (per branch)
```