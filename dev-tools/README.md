# dev-tools

Validation suite for POSYDON binary evolution. Evolves a fixed set of test binaries on a candidate branch and compares results against a stored baseline to catch regressions.

## Quick Start

```bash
# 1. Generate a baseline from the main branch (once)
./generate_baseline.sh main

# 2. Validate a candidate branch against that baseline
./validate_binaries.sh feature/my-branch
```

## Scripts

### `validate_binaries.sh`

Top-level entry point. Evolves test binaries on a candidate branch, then compares results against a baseline.

```bash
./validate_binaries.sh <candidate_branch> [baseline_branch] [metallicities] [--loose] [--rtol VALUE] [--atol VALUE]
```

By default, comparison is exact (rtol=0, atol=0). Use `--loose` for relaxed floating-point tolerances (rtol=1e-12, atol=1e-15), or set `--rtol`/`--atol` explicitly.

### `generate_baseline.sh`

Generates baseline HDF5 files from a designated branch or tag. Can also promote existing outputs to baseline with `--promote`.

```bash
./generate_baseline.sh <branch> [sha] [metallicities]
./generate_baseline.sh --promote <branch> [metallicities]
```

### `evolve_binaries.sh`

Clones a POSYDON branch, creates a conda environment, installs POSYDON, and runs the binary suite at all requested metallicities. Called by `validate_binaries.sh` and `generate_baseline.sh`; can also be run standalone.

```bash
./evolve_binaries.sh <branch> [sha] [metallicities]
```

### `binaries_suite.py`

Defines and evolves the set of 44 test binaries at a given metallicity. Each binary targets a specific edge case or past bug fix (e.g., matching failures, oRLO2 looping, SN type errors, NaN spins). Results are saved to an HDF5 file.

```bash
python binaries_suite.py --output results.h5 --metallicity 1
```

### `compare_runs.py`

Compares two HDF5 files produced by `binaries_suite.py`. Reports three categories of differences: structural (missing binaries, step count changes, new errors), qualitative (state/event/step name changes), and quantitative (numeric value changes).

```bash
python compare_runs.py baseline.h5 candidate.h5 [--loose] [--rtol VALUE] [--atol VALUE] [--verbose]
```

### `binaries_params.ini`

Configuration file for `SimulationProperties`. Defines the POSYDON evolution steps, supernova prescriptions, common envelope parameters, and output column selections. Metallicity is overridden at runtime by `binaries_suite.py`.

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

## Available Metallicities

The suite supports metallicities (in solar units): 2, 1, 0.45, 0.2, 0.1, 0.01, 0.001, 0.0001. All are run by default; pass a quoted subset to limit, e.g. `"1 0.45"`.

## Authors

Max Briel, Elizabeth Teng
