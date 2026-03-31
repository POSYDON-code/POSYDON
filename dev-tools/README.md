Validation suite for POSYDON binary evolution. Evolves a fixed set of test binaries on a candidate branch and compares results against a stored baseline to catch regressions. A baseline can be formed from any branch (`main` by default) and is represented by a set of results from `binary_suite.py`, saved HDF5 files, all stored in `dev-tools/baselines/<branch-name>`.

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

Top-level entry point. Evolves test binaries on a candidate branch, then compares results against an existing baseline. This script will look for baseline HDF5 files stored in `dev-tools/baseline/<branch-name>`, where `main` is the default `<branch-name>`.

```bash
./validate_binaries.sh <candidate_branch> [baseline_branch] [metallicities] [--loose] [--rtol VALUE] [--atol VALUE] [--skip-evolve]
```

By default, comparison is exact (rtol=0, atol=0). Use `--loose` for relaxed floating-point tolerances (rtol=1e-12, atol=1e-15), or set `--rtol`/`--atol` explicitly as per [np.allclose](https://numpy.org/devdocs/reference/generated/numpy.allclose.html). Use `--skip-evolve` to skip the evolution step and compare existing candidate outputs against the baseline.

### `generate_baseline.sh`

Generates baseline HDF5 files from a designated branch name and optionally a SHA to specify a commit. 

```bash
./generate_baseline.sh <branch> [sha] [metallicities]
```

If you already have results from prior runs of `evolve_binaries.sh` saved as HDF5 files in `outputs/<branch>/`, you can copy these directly into the baselines directory with the `--promote` option, skipping re-evolution:

```bash
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

The script also compares warning and error tables, reporting new, removed, or changed warnings per binary. The report header includes provenance metadata (branch, commit SHA, generation time, POSYDON data path) read from each file when available.

```bash
python compare_runs.py baseline.h5 candidate.h5 [--loose] [--rtol VALUE] [--atol VALUE] [--verbose]
```

### `binaries_params.ini`

Configuration file for `SimulationProperties`. Defines the POSYDON evolution steps, supernova prescriptions, common envelope parameters, and output column selections. Metallicity is overridden at runtime by `binaries_suite.py`.

## Running Scripts Manually

The shell scripts handle cloning, environment setup, orchestration, and execution. If you already have POSYDON installed in your current environment, you can execute the Python scripts directly.

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

## Interpreting Results

The comparison report groups differences into four categories. Here's how to read them:

**Structural** differences (missing/extra binaries, step count changes, newly failing/passing) almost always indicate a real change. A binary that newly fails or changes its number of evolution steps means the code is following a different evolutionary path. These warrant investigation regardless of tolerance settings.

**Qualitative** differences (state, event, step name, SN type) also represent real behavioral changes. Even a single qualitative diff means the binary is being classified differently, e.g. a different mass transfer history or a changed SN type. These are never tolerance-dependent.

**Quantitative** differences are more nuanced. With exact comparison (the default), any floating-point difference is reported. This is useful for detecting unintended changes, but expected after compiler/platform changes or numpy version bumps. If you see many quantitative diffs but zero structural/qualitative diffs, the evolution paths are the same and the differences are likely numerical noise — re-run with `--loose` or a custom `--rtol` to confirm. If quantitative diffs persist at `--rtol 1e-6` or larger, something meaningful has changed.

**Warning** differences are informational. New warnings may indicate a physics edge case being hit differently, or changes to warnings in the code, but are not failures on their own.

A healthy validation run after a non-physics code change should show zero structural and qualitative differences, and minimal quantitative changes. After an intentional physics change, expect significant diffs. If binaries unrelated to intentional physics changes show structural or qualitative diffs, that may indicate a problem with implementation.

## Tolerance Design

By default, comparison is exact (`rtol=0, atol=0`): any bitwise difference in a float is reported. The `--loose` flag sets `rtol=1e-12, atol=1e-15`, which is appropriate for filtering out platform-level floating-point noise while still catching meaningful changes.

For custom tolerances, `--rtol` and `--atol` follow the semantics of `np.allclose`: a value passes if `abs(baseline - candidate) <= atol + rtol * abs(baseline)`. In practice, `rtol` dominates for most columns (masses, periods, separations are all large numbers), while `atol` only matters near zero (e.g., eccentricity, certain hydrogen fractions).

Known limitation: when `baseline == 0` and `candidate != 0`, `rtol`-based comparison produces `0 + rtol * 0 = 0`, so any nonzero candidate value fails. This is correct behavior (a zero-to-nonzero change is meaningful), but be aware that the reverse (both values very small but nonzero) may pass even if the relative change is large, since `atol` provides a floor. For most POSYDON quantities this is not an issue, but it matters for quantities that are genuinely expected to be zero (e.g., eccentricity at ZAMS for circular binaries).

A single global tolerance works well for catching regressions but is a blunt instrument for columns spanning many orders of magnitude. Per-column or per-quantity scaling is a possible future improvement but is not currently implemented.

The `--loose` defaults (`rtol=1e-12, atol=1e-15`) were chosen just above float64 machine epsilon and may need to be adjusted if there are parts of the code that are non-deterministic. If parts of the POSYDON pipeline in the branches being tested introduce stochasticity (e.g. unseeded RNG), the irreducible noise floor may be higher. To calibrate, run the same branch against itself and check what tolerance is needed for a clean pass:

```bash
# Evolve the same branch twice under different output names
python binaries_suite.py --output /tmp/run_a.h5 --metallicity 1
python binaries_suite.py --output /tmp/run_b.h5 --metallicity 1

# Compare — any diffs here are the stochasticity floor
python compare_runs.py /tmp/run_a.h5 /tmp/run_b.h5

# Find the tolerance that absorbs the noise
python compare_runs.py /tmp/run_a.h5 /tmp/run_b.h5 --rtol 1e-10
```

The `--loose` defaults should sit just above whatever self-comparison noise you observe. If the self-comparison is clean at exact, the current defaults are fine.

**RNG reproducibility.** Several POSYDON evolution steps (Bondi-Hoyle accretion in `step_detached` and `MesaGridStep`, SN kicks in `step_SN`) use random number generation internally. Without a fixed seed, these produce nondeterministic results that appear as spurious `S1_lg_mdot` diffs in the validation suite. To ensure reproducibility, set the `entropy` parameter to a fixed integer in `binaries_params.ini`. This seeds the RNG passed to each step (see PR#826).

## Updating the Baseline

The baseline should be regenerated when the "expected correct" output changes. Typical triggers:

- **After a release or version tag.** Generate a baseline from the release tag so future development is compared against the release state: `./generate_baseline.sh v2.3`
- **After merging an intentional physics change.** If a PR deliberately changes evolution outcomes (e.g., a new SN prescription), validate the PR branch first to confirm only the expected binaries are affected, then regenerate the baseline from the updated main branch.
- **After updating POSYDON data grids.** Grid changes will alter interpolated values. Regenerate the baseline and record the new `PATH_TO_POSYDON_DATA` in `baseline_info.txt`.

Do not regenerate the baseline to silence unexpected diffs. If a validation run shows differences you don't understand, investigate them before updating the baseline.

The `--promote` flag on `generate_baseline.sh` is a convenience for skipping re-evolution when you've already run the suite and are satisfied with the outputs: `./generate_baseline.sh --promote main "1 0.45"`.

## Adding New Test Binaries

New binaries are added by appending entries to the `get_test_binaries()` function in `binaries_suite.py`. Each entry is a tuple of `(star1_kwargs, star2_kwargs, binary_kwargs, description)`.

When adding a binary:

- Choose initial conditions that reliably trigger the edge case or evolutionary pathway you want to test. Verify it does so at multiple metallicities if possible, since grid coverage varies.
- Use a descriptive string that references the PR or issue number if the binary guards a specific fix (e.g., `"PR574 - stepCE fix"`).
- After adding the binary, regenerate the baseline so it includes the new binary's expected output.
- The binary ID is assigned by list position. Appending to the end avoids changing IDs of existing binaries, which would invalidate old baselines against new code for no reason.
