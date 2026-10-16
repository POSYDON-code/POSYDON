# Grid setup rewrite: baseline tests

Temporary tests for the rewrite of `posydon-setup-grid`. They freeze what
today's script writes for every grid type, so the new implementation can be
checked against it. Delete this folder once the rewrite has landed.

They are kept out of `posydon/unit_tests/`, so the regular test run and its
coverage check do not pick them up. CI runs them through
`.github/workflows/grid_rewrite_tests.yml` (pushes to the rewrite branches and
pull requests touching the grid scripts or this folder); delete that workflow
together with this folder.

## Run

From the repository root:

    pytest grid_rewrite_tests

This takes about 1.5 minutes and needs `git` and `make`. No MESA, network or
cluster is used.

- `POSYDON_UPDATE_GOLDEN=1` rewrites `setup_grid_data/golden/`.
- `POSYDON_GOLDEN_MODE=semantic` ignores the sha256 lines, comparing only
  non-default namelist values, scripts, staged sources and folder trees.
- `POSYDON_YGGDRASIL_REFERENCE=<dir>` also replays real yggdrasil runs copied
  to `<dir>/grids/<type>/<Z>/<run>/` against today's setup (about 2.5
  minutes; needs a local `~/.posydon_mesa_inlists` clone).

## Contents

- `test_setup_grid_baseline.py`: harness, golden cases, folder-structure
  checks, scenario-vs-explicit-paths checks, the opt-in replay and strict
  xfail tests for known bugs.
- `setup_grid_data/inlists_main/`, `inlists_jeff_new_grid/`: vendored
  POSYDON-MESA-INLISTS snapshot (`a060b7b`, plus `HeMS-HMS/` from
  `jeff_new_grid`), with Fortran sources replaced by labelled stubs.
- `setup_grid_data/mesa_stub/`: stub `MESA_DIR` with no-op makefiles.
- `setup_grid_data/extra/`: extra fixture inlists and a psy-cris ini.
- `setup_grid_data/golden/`: one file per case.
