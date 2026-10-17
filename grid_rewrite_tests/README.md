# Grid setup rewrite: baseline tests

Temporary tests for the rewrite of `posydon-setup-grid`. They freeze what
today's script writes for every grid type, so the new implementation can be
checked against it. Delete this folder once the rewrite has landed.

They are kept out of `posydon/unit_tests/`, so the regular test run and its
coverage check do not pick them up. CI runs them in the `baseline` job of
`.github/workflows/grid_rewrite_tests.yml` on every push to the rewrite
branch; delete that job together with this folder. The same workflow also
runs `posydon/unit_tests` on the branch.

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
- `test_release_manifests.py`: the restructured inlist repo
  (`r11701/releases/v2.0` + reruns, branch `v2_sources`) against today's
  scenario rules on every archived production commit. Needs
  `POSYDON_INLIST_REPO=<clone at the commit pinned in the file>`; CI clones it.
- `test_new_ini.py`: `posydon-setup-grid` on a legacy ini and on its
  new-format equivalent (`[scenario] release = v2.0`) writes the same run
  directory, for every grid type. Needs `POSYDON_INLIST_REPO`.
- `test_new_vs_legacy.py`: each new module against the matching part of
  the pre-rewrite script.
- `legacy_posydon_setup_grid.py`: `bin/posydon-setup-grid` as it was before
  stage 3.6 replaced it with a thin script over `posydon.grids.setup`; the
  comparisons in `test_new_vs_legacy.py` run it. The goldens, folder checks,
  replay and xfails run the new `bin/posydon-setup-grid`.
- `setup_grid_data/inlists_main/`, `inlists_jeff_new_grid/`: vendored
  POSYDON-MESA-INLISTS snapshot (`a060b7b`, plus `HeMS-HMS/` from
  `jeff_new_grid`), with Fortran sources replaced by labelled stubs.
- `setup_grid_data/mesa_defaults/`: MESA r11701's `.defaults`, byte-identical
  to POSYDON-MESA-INLISTS `r11701/mesa/default/defaults/` (POSYDON no longer
  ships `grid_params/`).
- `setup_grid_data/mesa_stub/`: stub `MESA_DIR` with no-op makefiles.
- `setup_grid_data/extra/`: extra fixture inlists and a psy-cris ini.
- `setup_grid_data/golden/`: one file per case.
