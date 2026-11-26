# Grid setup notebooks

Hands-on notebooks for the rewritten grid setup (`posydon-grid`,
`posydon/grids/setup/`). Each takes one part of the setup apart, so you can
see how an ini file becomes a run directory, try changes, and reuse the
snippets in your own work.

| Notebook | What it covers | Main modules |
|---|---|---|
| [00_tour](00_tour.ipynb) | the whole thing from the command line: `types`, `init`, `explain`, `setup` | `posydon/CLI/grid/main.py` |
| [01_inlist_files](01_inlist_files.ipynb) | reading and writing MESA inlists | `inlist.py` |
| [02_layers](02_layers.ipynb) | stacking `mesa_inlists` < `posydon_inlists` < `user_inlists` < `generated` | `layers.py` |
| [03_ini_and_grid_types](03_ini_and_grid_types.ipynb) | the new ini format, the grid-type registry, errors, new vs legacy | `config.py`, `grid_types.py`, `template.py` |
| [04_releases_and_reruns](04_releases_and_reruns.ipynb) | the inlist repository: manifests, reruns, what each rerun changes | `manifest.py`, `resolve.py` |
| [05_building_inlists](05_building_inlists.ipynb) | MESA defaults, layers per target, ZAMS models, formation steps, generated keys | `defaults.py`, `build.py`, `explain.py` |
| [06_extras_and_columns](06_extras_and_columns.ipynb) | Fortran sources, makefiles, `mk`, column lists | `extras.py` |
| [07_submission_and_record](07_submission_and_record.ipynb) | the `posydon-run-grid` command, shell/SLURM scripts, `setup_record.ini` | `submit.py`, `record.py` |
| [08_step_by_step_and_legacy](08_step_by_step_and_legacy.ipynb) | `prepare` / `explain_grid` / `setup_grid`, comparing two run directories, legacy inis | `steps.py`, `legacy.py` |

Start with 00; the others stand on their own and can be read in any order.

## Running them

You need:

- this branch of POSYDON, and its Python environment (e.g.
  `posydon-debug`). The notebooks import `posydon` from this checkout, even
  if another POSYDON is installed in the environment;
- a checkout of
  [POSYDON-MESA-INLISTS](https://github.com/POSYDON-code/POSYDON-MESA-INLISTS)
  on the `v2_sources` branch. Point `POSYDON_INLIST_REPO` at it (default:
  `~/.posydon_mesa_inlists`). The notebooks only read it: they use
  `source = path`, so nothing is checked out.

```bash
export POSYDON_INLIST_REPO=~/POSYDON-MESA-INLISTS
cd grid_rewrite_tests/notebooks
jupyter lab
```

You don't need MESA. `nb_env.py` points `$MESA_DIR` at the test suite's stub
MESA (no-op makefiles, labelled stub sources), so a full setup takes seconds
and writes no executables. To use your real MESA instead, set
`POSYDON_NOTEBOOK_REAL_MESA=1`; the setup then really compiles, which takes
minutes.

Everything the notebooks write goes to `_work/` (git-ignored). Delete it
whenever you like.

## Adapting them

`nb_env.py` holds the small helpers the notebooks share: `env()` for the
environment, `fresh(name)` for an empty work folder, `write`, `show`, `tree`,
`grid_csv` and `cli` (runs this branch's `posydon-grid`). Copy a notebook,
change the ini text in it, and re-run.

The saved outputs were produced on one machine, so the paths in them are
that machine's. Re-run a notebook to see your own.
