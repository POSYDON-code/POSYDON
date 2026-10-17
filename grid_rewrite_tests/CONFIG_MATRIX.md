# What a grid setup can be: the configuration matrix

Everything `posydon-grid setup` accepts, and what each choice must do to
the run directory. `posydon/unit_tests/grids/test_setup_matrix.py` checks
every row against a synthetic inlist repository laid out like
POSYDON-MESA-INLISTS v2.0. Its expectations are written from this table and
the ini documentation (`docs/.../mesa_grids/inifile.rst`), not from the
code, so a failure means the code and the spec disagree.

## Axes

| Axis | Values |
|---|---|
| ini format | new (`[grid]`), legacy (`[run_parameters]`, `[mesa_inlists]`, ...) |
| grid type | HMS-HMS, CO-HMS_RLO, CO-HeMS, HeMS-HMS, single_HMS, single_HeMS |
| scenario | none; release (read from `inlist_repo` as it is); release + `rerun`; legacy `branch-commit` + system type |
| mode × backend | fixed: shell + job_array, slurm_array; dynamic: slurm_mpi, shell without job_array |
| ZAMS | none; `zams` (both); `zams_star1`; `zams_star2` |
| `[mesa_inlists]` | empty (release defaults); some targets; all targets |
| `[user_inlists]` | `binary`, `star1`, `star2`, a single namelist (`star1_controls`, `binary_job`, ...), `star1_formation`, `star1_formation_step<N>`, column lists, Fortran sources `run_star_extras`, `run_binary_extras`, `run_star1_extras`, `run_star2_extras` |
| `[output]` | per star: history, final_profile, final_model or none; binary history or none; history_interval |
| `[extras]` | MESA's work files only: `makefile_binary`, `makefile_star`, `binary_run`, `star_run` (each optional) |
| `[submission]` | account, partition, walltime, email (each optional); nodes, mpi_tasks, cpus_per_task, mem_per_cpu, work_dir, newgroup |
| `[grid]` extras | keep_profiles, keep_photos, run_directory |

## What each grid type writes

| grid type | star 1 | star 2 | binary inlists | star1 steps | star2 steps |
|---|---|---|---|---|---|
| HMS-HMS | zams | zams | yes | 0 | 0 |
| CO-HMS_RLO | zams | point mass | yes | 0 | 0 |
| CO-HeMS | formation | point mass | yes | 2 | 0 |
| HeMS-HMS | formation | zams | yes | 2 | 0 |
| single_HMS | zams (one step) | – | written, unused | 1 | 0 |
| single_HeMS | formation | – | written, unused | 3 (2 + the flag step) | 0 |

## Rules the run directory must follow

1. **Layers.** For every key, the value is the one from the highest layer
   that sets it: MESA defaults (or the ini's `[mesa_inlists]` file for that
   target, which replaces them) < release (< rerun overlay) <
   `[user_inlists]` < generated keys. There is no `[posydon_inlists]`:
   POSYDON's inlists come from the scenario. `test_every_layer_filled`
   fills every target of both ini layers at once and checks this for each
   written namelist.
2. **Where a user file goes.** `binary` → `binary/inlist_project`;
   `star1`/`star2` → `binary/inlist1`/`inlist2` in a binary grid, every
   step in a single-star grid; `star1_formation_step<N>` → only
   `star1/inlist_step<N>`; `star1_controls` → only `&controls`.
3. **Formation chain.** Step N saves `initial_star1_step<N>.mod`; step N>0
   loads step N−1's model; the binary star loads the last step's model. A
   single star's last step saves what `[output]` says instead (rule 5).
4. **ZAMS.** A star that starts from ZAMS gets `zams_filename` from
   `[grid]`. A formed star never does, whatever `[grid]` says. A single_HMS
   star gets it in its step.
5. **Output.** `[output] star<N>` sets `do_history_file`,
   `write_profile_when_terminate` and `save_model_when_terminate` in the
   inlist that runs the star; `binary = none` sets the binary
   `history_interval = -1`. A single star is run by its last step, so that
   step takes them (the old setup wrote them only into `binary/inlist1`,
   which a single-star grid never runs); its earlier steps keep saving
   `initial_star1_step<N>.mod`.
6. **Wiring.** In a binary grid, a star namelist that holds a grid column
   reads `inlist_grid_star<N>_binary_<part>`. Formation steps are not
   wired: `posydon-run-grid` wraps each step in an inlist that reads the
   step and then `inlist_grid_points`.
7. **Whole files.** Column lists: the highest layer wins
   (`[mesa_inlists]` < release < `[user_inlists]`). Fortran sources have
   one place, `[user_inlists]`: `run_star_extras` goes to every star's
   `src/`, `run_star1_extras`/`run_star2_extras` only to that star's and
   win over it there, `run_binary_extras` to `binary/src/`; all copied
   under the name MESA compiles. `[extras]` takes only MESA's work files;
   `mk` builds all three executables.
8. **Scripts.** One script per backend; job arrays cover every grid
   point; SLURM headers contain only the settings the ini gives.
9. **Refusals.** Fixed grids without one point per task, dynamic grids
   with a job array or without `psycris_inifile`, a scenario without a
   release, a rerun without a release, Fortran sources in `[mesa_inlists]`
   or `[extras]`, unknown sections, keys or values.
10. **Record.** `setup_record.ini` names the ini, grid type, mode, when,
    where and the command; the POSYDON, Python and MESA used; the release
    and the inlist repository as it was read (commit, dirty, fingerprint);
    the files stacked into every written namelist; and the source of every
    copied column list and Fortran source.
