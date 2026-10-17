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
| scenario | none; release + `ref`/`commit`; release + `source = path`; release + `rerun`; legacy `branch-commit` + system type |
| mode × backend | fixed: shell + job_array, slurm_array; dynamic: slurm_mpi, shell without job_array |
| ZAMS | none; `zams` (both); `zams_star1`; `zams_star2` |
| `[mesa_inlists]` | empty (release defaults); some targets; all targets |
| `[posydon_inlists]` | empty; extra files on top of the release |
| `[user_inlists]` | `binary`, `star1`, `star2`, a single namelist (`star1_controls`, `binary_job`, ...), `star1_formation`, `star1_formation_step<N>`, column lists, `run_binary_extras`, `run_star_extras` |
| `[output]` | per star: history, final_profile, final_model or none; binary history or none; history_interval |
| `[extras]` | empty (MESA's stock files); some keys; any key name |
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
   target, which replaces them) < `[posydon_inlists]` < release (< rerun
   overlay) < `[user_inlists]` < generated keys. The release comes after
   the ini's `[posydon_inlists]` in the same layer (tracker C7), so only a
   formation step's own `[posydon_inlists]` file beats the release's
   every-step file. `test_every_layer_filled` fills every target of all
   three ini layers at once and checks this for each written namelist.
2. **Where a user file goes.** `binary` → `binary/inlist_project`;
   `star1`/`star2` → `binary/inlist1`/`inlist2` in a binary grid, every
   step in a single-star grid; `star1_formation_step<N>` → only
   `star1/inlist_step<N>`; `star1_controls` → only `&controls`.
3. **Formation chain.** Step N saves `initial_star1_step<N>.mod`; step N>0
   loads step N−1's model; the binary star loads the last step's model.
4. **ZAMS.** A star that starts from ZAMS gets `zams_filename` from
   `[grid]`. A formed star never does, whatever `[grid]` says. A single_HMS
   star gets it in its step.
5. **Output.** `[output] star<N>` sets `do_history_file`,
   `write_profile_when_terminate` and `save_model_when_terminate` in the
   inlist that runs the star; `binary = none` sets the binary
   `history_interval = -1`. *Open:* single-star grids get these only in
   `binary/inlist1`, which they never run, as the old setup did (pinned as
   a strict xfail).
6. **Wiring.** In a binary grid, a star namelist that holds a grid column
   reads `inlist_grid_star<N>_binary_<part>`. Formation steps are not
   wired: `posydon-run-grid` wraps each step in an inlist that reads the
   step and then `inlist_grid_points`.
7. **Whole files.** Column lists and Fortran sources: the highest layer
   wins (`[mesa_inlists]` and `[posydon_inlists]` files lose to the
   release's, `[user_inlists]` files win), copied under the name MESA
   compiles; `mk` builds all three executables.
8. **Scripts.** One script per backend; job arrays cover every grid
   point; SLURM headers contain only the settings the ini gives.
9. **Refusals.** Fixed grids without one point per task, dynamic grids
   with a job array or without `psycris_inifile`, a scenario without a
   release, a rerun without a release, unknown sections, keys or values.
10. **Record.** `setup_record.ini` names the ini, grid type, mode, release
    and the commit used.
