"""Build a grid's inlists: stack the layers, add generated keys, write.

For every target the files of the four layers are stacked with
``layers.LayerStack`` (``mesa_inlists`` < ``posydon_inlists`` <
``user_inlists`` < ``generated``) and written into the run directory::

    binary/inlist_project   binary_controls, binary_job
    binary/inlist1          star1_controls, star1_job  (star 1 in the binary)
    binary/inlist2          star2_controls, star2_job
    star1/inlist_step<N>    star1_formation_step<N>_controls, _job
    star2/inlist_step<N>    star2_formation_step<N>_controls, _job

The *generated* keys are the ones setup writes itself, in the same order as
today's ``posydon-setup-grid``:

1. wiring: star inlists read the per-point grid inlists that
   ``posydon-run-grid`` writes (``read_extra_*_inlist1``) when a grid
   column is one of their keys; the binary points at its two star inlists
   (``inlist_names``);
2. formation: each step loads the previous step's model and saves its own;
   the binary star loads the last one;
3. output: history, final profile and final model, ``history_interval``;
4. ZAMS models: ``zams_filename``.

ZAMS models follow today's rules (``legacy_zams``), including two known
bugs (tracker S2, S3) that step 6 of the rewrite fixes. The files are
written in today's exact layout (``render_*``), so the output can be
compared byte for byte.
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

import os
from dataclasses import dataclass
from pathlib import Path

from posydon.grids.setup.grid_types import get_grid_type
from posydon.grids.setup.layers import (
    LAYER_ORDER,
    LayerStack,
    drop_wiring_keys,
    expand_every_step,
    fix_num_x_ctrls,
    formation_step_count,
)

STARS = ("star1", "star2")


###############################################################################
# which files go where
###############################################################################
def is_single_star(config):
    """Whether the grid evolves star 1 alone."""
    if config.grid.type is not None:
        return get_grid_type(config.grid.type).single_star
    mesa_inlists = config.legacy.get("mesa_inlists", {})
    return bool(mesa_inlists.get("single_star_grid", False))


def layer_targets(config, scenario=None):
    """Files per layer and target, with the scenario's files in place.

    The scenario's files go into ``posydon_inlists`` after any given in the
    ini, and its values (the single-star flag) after the files.

    A single-star grid runs star 1 as its formation steps, so files given
    for star 1 itself (``star1``, ``star1_controls``, ``star1_job``) apply
    to every step, after that layer's every-step files, as a legacy ini's
    ``star1_*_user`` keys do (before, a new-format ini's were dropped).

    Returns
    -------
    files : dict
        Layer to {target: [path, ...]}.
    flags : dict
        Target to values added to ``posydon_inlists`` after its files.
    """
    files = {layer: {target: list(paths) for target, paths in
                     config.inlists[layer].targets.items()}
             for layer in config.inlists}
    flags = {}
    if scenario is not None:
        posydon = files["posydon_inlists"]
        for target, paths in scenario.inlists.targets.items():
            posydon.setdefault(target, []).extend(paths)
        flags = {target: dict(values)
                 for target, values in scenario.flags.items()}
    if is_single_star(config):
        for targets in files.values():
            for part in ("controls", "job"):
                paths = targets.pop("star1_{}".format(part), [])
                if paths:
                    targets.setdefault("star1_formation_every_step_{}".format(
                        part), []).extend(paths)
    for star in STARS:
        n_steps = formation_step_count(star, list(files.values()) +
                                       [flags])
        for layer in files:
            files[layer] = expand_every_step(files[layer], star, n_steps)
        flags = expand_every_step(
            {target: [values] for target, values in flags.items()},
            star, n_steps)
        flags = {target: _merge_dicts(values)
                 for target, values in flags.items()}
    return files, flags


def _merge_dicts(dicts):
    merged = {}
    for values in dicts:
        merged.update(values)
    return merged


def _cleanup_for(target):
    """Today's clean-up of input keys: wiring dropped, and for star
    controls the num_x_ctrls placeholder replaced."""
    if target.startswith("star") and target.endswith("_controls"):
        return lambda values: fix_num_x_ctrls(drop_wiring_keys(values))
    return drop_wiring_keys


def stack_target(target, files, flags):
    """A LayerStack for one target, from ``layer_targets`` output."""
    stack = LayerStack(target)
    cleanup = _cleanup_for(target)
    for layer in LAYER_ORDER:
        for path in files.get(layer, {}).get(target, []):
            stack.add_file(layer, path, transform=cleanup)
        if layer == "posydon_inlists" and target in flags:
            stack.add(layer, flags[target], source="grid type flag")
    return stack


###############################################################################
# ZAMS models: today's rules
###############################################################################
def legacy_zams(config, scenario=None):
    """The ZAMS model of each star.

    In order: the ini's ``zams_filename_1/2`` (legacy); the ZAMS models the
    scenario names; the legacy ``zams_filename``, for both stars. A
    new-format ini's ``[grid] zams``/``zams_star1/2`` wins over the
    scenario, as the user's own setting (before, the release's model
    silently replaced it). Last, a star the scenario
    builds in formation steps has no ZAMS model, whatever the ini says
    (fixed S2: before, ``zams_filename`` came last, so a HeMS-HMS grid never
    formed star 1). A single-star grid whose star 1 is formed has no ZAMS
    model for star 2 either (fixed S3), as in the v2 single_HeMS grids.

    Returns
    -------
    dict
        ``{"star1": str or None, "star2": str or None}``, as written in the
        ini (paths are not normalised).
    """
    raw = config.legacy.get("mesa_inlists", {})
    zams = {"star1": raw.get("zams_filename_1"),
            "star2": raw.get("zams_filename_2")}
    if not raw:
        # a new-format config: what [grid] says
        zams = {star: None if path is None else str(path) for star, path in
                (("star1", config.grid.zams_star1),
                 ("star2", config.grid.zams_star2))}
    if scenario is not None:
        for star in STARS:
            found = getattr(scenario, "zams_" + star)
            if found is not None and (raw or zams[star] is None):
                zams[star] = str(found)
        if is_single_star(config) and raw:
            # legacy only; a new-format ini's [grid] zams_star2 is kept
            zams["star2"] = None
    if raw.get("zams_filename") is not None:
        zams["star1"] = zams["star2"] = raw["zams_filename"]
    if scenario is not None:
        for star in scenario.formation_stars:
            zams[star] = None  # S2
        if is_single_star(config) and "star1" in scenario.formation_stars:
            zams["star2"] = None  # S3
    return zams


###############################################################################
# writing, in today's exact layout
###############################################################################
def _lines(values):
    return "".join("\t{} = {}\n".format(key, value)
                   for key, value in values.items())


def render_binary_project(controls, job):
    """``binary/inlist_project`` text."""
    return ("&binary_controls\n\n" + _lines(controls) +
            "\n/ ! end of binary_controls namelist" + "\n\n" +
            "&binary_job\n\n" + _lines(job) +
            "\n/ ! end of binary_job namelist\n                ")


def render_binary_star(star, controls, job):
    """``binary/inlist1`` or ``binary/inlist2`` text."""
    return ("&controls\n\n" + _lines(controls) + "\n\n" +
            "\n/ ! end of {}_controls namelist\n\n".format(star) +
            "&star_job\n\n" + _lines(job) +
            "\n/ ! end of {}_job namelist\n\n".format(star))


def render_formation_step(star, controls, job):
    """``star<N>/inlist_step<M>`` text."""
    return ("&controls\n\n" + _lines(controls) + "\n\n" +
            "\n        / ! end of {}_controls namelist\n\n        ".format(
                star) +
            "&star_job\n\n" + _lines(job) +
            "\n        / ! end of {}_job namelist\n\n        ".format(star))


###############################################################################
# the whole build
###############################################################################
@dataclass
class BuiltInlists:
    """Paths of the written inlists, as posydon-run-grid needs them.

    ``star1_formation``/``star2_formation`` are today's command-line value:
    the step files, each preceded by a space, or ``None``.
    """

    binary_project: str
    binary_star1: str
    binary_star2: str
    star1_formation: str = None
    star2_formation: str = None


def build_inlists(config, run_directory, grid_parameters, scenario=None):
    """Stack, generate and write all inlists of a grid.

    Parameters
    ----------
    config : GridConfig
    run_directory : str or Path
        Where ``binary/``, ``star1/`` and ``star2/`` are written. Paths in
        the files are written exactly as given here.
    grid_parameters : iterable of str
        The grid's column names (or psy-cris's ``mesa_column_names``).
    scenario : ScenarioFiles or None

    Returns
    -------
    BuiltInlists
    """
    run = str(run_directory)
    grid_parameters = list(grid_parameters)
    single = is_single_star(config)
    raw = config.legacy.get("mesa_inlists", {})
    output = config.output
    files, flags = layer_targets(config, scenario)
    zams = legacy_zams(config, scenario)

    def merged(target):
        return stack_target(target, files, flags).merged()

    binary_project = os.path.join(run, "binary", "inlist_project")
    binary_star = {star: os.path.join(run, "binary", "inlist{}".format(
        star[-1])) for star in STARS}

    # 1. the binary and its two stars
    final = {}
    if not single:
        for target in ("binary_controls", "binary_job", "star1_controls",
                       "star1_job", "star2_controls", "star2_job"):
            final[target] = merged(target)
        for star in STARS:
            for part, extra in (("controls", "controls"),
                                ("job", "star_job")):
                target = "{}_{}".format(star, part)
                if any(param in final[target] for param in grid_parameters):
                    final[target]["read_extra_{}_inlist1".format(extra)] = \
                        ".true."
                    final[target]["extra_{}_inlist1_name".format(extra)] = \
                        "'inlist_grid_{}_binary_{}'".format(star, part)
        final["binary_job"]["inlist_names(1)"] = "'{}'".format(
            binary_star["star1"])
        final["binary_job"]["inlist_names(2)"] = "'{}'".format(
            binary_star["star2"])
    else:
        for target in ("binary_controls", "binary_job", "star1_controls",
                       "star1_job", "star2_controls", "star2_job"):
            final[target] = {}

    # 2. formation steps
    formation = {}
    for star in STARS:
        n_steps = formation_step_count(star, list(files.values()) + [flags])
        # today: a star with a ZAMS model is not formed, except star 1 of a
        # single-star grid, whose "formation" step is the whole run
        if star == "star1":
            builds = n_steps > 0 and (single or zams["star1"] is None)
        else:
            builds = n_steps > 0 and zams["star2"] is None
        if not builds:
            formation[star] = None
            continue
        paths = []
        for step in range(n_steps):
            controls = merged("{}_formation_step{}_controls".format(star,
                                                                    step))
            job = merged("{}_formation_step{}_job".format(star, step))
            model = "'initial_{}_step{}.mod'".format(star, step)
            binary_job = final["{}_job".format(star)]
            binary_job["create_pre_main_sequence_model"] = ".false."
            binary_job["load_saved_model"] = ".true."
            binary_job["saved_model_name"] = model
            if step > 0:
                job["create_pre_main_sequence_model"] = ".false."
                job["load_saved_model"] = ".true."
                job["saved_model_name"] = "'initial_{}_step{}.mod'".format(
                    star, step - 1)
            job["save_model_when_terminate"] = ".true."
            job["save_model_filename"] = model
            if star == "star1" and single:
                if zams["star1"] is not None:
                    controls["zams_filename"] = "'{}'".format(zams["star1"])
                else:
                    # formed without a ZAMS model: no zams_filename at all,
                    # not even MESA's default (S3)
                    controls.pop("zams_filename", None)
            path = os.path.join(run, star, "inlist_step{}".format(step))
            Path(path).parent.mkdir(parents=True, exist_ok=True)
            Path(path).write_text(render_formation_step(star, controls, job))
            paths.append(path)
        formation[star] = "".join(" {}".format(path) for path in paths)

    # 3. output options
    for star in STARS:
        job = final["{}_job".format(star)]
        controls = final["{}_controls".format(star)]
        if "final_profile" in getattr(output, star):
            job["write_profile_when_terminate"] = ".true."
            job["filename_for_profile_when_terminate"] = \
                "'final_profile_{}.data'".format(star)
        else:
            job["write_profile_when_terminate"] = ".false."
    for star in STARS:
        job = final["{}_job".format(star)]
        if "final_model" in getattr(output, star):
            job["save_model_when_terminate"] = ".true."
            job["save_model_filename"] = "'final_{}.mod'".format(star)
        else:
            job["save_model_when_terminate"] = ".false."
    for star in STARS:
        final["{}_controls".format(star)]["do_history_file"] = \
            ".true." if "history" in getattr(output, star) else ".false."
    history_interval = raw.get("history_interval", output.history_interval)
    for target in ("binary_controls", "star1_controls", "star2_controls"):
        final[target]["history_interval"] = history_interval
    if not output.binary_history:
        final["binary_controls"]["history_interval"] = "-1"

    # 4. ZAMS models
    for star in STARS:
        if zams[star] is not None:
            final["{}_controls".format(star)]["zams_filename"] = \
                "'{}'".format(zams[star])

    # write the binary inlists
    Path(binary_project).parent.mkdir(parents=True, exist_ok=True)
    Path(binary_project).write_text(render_binary_project(
        final["binary_controls"], final["binary_job"]))
    for star in STARS:
        Path(binary_star[star]).write_text(render_binary_star(
            star, final["{}_controls".format(star)],
            final["{}_job".format(star)]))
    return BuiltInlists(binary_project=binary_project,
                        binary_star1=binary_star["star1"],
                        binary_star2=binary_star["star2"],
                        star1_formation=formation["star1"],
                        star2_formation=formation["star2"])
