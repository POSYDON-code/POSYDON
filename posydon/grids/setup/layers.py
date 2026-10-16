"""Stack several inlists into one, in a fixed and named order.

A grid's inlists are built in layers. For every *target* (the binary's
``&binary_controls``, star 1's ``&controls``, a formation step, ...) the
values come from up to four sources, applied in this order::

    1. mesa_inlists      MESA's own defaults
    2. posydon_inlists   POSYDON's physics, picked by the scenario
    3. user_inlists      the user's own files
    4. generated         keys setup writes itself (output options, wiring)

A later layer overrides an earlier one key by key: it only changes the keys
it sets. A ``LayerStack`` collects the values of one target and merges them::

    >>> stack = LayerStack("star1_controls")
    >>> stack.add("mesa_inlists", {"mesh_delta_coeff": "1d0",
    ...                            "alpha_semiconvection": "0d0"})
    >>> stack.add("user_inlists", {"mesh_delta_coeff": "0.8d0"})
    >>> stack.merged()
    {'mesh_delta_coeff': '0.8d0', 'alpha_semiconvection': '0d0'}
    >>> [entry["layer"] for entry in stack.explain("mesh_delta_coeff")]
    ['mesa_inlists', 'user_inlists']

The order of ``add`` calls does not matter between layers, only within one.
Keys keep the position where they first appear, as in today's
``posydon-setup-grid``.

Column lists and Fortran sources are not merged: for those, the highest
layer that provides a file wins outright (``pick_whole_file``).
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

from posydon.grids.setup.inlist import LEGACY_SKIPPED_VALUES, read_inlist

LAYER_ORDER = ("mesa_inlists", "posydon_inlists", "user_inlists",
               "generated")


def section_for(target):
    """Name of the inlist section that feeds a target.

    Binary targets are read from the section with the same name; star
    targets (``star1_controls``, ``star2_job``, formation steps such as
    ``star1_formation_step0_controls``) from ``controls`` or ``star_job``.

    Examples
    --------
    >>> section_for("binary_job")
    'binary_job'
    >>> section_for("star1_controls")
    'controls'
    >>> section_for("star1_formation_step1_job")
    'star_job'
    """
    if target in ("binary_controls", "binary_job"):
        return target
    if target.startswith("star") and target.endswith("_controls"):
        return "controls"
    if target.startswith("star") and target.endswith("_job"):
        return "star_job"
    raise ValueError("unknown target {!r}: expected binary_controls, "
                     "binary_job, or a star target ending in _controls "
                     "or _job".format(target))


class LayerStack:
    """The values of one target, layer by layer.

    Parameters
    ----------
    target : str
        What the merged values are for, e.g. ``"star1_controls"``.
    layer_order : sequence of str
        Layer names, lowest first. Defaults to ``LAYER_ORDER``.
    """

    def __init__(self, target, layer_order=LAYER_ORDER):
        self.target = target
        self.section = section_for(target)
        self.layer_order = tuple(layer_order)
        # one list of (source, values) per layer, in the order they were added
        self._layers = {layer: [] for layer in self.layer_order}

    def add(self, layer, values, source=None):
        """Add values to a layer.

        Parameters
        ----------
        layer : str
            One of ``layer_order``.
        values : dict
            Key to value, as read by ``posydon.grids.setup.inlist``.
        source : str or None
            Where the values come from (usually a file path), for
            ``explain``.
        """
        if layer not in self._layers:
            raise ValueError("unknown layer {!r}; layers are {}".format(
                layer, ", ".join(self.layer_order)))
        self._layers[layer].append((source, dict(values)))

    def add_file(self, layer, path, skip_values=LEGACY_SKIPPED_VALUES,
                 transform=None):
        """Read this target's section from an inlist file and add it.

        Parameters
        ----------
        layer : str
            One of ``layer_order``.
        path : str or Path
            An inlist, or a MESA ``*.defaults`` file (which has no section
            lines; its keys are read as this target's section).
        skip_values : iterable of str
            Values to leave out; by default those today's setup drops.
        transform : callable or None
            Applied to the values before they are added, e.g.
            ``drop_wiring_keys``.

        Returns
        -------
        dict
            The values that were added (empty if the file has no such
            section).
        """
        sections = read_inlist(path, default_section=self.section,
                               skip_values=skip_values)
        values = sections.get(self.section, {})
        if transform is not None:
            values = transform(values)
        self.add(layer, values, source=str(path))
        return values

    def _in_order(self):
        """Every (layer, source, values) in the order they are applied."""
        for layer in self.layer_order:
            for source, values in self._layers[layer]:
                yield layer, source, values

    def merged(self):
        """The final values: later layers override earlier ones.

        Returns
        -------
        dict
            Key to value. A key keeps the position where it first appeared
            and the value of the last layer that set it.
        """
        result = {}
        for _, _, values in self._in_order():
            result.update(values)
        return result

    def explain(self, key):
        """Every layer that set a key, lowest first; the last one wins.

        Returns
        -------
        list of dict
            One ``{"layer", "source", "value"}`` per layer entry that sets
            ``key``; empty if no layer sets it.
        """
        return [{"layer": layer, "source": source, "value": values[key]}
                for layer, source, values in self._in_order()
                if key in values]


def pick_whole_file(candidates, layer_order=LAYER_ORDER):
    """Pick the file to use when files are not merged but replaced.

    Parameters
    ----------
    candidates : dict
        Layer name to a path, or to ``None`` if that layer has no file.
    layer_order : sequence of str
        Layer names, lowest first.

    Returns
    -------
    tuple or None
        ``(layer, path)`` of the highest layer with a file, or ``None``.

    Examples
    --------
    >>> pick_whole_file({"mesa_inlists": "a.list",
    ...                  "posydon_inlists": "b.list",
    ...                  "user_inlists": None})
    ('posydon_inlists', 'b.list')
    """
    unknown = set(candidates) - set(layer_order)
    if unknown:
        raise ValueError("unknown layer(s): {}".format(
            ", ".join(sorted(unknown))))
    for layer in reversed(tuple(layer_order)):
        if candidates.get(layer) is not None:
            return layer, candidates[layer]
    return None


def formation_step_count(star, targets_list):
    """How many formation steps a star has.

    Parameters
    ----------
    star : str
        ``"star1"`` or ``"star2"``.
    targets_list : iterable of dict
        Target-to-files maps (e.g. ``InlistFiles.targets`` of each layer).

    Returns
    -------
    int
        One more than the highest ``<star>_formation_step<N>_*`` target; 1
        if there are only ``<star>_formation_every_step_*`` targets; 0 if
        the star has no formation at all.

    Examples
    --------
    >>> formation_step_count("star1", [{"star1_formation_step1_job": []}])
    2
    """
    prefix = "{}_formation_".format(star)
    count = 0
    for targets in targets_list:
        for target in targets:
            if not target.startswith(prefix):
                continue
            rest = target[len(prefix):]
            if rest.startswith("every_step_"):
                count = max(count, 1)
            else:
                count = max(count, int(rest[len("step"):].split("_")[0]) + 1)
    return count


def expand_every_step(targets, star, n_steps):
    """Give each formation step its own targets, every-step files first.

    ``<star>_formation_every_step_<part>`` holds files that apply to all
    steps. They are put in front of each step's own files, as today's setup
    does with one file per step.

    Parameters
    ----------
    targets : dict
        Target to list of files, for one layer.
    star : str
        ``"star1"`` or ``"star2"``.
    n_steps : int
        From ``formation_step_count`` over all layers.

    Returns
    -------
    dict
        A new target map without ``every_step`` targets for ``star``.

    Examples
    --------
    >>> expanded = expand_every_step(
    ...     {"star1_formation_every_step_job": ["common"],
    ...      "star1_formation_step1_job": ["step2"]}, "star1", 2)
    >>> expanded["star1_formation_step0_job"]
    ['common']
    >>> expanded["star1_formation_step1_job"]
    ['common', 'step2']
    """
    every = "{}_formation_every_step_".format(star)
    result = {target: list(files) for target, files in targets.items()
              if not target.startswith(every)}
    for part in ("controls", "job"):
        shared = targets.get(every + part, [])
        if not shared:
            continue
        for step in range(n_steps):
            target = "{}_formation_step{}_{}".format(star, step, part)
            result[target] = list(shared) + result.get(target, [])
    return result


def drop_wiring_keys(values):
    """Remove keys that wire inlists together, as today's setup does.

    Input inlists may point MESA at further inlists (``read_extra_*`` and
    ``*inlist*`` keys). Setup writes its own wiring as generated keys, so it
    drops these from every input file.

    Examples
    --------
    >>> drop_wiring_keys({"read_extra_controls_inlist1": ".true.",
    ...                   "extra_controls_inlist1_name": "'x'",
    ...                   "mesh_delta_coeff": "1d0"})
    {'mesh_delta_coeff': '1d0'}
    """
    return {key: value for key, value in values.items()
            if "read_extra" not in key and "inlist" not in key}


def fix_num_x_ctrls(values):
    """Replace the ``num_x_ctrls`` placeholder in key names by ``1``.

    MESA's ``controls.defaults`` writes ``x_ctrl(1:num_x_ctrls) = 0d0``,
    which MESA itself cannot read back. Today's setup replaces the
    placeholder, so only element 1 gets the default (tracker item C2).

    Examples
    --------
    >>> fix_num_x_ctrls({"x_ctrl(1:num_x_ctrls)": "0d0"})
    {'x_ctrl(1:1)': '0d0'}
    """
    return {key.replace("num_x_ctrls", "1"): value
            for key, value in values.items()}
