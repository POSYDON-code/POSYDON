"""Say where each inlist value of a grid comes from (``--explain``).

For every target the setup writes, the layers are stacked as for the real
inlists (``build.layer_targets``, ``build.stack_target``), and every key
whose final value does not come from MESA's defaults is listed with the
layer and file that set it, and the values it overrode::

    star1_controls
      mesh_delta_coeff = 0.8d0  [user_inlists] my_inlist1
          over 1d0  [posydon_inlists] .../common/binary/inlist1

Nothing is written. Keys the setup generates itself (wiring, output
options, ZAMS models) come on top of these when the inlists are written.
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

from pathlib import Path

from posydon.grids.setup.build import (
    is_single_star,
    layer_targets,
    stack_target,
)

# targets a single-star grid does not write (it runs star 1 in steps)
BINARY_TARGETS = ("binary_controls", "binary_job", "star1_controls",
                  "star1_job", "star2_controls", "star2_job")


def _order(target):
    """binary, star1, star2, then formation steps by star and step."""
    if target.startswith("binary"):
        return (0, target)
    star, _, rest = target.partition("_")
    return (1 if not rest.startswith("formation") else 2, star, rest)


def explained_targets(config, scenario=None):
    """The targets the setup writes, each with its explained keys.

    Returns
    -------
    dict
        Target to a list of ``(key, history)``, where ``history`` is
        ``LayerStack.explain(key)`` (lowest layer first; the last entry
        wins), for every key not set by ``mesa_inlists`` alone.
    """
    files, flags = layer_targets(config, scenario)
    targets = {target for layer in files.values() for target in layer}
    targets |= set(flags)
    targets = {target for target in targets
               if "every_step" not in target}
    if is_single_star(config):
        targets -= set(BINARY_TARGETS)
    out = {}
    for target in sorted(targets, key=_order):
        stack = stack_target(target, files, flags)
        out[target] = [(key, stack.explain(key)) for key in stack.merged()
                       if stack.explain(key)[-1]["layer"] != "mesa_inlists"]
    return out


def _short(source, roots):
    """A source path relative to the first root that holds it."""
    for root in roots:
        try:
            return Path(source).relative_to(root).as_posix()
        except ValueError:
            continue
    return source


def explain_text(config, scenario=None, roots=()):
    """``explained_targets`` as text, for the terminal.

    Parameters
    ----------
    roots : sequence of str or Path
        Folders to show source paths relative to, e.g. the inlist
        repository.
    """
    lines = []
    for target, keys in explained_targets(config, scenario).items():
        lines.append(target)
        if not keys:
            lines.append("  (MESA defaults only)")
        for key, history in keys:
            final = history[-1]
            same = len(history) > 1 and \
                history[-2]["value"] == final["value"]
            lines.append("  {} = {}  [{}] {}{}".format(
                key, final["value"], final["layer"],
                _short(final["source"], roots),
                "  (same value as below)" if same else ""))
            for entry in reversed(history[:-1]):
                lines.append("      over {}  [{}] {}".format(
                    entry["value"], entry["layer"],
                    _short(entry["source"], roots)))
        lines.append("")
    lines.append("Keys the setup generates (wiring, output options, ZAMS "
                 "models) are added on top when the inlists are written.")
    return "\n".join(lines)
