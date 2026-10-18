"""Warn about inlist keys that MESA does not know.

MESA stops when a namelist holds a key it does not know, but only when the
run starts, on the cluster. A misspelled key (``max_modle_number``) or one
in the wrong namelist (``max_model_number`` in ``&star_job``) passes the
setup. ``check_keys`` finds these in the ini's own files and warns; it
does not stop the setup.

The keys of a namelist are those of MESA's ``.defaults`` file for it, of
the MESA version the grid uses:

- **with a release** (``[scenario] release``): the defaults of the
  release's MESA installation, ``<mesa_version>/mesa/<installation>/
  defaults/``, the files ``fill_defaults`` puts into ``mesa_inlists``
  (``manifest.mesa_installation``). This holds also for a target whose
  defaults ``[mesa_inlists]`` replaces: the release names the MESA
  version, so a replacement is checked against that version's keys, like
  any other file.
- **without a release** (a legacy ini, or a new-format ini without a
  scenario release): the ``.defaults`` files of the ``mesa_inlists`` layer
  for that namelist; a legacy ini names them in ``$MESA_DIR``, the MESA it
  compiles. A namelist with no such file is not checked.

Checked are the files of ``[user_inlists]`` and ``[mesa_inlists]`` (the
legacy ``*_user`` and ``*_mesa_defaults`` keys), other than the reference
files themselves; the release's own files are not (they are checked
against their installation when the inlist repository is tested).

A key is compared by name: Fortran names are case-insensitive, and an
array element (``x_logical_ctrl(5)``) is the array (``x_logical_ctrl``).
A ``.defaults`` file comments out some keys, e.g.
``! xa_central_lower_limit_species(1) = ''``, or only names them in a
heading, ``!### accretion_species_id``; those count as keys too.
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

import re
from dataclasses import dataclass
from pathlib import Path

from posydon.grids.setup.build import layer_targets
from posydon.grids.setup.inlist import read_inlist
from posydon.grids.setup.layers import drop_wiring_keys, section_for
from posydon.grids.setup.manifest import mesa_installation
from posydon.utils.posydonwarning import Pwarn

# manifest.MESA_DEFAULTS name -> the namelist its keys are for
DEFAULTS_NAMELISTS = {"binary_controls": "binary_controls",
                      "binary_job": "binary_job",
                      "star_controls": "controls",
                      "star_job": "star_job"}
CHECKED_LAYERS = ("mesa_inlists", "user_inlists")

# "key = value" or "key(1:3) = value", after a leading "!" is removed
_ASSIGNMENT = re.compile(r"^([a-z_][a-z0-9_]*)\s*(\([^)]*\))?\s*=",
                         re.IGNORECASE)
# a "!### key" heading of a .defaults file
_HEADING = re.compile(r"^###\s*([a-z_][a-z0-9_]*)\s*$", re.IGNORECASE)


def key_name(key):
    """The name MESA knows a key by: lower case, without an index.

    Examples
    --------
    >>> key_name("X_Logical_Ctrl(5)")
    'x_logical_ctrl'
    >>> key_name("accretion_species_id (1)")
    'accretion_species_id'
    """
    return key.split("(", 1)[0].strip().lower()


def defaults_keys(path):
    """The key names of a MESA ``.defaults`` file.

    Keys set in the file, keys in commented-out assignments and keys named
    in ``!###`` headings.

    Examples
    --------
    >>> import tempfile, os
    >>> with tempfile.NamedTemporaryFile("w", suffix=".defaults",
    ...                                  delete=False) as handle:
    ...     _ = handle.write("!### xa_central_lower_limit_species\\n"
    ...                      "  ! xa_central_lower_limit(1) = 0\\n"
    ...                      "  X_CTRL(1:num_x_ctrls) = 0d0  ! x\\n")
    >>> sorted(defaults_keys(handle.name))
    ['x_ctrl', 'xa_central_lower_limit', 'xa_central_lower_limit_species']
    >>> os.remove(handle.name)
    """
    keys = set()
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if line.startswith("!"):
            line = line[1:].strip()
            heading = _HEADING.match(line)
            if heading:
                keys.add(heading.group(1).lower())
                continue
        assignment = _ASSIGNMENT.match(line)
        if assignment:
            keys.add(assignment.group(1).lower())
    return keys


def reference_files(config, resolved=None):
    """The ``.defaults`` files that say each namelist's keys.

    Parameters
    ----------
    config : GridConfig
        After ``defaults.fill_defaults``.
    resolved : ResolvedScenario or None
        From ``resolve.resolve_scenario``; its release names the MESA
        installation.

    Returns
    -------
    dict
        Namelist (``controls``, ``star_job``, ``binary_controls``,
        ``binary_job``) to a list of existing files; a namelist without one
        is not checked.
    """
    release = resolved.release if resolved is not None else None
    if release is not None:
        return {DEFAULTS_NAMELISTS[name]: [path] for name, path in
                mesa_installation(release).items()}
    out = {}
    for target, paths in config.inlists["mesa_inlists"].targets.items():
        for path in paths:
            if Path(path).suffix == ".defaults" and Path(path).is_file():
                files = out.setdefault(section_for(target), [])
                if Path(path) not in files:
                    files.append(Path(path))
    return out


@dataclass
class UnknownKeys:
    """Keys of one file that are not keys of the namelist they land in."""

    layer: str
    path: Path
    namelist: str
    keys: list
    references: list

    def message(self):
        """What is wrong, naming the keys, the namelist and the file."""
        many = len(self.keys) > 1
        return ("[{}] {}: {} {} not {} of &{} in this MESA version "
                "({})".format(
                    self.layer, self.path, ", ".join(self.keys),
                    "are" if many else "is", "keys" if many else "a key",
                    self.namelist,
                    ", ".join(str(path) for path in self.references)))


def unknown_keys(config, resolved=None):
    """Keys of the ini's own files that MESA does not know.

    Parameters are those of ``reference_files``.

    Returns
    -------
    list of UnknownKeys
        One per file and namelist with unknown keys, in the order of the
        layers and targets. A file that feeds several targets of one
        namelist (star 1 and star 2, formation steps) is listed once.
    """
    references = reference_files(config, resolved)
    known = {namelist: set().union(*(defaults_keys(path) for path in paths))
             for namelist, paths in references.items()}
    scenario = resolved.files if resolved is not None else None
    files, _ = layer_targets(config, scenario)
    seen = set()
    out = []
    for layer in CHECKED_LAYERS:
        for target, paths in files[layer].items():
            namelist = section_for(target)
            if namelist not in known:
                continue
            for path in map(Path, paths):
                if (path, namelist) in seen or \
                        path in references[namelist] or not path.is_file():
                    # a missing file is reported by the step that reads it
                    continue
                seen.add((path, namelist))
                # wiring keys are dropped: setup writes its own
                values = drop_wiring_keys(read_inlist(
                    path, default_section=namelist).get(namelist, {}))
                keys = [key for key in values
                        if key_name(key) not in known[namelist]]
                if keys:
                    out.append(UnknownKeys(layer, path, namelist, keys,
                                           references[namelist]))
    return out


def check_keys(config, resolved=None):
    """Warn about every key MESA does not know; see ``unknown_keys``.

    One warning lists them all: ``Pwarn`` shows a warning once per line of
    code, so a warning per file would show only the first.

    Returns
    -------
    list of UnknownKeys
    """
    found = unknown_keys(config, resolved)
    if found:
        Pwarn("MESA stops at a key it does not know, when the run starts: "
              + "; ".join(unknown.message() for unknown in found),
              "InappropriateValueWarning")
    return found
