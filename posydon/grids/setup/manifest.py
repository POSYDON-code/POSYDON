"""Read a release of POSYDON-MESA-INLISTS from its manifests.

Since the restructure of the inlist repository, a grid type's files are no
longer found by rules in code (``scenario.scenario_files``) but listed in
the repository itself::

    r11701/                         one folder per MESA version
      mesa/<installation>/          MESA's stock defaults/, columns/, src/
      ZAMS_models/
      releases/<release>/
        release.ini                 [release] mesa_version, mesa, frozen
        common/ ...                 files shared by the grid types
        <grid type>/scenario.ini    which files each target uses, in order
        reruns/<rerun>/scenario.ini overlays and replacements

``release_files`` reads one grid type of one release, optionally with a
rerun, into the same ``ScenarioFiles`` that ``build.build_inlists`` takes,
so the rest of the setup does not change. A rerun works file by file:

- ``[overlays]``: an inlist holding only the keys the rerun changes, used
  right after the release file of the same path, wherever that file is
  used;
- ``[replace]``: a whole file (a Fortran source) used instead of the
  release's.

``mesa_installation`` gives the MESA defaults of the release's
installation, the ``mesa_inlists`` layer.
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

import configparser
from dataclasses import dataclass
from pathlib import Path

from posydon.grids.setup.config import ConfigError
from posydon.grids.setup.scenario import ScenarioFiles

# manifest target -> the two targets of the setup
PARTS = ("controls", "job")
STARTS = ("zams", "formation", "point_mass", "none")
# MESA's defaults per setup target, inside mesa/<installation>/defaults/
MESA_DEFAULTS = {"binary_controls": "binary/binary_controls.defaults",
                 "binary_job": "binary/binary_job.defaults",
                 "star_controls": "star/controls.defaults",
                 "star_job": "star/star_job.defaults"}


@dataclass
class Release:
    """A release folder and what its ``release.ini`` says."""

    name: str
    folder: Path
    mesa_version: str
    mesa: str
    frozen: bool


def _read_ini(path, what):
    if not Path(path).is_file():
        raise ConfigError("scenario", what, "{} not found".format(path))
    parser = configparser.ConfigParser(inline_comment_prefixes=(";",))
    parser.optionxform = str  # keys are case sensitive (x_logical_ctrl)
    parser.read(path)
    return parser


def _section(parser, name):
    return dict(parser[name]) if parser.has_section(name) else {}


def _paths(value):
    return [part.strip() for part in value.split(",") if part.strip()]


###############################################################################
# finding a release
###############################################################################
def find_release(repo, release, mesa_version=None):
    """The folder of a release, inferring the MESA version if not given.

    Parameters
    ----------
    repo : str or Path
        A checked-out POSYDON-MESA-INLISTS.
    release : str
        e.g. ``v2.0`` or ``development``.
    mesa_version : str or None
        e.g. ``r11701``; only needed if the release exists for several.

    Returns
    -------
    Release

    Raises
    ------
    ConfigError
        If the release is not found, or found for several MESA versions
        and ``mesa_version`` does not pick one.
    """
    repo = Path(repo)
    if mesa_version is not None:
        candidates = [repo / mesa_version / "releases" / release]
    else:
        candidates = sorted(repo.glob("*/releases/{}".format(release)))
    found = [folder for folder in candidates if folder.is_dir()]
    if not found:
        raise ConfigError("scenario", "release", "no release {!r} in {}"
                          .format(release, repo))
    if len(found) > 1:
        raise ConfigError("scenario", "mesa", "release {!r} exists for "
                          "several MESA versions ({}); set mesa".format(
                              release, ", ".join(sorted(
                                  f.parents[1].name for f in found))))
    folder = found[0]
    info = _section(_read_ini(folder / "release.ini", "release"), "release")
    return Release(name=release, folder=folder,
                   mesa_version=info.get("mesa_version",
                                         folder.parents[1].name),
                   mesa=info.get("mesa", "default"),
                   frozen=info.get("frozen", "false").lower() == "true")


def list_reruns(release):
    """The names of a release's reruns."""
    folder = release.folder / "reruns"
    return sorted(path.name for path in folder.iterdir()
                  if (path / "scenario.ini").is_file()) \
        if folder.is_dir() else []


###############################################################################
# a grid type, with a rerun
###############################################################################
def _rerun_files(release, rerun):
    """A rerun's overlays and replacements, by release-relative path."""
    folder = release.folder / "reruns" / rerun
    if not folder.is_dir():
        raise ConfigError("scenario", "rerun", "no rerun {!r} in release "
                          "{} (has: {})".format(rerun, release.name,
                                                ", ".join(list_reruns(
                                                    release)) or "none"))
    parser = _read_ini(folder / "scenario.ini", "rerun")
    overlays = {key: folder / value
                for key, value in _section(parser, "overlays").items()}
    replace = {key: folder / value
               for key, value in _section(parser, "replace").items()}
    return overlays, replace


def release_files(release, grid_type, rerun=None):
    """The files one grid type takes from a release, as ``ScenarioFiles``.

    Parameters
    ----------
    release : Release
        From ``find_release``.
    grid_type : str
        A grid-type name, e.g. ``CO-HeMS``; its folder holds the manifest.
    rerun : str or None
        A rerun of the release, applied file by file.

    Returns
    -------
    ScenarioFiles
        ``inlists`` per target (``binary`` -> ``binary_controls`` and
        ``binary_job``, and so on), column lists, extras, ZAMS models and
        the stars built from formation steps. Flags are files here, so
        ``flags`` stays empty.
    """
    manifest_path = release.folder / grid_type / "scenario.ini"
    parser = _read_ini(manifest_path, "grid type")
    overlays, replace = ({}, {}) if rerun is None else \
        _rerun_files(release, rerun)

    def resolve(rel):
        if rel in replace:
            return [replace[rel]]
        path = release.folder / rel
        if not path.is_file():
            raise ConfigError("scenario", "release", "{} names {}, which "
                              "does not exist".format(manifest_path, rel))
        return [path] + ([overlays[rel]] if rel in overlays else [])

    result = ScenarioFiles()
    for base, value in _section(parser, "inlists").items():
        paths = [path for rel in _paths(value) for path in resolve(rel)]
        for part in PARTS:
            for path in paths:
                result.inlists.add("{}_{}".format(base, part), path)
    for name, value in _section(parser, "columns").items():
        result.inlists.whole_files[name] = resolve(value)[0]
    for key, value in _section(parser, "extras").items():
        result.extras[key] = resolve(value)[0]
    zams = _section(parser, "zams")
    root = release.folder.parents[1]
    for star in ("star1", "star2"):
        if star in zams:
            setattr(result, "zams_" + star, root / zams[star])
    starts = _section(parser, "grid_type")
    for star in ("star1", "star2"):
        start = starts.get(star, "none")
        if start not in STARTS:
            raise ConfigError("scenario", "release", "{}: [grid_type] {} = "
                              "{!r}; expected one of {}".format(
                                  manifest_path, star, start,
                                  ", ".join(STARTS)))
        if start == "formation":
            result.formation_stars += (star,)
    return result


###############################################################################
# the MESA installation
###############################################################################
def mesa_installation(release):
    """MESA's default files for the release's installation.

    Returns
    -------
    dict
        ``binary_controls``, ``binary_job``, ``star_controls`` and
        ``star_job`` to the ``.defaults`` files in
        ``<mesa_version>/mesa/<installation>/defaults/``.
    """
    folder = release.folder.parents[1] / "mesa" / release.mesa / "defaults"
    paths = {name: folder / rel for name, rel in MESA_DEFAULTS.items()}
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise ConfigError("scenario", "release", "MESA installation {!r} "
                          "is missing {}".format(release.mesa,
                                                 ", ".join(missing)))
    return paths
