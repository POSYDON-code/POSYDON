"""The ``posydon-grid`` command line.

::

    posydon-grid types                    the grid types
    posydon-grid init CO-HeMS grid.ini    a starter ini for a grid type
    posydon-grid explain grid.ini         where each inlist value comes from;
                                          nothing is written
    posydon-grid setup grid.ini           write the run directory
    posydon-grid run ...                  run the grid: the arguments go to
                                          posydon-run-grid unchanged

A new-format ini is the only input of ``explain`` and ``setup``. A legacy
ini cannot say the grid type, run directory and submission type, so for
it they are given as ``--grid-type``, ``--run-directory`` and
``--submission-type``.

A mistake the user can fix (a ``SetupError``: a wrong ini value, a
missing file, an unset environment variable) is printed as one line,
``posydon-grid: error: ...``, with exit status 2, as argparse does for a
wrong option; anything else is a bug and keeps its traceback.

``posydon-setup-grid`` (``legacy_main``) still works with its old options,
says it is deprecated and prints the ``posydon-grid`` command to use.

Only the command line lives here; the steps of a setup are
``posydon.grids.setup.steps``.
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

import argparse
import os
import shlex
import sys

from posydon.grids.setup.config import GRID_MODES, SetupError
from posydon.grids.setup.grid_types import describe_grid_types
from posydon.grids.setup.steps import (
    SUBMISSION_TYPES,
    explain_grid,
    setup_grid,
)
from posydon.grids.setup.submit import find_run_grid
from posydon.grids.setup.template import ini_template
from posydon.utils.posydonwarning import Pwarn

LEGACY_OPTIONS = ("grid_type", "run_directory", "submission_type")


def _add_legacy_options(parser):
    group = parser.add_argument_group(
        "legacy ini files only", "a new-format ini says this in [grid] "
        "and [submission]")
    group.add_argument("--grid-type", choices=GRID_MODES,
                       help="fixed (a grid of points) or dynamic (points "
                       "sampled by psy-cris)")
    group.add_argument("--run-directory",
                       help="where the executables, inlists and output go "
                       "(default: the current directory)")
    group.add_argument("--submission-type", choices=SUBMISSION_TYPES,
                       help="a shell script or SLURM scripts (default: "
                       "shell)")


def build_parser():
    """The ``posydon-grid`` argument parser."""
    parser = argparse.ArgumentParser(
        prog="posydon-grid",
        description="Set up and run grids of MESA runs for POSYDON.")
    commands = parser.add_subparsers(title="commands", dest="command")
    commands.required = True

    commands.add_parser("types", help="list the grid types")

    init = commands.add_parser("init", help="write a starter ini file for "
                               "a grid type")
    init.add_argument("grid_type", metavar="GRID_TYPE",
                      help="see posydon-grid types")
    init.add_argument("inifile", metavar="INIFILE", help="the file to write")

    explain = commands.add_parser(
        "explain", help="show where each inlist value comes from; write "
        "nothing")
    explain.add_argument("inifile", metavar="INIFILE")
    _add_legacy_options(explain)

    setup = commands.add_parser("setup", help="write the run directory of "
                                "a grid")
    setup.add_argument("inifile", metavar="INIFILE")
    _add_legacy_options(setup)

    # listed for --help only: main() passes "run ..." on before parsing,
    # since argparse cannot pass on options it does not know
    commands.add_parser("run", help="run the grid (all arguments go to "
                        "posydon-run-grid unchanged)", add_help=False)
    return parser


def init_command(grid_type, inifile):
    """Write a starter ini; refuse to overwrite one."""
    if os.path.exists(inifile):
        raise SetupError("{} exists; not overwriting it".format(inifile))
    with open(inifile, "w") as handle:
        handle.write(ini_template(grid_type))
    print("Wrote {}; fill in the <...> values.".format(inifile))


def run_command(args):
    """Replace this process by ``posydon-run-grid`` with the same arguments.

    ``posydon-run-grid`` itself is not changed; this only gives it its
    place under ``posydon-grid``.
    """
    run_grid = find_run_grid()
    os.execv(sys.executable, [sys.executable, run_grid] + list(args))


def _options(args, overrides=False):
    options = {name: getattr(args, name) for name in LEGACY_OPTIONS}
    options["overrides"] = overrides
    return options


def fail(prog, error):
    """Print a user error as one line and exit with status 2.

    The status and the ``<prog>: error:`` prefix are those of argparse,
    so a wrong ini value reads like a wrong option.
    """
    print("{}: error: {}".format(prog, error), file=sys.stderr)
    sys.exit(2)


def main(argv=None):
    """Run ``posydon-grid`` with ``argv`` (default: the command line)."""
    argv = sys.argv[1:] if argv is None else list(argv)
    try:
        _main(argv)
    except SetupError as error:
        fail("posydon-grid", error)


def _main(argv):
    if argv[:1] == ["run"]:
        run_command(argv[1:])
        return
    args = build_parser().parse_args(argv)
    if args.command == "types":
        print(describe_grid_types())
    elif args.command == "init":
        init_command(args.grid_type, args.inifile)
    elif args.command == "explain":
        print(explain_grid(args.inifile, **_options(args)))
    else:
        setup_grid(args.inifile, command=" ".join(
            ["posydon-grid"] + [shlex.quote(arg) for arg in argv]),
            **_options(args))


###############################################################################
# posydon-setup-grid, deprecated
###############################################################################
def build_legacy_parser():
    """The argument parser of the deprecated ``posydon-setup-grid``."""
    parser = argparse.ArgumentParser(
        prog="posydon-setup-grid",
        description="Deprecated: use posydon-grid (types, init, explain, "
        "setup, run).")
    parser.add_argument("--inifile", help="Name of ini file of params")
    parser.add_argument("--grid-type", choices=GRID_MODES)
    parser.add_argument("--run-directory")
    parser.add_argument("--submission-type", choices=SUBMISSION_TYPES)
    parser.add_argument("-n", "--nproc", type=int, default=1,
                        help="not used")
    parser.add_argument("--verbose", action="store_true", help="not used")
    parser.add_argument("--list-grid-types", action="store_true")
    parser.add_argument("--init", metavar="GRID_TYPE")
    parser.add_argument("--explain", action="store_true")
    return parser


def equivalent_command(args):
    """The ``posydon-grid`` command doing what ``posydon-setup-grid`` does.

    Examples
    --------
    >>> args = build_legacy_parser().parse_args(
    ...     ["--inifile", "grid.ini", "--grid-type", "fixed"])
    >>> equivalent_command(args)
    'posydon-grid setup grid.ini --grid-type fixed'
    """
    if args.list_grid_types:
        words = ["types"]
    elif args.init is not None:
        words = ["init", args.init, args.inifile]
    else:
        words = ["explain" if args.explain else "setup", args.inifile]
        for name in LEGACY_OPTIONS:
            if getattr(args, name) is not None:
                words += ["--" + name.replace("_", "-"), getattr(args, name)]
    return "posydon-grid " + " ".join(shlex.quote(str(word))
                                      for word in words)


def legacy_main(argv=None):
    """Run the deprecated ``posydon-setup-grid`` with ``argv``."""
    parser = build_legacy_parser()
    args = parser.parse_args(argv)
    if args.inifile is None and not args.list_grid_types:
        parser.error("--inifile is required")
    try:
        _legacy_main(args, argv)
    except SetupError as error:
        fail("posydon-setup-grid", error)


def _legacy_main(args, argv):
    Pwarn("posydon-setup-grid is deprecated and will be removed; use "
          "'{}' instead (posydon-grid has the commands types, init, "
          "explain, setup and run)".format(equivalent_command(args)),
          "DeprecationWarning")
    if args.list_grid_types:
        print(describe_grid_types())
    elif args.init is not None:
        init_command(args.init, args.inifile)
    elif args.explain:
        print(explain_grid(args.inifile, **_options(args, overrides=True)))
    else:
        setup_grid(args.inifile, command=" ".join(
            ["posydon-setup-grid"] + [shlex.quote(arg) for arg in (
                sys.argv[1:] if argv is None else argv)]),
            **_options(args, overrides=True))
