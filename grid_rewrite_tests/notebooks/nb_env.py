"""Shared set-up for the grid-setup notebooks.

Every notebook starts with::

    from nb_env import env
    ENV = env()

which

- finds a checkout of POSYDON-MESA-INLISTS with the restructured layout
  (``r11701/releases/v2.0``): ``$POSYDON_INLIST_REPO`` if set, otherwise
  ``~/.posydon_mesa_inlists``. The notebooks read it with ``source = path``,
  so nothing in it is checked out or changed;
- points ``$MESA_DIR`` at a stub MESA (no-op makefiles, labelled stub
  Fortran sources) copied into ``_work/``, so a setup runs in seconds
  without MESA. Set ``POSYDON_NOTEBOOK_REAL_MESA=1`` to keep your own
  ``$MESA_DIR`` (the setup then really compiles, which takes minutes);
- makes Python import ``posydon`` from this checkout (first on
  ``sys.path`` and ``PYTHONPATH``), whatever POSYDON is installed in the
  environment, and puts this checkout's ``bin/`` first on the ``PATH``, so
  ``posydon-grid`` and ``posydon-run-grid`` are the ones of this branch.

Everything the notebooks write goes into ``_work/`` next to them, which git
ignores. Delete it whenever you like.
"""

import os
import shutil
import sys
from pathlib import Path
from types import SimpleNamespace

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]                      # the POSYDON checkout
WORK = HERE / "_work"
MESA_STUB = ROOT / "grid_rewrite_tests" / "setup_grid_data" / "mesa_stub"
EXTRA = ROOT / "grid_rewrite_tests" / "setup_grid_data" / "extra"

# this checkout's posydon, not the one installed in the environment
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if "posydon" in sys.modules and \
        not str(Path(sys.modules["posydon"].__file__)).startswith(str(ROOT)):
    raise RuntimeError("posydon was imported from {} before nb_env; restart "
                       "the kernel and import nb_env first".format(
                           sys.modules["posydon"].__file__))


def _inlist_repo():
    given = os.environ.get("POSYDON_INLIST_REPO")
    repo = Path(given or "~/.posydon_mesa_inlists").expanduser()
    if not (repo / "r11701" / "releases" / "v2.0").is_dir():
        raise RuntimeError(
            "{} is not a checkout of POSYDON-MESA-INLISTS with "
            "r11701/releases/v2.0. Clone it and check out v2_sources:\n"
            "  git clone https://github.com/POSYDON-code/"
            "POSYDON-MESA-INLISTS.git ~/inlists\n"
            "  git -C ~/inlists checkout v2_sources\n"
            "then set POSYDON_INLIST_REPO=~/inlists before starting "
            "Jupyter (or in the first cell, os.environ[...])."
            .format(repo))
    return repo


def env():
    """Prepare the environment; return the paths the notebooks use.

    Returns
    -------
    SimpleNamespace
        ``root`` (POSYDON), ``work`` (``_work/``), ``inlists`` (the inlist
        repository), ``release`` (its ``r11701/releases/v2.0``), ``mesa``
        (``$MESA_DIR``), ``extra`` (the test suite's example user inlists).
    """
    WORK.mkdir(exist_ok=True)
    repo = _inlist_repo()
    if os.environ.get("POSYDON_NOTEBOOK_REAL_MESA") == "1":
        mesa = Path(os.environ["MESA_DIR"])
    else:
        mesa = WORK / "mesa_stub"
        if not mesa.is_dir():
            shutil.copytree(MESA_STUB, mesa)
        os.environ["MESA_DIR"] = str(mesa)
    os.environ.setdefault("MESASDK_ROOT", "/opt/mesasdk-stub")
    pythonpath = os.environ.get("PYTHONPATH", "")
    if not pythonpath.startswith(str(ROOT)):
        os.environ["PYTHONPATH"] = os.pathsep.join(
            [str(ROOT)] + ([pythonpath] if pythonpath else []))
    path = os.environ.get("PATH", "")
    if not path.startswith(str(ROOT / "bin")):
        os.environ["PATH"] = os.pathsep.join([str(ROOT / "bin"), path])
    return SimpleNamespace(root=ROOT, work=WORK, inlists=repo,
                           release=repo / "r11701" / "releases" / "v2.0",
                           mesa=mesa, extra=EXTRA)


def fresh(name):
    """An empty folder ``_work/<name>`` (emptied if it exists)."""
    folder = WORK / name
    if folder.exists():
        shutil.rmtree(folder)
    folder.mkdir(parents=True)
    return folder


def write(path, text):
    """Write a text file (dedented, leading newline dropped); return it."""
    import textwrap
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(text).lstrip("\n"))
    return path


def show(path, head=None, grep=None):
    """Print a text file, or its first ``head`` lines, or only the lines
    containing ``grep``."""
    lines = Path(path).read_text().splitlines()
    if grep is not None:
        lines = [line for line in lines if grep in line]
    if head is not None and len(lines) > head:
        lines = lines[:head] + ["... ({} more lines)".format(
            len(lines) - head)]
    print("\n".join(lines))


def tree(folder, skip=("make", "src")):
    """Print the files under ``folder``, relative to it.

    Folders named in ``skip`` are shown with their file count only.
    """
    folder = Path(folder)
    for path in sorted(folder.rglob("*")):
        rel = path.relative_to(folder)
        if any(part in skip for part in rel.parts[:-1]):
            continue
        if path.is_dir() and path.name in skip:
            print("{}/  ({} files)".format(rel, sum(
                1 for p in path.rglob("*") if p.is_file())))
        elif path.is_file():
            print(rel)


def short(path, *roots):
    """A path relative to the first root that holds it, for printing."""
    for root in roots:
        try:
            return str(Path(path).relative_to(root))
        except ValueError:
            continue
    return str(path)


BINARY_POINTS = ("initial_z,Zbase,m1,m2,initial_period_in_days\n"
                 "0.0142,0.0142,30.0,20.0,10.0\n"
                 "0.0142,0.0142,50.0,25.0,3.0\n"
                 "0.0142,0.0142,12.5,10.0,100.0\n")
SINGLE_POINTS = ("initial_z,Zbase,initial_mass\n"
                 "0.0142,0.0142,10.0\n"
                 "0.0142,0.0142,20.0\n")


def grid_csv(folder, single=False):
    """Write a small ``grid.csv`` (three binaries or two single stars)."""
    return write(Path(folder) / "grid.csv",
                 SINGLE_POINTS if single else BINARY_POINTS)


def cli(*args, cwd=None, show_stderr=True):
    """Run this branch's ``posydon-grid`` with ``args`` and print its output.

    Uses the notebook's own Python, so it works whatever ``python`` is on
    the PATH. Returns the exit code.
    """
    import subprocess
    result = subprocess.run(
        [sys.executable, str(ROOT / "bin" / "posydon-grid"), *map(str, args)],
        cwd=cwd, capture_output=True, text=True)
    print(result.stdout, end="")
    if show_stderr and result.stderr.strip():
        print("--- stderr ---")
        print(result.stderr, end="")
    return result.returncode
