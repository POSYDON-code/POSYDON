"""Baseline (golden-file) tests of bin/posydon-setup-grid.

These tests freeze what today's posydon-setup-grid writes for every
supported grid type, with and without a scenario, so that the planned rewrite
can be checked against it. They run the script as a black box in a temporary
directory with:

- a stub MESA_DIR (no-op makefiles, labelled stub Fortran sources),
- a throwaway local git clone of a vendored POSYDON-MESA-INLISTS snapshot in
  a temporary HOME (so the scenario code never touches the network),
- a stub posydon-run-grid on PATH.

Each case is compared with grid_rewrite_tests/setup_grid_data/golden/.
Set POSYDON_UPDATE_GOLDEN=1 to regenerate the golden files. Set
POSYDON_GOLDEN_MODE=semantic to ignore the sha256 lines, i.e. to compare only
the non-default namelist values, the scripts and the staged files.
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

import hashlib
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SETUP_GRID = ROOT / "bin" / "posydon-setup-grid"
# today's script before the rewrite (stage 3.6), kept for the comparisons in
# test_new_vs_legacy.py
LEGACY_SETUP_GRID = Path(__file__).parent / "legacy_posydon_setup_grid.py"
DATA = Path(__file__).resolve().parent / "setup_grid_data"
GOLDEN = DATA / "golden"
# MESA r11701's defaults, from POSYDON-MESA-INLISTS r11701/mesa/default
MESA_DEFAULTS = Path(__file__).resolve().parent / "setup_grid_data" \
    / "mesa_defaults"

UPDATE_GOLDEN = os.environ.get("POSYDON_UPDATE_GOLDEN", "") == "1"
SEMANTIC = os.environ.get("POSYDON_GOLDEN_MODE", "strict") == "semantic"

# fixed identity and dates make the fixture commits reproducible
GIT_ENV = {
    "GIT_AUTHOR_NAME": "posydon-tests",
    "GIT_AUTHOR_EMAIL": "tests@posydon.invalid",
    "GIT_COMMITTER_NAME": "posydon-tests",
    "GIT_COMMITTER_EMAIL": "tests@posydon.invalid",
    "GIT_AUTHOR_DATE": "2026-01-01T00:00:00+0000",
    "GIT_COMMITTER_DATE": "2026-01-01T00:00:00+0000",
}

GRID_COLUMNS = {
    "binary": "initial_z,Zbase,m1,m2,initial_period_in_days",
    "binary_new_Z": "initial_z,Zbase,new_Z,m1,m2,initial_period_in_days",
    "single": "initial_z,Zbase,initial_mass",
    "single_new_Z": "initial_z,Zbase,new_Z,initial_mass",
}
GRID_ROWS = {
    "binary": ["0.0142,0.0142,30.0,20.0,10.0",
               "0.0142,0.0142,50.0,25.0,3.0",
               "0.0142,0.0142,12.5,10.0,100.0"],
    "binary_new_Z": ["0.0142,0.0142,0.0142,8.0,10.0,1.0",
                     "0.0142,0.0142,0.0142,15.0,20.0,5.0",
                     "0.0142,0.0142,0.0142,3.0,5.0,0.5"],
    "single": ["0.0142,0.0142,10.0",
               "0.0142,0.0142,30.0",
               "0.0142,0.0142,100.0"],
    "single_new_Z": ["0.0142,0.0142,0.0142,2.0",
                     "0.0142,0.0142,0.0142,8.0",
                     "0.0142,0.0142,0.0142,40.0"],
}

# output controls as in POSYDON-MESA-INLISTS r11701/running_scripts/*.ini
OUTPUT_ALL = dict(history_interval=1, binary_history=True,
                  history_star1=True, final_profile_star1=False,
                  final_model_star1=True, history_star2=True,
                  final_profile_star2=False, final_model_star2=True)
OUTPUT_CO = dict(OUTPUT_ALL, history_star2=False, final_model_star2=False)


@dataclass
class Case:
    """One posydon-setup-grid run."""

    name: str
    system_type: str          # scenario system type / inlist folder
    columns: str              # key of GRID_COLUMNS
    single_star: bool = False
    zams: str = "both"        # 'both', 'star2' or 'none'
    use_scenario: bool = True
    branch: str = "main"
    output: dict = field(default_factory=lambda: dict(OUTPUT_ALL))
    grid_type: str = "fixed"
    submission: str = "slurm"
    job_array: bool = True
    extra_inlists: list = field(default_factory=list)  # (key, value) pairs
    formation_mesa_defaults: bool = False


GRID_TYPES = {
    "HMS-HMS": Case("HMS-HMS", "HMS-HMS", "binary"),
    "CO-HMS_RLO": Case("CO-HMS_RLO", "CO-H_star", "binary",
                       output=dict(OUTPUT_CO)),
    "CO-HeMS": Case("CO-HeMS", "CO-He_star", "binary_new_Z", zams="none",
                    output=dict(OUTPUT_CO), formation_mesa_defaults=True),
    "single_HMS": Case("single_HMS", "HMS-HMS", "single", single_star=True),
    "single_HeMS": Case("single_HeMS", "CO-He_star", "single_new_Z",
                        single_star=True),
    "HeMS-HMS": Case("HeMS-HMS", "HeMS-HMS", "binary_new_Z", zams="star2",
                     branch="jeff_new_grid", formation_mesa_defaults=True),
}


def _variant(base, name, **changes):
    values = dict(base.__dict__)
    values.update(changes, name=name)
    values["output"] = dict(values["output"])
    values["extra_inlists"] = list(values["extra_inlists"])
    return Case(**values)


def _build_cases():
    cases = []
    for label, base in GRID_TYPES.items():
        cases.append(_variant(base, f"{label}__scenario"))
        # the v2 single_HeMS grids use no ZAMS model for either star
        zams = "none" if label == "single_HeMS" else base.zams
        cases.append(_variant(base, f"{label}__no_scenario",
                              use_scenario=False, zams=zams))
    hms = GRID_TYPES["HMS-HMS"]
    cases.append(_variant(hms, "HMS-HMS__scenario__shell", submission="shell"))
    cases.append(_variant(hms, "HMS-HMS__scenario__slurm_mpi",
                          job_array=False))
    cases.append(_variant(hms, "HMS-HMS__scenario__dynamic",
                          grid_type="dynamic", job_array=False))
    cases.append(_variant(hms, "HMS-HMS__scenario__dynamic_job_array",
                          grid_type="dynamic", job_array=True))
    cases.append(_variant(
        hms, "HMS-HMS__scenario__user_inlists",
        extra_inlists=[
            ("binary_controls_user", "@DATA@/extra/user_inlist_project"),
            ("star1_controls_user", "@DATA@/extra/user_inlist1"),
            ("star1_job_user", "@DATA@/extra/user_inlist1")]))
    cases.append(_variant(
        GRID_TYPES["CO-HMS_RLO"], "CO-HMS_RLO__scenario__user_inlists",
        extra_inlists=[
            ("binary_controls_user", "@DATA@/extra/user_inlist_project")]))
    return cases


CASES = _build_cases()
CASE_IDS = [case.name for case in CASES]


###############################################################################
# environment
###############################################################################
def _git(repo, *args):
    env = dict(os.environ, **GIT_ENV)
    return subprocess.run(["git", "-C", str(repo), *args], env=env,
                          check=True, capture_output=True,
                          text=True).stdout.strip()


@pytest.fixture(scope="module")
def inlist_origin(tmp_path_factory):
    """A local git repo standing in for POSYDON-MESA-INLISTS."""
    if shutil.which("git") is None or shutil.which("make") is None:
        pytest.skip("git and make are needed for the setup-grid baseline")
    origin = tmp_path_factory.mktemp("inlist_origin") / "POSYDON-MESA-INLISTS"
    shutil.copytree(DATA / "inlists_main", origin)
    _git(origin, "init", "-q", "-b", "main")
    _git(origin, "add", "-A")
    _git(origin, "commit", "-q", "-m", "vendored r11701 snapshot")
    main_hash = _git(origin, "rev-parse", "HEAD")
    _git(origin, "checkout", "-q", "-b", "jeff_new_grid")
    shutil.copytree(DATA / "inlists_jeff_new_grid" / "r11701" / "HeMS-HMS",
                    origin / "r11701" / "HeMS-HMS")
    _git(origin, "add", "-A")
    _git(origin, "commit", "-q", "-m", "Adding HeMS-HMS grid")
    jeff_hash = _git(origin, "rev-parse", "HEAD")
    _git(origin, "checkout", "-q", "main")
    return origin, {"main": main_hash, "jeff_new_grid": jeff_hash}


@pytest.fixture
def setup_env(tmp_path, inlist_origin):
    """Paths and environment for one posydon-setup-grid run."""
    return make_setup_env(tmp_path, inlist_origin)


def make_setup_env(tmp_path, inlist_origin):
    """Build a temporary HOME, MESA_DIR, PATH and run directory."""
    origin, hashes = inlist_origin
    home = tmp_path / "home"
    home.mkdir()
    clone = home / ".posydon_mesa_inlists"
    subprocess.run(["git", "clone", "-q", str(origin), str(clone)],
                   check=True, capture_output=True)
    mesa = tmp_path / "mesa"
    shutil.copytree(DATA / "mesa_stub", mesa)
    stub_bin = tmp_path / "bin"
    stub_bin.mkdir()
    run_grid = stub_bin / "posydon-run-grid"
    run_grid.write_text("#!/bin/sh\necho 'stub posydon-run-grid'\n")
    run_grid.chmod(0o755)
    run = tmp_path / "run"
    run.mkdir()
    env = {
        "HOME": str(home),
        "PATH": os.pathsep.join([str(stub_bin), "/usr/bin", "/bin"]),
        "MESA_DIR": str(mesa),
        "MESASDK_ROOT": "/opt/mesasdk-stub",
        "PYTHONPATH": str(ROOT),
        "LANG": "C.UTF-8",
        # importing posydon imports matplotlib, which would otherwise create
        # a platform-dependent config folder in HOME
        "MPLCONFIGDIR": str(tmp_path / "mplconfig"),
    }
    return dict(tmp=tmp_path, home=home, clone=clone, mesa=mesa,
                bin=stub_bin, run=run, env=env, hashes=hashes)


###############################################################################
# ini construction
###############################################################################
def _explicit_inlists(case):
    """The keys the scenario would add, spelled out as paths (no scenario).

    This mirrors find_inlist_from_scenario in bin/posydon-setup-grid.
    """
    common = "@INLISTS@/r11701/default_common_inlists"
    folder = f"@INLISTS@/r11701/{case.system_type}"
    keys = [
        ("star1_controls_posydon_defaults", f"{common}/binary/inlist1"),
        ("star1_job_posydon_defaults", f"{common}/binary/inlist1"),
        ("star2_controls_posydon_defaults", f"{common}/binary/inlist2"),
        ("star2_job_posydon_defaults", f"{common}/binary/inlist2"),
        ("binary_controls_posydon_defaults",
         f"{common}/binary/inlist_project"),
        ("binary_job_posydon_defaults", f"{common}/binary/inlist_project"),
    ]
    steps = f"{folder}/star1_formation/inlist_step1, " \
            f"{folder}/star1_formation/inlist_step2"
    two_common = f"{common}/binary/inlist1, {common}/binary/inlist1"
    if case.system_type in ("CO-H_star",) and not case.single_star:
        keys += [
            ("binary_controls_user", f"{folder}/binary/inlist_project"),
            ("binary_job_user", f"{folder}/binary/inlist_project"),
            ("star1_controls_user", f"{folder}/binary/inlist1"),
            ("star1_job_user", f"{folder}/binary/inlist1"),
            ("star2_controls_user", f"{folder}/binary/inlist2"),
            ("star2_job_user", f"{folder}/binary/inlist2"),
        ]
    if case.system_type == "CO-He_star" and not case.single_star:
        keys += [
            ("binary_controls_user", f"{folder}/binary/inlist_project"),
            ("binary_job_user", f"{folder}/binary/inlist_project"),
            ("star1_controls_user", f"{folder}/binary/inlist1"),
            ("star1_job_user", f"{folder}/binary/inlist1"),
            ("star1_formation_controls_posydon_defaults", two_common),
            ("star1_formation_job_posydon_defaults", two_common),
            ("star1_formation_controls_user", steps),
            ("star1_formation_job_user", steps),
        ]
    if case.system_type == "HeMS-HMS":
        keys += [
            ("star1_controls_user", f"{folder}/binary/inlist1"),
            ("star1_job_user", f"{folder}/binary/inlist1"),
            ("star1_formation_controls_posydon_defaults", two_common),
            ("star1_formation_job_posydon_defaults", two_common),
            ("star1_formation_controls_user", steps),
            ("star1_formation_job_user", steps),
        ]
    if case.single_star and case.system_type == "HMS-HMS":
        keys.append(("star1_controls_special",
                     "@DATA@/extra/special_single_star_user_inlist"))
    if case.single_star and case.system_type == "CO-He_star":
        # as the v2 single_HeMS grids on yggdrasil: the two He-star formation
        # steps, then a third step with the single-star flag
        special = "@DATA@/extra/special_single_He_star_user_inlist"
        keys += [
            ("star1_controls_user", f"{steps}, {special}"),
            ("star1_job_user", f"{steps}, {special}"),
        ]
    return keys


def _explicit_extras():
    common = "@INLISTS@/r11701/default_common_inlists"
    return [
        ("posydon_binary_extras", f"{common}/binary/src/run_binary_extras.f"),
        ("posydon_star_binary_extras",
         f"{common}/binary/src/run_star_extras.f"),
    ]


def build_ini(case, env):
    """Return the ini text for a case."""
    defaults = "@DATA@/mesa_defaults"
    columns_root = "@INLISTS@/r11701/default_common_inlists"
    lines = ["[slurm]",
             "number_of_nodes=1", "number_of_mpi_tasks=1",
             "number_of_cpus_per_task=4",
             f"job_array={case.job_array}",
             "user=tester", "partition=test-cpu", "account=test",
             "walltime='2-00:00:00'", "email=tester@posydon.invalid",
             "", "[mesa_inlists]",
             f"posydon_github_root={ROOT}"]
    if case.use_scenario:
        commit = env["hashes"][case.branch]
        lines.append(f"scenario = ['posydon', '{case.branch}-{commit}', "
                     f"'{case.system_type}']")
    # setup only writes this path into the inlists, it never opens the file
    zams = "/posydon-tests/ZAMS_models/zams_z1.42m2_y0.2703.data"
    if case.zams == "both":
        lines.append(f"zams_filename = {zams}")
    elif case.zams == "star2":
        if case.use_scenario:
            # the scenario unsets zams_filename_1 itself
            lines.append(f"zams_filename = {zams}")
        else:
            lines.append(f"zams_filename_2 = {zams}")
    lines.append(f"single_star_grid = {case.single_star}")
    if case.formation_mesa_defaults:
        lines += [
            f"star1_formation_job_mesa_defaults = "
            f"{defaults}/star/star_job.defaults",
            f"star1_formation_controls_mesa_defaults = "
            f"{defaults}/star/controls.defaults"]
    lines += [
        f"binary_controls_mesa_defaults = "
        f"{defaults}/binary/binary_controls.defaults",
        f"binary_job_mesa_defaults = {defaults}/binary/binary_job.defaults",
        f"star1_job_mesa_defaults = {defaults}/star/star_job.defaults",
        f"star1_controls_mesa_defaults = {defaults}/star/controls.defaults",
        f"star2_job_mesa_defaults = {defaults}/star/star_job.defaults",
        f"star2_controls_mesa_defaults = {defaults}/star/controls.defaults",
    ]
    if not case.use_scenario:
        lines += [f"{k} = {v}" for k, v in _explicit_inlists(case)]
    lines += [f"{k} = {v}" for k, v in case.extra_inlists]
    lines += [
        f"star_history_columns = {columns_root}/history_columns.list",
        f"binary_history_columns = "
        f"{columns_root}/binary_history_columns.list",
        f"profile_columns = {columns_root}/profile_columns.list",
    ]
    lines += [f"{k} = {v}" for k, v in case.output.items()]
    star1_extras = ("@INLISTS@/r11701/default_common_inlists/binary/src/"
                    "run_star_extras.f" if not case.use_scenario
                    else "${MESA_DIR}/star/work/src/run_star_extras.f")
    lines += ["", "[mesa_extras]",
              "makefile_binary = ${MESA_DIR}/binary/work/make/makefile",
              "makefile_star = ${MESA_DIR}/star/work/make/makefile",
              "mesa_binary_extras = "
              "${MESA_DIR}/binary/work/src/run_binary_extras.f",
              "mesa_star_binary_extras = "
              "${MESA_DIR}/binary/work/src/run_star_extras.f",
              f"mesa_star1_extras = {star1_extras}",
              "mesa_star2_extras = "
              "${MESA_DIR}/star/work/src/run_star_extras.f",
              "binary_run = ${MESA_DIR}/binary/work/src/binary_run.f",
              "star_run = ${MESA_DIR}/star/work/src/run.f"]
    if not case.use_scenario:
        lines += [f"{k} = {v}" for k, v in _explicit_extras()]
    lines += ["", "[run_parameters]", "grid = @RUN@/grid.csv"]
    if case.grid_type == "dynamic":
        lines.append("psycris_inifile = @DATA@/extra/psycris_dynamic.ini")
    text = "\n".join(lines) + "\n"
    inlists = env["clone"] if case.use_scenario else DATA / "inlists_main"
    if case.system_type == "HeMS-HMS" and not case.use_scenario:
        inlists = None
    replacements = {"@RUN@": str(env["run"]), "@DATA@": str(DATA)}
    if inlists is not None:
        replacements["@INLISTS@"] = str(inlists)
    else:
        # HeMS-HMS without scenario: common inlists from main, the grid
        # type's folder from the jeff_new_grid snapshot
        text = text.replace("@INLISTS@/r11701/HeMS-HMS",
                            str(DATA / "inlists_jeff_new_grid" / "r11701"
                                / "HeMS-HMS"))
        replacements["@INLISTS@"] = str(DATA / "inlists_main")
    for key, value in replacements.items():
        text = text.replace(key, value)
    return text


def run_setup(case, env, script=SETUP_GRID):
    """Write the grid and ini, then run posydon-setup-grid."""
    run = env["run"]
    rows = GRID_ROWS[case.columns]
    (run / "grid.csv").write_text(
        "\n".join([GRID_COLUMNS[case.columns]] + rows) + "\n")
    (run / "grid.ini").write_text(build_ini(case, env))
    return subprocess.run(
        [sys.executable, str(script), "--inifile", "grid.ini",
         "--grid-type", case.grid_type, "--run-directory", str(run),
         "--submission-type", case.submission],
        cwd=run, env=env["env"], capture_output=True, text=True, timeout=300)


###############################################################################
# snapshot of the outputs
###############################################################################
def read_namelists(path):
    """Parse a MESA namelist file into {section: {key: value}}.

    Files without '&section' lines (MESA .defaults) go under None. Values are
    kept as raw strings; later keys overwrite earlier ones.
    """
    sections = {}
    current = None
    for raw in Path(path).read_text(errors="replace").splitlines():
        line = raw.split("!")[0].strip()
        if not line:
            continue
        if line.startswith("&"):
            current = line[1:].strip()
            sections.setdefault(current, {})
        elif "=" in line:
            key, value = line.split("=", 1)
            sections.setdefault(current, {})[key.strip()] = value.strip()
    return sections


_DEFAULTS_CACHE = {}


def mesa_defaults(section):
    """MESA default values for a namelist section."""
    files = {"controls": "star/controls.defaults",
             "star_job": "star/star_job.defaults",
             "binary_controls": "binary/binary_controls.defaults",
             "binary_job": "binary/binary_job.defaults"}
    if section not in _DEFAULTS_CACHE:
        if section in files:
            parsed = read_namelists(MESA_DEFAULTS / files[section])
            _DEFAULTS_CACHE[section] = parsed.get(None, {})
        else:
            _DEFAULTS_CACHE[section] = {}
    return _DEFAULTS_CACHE[section]


def _normalizer(env):
    pairs = []
    for key, label in [("run", "<RUN>"), ("clone", "<INLIST_CLONE>"),
                       ("home", "<HOME>"), ("mesa", "<MESA_DIR>"),
                       ("bin", "<BIN>")]:
        path = env[key]
        pairs += [(str(path.resolve()), label), (str(path), label)]
    pairs += [(str(DATA), "<DATA>"), (str(ROOT), "<POSYDON>"),
              (sys.executable, "<PYTHON>")]
    for name, commit in env["hashes"].items():
        pairs.append((commit, f"<{name.upper()}_HASH>"))
    pairs.sort(key=lambda pair: -len(pair[0]))

    def normalize(text):
        for old, new in pairs:
            text = text.replace(old, new)
        return text
    return normalize


def _sha(text):
    return hashlib.sha256(text.encode()).hexdigest()


def snapshot(result, env):
    """A readable, path-independent summary of a setup run."""
    norm = _normalizer(env)
    run = env["run"]
    out = [f"returncode: {result.returncode}"]
    if result.returncode != 0:
        last = [line for line in result.stderr.splitlines() if line.strip()]
        out.append("error: " + norm(last[-1] if last else ""))
    inlists = sorted(
        [p for p in run.glob("binary/inlist*")]
        + [p for p in run.glob("star[12]/inlist_step*")])
    for path in inlists:
        rel = path.relative_to(run)
        text = norm(path.read_text(errors="replace"))
        out.append(f"\n## {rel}")
        out.append(f"sha256: {_sha(text)}")
        for section, values in read_namelists(path).items():
            defaults = mesa_defaults(section)
            changed = [(k, v) for k, v in values.items()
                       if defaults.get(k) != v]
            out.append(f"&{section}  ({len(values)} keys, "
                       f"{len(changed)} not MESA default)")
            out += [f"  {k} = {norm(v)}" for k, v in changed]
    for path in sorted(run.glob("column_lists/*")):
        text = path.read_text(errors="replace")
        out.append(f"\n## {path.relative_to(run)}")
        out.append(f"sha256: {_sha(text)}")
        out.append(f"lines: {len(text.splitlines())}")
    out.append("\n## staged sources")
    for path in sorted(run.glob("*/src/*")) + sorted(run.glob("*/make/*")):
        first = path.read_text(errors="replace").splitlines()[:1]
        label = norm(first[0]) if first else ""
        out.append(f"{path.relative_to(run)}: {label}")
    for name in ["mk", "grid_command.sh", "job_array_grid_submit.slurm",
                 "mpi_grid_submit.slurm", "cleanup.slurm", "run_grid.sh"]:
        path = run / name
        if path.exists():
            out.append(f"\n## {name}")
            out.append(norm(path.read_text()).rstrip())
    out.append("\n## run directory tree")
    out += tree(run, norm)
    out.append("\n## created outside the run directory")
    expected = {"home", "mesa", "bin", "run", "mplconfig"}
    outside = [f"<TMP>/{p.name}" for p in sorted(env["tmp"].iterdir())
               if p.name not in expected]
    outside += [f"<HOME>/{p.name}" for p in sorted(env["home"].iterdir())
                if p.name != ".posydon_mesa_inlists"]
    mesa_now = tree(env["mesa"], norm)
    if mesa_now != tree(DATA / "mesa_stub", norm):
        outside.append("<MESA_DIR> changed:")
        outside += mesa_now
    out += outside or ["(none)"]
    status = subprocess.run(["git", "-C", str(env["clone"]), "status",
                             "--porcelain"], capture_output=True, text=True)
    out.append("\n## untracked or changed files in the inlist clone")
    out += status.stdout.splitlines() or ["(none)"]
    return "\n".join(out) + "\n"


def tree(root, norm=str):
    """Every path under root, sorted, with its kind.

    d = directory, f = file, x = executable file, l = symlink (with target).
    """
    lines = []
    for path in sorted(Path(root).rglob("*")):
        rel = path.relative_to(root).as_posix()
        if path.is_symlink():
            lines.append(f"l {rel} -> {norm(os.readlink(path))}")
        elif path.is_dir():
            lines.append(f"d {rel}/")
        elif os.access(path, os.X_OK):
            lines.append(f"x {rel}")
        else:
            lines.append(f"f {rel}")
    return lines or ["(empty)"]


def _semantic(text):
    return "\n".join(line for line in text.splitlines()
                     if not line.startswith("sha256: "))


def _run_case(case, setup_env):
    result = run_setup(case, setup_env)
    return result, snapshot(result, setup_env)


###############################################################################
# tests
###############################################################################
@pytest.mark.parametrize("case", CASES, ids=CASE_IDS)
def test_matches_golden(case, setup_env):
    """Today's output for each case is frozen in a golden file."""
    _, snap = _run_case(case, setup_env)
    golden = GOLDEN / f"{case.name}.txt"
    if UPDATE_GOLDEN:
        golden.parent.mkdir(parents=True, exist_ok=True)
        golden.write_text(snap)
    assert golden.exists(), \
        f"no golden file for {case.name}; run with POSYDON_UPDATE_GOLDEN=1"
    expected = golden.read_text()
    if SEMANTIC:
        assert _semantic(snap) == _semantic(expected)
    else:
        assert snap == expected


def _tree_part(snap):
    """The folder-structure sections of a snapshot."""
    keep = []
    for block in snap.split("\n## ")[1:]:
        title = block.splitlines()[0]
        if title in ("run directory tree", "created outside the run "
                     "directory", "untracked or changed files in the "
                     "inlist clone"):
            keep.append(block.rstrip())
    return keep


@pytest.mark.parametrize("case", CASES, ids=CASE_IDS)
def test_folder_structure(case, setup_env):
    """The exact files and folders setup creates, and nothing elsewhere."""
    golden = GOLDEN / f"{case.name}.txt"
    if not golden.exists():
        pytest.skip("golden file missing; see test_matches_golden")
    _, snap = _run_case(case, setup_env)
    assert _tree_part(snap) == _tree_part(golden.read_text())


def _inlist_part(snap):
    """The inlist, column-list and source sections of a snapshot."""
    keep = []
    for block in snap.split("\n## ")[1:]:
        title = block.splitlines()[0]
        if title.startswith(("binary/", "star1/", "star2/", "column_lists/",
                             "staged sources")):
            keep.append(block)
    return "\n## ".join(keep)


# fixed in step 6: HeMS-HMS (S2)
_EQUIVALENCE_XFAIL = {}


@pytest.mark.parametrize("label", [
    pytest.param(label, marks=pytest.mark.xfail(
        strict=True, reason=_EQUIVALENCE_XFAIL[label]))
    if label in _EQUIVALENCE_XFAIL else label
    for label in GRID_TYPES])
def test_scenario_equals_explicit_paths(label, tmp_path, inlist_origin):
    """A scenario must give the same inlists as spelling its files out."""
    base = GRID_TYPES[label]
    snaps = {}
    for use_scenario in (True, False):
        sub = tmp_path / ("scenario" if use_scenario else "explicit")
        sub.mkdir()
        env = make_setup_env(sub, inlist_origin)
        zams = base.zams
        if label == "single_HeMS" and not use_scenario:
            # the scenario forms the He star without a ZAMS model and ignores
            # zams_filename (S3); spelled out, that is no ZAMS model at all
            zams = "none"
        case = _variant(base, base.name, use_scenario=use_scenario,
                        zams=zams)
        result, snap = _run_case(case, env)
        snaps[use_scenario] = (result.returncode, _inlist_part(snap))
    assert snaps[True] == snaps[False]


###############################################################################
# opt-in replay of real grids
###############################################################################
REFERENCE = os.environ.get("POSYDON_YGGDRASIL_REFERENCE", "")

# Run folders (<type>/<Z>/<run>) whose inlists were edited by hand after
# setup, so no ini reproduces them; not checked by the replay.
HAND_EDITED_RUNS = {
    # relax_initial_mass / new_mass set for a single-mass rerun
    "HMS-HMS/1e+00_Zsun/grid_low_res_rerun_LBV_wind_m2_283",
    "HMS-HMS/1e-03_Zsun/grid_low_res_rerun_113Msun",
    "HMS-HMS/1e-03_Zsun/grid_low_res_rerun_127Msun",
    "HMS-HMS/1e-03_Zsun/grid_low_res_rerun_142Msun",
    "HMS-HMS/1e-03_Zsun/grid_low_res_rerun_160Msun",
    "HMS-HMS/1e-03_Zsun/grid_low_res_rerun_202Msun",
    "HMS-HMS/1e-03_Zsun/grid_low_res_rerun_80Msun",
    "HMS-HMS/1e-03_Zsun/rerun_PISN_grid_low_res_rerun_202Msun",
    "HMS-HMS/2e+00_Zsun/grid_low_res_rerun_m1_227",
    "HMS-HMS/2e+00_Zsun/grid_low_res_rerun_m1_255",
    "HMS-HMS/2e+00_Zsun/grid_low_res_rerun_m1_286",
    "HMS-HMS/2e+00_Zsun/grid_low_res_rerun_m2_214",
    "HMS-HMS/2e+00_Zsun/grid_low_res_rerun_m2_229_q080",
    "HMS-HMS/2e+00_Zsun/grid_low_res_rerun_m2_243",
    "HMS-HMS/2e+00_Zsun/grid_low_res_rerun_m2_272",
    # other keys changed by hand (mesh, surface, max_tries1, do_jdot_ls, ...)
    "HMS-HMS/1e+00_Zsun/grid_low_res_0",
    "HMS-HMS/1e+00_Zsun/rerun_thermohaline_mixing_grid_random_combined",
    "HMS-HMS/1e-02_Zsun/grid_low_res_rerun_LBV_wind_fr_limit",
    "HMS-HMS/1e-02_Zsun/rerun_thermohaline_mixing_grid_random_combined",
    "HMS-HMS/1e-04_Zsun/rerun_thermohaline_mixing_grid_random_combined",
    "HMS-HMS/2e+00_Zsun/grid_low_res_rerun_LBV_wind_initsync",
    "HMS-HMS/2e+00_Zsun/grid_low_res_rerun_LBV_wind_initsync_extras",
}
# grid types whose setup is a pinned xfail bug today
KNOWN_FAILING_TYPES = set()


def _path_name(value):
    """Compare absolute paths by file name only (cluster roots differ)."""
    match = re.match(r"^'(.*/)?([^/']+)'$", value)
    return f"'{match.group(2)}'" if match and match.group(1) else value


@pytest.mark.skipif(not REFERENCE, reason="set POSYDON_YGGDRASIL_REFERENCE to "
                    "a local copy of yggdrasil run directories to replay them")
def test_replay_yggdrasil_runs(tmp_path):
    """Today's setup reproduces the inlists of real yggdrasil runs.

    POSYDON_YGGDRASIL_REFERENCE/grids/<type>/<Z>/<run>/ must hold the run's
    ini, its grid CSV and generated binary/ and star*/ inlists. The scenario
    refs are read from a copy of the local ~/.posydon_mesa_inlists clone.
    Skipped: hand-edited runs, known-failing grid types, ini files whose grid
    file is missing or empty (e.g. an unrelated ini copied into the folder),
    and scenario commits missing from the local clone.
    """
    source_clone = Path.home() / ".posydon_mesa_inlists"
    if not source_clone.is_dir():
        pytest.skip("needs a local ~/.posydon_mesa_inlists clone")
    home = tmp_path / "home"
    shutil.copytree(source_clone, home / ".posydon_mesa_inlists",
                    symlinks=True)
    subprocess.run(["git", "-C", str(home / ".posydon_mesa_inlists"),
                    "remote", "set-url", "origin", "file:///nonexistent"],
                   check=True)
    # the ini files' grid_params/ paths are read from POSYDON-MESA-INLISTS
    # at the pinned commit; without network, extract it from a local clone
    inlist_repo = os.environ.get("POSYDON_INLIST_REPO", "")
    if not inlist_repo:
        pytest.skip("needs POSYDON_INLIST_REPO for the grid_params/ files")
    from posydon.grids.setup.support import support_root
    support_root(home=home, clone=inlist_repo)
    mesa = tmp_path / "mesa"
    shutil.copytree(DATA / "mesa_stub", mesa)
    stub_bin = tmp_path / "bin"
    stub_bin.mkdir()
    (stub_bin / "posydon-run-grid").write_text("#!/bin/sh\n")
    (stub_bin / "posydon-run-grid").chmod(0o755)
    env = {"HOME": str(home), "PATH": f"{stub_bin}:/usr/bin:/bin",
           "MESA_DIR": str(mesa), "MESASDK_ROOT": "/opt/mesasdk-stub",
           "PYTHONPATH": str(ROOT), "LANG": "C.UTF-8",
           "MPLCONFIGDIR": str(tmp_path / "mplconfig")}
    failures, replayed = [], 0
    clone = home / ".posydon_mesa_inlists"
    for ini in sorted(Path(REFERENCE).glob("grids/*/*/*/*.ini")):
        ref = ini.parent
        rel = ref.relative_to(Path(REFERENCE) / "grids").as_posix()
        if not (ref / "binary" / "inlist1").exists() \
                or rel in HAND_EDITED_RUNS \
                or rel.split("/")[0] in KNOWN_FAILING_TYPES:
            continue
        text = ini.read_text()
        match = re.search(r"(?m)^grid\s*=\s*(\S+)", text)
        grid = ref / Path(match.group(1) if match else "grid.csv").name
        if not grid.exists() or grid.stat().st_size == 0:
            continue
        scenario = r"(?m)^scenario\s*=\s*\[[^,]*,\s*'[^']*-([0-9a-f]+)'"
        match = re.search(scenario, text)
        if match and subprocess.run(
                ["git", "-C", str(clone), "cat-file", "-e",
                 match.group(1) + "^{commit}"],
                capture_output=True).returncode:
            continue
        run = tmp_path / "runs" / rel / ini.stem
        run.mkdir(parents=True)
        shutil.copy(grid, run / "grid.csv")
        text = re.sub(r"(?m)^posydon_github_root\s*=.*$",
                      f"posydon_github_root={ROOT}", text)
        text = text.replace("POSYDON-MESA-INLISTS-public",
                            "POSYDON-MESA-INLISTS")
        text = re.sub(r"(?m)^grid\s*=.*$", f"grid = {run / 'grid.csv'}",
                      text)
        (run / "grid.ini").write_text(text)
        result = subprocess.run(
            [sys.executable, str(SETUP_GRID), "--inifile", "grid.ini",
             "--grid-type", "fixed", "--run-directory", str(run),
             "--submission-type", "slurm"], cwd=run, env=env,
            capture_output=True, text=True, timeout=600)
        replayed += 1
        if result.returncode != 0:
            failures.append(f"{ini}: setup failed")
            continue
        steps = "star[12]/inlist_step*"
        names = {p.relative_to(ref) for p in ref.glob("binary/inlist*")}
        names |= {p.relative_to(ref) for p in ref.glob(steps)}
        names |= {p.relative_to(run) for p in run.glob(steps)}
        for name in sorted(names):
            if not (ref / name).exists() or not (run / name).exists():
                failures.append(f"{ini}: {name} missing on one side")
                continue
            a, b = read_namelists(ref / name), read_namelists(run / name)
            for section in set(a) | set(b):
                sa, sb = a.get(section, {}), b.get(section, {})
                for key in set(sa) | set(sb):
                    va, vb = sa.get(key), sb.get(key)
                    if va is None or vb is None or \
                            _path_name(va) != _path_name(vb):
                        failures.append(f"{ini}: {name} &{section} {key}: "
                                        f"cluster={va} replay={vb}")
    assert replayed, "no replayable run directories found"
    assert not failures, "\n".join(failures[:50])


###############################################################################
# known bugs: the correct behaviour, expected to fail until fixed
###############################################################################
def test_single_HeMS_scenario_runs(setup_env):
    """As the v2 grids: 2 He-star formation steps plus a flagged third.

    Fixed in stage 3.6 (tracker S1); before, the scenario appended to
    mesa_inlists['star1_controls_special'], which was never created.
    """
    result = run_setup(_variant(GRID_TYPES["single_HeMS"], "single_HeMS"),
                       setup_env)
    assert result.returncode == 0, result.stderr[-2000:]
    steps = sorted((setup_env["run"] / "star1").glob("inlist_step*"))
    assert [p.name for p in steps] == ["inlist_step0", "inlist_step1",
                                       "inlist_step2"]
    last = read_namelists(steps[-1])["controls"]
    assert last.get("x_logical_ctrl(1)") == ".true."


def test_single_HeMS_scenario_uses_no_zams_model(setup_env):
    """As the v2 grids: no zams_filename anywhere. Fixed in stage 3.6 (S3)."""
    result = run_setup(_variant(GRID_TYPES["single_HeMS"], "single_HeMS"),
                       setup_env)
    assert result.returncode == 0, result.stderr[-2000:]
    steps = sorted((setup_env["run"] / "star1").glob("inlist_step*"))
    # the He star is relaxed from a model, not read from a ZAMS file
    for step in steps:
        assert "zams_filename" not in read_namelists(step)["controls"], step
    for name in ("inlist1", "inlist2"):
        assert "zams_filename" not in read_namelists(
            setup_env["run"] / "binary" / name)["controls"], name


def test_single_HMS_scenario_leaves_inlist_clone_clean(setup_env):
    """Fixed in stage 3.6 (tracker S7)."""
    run_setup(_variant(GRID_TYPES["single_HMS"], "single_HMS"), setup_env)
    status = subprocess.run(["git", "-C", str(setup_env["clone"]), "status",
                             "--porcelain"], capture_output=True, text=True)
    assert status.stdout.strip() == ""


def test_user_inlist_survives_scenario(setup_env):
    """Fixed in stage 3.6 (tracker S6): user_inlists above the scenario."""
    case = _variant(GRID_TYPES["CO-HMS_RLO"], "CO-HMS_RLO", extra_inlists=[
        ("binary_controls_user", "@DATA@/extra/user_inlist_project")])
    result = run_setup(case, setup_env)
    assert result.returncode == 0, result.stderr[-2000:]
    values = read_namelists(setup_env["run"] / "binary" / "inlist_project")
    assert values["binary_controls"].get("fm") == "0.002d0"


def test_single_HeMS_explicit_drops_zams_filename(setup_env):
    """Fixed in stage 3.6 (tracker S3); before, the branch meant to drop it
    checked for a key 'zams_filename_1' (legacy_posydon_setup_grid.py:608).
    """
    case = _variant(GRID_TYPES["single_HeMS"], "single_HeMS",
                    use_scenario=False, zams="none")
    result = run_setup(case, setup_env)
    assert result.returncode == 0, result.stderr[-2000:]
    for step in sorted((setup_env["run"] / "star1").glob("inlist_step*")):
        assert "zams_filename" not in read_namelists(step)["controls"], step


def test_HeMS_HMS_scenario_forms_star1(setup_env):
    """Fixed in step 6 (tracker S2); before, construct_static_inlist copied
    zams_filename into zams_filename_1 after the scenario set it to None
    (legacy_posydon_setup_grid.py:229 vs :395)."""
    result = run_setup(_variant(GRID_TYPES["HeMS-HMS"], "HeMS-HMS"),
                       setup_env)
    assert result.returncode == 0, result.stderr[-2000:]
    assert list((setup_env["run"] / "star1").glob("inlist_step*"))


@pytest.mark.xfail(strict=True, reason="a fixed grid without a job array "
                   "is refused: posydon-run-grid would run no points (R1); "
                   "before stage 3.6 it read run_parameters['psycris_inifile'] "
                   "(legacy_posydon_setup_grid.py:1125)")
def test_fixed_mpi_grid_needs_no_psycris_ini(setup_env):
    case = _variant(GRID_TYPES["HMS-HMS"], "HMS-HMS", job_array=False)
    result = run_setup(case, setup_env)
    assert result.returncode == 0, result.stderr[-2000:]


@pytest.mark.xfail(strict=True, reason="with job_array=True the command line "
                   "is built with grid type 'fixed', whatever --grid-type "
                   "says (legacy_posydon_setup_grid.py:1088-1107)")
def test_dynamic_grid_keeps_dynamic_with_job_array(setup_env):
    case = _variant(GRID_TYPES["HMS-HMS"], "HMS-HMS", grid_type="dynamic",
                    job_array=True)
    result = run_setup(case, setup_env)
    assert result.returncode == 0, result.stderr[-2000:]
    script = (setup_env["run"] / "job_array_grid_submit.slurm").read_text()
    assert "--grid-type dynamic" in script


def test_failed_setup_leaves_no_partial_run_directory(setup_env):
    """Fixed in stage 3.6 (tracker S10): all checks run before writing."""
    case = _variant(GRID_TYPES["HMS-HMS"], "HMS-HMS", job_array=False)
    result = run_setup(case, setup_env)
    if result.returncode == 0:
        pytest.skip("setup no longer fails for this case")
    created = {p.name for p in setup_env["run"].iterdir()}
    assert created <= {"grid.csv", "grid.ini"}, sorted(created)


@pytest.mark.xfail(strict=True, reason="the formation step's own src/ "
                   "(star1_formation/src/run_star_extras.f) is never staged")
def test_formation_extras_are_staged(setup_env):
    result = run_setup(_variant(GRID_TYPES["CO-HeMS"], "CO-HeMS"), setup_env)
    assert result.returncode == 0, result.stderr[-2000:]
    staged = (setup_env["run"] / "star1" / "src" / "run_star_extras.f")
    assert "CO-He_star/star1_formation/src" in staged.read_text()


def test_outputs_go_to_run_directory(setup_env):
    """Fixed in stage 3.6 (tracker S9)."""
    other = setup_env["tmp"] / "elsewhere"
    other.mkdir()
    case = _variant(GRID_TYPES["HMS-HMS"], "HMS-HMS")
    run = setup_env["run"]
    (run / "grid.csv").write_text(
        "\n".join([GRID_COLUMNS[case.columns]] + GRID_ROWS[case.columns])
        + "\n")
    (run / "grid.ini").write_text(build_ini(case, setup_env))
    subprocess.run([sys.executable, str(SETUP_GRID), "--inifile",
                    str(run / "grid.ini"), "--grid-type", "fixed",
                    "--run-directory", str(run), "--submission-type",
                    "slurm"], cwd=other, env=setup_env["env"],
                   capture_output=True, text=True, timeout=300)
    assert (run / "job_array_grid_submit.slurm").exists()
    assert not (other / "job_array_grid_submit.slurm").exists()


@pytest.mark.xfail(strict=True, reason="a shell script without a job array "
                   "has no newline between the command and compress-mesa "
                   "(legacy_posydon_setup_grid.py:1149-1153)")
def test_shell_script_without_job_array_ends_the_command(setup_env):
    case = _variant(GRID_TYPES["HMS-HMS"], "HMS-HMS", submission="shell",
                    grid_type="dynamic", job_array=False)
    result = run_setup(case, setup_env)
    assert result.returncode == 0, result.stderr[-2000:]
    script = (setup_env["run"] / "grid_command.sh").read_text()
    assert "\ncompress-mesa ." in script
