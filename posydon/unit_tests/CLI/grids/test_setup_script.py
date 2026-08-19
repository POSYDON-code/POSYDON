"""Unit tests locking the current behavior of bin/posydon-setup-grid."""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

import importlib.machinery
import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[4]


def _load_script(name, path):
    # The bin scripts have no .py extension, so spec_from_file_location returns
    # None; use an explicit SourceFileLoader instead.
    loader = importlib.machinery.SourceFileLoader(name, str(path))
    spec = importlib.util.spec_from_loader(name, loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    return mod


totest = _load_script("setup_grid_script", REPO_ROOT / "bin" / "posydon-setup-grid")


@pytest.fixture(autouse=True)
def _env(tmp_path, monkeypatch):
    """Always provide HOME/MESA_DIR and restore cwd after every test."""
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("MESA_DIR", str(tmp_path / "mesa"))
    original = os.getcwd()
    yield
    os.chdir(original)


@pytest.fixture
def inlist_files(tmp_path):
    """Write small MESA/posydon/user inlists for merge tests."""
    d = tmp_path / "inlists"
    d.mkdir()

    binary_controls = d / "binary_controls.inlist"
    binary_controls.write_text("&binary_controls\n  mdot_scheme = 'Kolb'\n/ ! end\n")

    binary_job = d / "binary_job.inlist"
    binary_job.write_text("&binary_job\n  evolve_both_stars = .true.\n/ ! end\n")

    star1_controls = d / "star1_controls.inlist"
    star1_controls.write_text(
        "&controls\n  m1 = 10.0\n  initial_mass = 10.0\n  num_x_ctrls = 3\n/ ! end\n"
    )

    star1_job = d / "star1_job.inlist"
    star1_job.write_text("&star_job\n  load_saved_model = .false.\n/ ! end\n")

    star2_controls = d / "star2_controls.inlist"
    star2_controls.write_text("&controls\n  m2 = 8.0\n/ ! end\n")

    star2_job = d / "star2_job.inlist"
    star2_job.write_text("&star_job\n/ ! end\n")

    return {
        "binary_controls": str(binary_controls),
        "binary_job": str(binary_job),
        "star1_controls": str(star1_controls),
        "star1_job": str(star1_job),
        "star2_controls": str(star2_controls),
        "star2_job": str(star2_job),
    }


def _binary_mesa_inlists(files):
    """Standard mesa_inlists dictionary for a binary-only grid."""
    return {
        "single_star_grid": False,
        "zams_filename_1": None,
        "zams_filename_2": None,
        "binary_controls_user": files["binary_controls"],
        "binary_job_user": files["binary_job"],
        "star1_controls_user": files["star1_controls"],
        "star1_job_user": files["star1_job"],
        "star2_controls_user": files["star2_controls"],
        "star2_job_user": files["star2_job"],
        "final_profile_star1": False,
        "final_profile_star2": False,
        "final_model_star1": False,
        "final_model_star2": False,
        "history_star1": True,
        "history_star2": True,
        "history_interval": 100,
        "binary_history": True,
    }


def _set_globals(mesa_inlists, mesa_extras):
    """Explicitly set the module globals mutated by find_inlist_from_scenario."""
    totest.mesa_inlists = mesa_inlists
    totest.mesa_extras = mesa_extras


class TestParseCommandline:
    """Tests for parse_commandline()."""

    @pytest.fixture(autouse=True)
    def _reset_argv(self, monkeypatch):
        monkeypatch.setattr(sys, "argv", ["posydon-setup-grid"])

    def test_valid_fixed_shell(self, monkeypatch):
        """Valid fixed grid with shell submission parses successfully."""
        monkeypatch.setattr(
            sys,
            "argv",
            ["posydon-setup-grid", "--inifile", "in.ini", "--grid-type", "fixed"],
        )
        args = totest.parse_commandline()
        assert args.inifile == "in.ini"
        assert args.grid_type == "fixed"

    def test_valid_dynamic_slurm(self, monkeypatch):
        """Valid dynamic grid with slurm submission parses successfully."""
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "posydon-setup-grid",
                "--inifile",
                "in.ini",
                "--grid-type",
                "dynamic",
                "--submission-type",
                "slurm",
                "-n",
                "4",
                "--run-directory",
                "/some/where",
                "--verbose",
            ],
        )
        args = totest.parse_commandline()
        assert args.grid_type == "dynamic"
        assert args.submission_type == "slurm"
        assert args.nproc == 4
        assert args.run_directory == "/some/where"
        assert args.verbose is True

    def test_invalid_grid_type_raises(self, monkeypatch):
        """An unknown grid type triggers parser.error (SystemExit)."""
        monkeypatch.setattr(
            sys,
            "argv",
            ["posydon-setup-grid", "--inifile", "in.ini", "--grid-type", "bogus"],
        )
        with pytest.raises(SystemExit):
            totest.parse_commandline()

    def test_invalid_submission_type_raises(self, monkeypatch):
        """An unknown submission type triggers parser.error (SystemExit)."""
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "posydon-setup-grid",
                "--inifile",
                "in.ini",
                "--grid-type",
                "fixed",
                "--submission-type",
                "bogus",
            ],
        )
        with pytest.raises(SystemExit):
            totest.parse_commandline()

    def test_defaults(self, monkeypatch, tmp_path):
        """run_directory defaults to cwd, submission_type to shell, nproc to 1."""
        os.chdir(tmp_path)
        monkeypatch.setattr(
            sys,
            "argv",
            ["posydon-setup-grid", "--inifile", "in.ini", "--grid-type", "fixed"],
        )
        args = totest.parse_commandline()
        assert args.run_directory == os.getcwd()
        assert args.submission_type == "shell"
        assert args.nproc == 1
        assert args.verbose is False


class TestFindInlistFromScenario:
    """Tests for find_inlist_from_scenario() with mocked git subprocesses."""

    @pytest.fixture
    def mock_git(self, monkeypatch):
        """Mock subprocess.Popen/call and record every call in order."""
        calls = []

        class FakePopen:
            def __init__(self, cmd, **kwargs):
                calls.append(("Popen", cmd))
                self.cmd = cmd

            def communicate(self):
                # Simulate the clone actually creating the destination repo.
                if self.cmd and self.cmd[1] == "clone":
                    os.makedirs(self.cmd[-1], exist_ok=True)
                return b"output", b"error"

            def wait(self):
                return 0

        def fake_popen(cmd, **kwargs):
            return FakePopen(cmd, **kwargs)

        def fake_call(cmd, **kwargs):
            calls.append(("call", cmd))
            return 0

        monkeypatch.setattr(subprocess, "Popen", fake_popen)
        monkeypatch.setattr(subprocess, "call", fake_call)
        return calls

    def _posydon_repo(self, home, system_type=None):
        """Create the ~/.posydon_mesa_inlists clone layout for a scenario."""
        repo = Path(home) / ".posydon_mesa_inlists"
        if not repo.is_dir():
            repo.mkdir()
        return repo

    def test_posydon_source_repo_present_only_checkout_pull(self, mock_git, tmp_path):
        """Repo already cloned: clone is skipped, checkout/pull order preserved."""
        repo = self._posydon_repo(tmp_path)
        _set_globals({"single_star_grid": False}, {})
        totest.find_inlist_from_scenario("posydon", "master-abc123", "HeMS_HeMS")

        commands = [kind for kind, _ in mock_git]
        assert commands == ["Popen", "call", "Popen"]
        assert mock_git[0][1] == ["git", "checkout", "master"]
        assert mock_git[1][1] == ["git", "pull"]
        assert mock_git[2][1] == ["git", "checkout", "abc123"]

    def test_posydon_source_repo_absent_clone(self, mock_git, tmp_path):
        """Repo missing: clone runs first, then checkout/pull/checkout order."""
        _set_globals({"single_star_grid": False}, {})
        totest.find_inlist_from_scenario("posydon", "master-abc123", "HeMS_HeMS")

        commands = [kind for kind, _ in mock_git]
        assert commands == ["Popen", "Popen", "call", "Popen"]
        clone = mock_git[0][1]
        assert clone[0] == "git"
        assert clone[1] == "clone"
        assert clone[2] == "https://github.com/POSYDON-code/POSYDON-MESA-INLISTS.git"
        assert clone[3] == os.path.join(tmp_path, ".posydon_mesa_inlists")
        assert mock_git[1][1] == ["git", "checkout", "master"]
        assert mock_git[2][1] == ["git", "pull"]
        assert mock_git[3][1] == ["git", "checkout", "abc123"]

    def test_user_source_invalid_gitcommit(self, mock_git, tmp_path):
        """Invalid user gitcommit (no '-') raises ValueError with exact message."""
        _set_globals({}, {})
        with pytest.raises(
            ValueError,
            match="You have supplied an invalid user gitcommit format, must be of format 'branch-githash'",
        ):
            totest.find_inlist_from_scenario("user", "notagithash", "HeMS-HMS")

    def test_invalid_source(self, mock_git, tmp_path):
        """Unknown source raises ValueError with the exact message."""
        _set_globals({}, {})
        with pytest.raises(
            ValueError,
            match="supplied source is not valid/understood. Valid sources are user and posydon",
        ):
            totest.find_inlist_from_scenario("other", "branch-hash", "HeMS-HMS")

    def test_zams_filename_detected_in_inlist1(self, mock_git, tmp_path):
        """A zams_filename line in binary/inlist1 sets zams_filename_1."""
        repo = self._posydon_repo(tmp_path)
        common = repo / "r11701" / "default_common_inlists"
        (common / "binary").mkdir(parents=True)
        (common / "binary" / "inlist1").write_text(
            "&controls\n  zams_filename = 'zams1.data'\n/ ! end\n"
        )
        zams = repo / "r11701" / "ZAMS_models" / "zams1.data"
        zams.parent.mkdir(parents=True)
        zams.write_text("zams data")

        _set_globals({"single_star_grid": False}, {})
        totest.find_inlist_from_scenario("posydon", "master-abc123", "HeMS_HeMS")

        assert totest.mesa_inlists["zams_filename_1"] == str(zams)
        assert "zams_filename_2" not in totest.mesa_inlists

    def test_hems_hms_sets_star1_formation_and_unsets_zams(self, mock_git, tmp_path):
        """HeMS-HMS: star1 formation arrays populated and zams_filename_1 set None."""
        repo = self._posydon_repo(tmp_path)
        common = repo / "r11701" / "default_common_inlists"
        (common / "binary").mkdir(parents=True)
        (common / "binary" / "inlist1").write_text(
            "&controls\n  zams_filename = 'zams1.data'\n/ ! end\n"
        )
        (repo / "r11701" / "ZAMS_models").mkdir(parents=True)
        (repo / "r11701" / "ZAMS_models" / "zams1.data").write_text("zams data")

        system = repo / "r11701" / "HeMS-HMS"
        (system / "star1_formation").mkdir(parents=True)
        step0 = system / "star1_formation" / "step0"
        step1 = system / "star1_formation" / "step1"
        step0.write_text("step0")
        step1.write_text("step1")

        _set_globals({"single_star_grid": False}, {})
        totest.find_inlist_from_scenario("posydon", "master-abc123", "HeMS-HMS")

        # the script layers '{inlists_location_common}/binary/inlist1' where
        # inlists_location_common has a trailing slash, giving a double slash
        common_inlist1 = "{0}/r11701/default_common_inlists//binary/inlist1".format(
            repo
        )
        assert totest.mesa_inlists["zams_filename_1"] is None
        assert totest.mesa_inlists["star1_formation_controls_user"] == [
            str(step0),
            str(step1),
        ]
        assert totest.mesa_inlists["star1_formation_job_user"] == [
            str(step0),
            str(step1),
        ]
        assert totest.mesa_inlists["star1_formation_controls_posydon_defaults"] == [
            common_inlist1,
            common_inlist1,
        ]
        assert totest.mesa_inlists["star1_formation_job_posydon_defaults"] == [
            common_inlist1,
            common_inlist1,
        ]

    def test_hms_hms_single_star_special_inlist(self, mock_git, tmp_path):
        """HMS-HMS single star: special inlist written with x_logical_ctrl(1)=.true."""
        repo = self._posydon_repo(tmp_path)
        _set_globals({"single_star_grid": True}, {})
        totest.find_inlist_from_scenario("posydon", "master-abc123", "HMS-HMS")

        special = repo / "special_single_star_user_inlist"
        assert totest.mesa_inlists["star1_controls_special"] == str(special)
        content = special.read_bytes()
        assert b"&controls\n\n" in content
        assert b"\tx_logical_ctrl(1) = .true.\n" in content
        assert b"/ ! end of star1_controls namelist" in content
        assert b"&star_job" not in content

    def test_co_he_star_single_star_special_inlist(self, mock_git, tmp_path):
        """CO-He_star single star: special inlist gets an extra &star_job section."""
        repo = self._posydon_repo(tmp_path)
        system = repo / "r11701" / "CO-He_star"
        (system / "star1_formation").mkdir(parents=True)
        step0 = system / "star1_formation" / "step0"
        step0.write_text("step0")

        _set_globals({"single_star_grid": True, "star1_controls_special": []}, {})
        totest.find_inlist_from_scenario("posydon", "master-abc123", "CO-He_star")

        special = repo / "special_single_star_user_inlist"
        assert totest.mesa_inlists["zams_filename_2"] is None
        assert totest.mesa_inlists["star1_controls_user"] == [str(step0)]
        assert totest.mesa_inlists["star1_controls_special"] == [str(special)]
        content = special.read_bytes()
        assert b"\tx_logical_ctrl(1) = .true.\n" in content
        assert b"&star_job\n\n" in content
        assert b"/ / ! end of star_job namelist" in content

    def test_generic_binary_star2_formation(self, mock_git, tmp_path):
        """Generic binary: star2_formation dir unsets zams_filename_2 and populates arrays."""
        repo = self._posydon_repo(tmp_path)
        common = repo / "r11701" / "default_common_inlists"
        (common / "binary").mkdir(parents=True)
        (common / "binary" / "inlist2").write_text(
            "&controls\n  zams_filename = 'zams2.data'\n/ ! end\n"
        )
        (repo / "r11701" / "ZAMS_models").mkdir(parents=True)
        (repo / "r11701" / "ZAMS_models" / "zams2.data").write_text("zams data")

        system = repo / "r11701" / "HeMS_HeMS"
        (system / "binary").mkdir(parents=True)
        (system / "binary" / "inlist_project").write_text("project")
        (system / "star2_formation").mkdir(parents=True)
        step0 = system / "star2_formation" / "step0"
        step1 = system / "star2_formation" / "step1"
        step0.write_text("step0")
        step1.write_text("step1")

        _set_globals({"single_star_grid": False}, {})
        totest.find_inlist_from_scenario("posydon", "master-abc123", "HeMS_HeMS")

        assert totest.mesa_inlists["zams_filename_2"] is None
        assert totest.mesa_inlists["binary_controls_user"] == str(
            system / "binary" / "inlist_project"
        )
        assert totest.mesa_inlists["binary_job_user"] == str(
            system / "binary" / "inlist_project"
        )
        assert totest.mesa_inlists["star2_formation_controls_user"] == [
            str(step0),
            str(step1),
        ]
        assert totest.mesa_inlists["star2_formation_job_user"] == [
            str(step0),
            str(step1),
        ]

    def test_cwd_restored_after_return(self, mock_git, tmp_path):
        """find_inlist_from_scenario restores the original cwd before returning."""
        original = os.getcwd()
        _set_globals({"single_star_grid": False}, {})
        totest.find_inlist_from_scenario("posydon", "master-abc123", "HeMS-HMS")
        assert os.getcwd() == original


class TestConstructStaticInlist:
    """Tests for construct_static_inlist()."""

    def _ensure_grid_dirs(self, working_directory):
        """The binary/star1/star2 dirs are normally created by make_executables."""
        for sub in ["binary", "star1", "star2"]:
            os.makedirs(os.path.join(working_directory, sub), exist_ok=True)

    def test_binary_only_grid(self, tmp_path, inlist_files, capsys):
        """Binary-only grid: writes 3 inlists, merges params, returns 5-tuple."""
        self._ensure_grid_dirs(tmp_path)
        mesa_inlists = _binary_mesa_inlists(inlist_files)
        grid_parameters = ["mdot_scheme", "evolve_both_stars", "m1", "m2"]
        result = totest.construct_static_inlist(
            mesa_inlists, grid_parameters, working_directory=str(tmp_path)
        )
        (
            inlist_star1_formation,
            inlist_star2_formation,
            inlist_binary_project,
            inlist_star1_binary,
            inlist_star2_binary,
        ) = result

        assert inlist_star1_formation is None
        assert inlist_star2_formation is None
        assert inlist_binary_project == os.path.join(
            tmp_path, "binary", "inlist_project"
        )
        assert inlist_star1_binary == os.path.join(tmp_path, "binary", "inlist1")
        assert inlist_star2_binary == os.path.join(tmp_path, "binary", "inlist2")

        project = Path(inlist_binary_project).read_text()
        assert "&binary_controls" in project
        assert "\tmdot_scheme = 'Kolb'" in project
        assert "\thistory_interval = 100" in project
        assert "&binary_job" in project
        assert "\tevolve_both_stars = .true." in project
        assert "\tinlist_names(1) = '{0}'".format(inlist_star1_binary) in project
        assert "\tinlist_names(2) = '{0}'".format(inlist_star2_binary) in project

        inlist1 = Path(inlist_star1_binary).read_text()
        assert "&controls" in inlist1
        assert "\tm1 = 10.0" in inlist1
        assert "\t1 = 3" in inlist1  # num_x_ctrls renamed
        assert "\tdo_history_file = .true." in inlist1
        assert "&star_job" in inlist1
        assert "\tload_saved_model = .false." in inlist1
        assert "\twrite_profile_when_terminate = .false." in inlist1

        inlist2 = Path(inlist_star2_binary).read_text()
        assert "&controls" in inlist2
        assert "\tm2 = 8.0" in inlist2

        captured = capsys.readouterr()
        assert (
            "Grid parameters that effect binary_controls: mdot_scheme" in captured.out
        )
        assert (
            "Grid parameters that effect binary_job: evolve_both_stars" in captured.out
        )
        assert "Grid parameters that effect star1_binary_controls: m1" in captured.out
        assert "Grid parameters that effect star2_binary_controls: m2" in captured.out

    def test_read_extra_wiring(self, tmp_path, inlist_files, capsys):
        """Grid params affecting a section trigger read_extra_* lines."""
        mesa_inlists = _binary_mesa_inlists(inlist_files)
        self._ensure_grid_dirs(tmp_path)
        totest.construct_static_inlist(
            mesa_inlists, ["m1", "m2"], working_directory=str(tmp_path)
        )
        inlist1 = Path(tmp_path, "binary", "inlist1").read_text()
        inlist2 = Path(tmp_path, "binary", "inlist2").read_text()
        assert "\tread_extra_controls_inlist1 = .true." in inlist1
        assert (
            "\textra_controls_inlist1_name = 'inlist_grid_star1_binary_controls'"
            in inlist1
        )
        assert "\tread_extra_controls_inlist1 = .true." in inlist2
        assert (
            "\textra_controls_inlist1_name = 'inlist_grid_star2_binary_controls'"
            in inlist2
        )

    def test_zams_filename_suppresses_star1_formation(self, tmp_path, inlist_files):
        """A provided zams_filename_1 skips formation and writes zams_filename."""
        mesa_inlists = _binary_mesa_inlists(inlist_files)
        mesa_inlists["zams_filename_1"] = "/some/path/zams1.data"
        self._ensure_grid_dirs(tmp_path)
        result = totest.construct_static_inlist(
            mesa_inlists, [], working_directory=str(tmp_path)
        )
        assert result[0] is None
        inlist1 = Path(tmp_path, "binary", "inlist1").read_text()
        assert "\tzams_filename = '/some/path/zams1.data'" in inlist1

    def test_star1_formation_two_steps(self, tmp_path, inlist_files):
        """2-step star1 formation produces step files with save/load chaining."""
        step0 = tmp_path / "formation" / "step0"
        step1 = tmp_path / "formation" / "step1"
        step0.parent.mkdir(parents=True)
        step0.write_text(
            "&controls\n  initial_mass = 10.0\n/ ! end\n&star_job\n/ ! end\n"
        )
        step1.write_text(
            "&controls\n  initial_mass = 10.0\n/ ! end\n&star_job\n/ ! end\n"
        )

        mesa_inlists = _binary_mesa_inlists(inlist_files)
        mesa_inlists.pop("star1_controls_user")
        mesa_inlists.pop("star1_job_user")
        mesa_inlists["star1_formation_controls_user"] = [str(step0), str(step1)]
        mesa_inlists["star1_formation_job_user"] = [str(step0), str(step1)]
        mesa_inlists["final_model_star1"] = True

        self._ensure_grid_dirs(tmp_path)
        result = totest.construct_static_inlist(
            mesa_inlists, [], working_directory=str(tmp_path)
        )
        step_file0 = os.path.join(tmp_path, "star1", "inlist_step0")
        step_file1 = os.path.join(tmp_path, "star1", "inlist_step1")
        assert result[0] == " {0} {1}".format(step_file0, step_file1)

        s0 = Path(step_file0).read_text()
        assert "&controls" in s0
        assert "\tsave_model_when_terminate = .true." in s0
        assert "\tsave_model_filename = 'initial_star1_step0.mod'" in s0

        s1 = Path(step_file1).read_text()
        assert "\tcreate_pre_main_sequence_model = .false." in s1
        assert "\tload_saved_model = .true." in s1
        assert "\tsaved_model_name = 'initial_star1_step0.mod'" in s1
        assert "\tsave_model_when_terminate = .true." in s1
        assert "\tsave_model_filename = 'initial_star1_step1.mod'" in s1

        # binary inlist1 is wired to load the final (step1) formation model
        inlist1 = Path(tmp_path, "binary", "inlist1").read_text()
        assert "\tload_saved_model = .true." in inlist1
        assert "\tsaved_model_name = 'initial_star1_step1.mod'" in inlist1
        assert "\tsave_model_filename = 'final_star1.mod'" in inlist1

    def test_single_star_grid_gathers_star1_keys(self, tmp_path, inlist_files):
        """single_star_grid=True gathers star1_* keys into the formation steps."""
        special = tmp_path / "special.inlist"
        special.write_text(
            "&controls\n  initial_mass = 10.0\n/ ! end\n&star_job\n/ ! end\n"
        )

        mesa_inlists = _binary_mesa_inlists(inlist_files)
        mesa_inlists["single_star_grid"] = True
        mesa_inlists["star1_controls_special"] = str(special)
        mesa_inlists["star1_job_special"] = str(special)

        self._ensure_grid_dirs(tmp_path)
        result = totest.construct_static_inlist(
            mesa_inlists, [], working_directory=str(tmp_path)
        )
        assert result[0] is not None
        step_file = os.path.join(tmp_path, "star1", "inlist_step0")
        assert result[0] == " {0}".format(step_file)
        content = Path(step_file).read_text()
        assert "\tinitial_mass = 10.0" in content

    def test_output_flags(self, tmp_path, inlist_files):
        """Output control flags produce the corresponding inlist lines."""
        mesa_inlists = _binary_mesa_inlists(inlist_files)
        mesa_inlists["final_profile_star1"] = True
        mesa_inlists["final_model_star2"] = True
        mesa_inlists["binary_history"] = False

        self._ensure_grid_dirs(tmp_path)
        totest.construct_static_inlist(
            mesa_inlists, [], working_directory=str(tmp_path)
        )
        inlist1 = Path(tmp_path, "binary", "inlist1").read_text()
        assert "\twrite_profile_when_terminate = .true." in inlist1
        assert (
            "\tfilename_for_profile_when_terminate = 'final_profile_star1.data'"
            in inlist1
        )

        inlist2 = Path(tmp_path, "binary", "inlist2").read_text()
        assert "\tsave_model_when_terminate = .true." in inlist2
        assert "\tsave_model_filename = 'final_star2.mod'" in inlist2

        project = Path(tmp_path, "binary", "inlist_project").read_text()
        assert "\thistory_interval = -1" in project

    def test_num_x_ctrls_renamed(self, tmp_path, inlist_files):
        """num_x_ctrls params are renamed by replacing the substring with '1'."""
        ctrl = tmp_path / "ctrl.inlist"
        ctrl.write_text("&controls\n  num_x_ctrls = 5\n/ ! end\n")
        mesa_inlists = _binary_mesa_inlists(inlist_files)
        mesa_inlists["star1_controls_user"] = str(ctrl)
        self._ensure_grid_dirs(tmp_path)
        totest.construct_static_inlist(
            mesa_inlists, [], working_directory=str(tmp_path)
        )
        inlist1 = Path(tmp_path, "binary", "inlist1").read_text()
        assert "\t1 = 5" in inlist1


class TestMakeExecutables:
    """Tests for make_executables()."""

    def test_folder_layout_and_mk_script(self, tmp_path, monkeypatch, capsys):
        """Creates the folder layout, writes mk, routes extras, returns exe paths."""
        os.chdir(tmp_path)
        extras = tmp_path / "extras"
        extras.mkdir()
        binary_extras = extras / "run_binary_extras.f"
        binary_extras.write_text("! binary extras")
        binary_run = extras / "run_binary.f"
        binary_run.write_text("! binary run")
        star_run = extras / "run.f"
        star_run.write_text("! star run")
        star1_extras = extras / "star1_extras.f"
        star1_extras.write_text("! star1 extras")
        star2_extras = extras / "star2_extras.f"
        star2_extras.write_text("! star2 extras")
        makefile_binary = extras / "makefile_binary"
        makefile_binary.write_text("! binary makefile")
        makefile_star = extras / "makefile_star"
        makefile_star.write_text("! star makefile")
        other = extras / "other.txt"
        other.write_text("other")

        mesa_extras = {
            "mesa_binary_extras": str(binary_extras),
            "mesa_binary_run": str(binary_run),
            "mesa_star_run": str(star_run),
            "star1_extras": str(star1_extras),
            "star2_extras": str(star2_extras),
            "makefile_binary": str(makefile_binary),
            "makefile_star": str(makefile_star),
            "mesa_dir": "/some/mesa",
            "misc_extra": str(other),
        }

        system_calls = []
        monkeypatch.setattr(totest.os, "system", system_calls.append)

        work = tmp_path / "work"
        result = totest.make_executables(mesa_extras, working_directory=str(work))

        # folder layout
        for sub in [
            "star1/src",
            "star1/make",
            "star2/src",
            "star2/make",
            "binary/src",
            "binary/make",
        ]:
            assert os.path.isdir(os.path.join(work, sub))

        # extras routing
        assert os.path.isfile(
            os.path.join(work, "binary", "src", "run_binary_extras.f")
        )
        assert os.path.isfile(os.path.join(work, "binary", "src", "run_binary.f"))
        assert os.path.isfile(os.path.join(work, "star1", "src", "run.f"))
        assert os.path.isfile(os.path.join(work, "star2", "src", "run.f"))
        assert os.path.isfile(os.path.join(work, "star1", "src", "star1_extras.f"))
        assert os.path.isfile(os.path.join(work, "star2", "src", "star2_extras.f"))
        assert os.path.isfile(os.path.join(work, "binary", "make", "makefile_binary"))
        assert os.path.isfile(os.path.join(work, "star1", "make", "makefile_star"))
        assert os.path.isfile(os.path.join(work, "star2", "make", "makefile_star"))
        assert os.path.isfile(os.path.join(work, "other.txt"))

        # mk script content (note: mk is written to cwd, not working_directory)
        mk = Path(tmp_path, "mk").read_text()
        assert "cd {0}\n".format(os.path.join(work, "binary", "make")) in mk
        assert "make -f makefile_binary\n" in mk
        assert "cd {0}\n".format(os.path.join(work, "star1", "make")) in mk
        assert "cd {0}\n".format(os.path.join(work, "star2", "make")) in mk
        assert "make -f makefile_star\n" in mk

        # os.system invocations
        assert system_calls == ["chmod 755 mk", "./mk"]

        # returned executable paths
        assert result == (
            os.path.join(work, "binary", "binary"),
            os.path.join(work, "star1", "star"),
            os.path.join(work, "star2", "star"),
        )


class TestConstructCommandLine:
    """Tests for construct_command_line()."""

    def test_fixed_grid(self):
        """Fixed grid command line uses the 'python {15} ...' template."""
        command = totest.construct_command_line(
            1,
            "/grid/points.csv",
            "/work/binary/binary",
            "/work/star1/star",
            "/work/star2/star",
            "/work/binary/inlist_project",
            "/work/binary/inlist1",
            "/work/binary/inlist2",
            "/work/star1/inlist_step0",
            None,
            "/work/column_lists/history_columns.list",
            "/work/column_lists/binary_history_columns.list",
            "/work/column_lists/profile_columns.list",
            "/work",
            "fixed",
            "/path/to/posydon-run-grid",
            None,
        )
        expected = (
            "python /path/to/posydon-run-grid --mesa-grid /grid/points.csv "
            "--mesa-binary-executable /work/binary/binary "
            "--mesa-star1-executable /work/star1/star "
            "--mesa-star2-executable /work/star2/star "
            "--mesa-binary-inlist-project /work/binary/inlist_project "
            "--mesa-binary-inlist1 /work/binary/inlist1 "
            "--mesa-binary-inlist2 /work/binary/inlist2 "
            "--mesa-star1-inlist-project /work/star1/inlist_step0 "
            "--mesa-star2-inlist-project None "
            "--mesa-star-history-columns /work/column_lists/history_columns.list "
            "--mesa-binary-history-columns /work/column_lists/binary_history_columns.list "
            "--mesa-profile-columns /work/column_lists/profile_columns.list "
            "--output-directory /work --grid-type fixed --psycris-inifile None"
        )
        assert command == expected

    def test_dynamic_grid(self):
        """Dynamic grid command line wraps the run-grid call in mpirun/mpi4py."""
        command = totest.construct_command_line(
            4,
            "/grid/points.h5",
            "/work/binary/binary",
            "/work/star1/star",
            "/work/star2/star",
            "/work/binary/inlist_project",
            "/work/binary/inlist1",
            "/work/binary/inlist2",
            "/work/star1/inlist_step0",
            None,
            "/work/history.list",
            "/work/binary_history.list",
            "/work/profile.list",
            "/work",
            "dynamic",
            "/path/to/posydon-run-grid",
            "/work/psycris.ini",
        )
        assert command.startswith(
            "mpirun --bind-to none -np 4 python -m mpi4py /path/to/posydon-run-grid "
            "--mesa-grid /grid/points.h5 --mesa-binary-executable /work/binary/binary "
        )
        assert " --grid-type dynamic " in command
        assert " --psycris-inifile /work/psycris.ini" in command

    def test_keep_profiles_and_photos(self):
        """keep_profiles/keep_photos append the corresponding flags."""
        command = totest.construct_command_line(
            1,
            "/grid.csv",
            "/b",
            "/s1",
            "/s2",
            "/p",
            "/i1",
            "/i2",
            None,
            None,
            "/h",
            "/bh",
            "/pr",
            "/out",
            "fixed",
            "/path/to/posydon-run-grid",
            None,
            keep_profiles=True,
            keep_photos=True,
        )
        assert command.endswith(" --keep_profiles --keep_photos")
