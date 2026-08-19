"""Unit tests for posydon.CLI.grids.scenario."""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

import os
import subprocess
from pathlib import Path

import pytest

from posydon.CLI.grids import scenario

totest = scenario


@pytest.fixture(autouse=True)
def _env(tmp_path, monkeypatch):
    """Always provide HOME/MESA_DIR and restore cwd after every test."""
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("MESA_DIR", str(tmp_path / "mesa"))
    original = os.getcwd()
    yield
    os.chdir(original)


class TestSetupInlistsFromScenario:
    """Tests for setup_inlists_from_scenario() with a mocked _run_git."""

    @pytest.fixture
    def mock_git(self, monkeypatch):
        """Mock _run_git and record every call in order."""
        calls = []

        def fake_run_git(args, cwd=None):
            calls.append(list(args))
            if args[0] == "git" and args[1] == "clone":
                os.makedirs(args[-1], exist_ok=True)
            return subprocess.CompletedProcess(args, 0, b"", b"")

        monkeypatch.setattr(totest, "_run_git", fake_run_git)
        return calls

    def _posydon_repo(self, home):
        """Create the ~/.posydon_mesa_inlists clone layout."""
        repo = Path(home) / ".posydon_mesa_inlists"
        if not repo.is_dir():
            repo.mkdir()
        return repo

    def test_posydon_source_repo_present_only_checkout_pull(self, mock_git, tmp_path):
        """Repo already cloned: clone is skipped, checkout/pull order preserved."""
        repo = self._posydon_repo(tmp_path)
        mesa_inlists = {"single_star_grid": False}
        totest.setup_inlists_from_scenario(
            "posydon", "master-abc123", "HeMS_HeMS", mesa_inlists, {}
        )

        assert mock_git == [
            ["git", "checkout", "master"],
            ["git", "pull"],
            ["git", "checkout", "abc123"],
        ]
        assert repo.is_dir()

    def test_posydon_source_repo_absent_clone(self, mock_git, tmp_path):
        """Repo missing: clone runs first, then checkout/pull/checkout order."""
        mesa_inlists = {"single_star_grid": False}
        totest.setup_inlists_from_scenario(
            "posydon", "master-abc123", "HeMS_HeMS", mesa_inlists, {}
        )

        assert mock_git == [
            [
                "git",
                "clone",
                "https://github.com/POSYDON-code/POSYDON-MESA-INLISTS.git",
                os.path.join(tmp_path, ".posydon_mesa_inlists"),
            ],
            ["git", "checkout", "master"],
            ["git", "pull"],
            ["git", "checkout", "abc123"],
        ]

    def test_user_source_invalid_gitcommit(self, mock_git, tmp_path):
        """Invalid user gitcommit (no '-') raises ValueError with exact message."""
        with pytest.raises(
            ValueError,
            match="You have supplied an invalid user gitcommit format, must be of format 'branch-githash'",
        ):
            totest.setup_inlists_from_scenario(
                "user", "notagithash", "HeMS-HMS", {}, {}
            )

    def test_user_source_repo_absent_clone(self, mock_git, tmp_path):
        """User repo missing: USER-MESA-INLISTS is cloned first."""
        mesa_inlists = {"single_star_grid": False}
        totest.setup_inlists_from_scenario(
            "user", "master-abc123", "HeMS-HMS", mesa_inlists, {}
        )

        assert mock_git == [
            [
                "git",
                "clone",
                "https://github.com/POSYDON-code/USER-MESA-INLISTS.git",
                os.path.join(tmp_path, ".user_mesa_inlists"),
            ],
            ["git", "checkout", "master"],
            ["git", "pull"],
            ["git", "checkout", "abc123"],
        ]

    def test_user_source_repo_present_pwarn(self, mock_git, tmp_path, monkeypatch):
        """User repo already present: OverwriteWarning is issued."""
        repo = Path(tmp_path) / ".user_mesa_inlists"
        repo.mkdir()
        warnings_calls = []
        monkeypatch.setattr(
            totest,
            "Pwarn",
            lambda message, category: warnings_calls.append((message, category)),
        )
        mesa_inlists = {"single_star_grid": False}
        totest.setup_inlists_from_scenario(
            "user", "master-abc123", "HeMS-HMS", mesa_inlists, {}
        )

        assert warnings_calls == [
            ("git repository is already there, using that", "OverwriteWarning")
        ]
        assert mock_git == [
            ["git", "checkout", "master"],
            ["git", "pull"],
            ["git", "checkout", "abc123"],
        ]

    def test_invalid_source(self, mock_git, tmp_path):
        """Unknown source raises ValueError with the exact message."""
        with pytest.raises(
            ValueError,
            match="supplied source is not valid/understood. Valid sources are user and posydon",
        ):
            totest.setup_inlists_from_scenario(
                "other", "branch-hash", "HeMS-HMS", {}, {}
            )

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

        mesa_inlists = {"single_star_grid": False}
        mesa_extras = {}
        totest.setup_inlists_from_scenario(
            "posydon", "master-abc123", "HeMS_HeMS", mesa_inlists, mesa_extras
        )

        assert mesa_inlists["zams_filename_1"] == str(zams)
        assert "zams_filename_2" not in mesa_inlists

        common_location = "{0}/r11701/default_common_inlists/".format(repo)
        assert (
            mesa_inlists["star1_controls_posydon_defaults"]
            == "{0}/binary/inlist1".format(common_location)
        )
        assert (
            mesa_inlists["star1_job_posydon_defaults"]
            == "{0}/binary/inlist1".format(common_location)
        )
        assert (
            mesa_inlists["star2_controls_posydon_defaults"]
            == "{0}/binary/inlist2".format(common_location)
        )
        assert (
            mesa_inlists["star2_job_posydon_defaults"]
            == "{0}/binary/inlist2".format(common_location)
        )
        assert (
            mesa_inlists["binary_controls_posydon_defaults"]
            == "{0}/binary/inlist_project".format(common_location)
        )
        assert (
            mesa_inlists["binary_job_posydon_defaults"]
            == "{0}/binary/inlist_project".format(common_location)
        )
        assert (
            mesa_inlists["star_history_columns"]
            == "{0}/history_columns.list".format(common_location)
        )
        assert (
            mesa_inlists["binary_history_columns"]
            == "{0}/binary_history_columns.list".format(common_location)
        )
        assert (
            mesa_inlists["profile_columns"]
            == "{0}/profile_columns.list".format(common_location)
        )
        assert (
            mesa_extras["posydon_binary_extras"]
            == "{0}/binary/src/run_binary_extras.f".format(common_location)
        )
        assert (
            mesa_extras["posydon_star_binary_extras"]
            == "{0}/binary/src/run_star_extras.f".format(common_location)
        )
        assert (
            mesa_extras["mesa_star1_extras"]
            == "{0}/binary/src/run_star_extras.f".format(common_location)
        )

    def test_zams_line_present_but_file_missing(self, mock_git, tmp_path):
        """A zams line whose data file is missing does not set zams_filename_1."""
        repo = self._posydon_repo(tmp_path)
        common = repo / "r11701" / "default_common_inlists"
        (common / "binary").mkdir(parents=True)
        (common / "binary" / "inlist1").write_text(
            "&controls\n  zams_filename = 'missing.data'\n/ ! end\n"
        )

        mesa_inlists = {"single_star_grid": False}
        totest.setup_inlists_from_scenario(
            "posydon", "master-abc123", "HeMS_HeMS", mesa_inlists, {}
        )

        assert "zams_filename_1" not in mesa_inlists

    def test_inlist1_no_zams_line(self, mock_git, tmp_path):
        """An inlist1 without a zams_filename line sets nothing."""
        repo = self._posydon_repo(tmp_path)
        common = repo / "r11701" / "default_common_inlists"
        (common / "binary").mkdir(parents=True)
        (common / "binary" / "inlist1").write_text(
            "&controls\n  initial_mass = 10.0\n/ ! end\n"
        )

        mesa_inlists = {"single_star_grid": False}
        totest.setup_inlists_from_scenario(
            "posydon", "master-abc123", "HeMS_HeMS", mesa_inlists, {}
        )

        assert "zams_filename_1" not in mesa_inlists

    def test_inlist2_zams_line_present_but_file_missing(self, mock_git, tmp_path):
        """A zams line in inlist2 whose data file is missing sets nothing."""
        repo = self._posydon_repo(tmp_path)
        common = repo / "r11701" / "default_common_inlists"
        (common / "binary").mkdir(parents=True)
        (common / "binary" / "inlist2").write_text(
            "&controls\n  zams_filename = 'missing.data'\n/ ! end\n"
        )

        mesa_inlists = {"single_star_grid": False}
        totest.setup_inlists_from_scenario(
            "posydon", "master-abc123", "HeMS_HeMS", mesa_inlists, {}
        )

        assert "zams_filename_2" not in mesa_inlists

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

        mesa_inlists = {"single_star_grid": False}
        totest.setup_inlists_from_scenario(
            "posydon", "master-abc123", "HeMS-HMS", mesa_inlists, {}
        )

        common_inlist1 = "{0}/r11701/default_common_inlists//binary/inlist1".format(
            repo
        )
        assert mesa_inlists["zams_filename_1"] is None
        assert mesa_inlists["star1_formation_controls_user"] == [
            str(step0),
            str(step1),
        ]
        assert mesa_inlists["star1_formation_job_user"] == [
            str(step0),
            str(step1),
        ]
        assert mesa_inlists["star1_formation_controls_posydon_defaults"] == [
            common_inlist1,
            common_inlist1,
        ]
        assert mesa_inlists["star1_formation_job_posydon_defaults"] == [
            common_inlist1,
            common_inlist1,
        ]

    def test_hems_hms_no_star1_formation_dir(self, mock_git, tmp_path):
        """HeMS-HMS without a star1_formation dir leaves formation arrays unset."""
        repo = self._posydon_repo(tmp_path)
        system = repo / "r11701" / "HeMS-HMS"
        (system / "binary").mkdir(parents=True)
        inlist1 = system / "binary" / "inlist1"
        inlist2 = system / "binary" / "inlist2"
        inlist1.write_text("inlist1")
        inlist2.write_text("inlist2")

        mesa_inlists = {"single_star_grid": False}
        totest.setup_inlists_from_scenario(
            "posydon", "master-abc123", "HeMS-HMS", mesa_inlists, {}
        )

        assert mesa_inlists["star1_controls_user"] == str(inlist1)
        assert mesa_inlists["star1_job_user"] == str(inlist1)
        assert mesa_inlists["star2_controls_user"] == str(inlist2)
        assert mesa_inlists["star2_job_user"] == str(inlist2)
        assert "star1_formation_controls_user" not in mesa_inlists
        assert "star1_formation_job_user" not in mesa_inlists

    def test_hms_hms_single_star_special_inlist(self, mock_git, tmp_path):
        """HMS-HMS single star: special inlist written with x_logical_ctrl(1)=.true."""
        repo = self._posydon_repo(tmp_path)
        mesa_inlists = {"single_star_grid": True}
        totest.setup_inlists_from_scenario(
            "posydon", "master-abc123", "HMS-HMS", mesa_inlists, {}
        )

        special = repo / "special_single_star_user_inlist"
        assert mesa_inlists["star1_controls_special"] == str(special)
        content = special.read_bytes()
        assert b"&controls\n\n" in content
        assert b"\tx_logical_ctrl(1) = .true.\n" in content
        assert b"/ ! end of star1_controls namelist" in content
        assert b"&star_job" not in content

    def test_hms_hms_single_star_special_inlist_exists(
        self, mock_git, tmp_path, monkeypatch
    ):
        """HMS-HMS single star: an existing special inlist is overwritten."""
        repo = self._posydon_repo(tmp_path)
        special = repo / "special_single_star_user_inlist"
        special.write_bytes(b"old content")
        warnings_calls = []
        monkeypatch.setattr(
            totest,
            "Pwarn",
            lambda message, category: warnings_calls.append((message, category)),
        )
        mesa_inlists = {"single_star_grid": True}
        totest.setup_inlists_from_scenario(
            "posydon", "master-abc123", "HMS-HMS", mesa_inlists, {}
        )

        assert warnings_calls == [
            ("git repository is already there, using that", "OverwriteWarning"),
            ("Replace " + str(special), "OverwriteWarning"),
        ]
        assert mesa_inlists["star1_controls_special"] == str(special)
        content = special.read_bytes()
        assert b"&controls\n\n" in content
        assert b"\tx_logical_ctrl(1) = .true.\n" in content
        assert b"/ ! end of star1_controls namelist" in content
        assert b"&star_job" not in content

    def test_hms_hms_single_star_grid_false_falls_through(self, mock_git, tmp_path):
        """HMS-HMS with single_star_grid=False does not write a special inlist."""
        repo = self._posydon_repo(tmp_path)
        mesa_inlists = {"single_star_grid": False}
        totest.setup_inlists_from_scenario(
            "posydon", "master-abc123", "HMS-HMS", mesa_inlists, {}
        )

        assert "star1_controls_special" not in mesa_inlists
        assert not (repo / "special_single_star_user_inlist").exists()

    def test_co_he_star_single_star_special_inlist(self, mock_git, tmp_path):
        """CO-He_star single star: special inlist gets an extra &star_job section."""
        repo = self._posydon_repo(tmp_path)
        system = repo / "r11701" / "CO-He_star"
        (system / "star1_formation").mkdir(parents=True)
        step0 = system / "star1_formation" / "step0"
        step0.write_text("step0")

        mesa_inlists = {"single_star_grid": True, "star1_controls_special": []}
        totest.setup_inlists_from_scenario(
            "posydon", "master-abc123", "CO-He_star", mesa_inlists, {}
        )

        special = repo / "special_single_star_user_inlist"
        assert mesa_inlists["zams_filename_2"] is None
        assert mesa_inlists["star1_controls_user"] == [str(step0)]
        assert mesa_inlists["star1_job_user"] == [str(step0)]
        assert mesa_inlists["star1_controls_special"] == [str(special)]
        content = special.read_bytes()
        assert b"\tx_logical_ctrl(1) = .true.\n" in content
        assert b"&star_job\n\n" in content
        assert b"/ / ! end of star_job namelist" in content

    def test_co_he_star_special_inlist_exists(self, mock_git, tmp_path, monkeypatch):
        """CO-He_star single star: an existing special inlist is overwritten."""
        repo = self._posydon_repo(tmp_path)
        system = repo / "r11701" / "CO-He_star"
        (system / "star1_formation").mkdir(parents=True)
        step0 = system / "star1_formation" / "step0"
        step0.write_text("step0")
        special = repo / "special_single_star_user_inlist"
        special.write_bytes(b"old content")
        warnings_calls = []
        monkeypatch.setattr(
            totest,
            "Pwarn",
            lambda message, category: warnings_calls.append((message, category)),
        )
        mesa_inlists = {"single_star_grid": True, "star1_controls_special": []}
        totest.setup_inlists_from_scenario(
            "posydon", "master-abc123", "CO-He_star", mesa_inlists, {}
        )

        assert warnings_calls == [
            ("git repository is already there, using that", "OverwriteWarning"),
            ("Replace " + str(special), "OverwriteWarning"),
        ]
        assert mesa_inlists["star1_controls_special"] == [str(special)]
        content = special.read_bytes()
        assert b"&controls\n\n" in content
        assert b"\tx_logical_ctrl(1) = .true.\n" in content
        assert b"/ ! end of star1_controls namelist" in content
        assert b"&star_job\n\n" in content
        assert b"/ / ! end of star_job namelist" in content

    def test_co_he_star_single_star_grid_false_no_special_inlist(
        self, mock_git, tmp_path
    ):
        """CO-He_star with single_star_grid=False does not write a special inlist."""
        repo = self._posydon_repo(tmp_path)
        mesa_inlists = {"single_star_grid": False}
        totest.setup_inlists_from_scenario(
            "posydon", "master-abc123", "CO-He_star", mesa_inlists, {}
        )

        assert "star1_controls_special" not in mesa_inlists
        assert not (repo / "special_single_star_user_inlist").exists()

    def test_generic_binary_star1_formation_unsets_zams(self, mock_git, tmp_path):
        """Generic binary: star1_formation dir unsets zams_filename_1."""
        repo = self._posydon_repo(tmp_path)
        common = repo / "r11701" / "default_common_inlists"
        (common / "binary").mkdir(parents=True)
        (common / "binary" / "inlist1").write_text(
            "&controls\n  zams_filename = 'zams1.data'\n/ ! end\n"
        )
        (repo / "r11701" / "ZAMS_models").mkdir(parents=True)
        (repo / "r11701" / "ZAMS_models" / "zams1.data").write_text("zams data")

        system = repo / "r11701" / "HeMS_HeMS"
        (system / "star1_formation").mkdir(parents=True)
        step0 = system / "star1_formation" / "step0"
        step1 = system / "star1_formation" / "step1"
        step0.write_text("step0")
        step1.write_text("step1")

        mesa_inlists = {"single_star_grid": False}
        totest.setup_inlists_from_scenario(
            "posydon", "master-abc123", "HeMS_HeMS", mesa_inlists, {}
        )

        common_inlist1 = "{0}/r11701/default_common_inlists//binary/inlist1".format(
            repo
        )
        assert mesa_inlists["zams_filename_1"] is None
        assert mesa_inlists["star1_formation_controls_user"] == [
            str(step0),
            str(step1),
        ]
        assert mesa_inlists["star1_formation_job_user"] == [
            str(step0),
            str(step1),
        ]
        assert mesa_inlists["star1_formation_controls_posydon_defaults"] == [
            common_inlist1,
            common_inlist1,
        ]
        assert mesa_inlists["star1_formation_job_posydon_defaults"] == [
            common_inlist1,
            common_inlist1,
        ]

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

        mesa_inlists = {"single_star_grid": False}
        totest.setup_inlists_from_scenario(
            "posydon", "master-abc123", "HeMS_HeMS", mesa_inlists, {}
        )

        common_inlist1 = "{0}/r11701/default_common_inlists//binary/inlist1".format(
            repo
        )
        assert mesa_inlists["zams_filename_2"] is None
        assert mesa_inlists["binary_controls_user"] == str(
            system / "binary" / "inlist_project"
        )
        assert mesa_inlists["binary_job_user"] == str(
            system / "binary" / "inlist_project"
        )
        assert mesa_inlists["star2_formation_controls_user"] == [
            str(step0),
            str(step1),
        ]
        assert mesa_inlists["star2_formation_job_user"] == [
            str(step0),
            str(step1),
        ]
        assert mesa_inlists["star2_formation_controls_posydon_defaults"] == [
            common_inlist1,
            common_inlist1,
        ]
        assert mesa_inlists["star2_formation_job_posydon_defaults"] == [
            common_inlist1,
            common_inlist1,
        ]

    def test_generic_binary_columns_and_extras(self, mock_git, tmp_path):
        """Generic binary: column lists and src extras get picked up."""
        repo = self._posydon_repo(tmp_path)
        system = repo / "r11701" / "HeMS_HeMS"
        (system / "binary").mkdir(parents=True)
        (system / "binary" / "inlist1").write_text("i1")
        (system / "binary" / "inlist2").write_text("i2")
        (system / "history_columns.list").write_text("h")
        (system / "binary_history_columns.list").write_text("bh")
        (system / "profile_columns.list").write_text("p")
        (system / "src").mkdir(parents=True)
        (system / "src" / "run_binary_extras.f").write_text("be")
        (system / "src" / "run_star_extras.f").write_text("se")

        mesa_inlists = {"single_star_grid": False}
        mesa_extras = {}
        totest.setup_inlists_from_scenario(
            "posydon", "master-abc123", "HeMS_HeMS", mesa_inlists, mesa_extras
        )

        assert mesa_inlists["star1_controls_user"] == str(system / "binary" / "inlist1")
        assert mesa_inlists["star1_job_user"] == str(system / "binary" / "inlist1")
        assert mesa_inlists["star2_controls_user"] == str(system / "binary" / "inlist2")
        assert mesa_inlists["star2_job_user"] == str(system / "binary" / "inlist2")
        assert mesa_inlists["star_history_columns"] == str(
            system / "history_columns.list"
        )
        assert mesa_inlists["binary_history_columns"] == str(
            system / "binary_history_columns.list"
        )
        assert mesa_inlists["profile_columns"] == str(system / "profile_columns.list")
        assert mesa_extras["user_binary_extras"] == str(
            system / "src" / "run_binary_extras.f"
        )
        assert mesa_extras["user_star_binary_extras"] == str(
            system / "src" / "run_star_extras.f"
        )
        assert "star1_formation_controls_user" not in mesa_inlists
        assert "star2_formation_controls_user" not in mesa_inlists

    def test_cwd_restored_after_return(self, mock_git, tmp_path):
        """setup_inlists_from_scenario restores the original cwd before returning."""
        original = os.getcwd()
        mesa_inlists = {"single_star_grid": False}
        totest.setup_inlists_from_scenario(
            "posydon", "master-abc123", "HeMS-HMS", mesa_inlists, {}
        )
        assert os.getcwd() == original

    def test_cwd_restored_after_exception(self, mock_git, tmp_path):
        """setup_inlists_from_scenario restores the cwd even when it raises."""
        original = os.getcwd()
        with pytest.raises(ValueError):
            totest.setup_inlists_from_scenario(
                "other", "branch-hash", "HeMS-HMS", {}, {}
            )
        assert os.getcwd() == original


class TestRunGit:
    """Tests for the _run_git helper."""

    def test_run_git_returns_completed_process(self):
        """_run_git runs a command and returns its CompletedProcess."""
        result = totest._run_git(["echo", "hello"])
        assert result.returncode == 0
        assert result.stdout == b"hello\n"
        assert result.stderr == b""

    def test_run_git_with_cwd(self, tmp_path):
        """_run_git passes the cwd argument to subprocess.run."""
        result = totest._run_git(["pwd"], cwd=str(tmp_path))
        assert result.returncode == 0
        assert result.stdout == (os.path.realpath(str(tmp_path)) + "\n").encode()
