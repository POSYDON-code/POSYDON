"""Unit tests locking the current behavior of bin/posydon-run-grid."""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

import os
import sys
from types import SimpleNamespace

import numpy as np
import pandas
import pytest

from posydon.CLI.grids import run
from posydon.utils.posydonwarning import OverwriteWarning

totest = run


@pytest.fixture(autouse=True)
def _env(tmp_path, monkeypatch):
    """Provide HOME/MESA_DIR and restore cwd after every test."""
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("MESA_DIR", str(tmp_path / "mesa"))
    original = os.getcwd()
    yield
    os.chdir(original)


def _make_args(tmp_path, **overrides):
    """Build an argparse-like namespace with the fields the helpers read."""
    args = SimpleNamespace(
        temporary_directory=str(tmp_path / "temp"),
        output_directory=str(tmp_path / "output"),
        grid_point_index=None,
        grid_type="fixed",
        mesa_star_history_columns=str(tmp_path / "history_columns.list"),
        mesa_binary_history_columns=str(tmp_path / "binary_history_columns.list"),
        mesa_profile_columns=str(tmp_path / "profile_columns.list"),
    )
    for key, value in overrides.items():
        setattr(args, key, value)
    return args


def _make_main_args(tmp_path, **overrides):
    """Build an argparse-like namespace with all fields main() reads."""
    args = _make_args(
        tmp_path,
        mesa_grid=str(tmp_path / "grid.csv"),
        grid_type="fixed",
        grid_point_index=0,
        mesa_binary_executable="/exe/binary",
        mesa_star1_executable="/exe/star1",
        mesa_star2_executable="/exe/star2",
        mesa_binary_inlist_project="/inlist/project",
        mesa_binary_inlist1="/inlist/1",
        mesa_binary_inlist2="/inlist/2",
        mesa_star1_inlist_project=["/inlist/star1"],
        mesa_star2_inlist_project=["/inlist/star2"],
        psycris_inifile="/path/psycris.ini",
        keep_profiles=False,
        keep_photos=False,
        job_end=1000,
        job_before_end=300,
    )
    for key, value in overrides.items():
        setattr(args, key, value)
    return args


def _fake_clean_inlist(path):
    """Return a fixed-grid inlist structure used by most main() tests."""
    return {
        "&binary_controls": {"m1": None, "mdot_scheme": None},
        "&binary_job": {},
        "&controls": {"m2": None},
        "&star_job": {},
    }


class TestIdGenerator:
    """Tests for id_generator()."""

    def test_length(self):
        """id_generator returns a string of the requested length."""
        assert len(totest.id_generator()) == 10
        assert len(totest.id_generator(size=5)) == 5

    def test_charset(self):
        """id_generator only uses the configured alphanumeric charset."""
        value = totest.id_generator()
        assert all(c.isalnum() for c in value)


class TestCreateWorkingDirectory:
    """Tests for create_working_directory()."""

    def test_directory_name_formatting(self, tmp_path):
        """Non-log params use .4f, log params use .4e formatting."""
        os.makedirs(tmp_path / "temp")
        os.makedirs(tmp_path / "output")
        hist = tmp_path / "history_columns.list"
        hist.write_text("history")
        bhist = tmp_path / "binary_history_columns.list"
        bhist.write_text("bhistory")
        prof = tmp_path / "profile_columns.list"
        prof.write_text("profile")
        args = _make_args(tmp_path)

        grid_param_dict = {
            "initial_z": 0.0142,
            "m1": 10.0,
            "initial_period_in_days": 3.5,
        }
        work_dir, final_dir = totest.create_working_directory(grid_param_dict, args)

        name = "m1_10.0000_initial_z_1.4200e-02_initial_period_in_days_3.5000e+00"
        assert work_dir == os.path.join(args.temporary_directory, name)
        assert final_dir == os.path.join(args.output_directory, name)
        assert os.path.isdir(work_dir)
        assert os.path.isfile(os.path.join(work_dir, "history_columns.list"))
        assert os.path.isfile(os.path.join(work_dir, "binary_history_columns.list"))
        assert os.path.isfile(os.path.join(work_dir, "profile_columns.list"))

    def test_parentheses_dropped_from_name(self, tmp_path):
        """Parentheses in parameter names are dropped from the directory name."""
        os.makedirs(tmp_path / "temp")
        os.makedirs(tmp_path / "output")
        hist = tmp_path / "history_columns.list"
        hist.write_text("h")
        bhist = tmp_path / "binary_history_columns.list"
        bhist.write_text("b")
        prof = tmp_path / "profile_columns.list"
        prof.write_text("p")
        args = _make_args(tmp_path)

        work_dir, _ = totest.create_working_directory({"m1(msun)": 5.0}, args)
        assert "m1msun_5.0000_" in os.path.basename(work_dir)

    def test_grid_point_index_suffix(self, tmp_path):
        """grid_point_index appends _grid_index_{idx} to the name."""
        os.makedirs(tmp_path / "temp")
        os.makedirs(tmp_path / "output")
        hist = tmp_path / "history_columns.list"
        hist.write_text("h")
        bhist = tmp_path / "binary_history_columns.list"
        bhist.write_text("b")
        prof = tmp_path / "profile_columns.list"
        prof.write_text("p")
        args = _make_args(tmp_path, grid_point_index=7)

        work_dir, _ = totest.create_working_directory({"m1": 10.0}, args)
        assert work_dir.endswith("_grid_index_7")

    def test_dynamic_appends_random_id(self, tmp_path, monkeypatch):
        """Dynamic grids append a random id generated by id_generator."""
        os.makedirs(tmp_path / "temp")
        os.makedirs(tmp_path / "output")
        hist = tmp_path / "history_columns.list"
        hist.write_text("h")
        bhist = tmp_path / "binary_history_columns.list"
        bhist.write_text("b")
        prof = tmp_path / "profile_columns.list"
        prof.write_text("p")
        monkeypatch.setattr(totest, "id_generator", lambda: "ABC123XYZ")
        args = _make_args(tmp_path, grid_type="dynamic")

        work_dir, _ = totest.create_working_directory({"m1": 10.0}, args)
        assert work_dir.endswith("_ABC123XYZ")

    def test_existing_work_dir_is_removed(self, tmp_path, capsys):
        """An existing work dir is warned about and removed."""
        os.makedirs(tmp_path / "temp")
        os.makedirs(tmp_path / "output")
        hist = tmp_path / "history_columns.list"
        hist.write_text("h")
        bhist = tmp_path / "binary_history_columns.list"
        bhist.write_text("b")
        prof = tmp_path / "profile_columns.list"
        prof.write_text("p")
        args = _make_args(tmp_path)

        existing = tmp_path / "temp" / "m1_10.0000_"
        existing.mkdir()
        (existing / "old.dat").write_text("data")

        work_dir, _ = totest.create_working_directory({"m1": 10.0}, args)
        assert work_dir == str(existing)
        assert not (existing / "old.dat").exists()


class TestCreateBinaryInlists:
    """Tests for create_binary_inlists()."""

    def test_file_contents_and_fortran_formatting(self, tmp_path):
        """Grid-point inlists format floats with d-notation, bools, and ints."""
        work_dir = tmp_path
        binary_controls = {"mass_transfer_alpha": 0.5}
        binary_job = {"evolve_both_stars": True, "some_flag": False, "num_steps": 2.0}
        star1_binary_controls = {"m1": 10.0}
        star1_binary_job = {}
        star2_binary_controls = {}
        star2_binary_job = {"saved_model_name": "initial.mod"}

        totest.create_binary_inlists(
            binary_controls,
            binary_job,
            star1_binary_controls,
            star1_binary_job,
            star2_binary_controls,
            star2_binary_job,
            str(work_dir),
            "/path/to/inlist_project",
        )

        inlist = (work_dir / "inlist").read_text()
        assert (
            "extra_binary_controls_inlist1_name = '/path/to/inlist_project'" in inlist
        )
        assert "extra_binary_job_inlist1_name = '/path/to/inlist_project'" in inlist

        points = (work_dir / "inlist_grid_points").read_text()
        assert "&binary_controls\n" in points
        assert "\tmass_transfer_alpha = 5.0000000000d-01" in points
        assert "&binary_job\n" in points
        assert "\tevolve_both_stars = .true." in points
        assert "\tsome_flag = .false." in points
        assert "\tnum_steps = 2" in points

        star1 = (work_dir / "inlist_grid_star1_binary_controls").read_text()
        assert "&controls\n" in star1
        assert "\tm1 = 10" in star1

        # empty dicts must not produce files
        assert not (work_dir / "inlist_grid_star2_binary_controls").exists()
        assert not (work_dir / "inlist_grid_star1_binary_job").exists()

        star2 = (work_dir / "inlist_grid_star2_binary_job").read_text()
        assert "&star_job\n" in star2
        assert "saved_model_name = initial.mod" in star2

    def test_all_value_types_and_overwrite_warning(self, tmp_path):
        """Every value type is formatted in all six sections and existing
        files are warned about and replaced."""
        work_dir = tmp_path
        for name in [
            "inlist",
            "inlist_grid_points",
            "inlist_grid_star1_binary_controls",
            "inlist_grid_star2_binary_controls",
            "inlist_grid_star1_binary_job",
            "inlist_grid_star2_binary_job",
        ]:
            (work_dir / name).write_text("old")

        values = {
            "float_non_int": 0.5,
            "float_int": 2.0,
            "true_flag": True,
            "false_flag": False,
            "str_value": "initial.mod",
            "int_value": 7,
        }

        with pytest.warns(OverwriteWarning):
            totest.create_binary_inlists(
                values, values, values, values, values, values,
                str(work_dir), "/path/to/inlist_project",
            )

        points = (work_dir / "inlist_grid_points").read_text()
        assert "\tfloat_non_int = 5.0000000000d-01" in points
        assert "\tfloat_int = 2" in points
        assert "\ttrue_flag = .true." in points
        assert "\tfalse_flag = .false." in points
        assert "str_value = initial.mod" in points
        assert "int_value = 7" in points

        for name in [
            "inlist_grid_star1_binary_controls",
            "inlist_grid_star2_binary_controls",
            "inlist_grid_star1_binary_job",
            "inlist_grid_star2_binary_job",
        ]:
            content = (work_dir / name).read_text()
            assert "\tfloat_non_int = 5.0000000000d-01" in content
            assert "\tfloat_int = 2" in content
            assert "\ttrue_flag = .true." in content
            assert "\tfalse_flag = .false." in content
            assert "str_value = initial.mod" in content
            assert "int_value = 7" in content

    def test_star2_controls_and_star1_job_without_preexisting(self, tmp_path):
        """star2 controls and star1 job sections are written when the target
        files do not already exist."""
        work_dir = tmp_path
        totest.create_binary_inlists(
            {},
            {},
            {},
            {"star1_job_flag": True},
            {"star2_control": 0.5},
            {},
            str(work_dir),
            "/path/to/inlist_project",
        )

        star1 = (work_dir / "inlist_grid_star1_binary_job").read_text()
        assert "&star_job\n" in star1
        assert "\tstar1_job_flag = .true." in star1

        star2 = (work_dir / "inlist_grid_star2_binary_controls").read_text()
        assert "&controls\n" in star2
        assert "\tstar2_control = 5.0000000000d-01" in star2

    def test_all_empty_dicts_writes_only_inlist(self, tmp_path):
        """Empty dicts skip all grid-point sections and write no per-star files."""
        work_dir = tmp_path
        totest.create_binary_inlists(
            {}, {}, {}, {}, {}, {}, str(work_dir), "/path/to/inlist_project"
        )

        inlist = (work_dir / "inlist").read_text()
        assert "&binary_controls\n" in inlist
        assert "&binary_job\n" in inlist

        points = (work_dir / "inlist_grid_points").read_text()
        assert points.strip() == ""

        for name in [
            "inlist_grid_star1_binary_controls",
            "inlist_grid_star2_binary_controls",
            "inlist_grid_star1_binary_job",
            "inlist_grid_star2_binary_job",
        ]:
            assert not (work_dir / name).exists()


class TestCreateStarFormation:
    """Tests for create_star_formation()."""

    def test_with_initial_z_and_new_Z(self, tmp_path):
        """initial_mass, initial_z, Zbase, and new_Z lines are written."""
        work_dir = tmp_path
        totest.create_star_formation(
            10.0, str(work_dir), "/path/inlist", initial_z=0.0142, new_Z=True
        )

        inlist = (work_dir / "inlist").read_text()
        assert "extra_controls_inlist1_name = '/path/inlist'" in inlist
        assert "extra_star_job_inlist1_name = '/path/inlist'" in inlist

        points = (work_dir / "inlist_grid_points").read_text()
        assert "&controls\n" in points
        assert "initial_mass = 1.0000000000d+01" in points
        assert "initial_z = 1.4200000000d-02" in points
        assert "Zbase = 1.4200000000d-02" in points
        assert "&star_job\n" in points
        assert "new_Z = 1.4200000000d-02" in points

    def test_without_initial_z(self, tmp_path):
        """Without initial_z no Zbase/new_Z lines are written."""
        work_dir = tmp_path
        totest.create_star_formation(10.0, str(work_dir), "/path/inlist")

        points = (work_dir / "inlist_grid_points").read_text()
        assert "initial_mass = 1.0000000000d+01" in points
        assert "initial_z" not in points
        assert "Zbase" not in points
        assert "new_Z" not in points

    def test_existing_inlists_warned_and_overwritten(self, tmp_path):
        """Existing inlist files are warned about and replaced."""
        work_dir = tmp_path
        (work_dir / "inlist").write_text("old")
        (work_dir / "inlist_grid_points").write_text("old")

        with pytest.warns(OverwriteWarning):
            totest.create_star_formation(
                10.0, str(work_dir), "/path/inlist", initial_z=0.0142, new_Z=True
            )

        inlist = (work_dir / "inlist").read_text()
        assert "extra_controls_inlist1_name = '/path/inlist'" in inlist
        points = (work_dir / "inlist_grid_points").read_text()
        assert "initial_mass = 1.0000000000d+01" in points
        assert "new_Z = 1.4200000000d-02" in points


class TestMoveMesaOutput:
    """Tests for move_mesa_output()."""

    def test_move(self, tmp_path):
        """Default behavior moves work_dir into final_dir."""
        work = tmp_path / "work"
        work.mkdir()
        (work / "out.txt").write_text("done")
        final = tmp_path / "final"

        totest.move_mesa_output(str(work), str(final))
        assert not work.exists()
        assert (final / "out.txt").read_text() == "done"

    def test_copy(self, tmp_path):
        """copyinstead=True leaves work_dir in place and copies it."""
        work = tmp_path / "work"
        work.mkdir()
        (work / "out.txt").write_text("done")
        final = tmp_path / "final"

        totest.move_mesa_output(str(work), str(final), copyinstead=True)
        assert work.exists()
        assert (final / "out.txt").read_text() == "done"

    def test_overwrite_existing_final(self, tmp_path, capsys):
        """An existing final_dir is removed before moving."""
        work = tmp_path / "work"
        work.mkdir()
        (work / "out.txt").write_text("new")
        final = tmp_path / "final"
        final.mkdir()
        (final / "old.txt").write_text("old")

        totest.move_mesa_output(str(work), str(final))
        assert not work.exists()
        assert (final / "out.txt").read_text() == "new"
        assert not (final / "old.txt").exists()

    def test_equal_dirs_is_noop(self, tmp_path):
        """Moving to the same directory does nothing."""
        work = tmp_path / "work"
        work.mkdir()
        (work / "out.txt").write_text("done")
        totest.move_mesa_output(str(work), str(work))
        assert (work / "out.txt").read_text() == "done"


class TestRunMesa:
    """Tests for run_mesa()."""

    def _build_work_dir(self, tmp_path):
        work = tmp_path / "work"
        work.mkdir()
        (work / "LOGS").mkdir()
        (work / "LOGS" / "profile1.data").write_text("p")
        (work / "LOGS1").mkdir()
        (work / "LOGS1" / "profile2.data").write_text("p")
        (work / "LOGS2").mkdir()
        (work / "LOGS2" / "profile3.data").write_text("p")
        (work / ".mesa_temp_cache").mkdir()
        (work / "photos1").write_text("photo")
        return work

    def test_parent_process_cleans_up(self, tmp_path, monkeypatch):
        """Parent process runs MESA and removes profiles/photos/cache."""
        work = self._build_work_dir(tmp_path)
        final = tmp_path / "final"

        monkeypatch.setattr(totest.os, "fork", lambda: 1)
        system_calls = []
        monkeypatch.setattr(totest.os, "system", system_calls.append)
        killed = []
        monkeypatch.setattr(totest.os, "kill", lambda pid, sig: killed.append(pid))

        totest.run_mesa("/path/mesa_binary", str(work), str(final))

        assert system_calls[0] == "/path/mesa_binary &> {0}".format(
            os.path.join(work, "out.txt")
        )
        assert "rm LOGS/profile*" in system_calls
        assert "rm LOGS1/profile*" in system_calls
        assert "rm LOGS2/profile*" in system_calls
        assert "rm -rf .mesa_temp_cache" in system_calls
        assert "rm -rf photos*" in system_calls
        assert killed == [1]
        # output moved to final_dir
        assert os.path.isdir(os.path.join(final, "LOGS"))
        assert not work.exists()

    def test_keep_profiles_and_photos(self, tmp_path, monkeypatch):
        """keep_profiles/keep_photos skip the corresponding removals."""
        work = self._build_work_dir(tmp_path)
        final = tmp_path / "final"

        monkeypatch.setattr(totest.os, "fork", lambda: 1)
        system_calls = []
        monkeypatch.setattr(totest.os, "system", system_calls.append)
        monkeypatch.setattr(totest.os, "kill", lambda pid, sig: None)

        totest.run_mesa(
            "/path/mesa_binary",
            str(work),
            str(final),
            keep_profiles=True,
            keep_photos=True,
        )

        assert not any("rm LOGS" in call for call in system_calls)
        assert not any("rm -rf photos*" in call for call in system_calls)
        assert "rm -rf .mesa_temp_cache" in system_calls

    def test_child_process_sleeps_then_copies(self, tmp_path, monkeypatch):
        """Child process waits then copies output before exiting."""
        work = self._build_work_dir(tmp_path)
        final = tmp_path / "final"

        monkeypatch.setattr(totest.os, "fork", lambda: 0)
        monkeypatch.setattr(totest.time, "time", lambda: 5)
        sleeps = []
        monkeypatch.setattr(totest.time, "sleep", sleeps.append)

        with pytest.raises(SystemExit):
            totest.run_mesa(
                "/path/mesa_binary", str(work), str(final), job_end=10, job_before_end=2
            )

        assert sleeps == [3]
        # copyinstead leaves the work dir in place
        assert work.exists()
        assert os.path.isdir(os.path.join(final, "LOGS"))

    def test_equal_dirs_does_not_fork(self, tmp_path, monkeypatch):
        """When work_dir equals final_dir no fork happens and nothing is killed."""
        work = tmp_path / "work"
        work.mkdir()

        fork_calls = []
        monkeypatch.setattr(totest.os, "fork", lambda: fork_calls.append(1) or 999)
        system_calls = []
        monkeypatch.setattr(totest.os, "system", system_calls.append)
        killed = []
        monkeypatch.setattr(totest.os, "kill", lambda pid, sig: killed.append(pid))

        totest.run_mesa("/path/mesa_binary", str(work), str(work))

        assert fork_calls == []
        assert system_calls == [
            "/path/mesa_binary &> {0}".format(os.path.join(work, "out.txt")),
            "rm -rf photos*",
        ]
        assert killed == []

    def test_child_waits_minimum_when_past_job_end(self, tmp_path, monkeypatch):
        """A child with a job_end already in the past waits the minimum 60s."""
        work = tmp_path / "work"
        work.mkdir()
        final = tmp_path / "final"

        monkeypatch.setattr(totest.os, "fork", lambda: 0)
        monkeypatch.setattr(totest.time, "time", lambda: 100)
        sleeps = []
        monkeypatch.setattr(totest.time, "sleep", sleeps.append)

        with pytest.raises(SystemExit):
            totest.run_mesa(
                "/path/mesa_binary",
                str(work),
                str(final),
                job_end=50,
                job_before_end=300,
            )

        assert sleeps == [60]

    def test_negative_fork_result_skips_both_paths(self, tmp_path, monkeypatch):
        """A negative fork result (an error) skips parent and child logic."""
        work = tmp_path / "work"
        work.mkdir()
        final = tmp_path / "final"
        final.mkdir()

        monkeypatch.setattr(totest.os, "fork", lambda: -1)
        system_calls = []
        monkeypatch.setattr(totest.os, "system", system_calls.append)
        killed = []
        monkeypatch.setattr(totest.os, "kill", lambda pid, sig: killed.append(pid))

        totest.run_mesa("/path/mesa_binary", str(work), str(final))

        assert system_calls == []
        assert killed == []


class TestConvertInputColsToMesaCols:
    """Tests for convert_input_cols_to_mesa_cols()."""

    def test_column_renames(self):
        """Convenience columns are derived from the standard names."""
        grid = pandas.DataFrame(
            {
                "initial_period_days": [3.0, 4.0],
                "initial_star_1_mass": [10.0, 12.0],
                "initial_star_2_mass": [8.0, 9.0],
            }
        )
        out = totest.convert_input_cols_to_mesa_cols(grid)
        assert list(out["initial_period_in_days"]) == [3.0, 4.0]
        assert list(out["m1"]) == [10.0, 12.0]
        assert list(out["m2"]) == [8.0, 9.0]

    def test_log_p_and_q_derivations(self):
        """log_p maps to 10**log_p and q*m1 maps to m2."""
        grid = pandas.DataFrame(
            {
                "log_p": [0.5],
                "initial_star_1_mass": [10.0],
                "q": [0.8],
            }
        )
        out = totest.convert_input_cols_to_mesa_cols(grid)
        assert out["initial_period_in_days"].iloc[0] == pytest.approx(10**0.5)
        assert out["m2"].iloc[0] == pytest.approx(8.0)

    def test_without_initial_star_1_mass(self):
        """m2 is still derived when only initial_star_2_mass is present."""
        grid = pandas.DataFrame({"initial_star_2_mass": [8.0], "q": [0.8]})
        out = totest.convert_input_cols_to_mesa_cols(grid)
        assert "m1" not in out.columns
        assert list(out["m2"]) == [8.0]


class TestExtractMesaResults:
    """Tests for extract_mesa_results()."""

    def test_success_returns_dataframe(self, tmp_path, monkeypatch):
        """A successful PSyGrid extraction returns the resulting DataFrame."""
        os.chdir(tmp_path)

        class FakeGrid:
            def create(self, *args, **kwargs):
                return self

            def load(self, path):
                return self

            def get_pandas_initial_final(self):
                return pandas.DataFrame({"m1": [10.0]})

            def close(self):
                return None

        monkeypatch.setattr(totest, "PSyGrid", lambda: FakeGrid())
        df = totest.extract_mesa_results(str(tmp_path / "final"))
        assert not df.empty
        assert list(df["m1"]) == [10.0]

    def test_failure_returns_empty_dataframe(self, tmp_path, monkeypatch):
        """Any exception during extraction yields an empty DataFrame."""
        os.chdir(tmp_path)

        class BrokenGrid:
            def create(self, *args, **kwargs):
                raise RuntimeError("boom")

        monkeypatch.setattr(totest, "PSyGrid", lambda: BrokenGrid())
        df = totest.extract_mesa_results(str(tmp_path / "final"))
        assert df.empty


class TestGetNextGridPoint:
    """Tests for get_next_grid_point()."""

    @pytest.fixture
    def mesa_kwargs(self):
        return {
            "TableData_kwargs": {
                "input_cols": ["m1", "m2"],
                "output_cols": ["period_days"],
                "class_col_name": "class",
            },
            "posydon_dynamic_sampling_kwargs": {
                "in_scaling": ["log", "log"],
                "out_scaling": ["log"],
            },
        }

    @pytest.fixture
    def mock_psy_cris(self, monkeypatch, mesa_kwargs):
        monkeypatch.setattr(
            "posydon.active_learning.psy_cris.utils.parse_inifile",
            lambda path: mesa_kwargs,
        )
        monkeypatch.setattr(
            "posydon.active_learning.psy_cris.utils.get_new_query_points",
            lambda N_new_points, **kwargs: (np.array([[10.5, 8.5]]), np.array([1])),
        )

        class FakeScaler:
            def fit_and_transform(self, data, method=None):
                return data

            def inv_transform(self, data):
                return data

        monkeypatch.setattr("posydon.interpolation.data_scaling.DataScaler", FakeScaler)

    def test_returns_query_points(self, mock_psy_cris, mesa_kwargs):
        """A valid grid returns the inv-transformed query points and classes."""
        grid = pandas.DataFrame(
            {
                "m1": [10.0, 11.0],
                "m2": [8.0, 9.0],
                "period_days": [1.0, 2.0],
            }
        )
        query_points, classes = totest.get_next_grid_point(
            "/path/psycris.ini", grid, n_new_points=1
        )
        assert np.array_equal(query_points, np.array([[10.5, 8.5]]))
        assert np.array_equal(classes, np.array([1]))

    def test_bad_input_columns(self, mock_psy_cris, mesa_kwargs):
        """Missing input columns raise ValueError mentioning the bad columns."""
        grid = pandas.DataFrame({"m1": [10.0], "period_days": [1.0]})
        with pytest.raises(
            ValueError,
            match=r"Bad input columns given: \['m1', 'm2'\],\nMust be in:",
        ):
            totest.get_next_grid_point("/path/psycris.ini", grid, n_new_points=1)

    def test_bad_output_columns(self, mock_psy_cris, mesa_kwargs):
        """Missing output columns raise ValueError mentioning the bad columns."""
        grid = pandas.DataFrame({"m1": [10.0], "m2": [8.0]})
        with pytest.raises(
            ValueError,
            match=r"Bad output columns given: \['period_days'\]\nMust be in:",
        ):
            totest.get_next_grid_point("/path/psycris.ini", grid, n_new_points=1)

    def test_class_output_column_skipped(self, mock_psy_cris, mesa_kwargs, monkeypatch):
        """The class column is excluded from output scaling."""
        mesa_kwargs["TableData_kwargs"]["output_cols"] = ["period_days", "class"]

        class FakeScaler:
            transformed = []

            def fit_and_transform(self, data, method=None):
                FakeScaler.transformed.append(np.asarray(data).tolist())
                return data

            def inv_transform(self, data):
                return data

        monkeypatch.setattr(
            "posydon.interpolation.data_scaling.DataScaler", FakeScaler
        )

        grid = pandas.DataFrame(
            {
                "m1": [10.0, 11.0],
                "m2": [8.0, 9.0],
                "period_days": [1.0, 2.0],
                "class": [0, 1],
            }
        )
        query_points, classes = totest.get_next_grid_point(
            "/path/psycris.ini", grid, n_new_points=1
        )

        assert [0, 1] not in FakeScaler.transformed
        assert [1.0, 2.0] in FakeScaler.transformed
        assert np.array_equal(query_points, np.array([[10.5, 8.5]]))
        assert np.array_equal(classes, np.array([1]))


class TestMainFlow:
    """Tests for main()."""

    def test_main_mpi_not_imported(self, monkeypatch):
        """main raises ImportError when mpi4py was not imported."""
        monkeypatch.setattr(totest, "MPI_IMPORTED", False)
        monkeypatch.setattr(sys, "argv", ["posydon-run-grid"])
        with pytest.raises(ImportError, match="MPI module mpi4py not installed"):
            totest.main()

    def test_main_unrecognized_grid_format(self, monkeypatch):
        """main raises ValueError when the grid format is not recognized."""
        monkeypatch.setattr(totest, "MPI_IMPORTED", True)
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "posydon-run-grid",
                "--mesa-grid", "/path/to/grid.unknown",
                "--grid-type", "fixed",
                "--output-directory", "/out",
                "--temporary-directory", "/tmp",
                "--mesa-binary-executable", "/exe/binary",
                "--mesa-binary-inlist-project", "/inlist/project",
                "--mesa-binary-inlist1", "/inlist/1",
                "--mesa-binary-inlist2", "/inlist/2",
                "--mesa-star1-inlist-project", "/inlist/star1",
                "--mesa-star2-inlist-project", "/inlist/star2",
                "--mesa-star-history-columns", "/cols/history",
                "--mesa-binary-history-columns", "/cols/binary_history",
                "--mesa-profile-columns", "/cols/profile",
            ],
        )
        with pytest.raises(ValueError, match="Grid format not recognized"):
            totest.main()


class TestParseCommandline:
    """Tests for parse_commandline()."""

    def _argv(self, grid_type="fixed", extra=None):
        argv = [
            "posydon-run-grid",
            "--mesa-grid", "/path/grid.csv",
            "--grid-type", grid_type,
            "--output-directory", "/out",
            "--mesa-binary-executable", "/exe/binary",
            "--mesa-binary-inlist-project", "/inlist/project",
            "--mesa-binary-inlist1", "/inlist/1",
            "--mesa-binary-inlist2", "/inlist/2",
            "--mesa-star1-inlist-project", "/inlist/star1",
            "--mesa-star2-inlist-project", "/inlist/star2",
            "--mesa-star-history-columns", "/cols/history",
            "--mesa-binary-history-columns", "/cols/binary_history",
            "--mesa-profile-columns", "/cols/profile",
        ]
        if extra:
            argv.extend(extra)
        return argv

    def test_valid_fixed(self, monkeypatch):
        """A valid fixed-grid invocation returns the parsed arguments."""
        monkeypatch.setattr(sys, "argv", self._argv())
        args = totest.parse_commandline()
        assert args.grid_type == "fixed"
        assert args.grid_point_index is None
        assert args.Niter == 200
        assert args.job_before_end == 300

    def test_valid_dynamic(self, monkeypatch):
        """A valid dynamic-grid invocation requires a psycris inifile."""
        monkeypatch.setattr(
            sys, "argv", self._argv("dynamic", ["--psycris-inifile", "/path/ini"])
        )
        args = totest.parse_commandline()
        assert args.grid_type == "dynamic"
        assert args.psycris_inifile == "/path/ini"

    def test_invalid_grid_type(self, monkeypatch):
        """An unknown grid type makes argparse exit."""
        monkeypatch.setattr(sys, "argv", self._argv("unknown"))
        with pytest.raises(SystemExit):
            totest.parse_commandline()

    def test_dynamic_without_psycris_inifile(self, monkeypatch):
        """A dynamic grid without --psycris-inifile makes argparse exit."""
        monkeypatch.setattr(sys, "argv", self._argv("dynamic"))
        with pytest.raises(SystemExit):
            totest.parse_commandline()


class TestDoRootProcessLogic:
    """Tests for do_root_process_logic()."""

    def _setup(self, tmp_path, monkeypatch):
        monkeypatch.setattr(totest, "size", 2, raising=False)
        sent = []
        received = []

        class FakeReq:
            def wait(self):
                return None

        class FakeComm:
            def isend(self, obj, dest):
                sent.append((obj, dest))
                return FakeReq()

            def recv(self, source, tag, status):
                received.append((source, tag))
                return pandas.DataFrame(
                    {
                        "initial_star_1_mass": [11.0],
                        "initial_star_2_mass": [9.0],
                        "initial_period_days": [4.0],
                    }
                )

        comm = FakeComm()

        class FakeStatus:
            def Get_source(self):
                return 1

        class FakeMPI:
            ANY_SOURCE = "ANY_SOURCE"
            ANY_TAG = "ANY_TAG"

            @staticmethod
            def Status():
                return FakeStatus()

        monkeypatch.setattr(totest, "MPI", FakeMPI, raising=False)
        monkeypatch.setattr(
            "posydon.active_learning.psy_cris.utils.parse_inifile",
            lambda path: {"TableData_kwargs": {"input_cols": ["m1", "m2"]}},
        )
        sampled = []
        monkeypatch.setattr(
            totest,
            "get_next_grid_point",
            lambda *a, **k: sampled.append(k)
            or (np.array([[10.5, 8.5]]), np.array([1])),
        )
        grid = pandas.DataFrame(
            {
                "initial_star_1_mass": [10.0],
                "initial_star_2_mass": [8.0],
                "initial_period_days": [3.0],
            }
        )
        return comm, grid, sent, received, sampled

    def test_root_loop(self, tmp_path, monkeypatch):
        """Root sends points, receives a result, saves it, and samples again."""
        (tmp_path / "output").mkdir()
        args = _make_args(tmp_path, psycris_inifile="/path/ini", Niter=1)
        comm, grid, sent, received, sampled = self._setup(tmp_path, monkeypatch)
        grid.append = lambda other, sort=False: pandas.concat([grid, other], sort=sort)

        totest.do_root_process_logic(
            comm, grid, args, [], [], [], [], [], []
        )

        assert len(received) == 1
        assert len(sent) == 2
        assert len(sampled) == 2
        results = tmp_path / "output" / "grid_results.csv"
        assert results.exists()
        content = results.read_text()
        assert "initial_star_1_mass" in content
        assert "initial_star_2_mass" in content

    def test_root_overwrites_existing_results(self, tmp_path, monkeypatch):
        """An existing grid_results.csv is warned about and replaced."""
        output = tmp_path / "output"
        output.mkdir()
        (output / "grid_results.csv").write_text("old")
        args = _make_args(tmp_path, psycris_inifile="/path/ini", Niter=1)
        comm, grid, sent, received, sampled = self._setup(tmp_path, monkeypatch)
        grid.append = lambda other, sort=False: pandas.concat([grid, other], sort=sort)

        with pytest.warns(OverwriteWarning):
            totest.do_root_process_logic(
                comm, grid, args, [], [], [], [], [], []
            )

        content = (output / "grid_results.csv").read_text()
        assert "old" not in content
        assert "initial_star_1_mass" in content

    def test_root_result_without_optional_columns(self, tmp_path, monkeypatch):
        """Root skips q/log_p derivation when the result lacks those columns."""
        (tmp_path / "output").mkdir()
        args = _make_args(tmp_path, psycris_inifile="/path/ini", Niter=1)
        monkeypatch.setattr(totest, "size", 2, raising=False)
        sent = []

        class FakeReq:
            def wait(self):
                return None

        class FakeComm:
            def isend(self, obj, dest):
                sent.append((obj, dest))
                return FakeReq()

            def recv(self, source, tag, status):
                return pandas.DataFrame({"initial_star_1_mass": [11.0]})

        comm = FakeComm()

        class FakeStatus:
            def Get_source(self):
                return 1

        class FakeMPI:
            ANY_SOURCE = "ANY_SOURCE"
            ANY_TAG = "ANY_TAG"

            @staticmethod
            def Status():
                return FakeStatus()

        monkeypatch.setattr(totest, "MPI", FakeMPI, raising=False)
        monkeypatch.setattr(
            "posydon.active_learning.psy_cris.utils.parse_inifile",
            lambda path: {"TableData_kwargs": {"input_cols": ["m1"]}},
        )
        sampled = []
        monkeypatch.setattr(
            totest,
            "get_next_grid_point",
            lambda *a, **k: sampled.append(k)
            or (np.array([[10.5]]), np.array([1])),
        )
        grid = pandas.DataFrame({"initial_star_1_mass": [10.0]})
        grid.append = lambda other, sort=False: pandas.concat([grid, other], sort=sort)

        totest.do_root_process_logic(
            comm, grid, args, [], [], [], [], [], []
        )

        assert len(sent) == 2
        results = tmp_path / "output" / "grid_results.csv"
        assert results.exists()

    def test_root_append_failure_raises(self, tmp_path, monkeypatch):
        """A failed grid append raises the version-mismatch ValueError."""
        (tmp_path / "output").mkdir()
        args = _make_args(tmp_path, psycris_inifile="/path/ini", Niter=1)
        comm, grid, sent, received, sampled = self._setup(tmp_path, monkeypatch)

        with pytest.raises(
            ValueError, match="version of POSYDON used to make the fixed grid"
        ):
            totest.do_root_process_logic(
                comm, grid, args, [], [], [], [], [], []
            )


class TestDoChildProcessLogic:
    """Tests for do_child_process_logic()."""

    def test_receives_runs_and_sends_result(self, tmp_path, monkeypatch, capsys):
        """Child receives a grid dict, runs the point, and sends the result back."""
        monkeypatch.setattr(totest, "grid_params_binary_controls", ["mdot_scheme"], raising=False)
        monkeypatch.setattr(totest, "grid_params_binary_job", [], raising=False)
        monkeypatch.setattr(totest, "grid_params_star1_binary_controls", [], raising=False)
        monkeypatch.setattr(totest, "grid_params_star1_binary_job", [], raising=False)
        monkeypatch.setattr(totest, "grid_params_star2_binary_controls", [], raising=False)
        monkeypatch.setattr(totest, "grid_params_star2_binary_job", [], raising=False)
        monkeypatch.setattr(totest, "rank", 1, raising=False)

        grid_dict = {"mdot_scheme": 0.5, "m1": 10.0}
        sent = []

        class FakeReq:
            def wait(self):
                return grid_dict

        class FakeComm:
            def irecv(self, source):
                return FakeReq()

            def send(self, obj, dest):
                sent.append((obj, dest))

        result = pandas.DataFrame({"m1": [10.0]})
        monkeypatch.setattr(totest, "run_grid_point", lambda *a, **k: result)
        args = _make_args(tmp_path)

        totest.do_child_process_logic(
            FakeComm(), pandas.DataFrame({"x": [1]}), args,
            star1_formation=True, star2_formation=False,
        )

        assert sent == [(result, 0)]
        out = capsys.readouterr().out
        assert "Process 1 is running with mdot_scheme : 0.5000" in out
        assert "Job completed: From process 1" in out


class TestRunGridPoint:
    """Tests for run_grid_point()."""

    def test_star_formation_steps_and_binary(self, tmp_path, monkeypatch):
        """Both star-formation loops run in order, then the binary."""
        calls = []
        monkeypatch.setattr(
            totest, "create_working_directory",
            lambda *a, **k: ("/tmp/work", "/tmp/final"),
        )
        monkeypatch.setattr(
            totest, "create_star_formation",
            lambda *a, **k: calls.append(("create_star_formation", a, k)),
        )
        monkeypatch.setattr(
            totest, "run_mesa",
            lambda *a, **k: calls.append(("run_mesa", a, k)),
        )
        monkeypatch.setattr(
            totest, "create_binary_inlists",
            lambda *a, **k: calls.append(("create_binary_inlists", a, k)),
        )
        result = pandas.DataFrame({"m1": [10.0]})
        monkeypatch.setattr(totest, "extract_mesa_results", lambda final_dir: result)

        args = _make_main_args(
            tmp_path,
            mesa_star1_inlist_project=["step0", "step1"],
            mesa_star2_inlist_project=["step2"],
            keep_profiles=True,
            keep_photos=False,
            job_end=500,
            job_before_end=50,
        )
        grid_param_dict = {"m1": 10.0, "m2": 8.0, "initial_z": 0.0142}
        binary_controls = {"mdot_scheme": 0.5}
        binary_job = {}
        s1c, s1j, s2c, s2j = {}, {}, {}, {}

        out = totest.run_grid_point(
            pandas.DataFrame({"x": [1]}), True, True, grid_param_dict,
            binary_controls, binary_job, s1c, s1j, s2c, s2j, args,
        )

        names = [c[0] for c in calls]
        assert names == [
            "create_star_formation", "run_mesa",
            "create_star_formation", "run_mesa",
            "create_star_formation", "run_mesa",
            "create_binary_inlists", "run_mesa",
        ]
        assert calls[0][1] == (10.0, "/tmp/work", "step0")
        assert calls[1][2] == {
            "keep_work_dir": True,
            "outfile_name": "out_star1_formation_step0.txt",
        }
        assert calls[2][1] == (10.0, "/tmp/work", "step1", 0.0142, True)
        assert calls[3][2]["outfile_name"] == "out_star1_formation_step1.txt"
        assert calls[4][1] == (8.0, "/tmp/work", "step2")
        assert calls[5][2]["outfile_name"] == "out_star2_formation_step0.txt"
        assert calls[6][1][:6] == (binary_controls, binary_job, s1c, s1j, s2c, s2j)
        assert calls[6][1][6:] == ("/tmp/work", "/inlist/project")
        assert calls[7][1] == ("/exe/binary", "/tmp/work", "/tmp/final")
        assert calls[7][2] == {
            "keep_profiles": True,
            "keep_photos": False,
            "job_end": 500,
            "job_before_end": 50,
        }
        assert out is result

    def test_no_formation_skips_star_loops(self, tmp_path, monkeypatch):
        """With both formations disabled only the binary is run."""
        calls = []
        monkeypatch.setattr(
            totest, "create_working_directory",
            lambda *a, **k: ("/tmp/work", "/tmp/final"),
        )
        monkeypatch.setattr(
            totest, "create_star_formation",
            lambda *a, **k: calls.append("create_star_formation"),
        )
        monkeypatch.setattr(
            totest, "run_mesa",
            lambda *a, **k: calls.append("run_mesa"),
        )
        monkeypatch.setattr(
            totest, "create_binary_inlists",
            lambda *a, **k: calls.append("create_binary_inlists"),
        )
        monkeypatch.setattr(
            totest, "extract_mesa_results", lambda final_dir: pandas.DataFrame()
        )
        args = _make_main_args(tmp_path)

        out = totest.run_grid_point(
            pandas.DataFrame({"x": [1]}), False, False, {"m1": 10.0},
            {}, {}, {}, {}, {}, {}, args,
        )

        assert calls == ["create_binary_inlists", "run_mesa"]
        assert out.empty

    def test_three_step_star1_and_two_step_star2(self, tmp_path, monkeypatch):
        """Star1 formation indexes beyond 1 and star2 index 1 are handled."""
        calls = []
        monkeypatch.setattr(
            totest, "create_working_directory",
            lambda *a, **k: ("/tmp/work", "/tmp/final"),
        )
        monkeypatch.setattr(
            totest, "create_star_formation",
            lambda *a, **k: calls.append(("create_star_formation", a, k)),
        )
        monkeypatch.setattr(
            totest, "run_mesa",
            lambda *a, **k: calls.append(("run_mesa", a, k)),
        )
        monkeypatch.setattr(
            totest, "create_binary_inlists",
            lambda *a, **k: calls.append(("create_binary_inlists", a, k)),
        )
        monkeypatch.setattr(
            totest, "extract_mesa_results", lambda final_dir: pandas.DataFrame()
        )

        args = _make_main_args(
            tmp_path,
            mesa_star1_inlist_project=["s0", "s1", "s2"],
            mesa_star2_inlist_project=["s0", "s1", "s2"],
        )
        grid_param_dict = {"m1": 10.0, "m2": 8.0, "initial_z": 0.0142}

        totest.run_grid_point(
            pandas.DataFrame({"x": [1]}), True, True, grid_param_dict,
            {}, {}, {}, {}, {}, {}, args,
        )

        sf = [c for c in calls if c[0] == "create_star_formation"]
        assert sf[0][1] == (10.0, "/tmp/work", "s0")
        assert sf[1][1] == (10.0, "/tmp/work", "s1", 0.0142, True)
        assert sf[2][1] == (8.0, "/tmp/work", "s0")
        assert sf[3][1] == (8.0, "/tmp/work", "s1", 0.0142, True)
        assert len(sf) == 4
        rm = [c for c in calls if c[0] == "run_mesa"]
        assert len(rm) == 7
        assert rm[2][2]["outfile_name"] == "out_star1_formation_step2.txt"
        assert rm[5][2]["outfile_name"] == "out_star2_formation_step2.txt"
        assert "outfile_name" not in rm[6][2]


class TestMain:
    """Tests for main() beyond the error paths in TestMainFlow."""

    def test_main_fixed_pass(self, tmp_path, monkeypatch, capsys):
        """A fixed grid without a grid-point index only prints the parameters."""
        grid_path = tmp_path / "grid.csv"
        pandas.DataFrame({"initial_mass": [10.0]}).to_csv(grid_path, index=False)
        args = _make_main_args(
            tmp_path,
            mesa_grid=str(grid_path),
            grid_point_index=None,
            mesa_star1_inlist_project=["None"],
            mesa_star2_inlist_project=["None"],
        )
        monkeypatch.setattr(totest, "MPI_IMPORTED", True)
        monkeypatch.setattr(totest, "parse_commandline", lambda: args)
        monkeypatch.setattr(totest.utils, "clean_inlist_file", _fake_clean_inlist)
        run_calls = []
        monkeypatch.setattr(
            totest, "run_grid_point",
            lambda *a, **k: run_calls.append(a) or pandas.DataFrame(),
        )

        totest.main()

        assert run_calls == []
        out = capsys.readouterr().out
        assert "Grid parameters that effect binary_controls:" in out
        assert "Grid parameters that effect star2_binary_job:" in out

    def test_main_fixed_grid_point_writes_and_appends(self, tmp_path, monkeypatch):
        """A fixed grid-point index runs the point and writes/extends results."""
        grid_path = tmp_path / "grid.csv"
        row = {
            "initial_mass": 10.0,
            "initial_z": 0.0142,
            "m1": 10.0,
            "m2": 8.0,
            "initial_star_1_mass": 10.0,
            "initial_star_2_mass": 8.0,
            "initial_period_days": 3.0,
            "initial_period_in_days": 3.0,
            "mdot_scheme": "Kolb",
        }
        pandas.DataFrame([row]).to_csv(grid_path, index=False)
        args = _make_main_args(tmp_path, mesa_grid=str(grid_path))
        monkeypatch.setattr(totest, "MPI_IMPORTED", True)
        monkeypatch.setattr(totest, "parse_commandline", lambda: args)
        monkeypatch.setattr(totest.utils, "clean_inlist_file", _fake_clean_inlist)
        results = [
            pandas.DataFrame({"m1": [10.0], "log_p": [0.5]}),
            pandas.DataFrame({"m1": [11.0], "log_p": [0.6]}),
        ]
        run_calls = []
        monkeypatch.setattr(
            totest, "run_grid_point",
            lambda *a, **k: run_calls.append(a) or results.pop(0),
        )
        (tmp_path / "output").mkdir()

        totest.main()
        totest.main()

        assert len(run_calls) == 2
        csv = tmp_path / "output" / "grid_results.csv"
        assert csv.exists()
        lines = csv.read_text().strip().splitlines()
        assert len(lines) == 3
        assert "10.0" in lines[1]
        assert "11.0" in lines[2]

    def test_main_fixed_empty_mesa_result(self, tmp_path, monkeypatch):
        """An empty mesa result writes no grid_results.csv."""
        grid_path = tmp_path / "grid.csv"
        pandas.DataFrame(
            [{"m1": 10.0, "m2": 8.0, "mdot_scheme": "Kolb"}]
        ).to_csv(grid_path, index=False)
        args = _make_main_args(tmp_path, mesa_grid=str(grid_path))
        monkeypatch.setattr(totest, "MPI_IMPORTED", True)
        monkeypatch.setattr(totest, "parse_commandline", lambda: args)
        monkeypatch.setattr(totest.utils, "clean_inlist_file", _fake_clean_inlist)
        monkeypatch.setattr(totest, "run_grid_point", lambda *a, **k: pandas.DataFrame())
        (tmp_path / "output").mkdir()

        totest.main()

        assert not (tmp_path / "output" / "grid_results.csv").exists()

    def test_main_h5_grid(self, tmp_path, monkeypatch):
        """An .h5 grid is loaded through PSyGrid."""
        args = _make_main_args(tmp_path, mesa_grid="/path/grid.h5", grid_point_index=None)
        monkeypatch.setattr(totest, "MPI_IMPORTED", True)
        monkeypatch.setattr(totest, "parse_commandline", lambda: args)
        monkeypatch.setattr(totest.utils, "clean_inlist_file", _fake_clean_inlist)

        class FakePSyGrid:
            def __init__(self):
                self.loaded = None

            def load(self, path):
                self.loaded = path
                return self

            def get_pandas_initial_final(self):
                return pandas.DataFrame({"m1": [10.0]})

            def close(self):
                return None

        grid = FakePSyGrid()
        monkeypatch.setattr(totest, "PSyGrid", lambda: grid)

        totest.main()

        assert grid.loaded == "/path/grid.h5"

    def test_main_single_star_grid(self, tmp_path, monkeypatch, capsys):
        """A single-star grid runs its formation steps and exits."""
        grid_path = tmp_path / "grid.csv"
        pandas.DataFrame(
            [{"initial_mass": 10.0, "initial_z": 0.0142}]
        ).to_csv(grid_path, index=False)
        args = _make_main_args(
            tmp_path,
            mesa_grid=str(grid_path),
            grid_point_index=0,
            mesa_star1_inlist_project=["step0", "step1"],
            mesa_star2_executable=None,
        )
        monkeypatch.setattr(totest, "MPI_IMPORTED", True)
        monkeypatch.setattr(totest, "parse_commandline", lambda: args)
        monkeypatch.setattr(
            totest.utils,
            "clean_inlist_file",
            lambda path: {
                "&binary_controls": {"m1": None},
                "&binary_job": {},
                "&controls": {"initial_mass": None, "initial_z": None},
                "&star_job": {},
            },
        )
        monkeypatch.setattr(
            totest, "create_working_directory",
            lambda *a, **k: ("/tmp/work", "/tmp/final"),
        )
        formation_calls = []
        run_calls = []
        monkeypatch.setattr(
            totest, "create_star_formation",
            lambda *a, **k: formation_calls.append((a, k)),
        )
        monkeypatch.setattr(
            totest, "run_mesa",
            lambda *a, **k: run_calls.append((a, k)),
        )

        with pytest.raises(SystemExit):
            totest.main()

        assert formation_calls[0][0] == (10.0, "/tmp/work", "step0")
        assert formation_calls[1][0] == (10.0, "/tmp/work", "step1", 0.0142, False)
        assert run_calls[0][1]["keep_work_dir"] is True
        assert run_calls[1][1]["keep_work_dir"] is False
        assert run_calls[0][1]["outfile_name"] == "out_star1_formation_step0.txt"
        assert run_calls[1][1]["outfile_name"] == "out_star1_formation_step1.txt"
        out = capsys.readouterr().out
        assert "You are running a single star grid" in out

    def test_main_single_star_three_steps(self, tmp_path, monkeypatch, capsys):
        """A three-step single-star grid hits the index==1 formation branch."""
        grid_path = tmp_path / "grid.csv"
        pandas.DataFrame(
            [{"initial_mass": 10.0, "initial_z": 0.0142}]
        ).to_csv(grid_path, index=False)
        args = _make_main_args(
            tmp_path,
            mesa_grid=str(grid_path),
            grid_point_index=0,
            mesa_star1_inlist_project=["step0", "step1", "step2"],
            mesa_star2_executable=None,
        )
        monkeypatch.setattr(totest, "MPI_IMPORTED", True)
        monkeypatch.setattr(totest, "parse_commandline", lambda: args)
        monkeypatch.setattr(
            totest.utils,
            "clean_inlist_file",
            lambda path: {
                "&binary_controls": {"m1": None},
                "&binary_job": {},
                "&controls": {"initial_mass": None, "initial_z": None},
                "&star_job": {},
            },
        )
        monkeypatch.setattr(
            totest, "create_working_directory",
            lambda *a, **k: ("/tmp/work", "/tmp/final"),
        )
        formation_calls = []
        run_calls = []
        monkeypatch.setattr(
            totest, "create_star_formation",
            lambda *a, **k: formation_calls.append((a, k)),
        )
        monkeypatch.setattr(
            totest, "run_mesa",
            lambda *a, **k: run_calls.append((a, k)),
        )

        with pytest.raises(SystemExit):
            totest.main()

        assert formation_calls[0][0] == (10.0, "/tmp/work", "step0")
        assert formation_calls[1][0] == (10.0, "/tmp/work", "step1", 0.0142, True)
        assert formation_calls[2][0] == (10.0, "/tmp/work", "step2", 0.0142, False)
        assert run_calls[1][1]["keep_work_dir"] is True
        assert run_calls[2][1]["keep_work_dir"] is False

    def test_main_unknown_grid_type_prints_and_returns(self, tmp_path, monkeypatch, capsys):
        """A grid type other than fixed/dynamic prints and returns without running."""
        grid_path = tmp_path / "grid.csv"
        pandas.DataFrame({"initial_mass": [10.0]}).to_csv(grid_path, index=False)
        args = _make_main_args(
            tmp_path,
            mesa_grid=str(grid_path),
            grid_type="other",
            grid_point_index=None,
            mesa_star1_inlist_project=["None"],
            mesa_star2_inlist_project=["None"],
        )
        monkeypatch.setattr(totest, "MPI_IMPORTED", True)
        monkeypatch.setattr(totest, "parse_commandline", lambda: args)
        monkeypatch.setattr(totest.utils, "clean_inlist_file", _fake_clean_inlist)

        totest.main()

        out = capsys.readouterr().out
        assert "Grid parameters that effect binary_controls:" in out

    def test_main_dynamic_root(self, tmp_path, monkeypatch, capsys):
        """The dynamic root process prints parameters and starts the root logic."""
        grid_path = tmp_path / "grid.csv"
        pandas.DataFrame([{"m1": 10.0, "initial_star_1_mass": 10.0}]).to_csv(
            grid_path, index=False
        )
        args = _make_main_args(
            tmp_path,
            mesa_grid=str(grid_path),
            grid_type="dynamic",
            grid_point_index=None,
            mesa_star1_inlist_project=["None"],
            mesa_star2_inlist_project=["None"],
        )
        monkeypatch.setattr(totest, "MPI_IMPORTED", True)
        monkeypatch.setattr(totest, "parse_commandline", lambda: args)
        monkeypatch.setattr(totest.utils, "clean_inlist_file", _fake_clean_inlist)
        monkeypatch.setattr(
            totest,
            "parse_inifile",
            lambda path: {
                "posydon_dynamic_sampling_kwargs": {"mesa_column_names": ["m1"]}
            },
        )

        class FakeMPI:
            ANY_SOURCE = "ANY_SOURCE"
            ANY_TAG = "ANY_TAG"

            class COMM_WORLD:
                @staticmethod
                def Get_size():
                    return 2

                @staticmethod
                def Get_rank():
                    return 0

            @staticmethod
            def Get_processor_name():
                return "node"

        monkeypatch.setattr(totest, "MPI", FakeMPI, raising=False)
        root_calls = []
        monkeypatch.setattr(totest, "do_root_process_logic", lambda *a: root_calls.append(a))

        totest.main()

        assert len(root_calls) == 1
        call = root_calls[0]
        assert call[0] is FakeMPI.COMM_WORLD
        assert call[2] is args
        assert call[3] == ["m1"]
        assert all(part == [] for part in call[4:])
        out = capsys.readouterr().out
        for msg in [
            "binary_controls",
            "binary_job",
            "star1_binary_controls",
            "star1_binary_job",
            "star2_binary_controls",
            "star2_binary_job",
        ]:
            assert "Grid parameters that effect {0}:".format(msg) in out

    def test_main_dynamic_child(self, tmp_path, monkeypatch):
        """The dynamic child loops on do_child_process_logic until exit."""
        grid_path = tmp_path / "grid.csv"
        pandas.DataFrame([{"m1": 10.0}]).to_csv(grid_path, index=False)
        args = _make_main_args(
            tmp_path,
            mesa_grid=str(grid_path),
            grid_type="dynamic",
            grid_point_index=None,
            mesa_star1_inlist_project=["None"],
            mesa_star2_inlist_project=["None"],
        )
        monkeypatch.setattr(totest, "MPI_IMPORTED", True)
        monkeypatch.setattr(totest, "parse_commandline", lambda: args)
        monkeypatch.setattr(totest.utils, "clean_inlist_file", _fake_clean_inlist)
        monkeypatch.setattr(
            totest,
            "parse_inifile",
            lambda path: {
                "posydon_dynamic_sampling_kwargs": {"mesa_column_names": ["m1"]}
            },
        )

        class FakeMPI:
            ANY_SOURCE = "ANY_SOURCE"
            ANY_TAG = "ANY_TAG"

            class COMM_WORLD:
                @staticmethod
                def Get_size():
                    return 2

                @staticmethod
                def Get_rank():
                    return 1

            @staticmethod
            def Get_processor_name():
                return "node"

        monkeypatch.setattr(totest, "MPI", FakeMPI, raising=False)
        child_calls = []

        def fake_child(*a, **k):
            child_calls.append((a, k))
            raise SystemExit(0)

        monkeypatch.setattr(totest, "do_child_process_logic", fake_child)

        with pytest.raises(SystemExit):
            totest.main()

        assert len(child_calls) == 1
        assert child_calls[0][1] == {
            "star1_formation": False,
            "star2_formation": False,
        }

    def test_main_dynamic_grid_point_index_error(self, tmp_path, monkeypatch):
        """Requesting a specific point with multiple MPI tasks raises ValueError."""
        grid_path = tmp_path / "grid.csv"
        pandas.DataFrame([{"m1": 10.0}]).to_csv(grid_path, index=False)
        args = _make_main_args(
            tmp_path, mesa_grid=str(grid_path), grid_type="dynamic", grid_point_index=0
        )
        monkeypatch.setattr(totest, "MPI_IMPORTED", True)
        monkeypatch.setattr(totest, "parse_commandline", lambda: args)

        class FakeMPI:
            ANY_SOURCE = "ANY_SOURCE"
            ANY_TAG = "ANY_TAG"

            class COMM_WORLD:
                @staticmethod
                def Get_size():
                    return 2

                @staticmethod
                def Get_rank():
                    return 1

            @staticmethod
            def Get_processor_name():
                return "node"

        monkeypatch.setattr(totest, "MPI", FakeMPI, raising=False)

        with pytest.raises(
            ValueError,
            match="You can not have multiple MPI tasks and request to run a "
            "specific grid point.Please make sure you cal with mpirun -np 1.",
        ):
            totest.main()
