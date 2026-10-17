"""Unit tests of posydon/grids/setup/steps.py
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

# import other needed code for the tests, which is not already imported in the
# module you like to test
import configparser
from inspect import isroutine
from pathlib import Path

import pandas
import pytest

# import the module which will be tested
import posydon.grids.setup.steps as totest
from posydon.grids.setup.config import ConfigError, GridConfig
from posydon.grids.setup.resolve import ResolvedScenario
from posydon.grids.setup.scenario import ScenarioFiles

LEGACY_INI = """\
[slurm]
job_array=True
partition=p
account=a
walltime='1:00:00'
email=e@e

[mesa_inlists]
binary_controls_mesa_defaults = {d}/binary_controls.defaults

[mesa_extras]
binary_run = {d}/binary_run.f

[run_parameters]
grid = {d}/grid.csv
"""


def write(path, text=""):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


@pytest.fixture
def grid(tmp_path, monkeypatch):
    """A complete new-format grid without a scenario."""
    d = tmp_path / "in"
    for name in ("binary_controls", "binary_job"):
        write(d / (name + ".defaults"), "history_interval = 5\n")
    write(d / "controls.defaults", "mesh_delta_coeff = 1d0\n")
    write(d / "star_job.defaults", "pgstar_flag = .false.\n")
    for name in ("history", "binary_history", "profile"):
        write(d / (name + "_columns.list"), name)
    write(d / "binary_run.f", "! binary_run")
    write(d / "grid.csv", "m1,m2,initial_period_in_days\n10,8,3\n12,9,4\n")
    ini = write(tmp_path / "grid.ini", """\
[grid]
type = HMS-HMS
points = {d}/grid.csv
run_directory = {run}

[mesa_inlists]
binary_controls = {d}/binary_controls.defaults
binary_job = {d}/binary_job.defaults
star1_controls = {d}/controls.defaults
star1_job = {d}/star_job.defaults
star2_controls = {d}/controls.defaults
star2_job = {d}/star_job.defaults

[posydon_inlists]
history_columns = {d}/history_columns.list
binary_history_columns = {d}/binary_history_columns.list
profile_columns = {d}/profile_columns.list

[extras]
binary_run = {d}/binary_run.f
""".format(d=d, run=tmp_path / "run"))
    # a stub MESA: the files fill_defaults takes for what [extras] leaves out
    from posydon.grids.setup.defaults import MESA_EXTRAS
    for rel in MESA_EXTRAS.values():
        write(tmp_path / "mesa" / rel, "! MESA " + rel + "\n")
    monkeypatch.setenv("MESA_DIR", str(tmp_path / "mesa"))
    monkeypatch.setenv("MESASDK_ROOT", "/sdk")
    monkeypatch.setattr(totest, "find_run_grid", lambda: "/bin/run-grid")
    return {"ini": ini, "dir": d, "run": tmp_path / "run"}


class TestElements:
    def test_dir(self):
        assert {'SUBMISSION_TYPES', 'Prepared',
                'check_config', 'chosen_extras', 'explain_grid',
                'grid_parameters', 'is_new_format', 'prepare',
                'read_any_config', 'read_grid_points',
                'setup_grid'} <= set(dir(totest))

    @pytest.mark.parametrize("name", [
        'is_new_format', 'read_any_config', 'check_config',
        'read_grid_points', 'grid_parameters', 'chosen_extras', 'prepare',
        'explain_grid', 'setup_grid'])
    def test_instance_functions(self, name):
        assert isroutine(getattr(totest, name))


class TestSetupGrid:
    def test_writes_the_run_directory(self, grid):
        run = totest.setup_grid(grid["ini"])
        assert run == grid["run"]
        names = sorted(p.relative_to(run).as_posix()
                       for p in run.rglob("*") if p.is_file())
        assert "binary/inlist_project" in names
        assert "binary/src/binary_run.f" in names
        # the ini names one extra; MESA's makefiles are still filled in
        assert "binary/make/makefile_binary" in names
        assert (run / "mk").read_text().count("make -f") == 3
        assert (run / "binary" / "src" / "binary_run.f").read_text() == \
            "! binary_run"
        assert "column_lists/profile_columns.list" in names
        assert "grid_command.sh" in names
        record = configparser.ConfigParser(interpolation=None)
        record.read(run / "setup_record.ini")
        assert record["setup"]["grid_type"] == "HMS-HMS"
        assert "for SLURM_ARRAY_TASK_ID in 0 1 ; do python /bin/run-grid" \
            in (run / "grid_command.sh").read_text()

    def test_needs_the_environment(self, grid, monkeypatch):
        monkeypatch.delenv("MESASDK_ROOT")
        with pytest.raises(ValueError, match="MESASDK_ROOT must be defined"):
            totest.setup_grid(grid["ini"])

    def test_explain(self, grid, capsys):
        text = totest.explain_grid(grid["ini"])
        assert text.startswith("binary_controls\n")
        assert not grid["run"].exists()

    def test_explain_with_a_scenario(self, grid, monkeypatch, tmp_path):
        d = grid["dir"]
        files = ScenarioFiles()
        files.inlists.add("star1_controls", write(
            tmp_path / "repo" / "inlist1", "&controls\n x = 1\n/\n"))
        monkeypatch.setattr(totest, "resolve_scenario", lambda config:
                            ResolvedScenario(files=files,
                                             repo=tmp_path / "repo"))
        text = totest.explain_grid(grid["ini"])
        assert "  x = 1  [posydon_inlists] inlist1" in text
        assert d.is_dir()


class TestPrepare:
    def test_replaced_column_lists_warn(self, grid, monkeypatch, tmp_path):
        files = ScenarioFiles()
        files.inlists.whole_files["history_columns"] = write(
            tmp_path / "other.list", "other")
        monkeypatch.setattr(totest, "resolve_scenario", lambda config:
                            ResolvedScenario(files=files, repo=tmp_path))
        with pytest.warns(match="history_columns in the ini file is not "
                          "used"):
            ready = totest.prepare(grid["ini"])
        assert ready.columns["history_columns"] == tmp_path / "other.list"
        assert ready.scenario is files

    def test_missing_column_lists(self, grid):
        text = grid["ini"].read_text()
        grid["ini"].write_text(text.split("[posydon_inlists]")[0] +
                               "[extras]" + text.split("[extras]")[1])
        with pytest.raises(ConfigError, match="no file for the column "
                           "lists history_columns, binary_history_columns, "
                           "profile_columns"):
            totest.prepare(grid["ini"])

    def test_mpi_grid_path(self, grid, monkeypatch):
        grid["ini"].write_text(grid["ini"].read_text().replace(
            "type = HMS-HMS", "type = HMS-HMS\nmode = dynamic\n"
            "psycris_inifile = psycris.ini") +
            "\n[submission]\nbackend = slurm_mpi\n")
        monkeypatch.setattr(totest, "grid_parameters",
                            lambda config, df: ["m1"])
        ready = totest.prepare(grid["ini"])
        assert ready.config.grid.points == grid["dir"] / "grid.csv"
        assert ready.parameters == ["m1"]

    def test_needs_mesa_dir(self, grid, monkeypatch):
        monkeypatch.delenv("MESA_DIR")
        with pytest.raises(ValueError, match="MESA_DIR must be defined"):
            totest.prepare(grid["ini"])


class TestReadAnyConfig:
    def test_is_new_format(self, grid, tmp_path):
        assert totest.is_new_format(grid["ini"])
        assert not totest.is_new_format(write(tmp_path / "old.ini",
                                              "[slurm]\n"))

    def test_new_format_refuses_options(self, grid):
        with pytest.raises(ValueError, match="--grid-type, "
                           "--submission-type only apply to legacy"):
            totest.read_any_config(grid["ini"], grid_type="fixed",
                                   submission_type="shell")

    def test_new_format_refuses_one_option(self, grid):
        with pytest.raises(ValueError, match="--grid-type only applies to "
                           "legacy"):
            totest.read_any_config(grid["ini"], grid_type="fixed")

    def test_new_format_overrides(self, grid):
        config = totest.read_any_config(
            grid["ini"], grid_type="dynamic", run_directory="/elsewhere",
            submission_type="shell", overrides=True)
        assert (config.grid.mode, config.grid.run_directory,
                config.submission.backend) == ("dynamic", Path("/elsewhere"),
                                               "shell")

    @pytest.mark.parametrize("job_array, backend", [
        (True, "slurm_array"), (False, "slurm_mpi")])
    def test_slurm_override(self, grid, job_array, backend):
        grid["ini"].write_text(grid["ini"].read_text() + "\n[submission]\n"
                               "job_array = {}\n".format(job_array))
        config = totest.read_any_config(grid["ini"],
                                        submission_type="slurm",
                                        overrides=True)
        assert config.submission.backend == backend

    def test_legacy(self, grid, tmp_path):
        ini = write(tmp_path / "old.ini", LEGACY_INI.format(d=grid["dir"]))
        with pytest.raises(ValueError, match="--grid-type is needed"):
            totest.read_any_config(ini)
        config = totest.read_any_config(ini, grid_type="fixed")
        assert config.source_format == "legacy"
        assert config.submission.backend == "shell"
        assert config.grid.run_directory == Path.cwd()


class TestCheckConfig:
    def config(self, tmp_path, **grid):
        config = GridConfig()
        config.grid.points = write(tmp_path / "grid.csv", "m1\n1\n")
        for key, value in grid.items():
            setattr(config.grid, key, value)
        return config

    def test_ok_and_notes(self, tmp_path):
        config = self.config(tmp_path)
        config.notes = ["legacy ini format; see the docs", "MS-MS is old"]
        with pytest.warns(match="MS-MS is old") as caught:
            totest.check_config(config)
        assert len(caught) == 1

    def test_missing_points(self, tmp_path):
        with pytest.raises(ValueError, match="Supplied grid does not exist"):
            totest.check_config(self.config(tmp_path,
                                            points=tmp_path / "nope"))

    def test_dynamic_needs_psycris(self, tmp_path):
        with pytest.raises(ValueError, match="psycris inifile"):
            totest.check_config(self.config(tmp_path, mode="dynamic"))

    def test_fixed_needs_one_point_per_task(self, tmp_path):
        config = self.config(tmp_path)
        config.submission.backend = "slurm_mpi"
        with pytest.raises(ConfigError, match=r"^\[submission\] backend: a "
                           "fixed grid runs one point per task: use "
                           "slurm_array"):
            totest.check_config(config)
        config.source_format = "legacy"
        with pytest.raises(ConfigError, match=r"^\[slurm\] job_array: a "
                           "fixed grid needs job_array = True"):
            totest.check_config(config)

    def test_dynamic_refuses_a_job_array(self, tmp_path):
        config = self.config(tmp_path, mode="dynamic",
                             psycris_inifile=tmp_path / "p.ini")
        with pytest.raises(ConfigError, match=r"^\[submission\] backend: a "
                           "dynamic grid runs as one MPI job: use slurm_mpi"):
            totest.check_config(config)
        config.source_format = "legacy"
        with pytest.raises(ConfigError, match=r"^\[slurm\] job_array: a "
                           "dynamic grid runs as one MPI job: set job_array "
                           "= False"):
            totest.check_config(config)

    def test_unknown_system_type(self, tmp_path):
        config = self.config(tmp_path)
        config.scenario.system_type = "XX"
        with pytest.raises(ConfigError, match="unknown system type 'XX'"):
            totest.check_config(config)


class FakePSyGrid:
    created = []

    def create(self, folder, h5_path, slim):
        self.created.append((folder, h5_path, slim))

    def load(self, path):
        self.path = path

    def get_pandas_initial_final(self):
        return pandas.DataFrame({"m1": [1.0]})

    def close(self):
        pass


class TestReadGridPoints:
    def grid(self, points):
        config = GridConfig()
        config.grid.points = points
        return config.grid

    def test_csv(self, tmp_path):
        path = write(tmp_path / "grid.csv", "m1\n1\n2\n")
        df, mpi = totest.read_grid_points(self.grid(path))
        assert list(df["m1"]) == [1, 2] and mpi == str(path)

    def test_psygrid_and_folder(self, tmp_path, monkeypatch):
        import posydon.grids.psygrid as psygrid
        monkeypatch.setattr(psygrid, "PSyGrid", FakePSyGrid)
        monkeypatch.chdir(tmp_path)
        df, mpi = totest.read_grid_points(self.grid(tmp_path / "g.h5"))
        assert mpi == str(tmp_path / "g.h5") and len(df) == 1
        folder = tmp_path / "runs"
        folder.mkdir()
        df, mpi = totest.read_grid_points(self.grid(folder))
        assert mpi == str(tmp_path / "fixed_grid_results.h5")
        assert FakePSyGrid.created[-1] == (str(folder), mpi, True)

    def test_unknown(self, tmp_path, monkeypatch):
        import posydon.grids.psygrid as psygrid
        monkeypatch.setattr(psygrid, "PSyGrid", FakePSyGrid)
        with pytest.raises(ValueError, match="Grid format not recognized"):
            totest.read_grid_points(self.grid(tmp_path / "grid.txt"))


class TestGridParameters:
    def test_fixed(self):
        config = GridConfig()
        df = pandas.DataFrame({"m1": [1], "m2": [2]})
        assert list(totest.grid_parameters(config, df)) == ["m1", "m2"]

    def test_dynamic(self, monkeypatch):
        import posydon.active_learning.psy_cris.utils as utils
        monkeypatch.setattr(utils, "parse_inifile", lambda path: {
            "posydon_dynamic_sampling_kwargs": {
                "mesa_column_names": ["m1", "initial_period_in_days"]}})
        config = GridConfig()
        config.grid.mode = "dynamic"
        config.grid.psycris_inifile = Path("p.ini")
        assert totest.grid_parameters(config, None) == [
            "m1", "initial_period_in_days"]


class TestChosenExtras:
    def test_lower_layers_dropped_with_a_warning(self):
        config = GridConfig(source_format="legacy")
        config.extras = {"mesa_binary_extras": Path("m/run_binary_extras.f")}
        scenario = ScenarioFiles()
        scenario.extras = {"posydon_binary_extras":
                           Path("p/run_binary_extras.f")}
        with pytest.warns(match="mesa_binary_extras is being set to None"):
            kept = totest.chosen_extras(config, scenario)
        assert kept == {"posydon_binary_extras":
                        Path("p/run_binary_extras.f")}

    def test_new_format_names_the_extras_key(self):
        config = GridConfig()
        config.extras = {"mesa_binary_extras": Path("m/run_binary_extras.f")}
        scenario = ScenarioFiles()
        scenario.extras = {"posydon_binary_extras":
                           Path("p/run_binary_extras.f")}
        with pytest.warns(match=r"\[extras\] mesa_binary_extras is not "
                          "used: a higher layer's file goes to the same "
                          "place"):
            totest.chosen_extras(config, scenario)

    def test_filled_in_mesa_extras_drop_silently(self, recwarn):
        config = GridConfig(default_extras=("mesa_binary_extras",))
        config.extras = {"mesa_binary_extras": Path("m/run_binary_extras.f")}
        scenario = ScenarioFiles()
        scenario.extras = {"posydon_binary_extras":
                           Path("p/run_binary_extras.f")}
        kept = totest.chosen_extras(config, scenario)
        assert list(kept) == ["posydon_binary_extras"]
        assert len(recwarn) == 0

    def test_scenario_files_a_user_file_replaces_drop_silently(
            self, recwarn):
        config = GridConfig()
        config.inlists["user_inlists"].whole_files["run_binary_extras"] = \
            Path("u/mine.f")
        scenario = ScenarioFiles()
        scenario.extras = {"posydon_binary_extras":
                           Path("p/run_binary_extras.f")}
        kept = totest.chosen_extras(config, scenario)
        assert kept == {"user_binary_extras": Path("u/mine.f")}
        assert len(recwarn) == 0

    def test_whole_file_star_extras_go_everywhere_a_star_runs(self):
        config = GridConfig()
        config.inlists["user_inlists"].whole_files["run_star_extras"] = \
            Path("u/mine.f")
        scenario = ScenarioFiles()
        scenario.extras = {
            "posydon_star_binary_extras": Path("p/run_star_extras.f"),
            "mesa_star1_extras": Path("p/run_star_extras.f")}
        kept = totest.chosen_extras(config, scenario)
        assert kept == {"user_star_binary_extras": Path("u/mine.f"),
                        "user_star1_extras": Path("u/mine.f"),
                        "user_star2_extras": Path("u/mine.f")}

    def test_no_scenario(self):
        config = GridConfig()
        config.extras = {"binary_run": Path("b.f")}
        assert totest.chosen_extras(config, None) == config.extras
