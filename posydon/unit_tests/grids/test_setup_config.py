"""Unit tests of posydon/grids/setup/config.py
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

# import other needed code for the tests, which is not already imported in the
# module you like to test
from dataclasses import is_dataclass
from inspect import isclass, isroutine

import pytest

# import the module which will be tested
import posydon.grids.setup.config as totest

SPEC_EXAMPLE = """
[grid]
type = CO-HeMS            ; key in the grid-type registry
mode = fixed              ; fixed | dynamic
points = ./grid.csv
run_directory = ./run

[scenario]                ; selects POSYDON's inlists
release = v2.0
rerun = dedt_hepulse
inlist_repo = ~/inlists   ; read as it is

[user_inlists]            ; your MESA inlists, on top of the release
star1_controls = ./my_inlist1
binary = ./my_inlist_project

[output]
history_interval = 1
star1 = history, final_model
star2 = none
binary = history

[submission]
backend = slurm_array     ; shell | slurm_array | slurm_mpi
partition = public-cpu
walltime = 2-00:00:00
mem_per_cpu = 4G
"""

MINIMAL = "[grid]\ntype = HMS-HMS\npoints = grid.csv\n"


@pytest.fixture
def write_ini(tmp_path):
    def write(text):
        path = tmp_path / "grid.ini"
        path.write_text(text)
        return path
    return write


class TestElements:
    def test_dir(self):
        elements = {'BACKENDS', 'ConfigError', 'GRID_MODES',
                    'GridConfig', 'GridSection', 'INLIST_SHORTHANDS',
                    'INLIST_TARGETS', 'InlistFiles', 'OUTPUT_OPTIONS',
                    'OutputSection', 'Path', 'GIT_KEYS', 'INI_LAYERS', 'SECTIONS',
                    'ScenarioSection', 'SetupError', 'SubmissionSection',
                    'WHOLE_FILES',
                    'COLUMN_FILES', 'SOURCE_FILES', 'EXTRAS_KEYS',
                    '_read_extras',
                    '__authors__', '__builtins__', '__cached__', '__doc__',
                    '__file__', '__loader__', '__name__', '__package__',
                    '__spec__', '_bool', '_choice', '_int', '_list', '_path',
                    '_read_grid', '_read_inlist_files', '_read_output',
                    '_read_scenario', '_read_submission', 'configparser',
                    'dataclass', 'field', 're', 'read_config'}
        assert set(dir(totest)) == elements

    def test_instance_constants(self):
        assert totest.GRID_MODES == ("fixed", "dynamic")
        assert totest.BACKENDS == ("shell", "slurm_array", "slurm_mpi")
        assert set(totest.INLIST_SHORTHANDS) == {"binary", "star1", "star2"}
        assert "profile_columns" in totest.WHOLE_FILES

    @pytest.mark.parametrize("name", ['GridSection', 'ScenarioSection',
                                      'InlistFiles', 'OutputSection',
                                      'SubmissionSection', 'GridConfig'])
    def test_instance_dataclasses(self, name):
        assert is_dataclass(getattr(totest, name))

    def test_instance_ConfigError(self):
        assert isclass(totest.ConfigError)
        assert issubclass(totest.ConfigError, ValueError)

    def test_instance_read_config(self):
        assert isroutine(totest.read_config)


class TestConfigError:
    def test_message_with_key(self):
        error = totest.ConfigError("grid", "type", "missing")
        assert str(error) == "[grid] type: missing"
        assert (error.section, error.key) == ("grid", "type")

    def test_message_without_key(self):
        assert str(totest.ConfigError("grid", None, "missing")) == \
            "[grid]: missing"


class TestDataclasses:
    def test_defaults(self):
        config = totest.GridConfig()
        assert config.grid.mode == "fixed"
        assert config.scenario.is_set is False
        assert set(config.inlists) == {"mesa_inlists", "posydon_inlists",
                                       "user_inlists"}
        assert config.output.star1 == ("history", "final_model")
        assert config.submission.backend == "shell"
        assert config.source_format == "new"
        assert config.notes == []
        assert config.legacy == {}

    def test_scenario_is_set(self):
        assert totest.ScenarioSection(release="v2.0").is_set
        # a legacy scenario names a commit
        assert totest.ScenarioSection(commit="18710e9").is_set
        assert not totest.ScenarioSection(source="user").is_set
        assert not totest.ScenarioSection(inlist_repo="repo").is_set

    def test_inlist_files_add(self):
        files = totest.InlistFiles()
        files.add("star1_controls", "a")
        files.add("star1_controls", "b")
        assert files.targets == {"star1_controls": ["a", "b"]}


class TestReadConfig:
    def test_spec_example(self, write_ini, tmp_path):
        config = totest.read_config(write_ini(SPEC_EXAMPLE))
        assert config.grid.type == "CO-HeMS"
        assert config.grid.points == tmp_path / "grid.csv"
        assert config.grid.run_directory == tmp_path / "run"
        assert (config.scenario.release, config.scenario.rerun) == \
            ("v2.0", "dedt_hepulse")
        assert config.scenario.is_set
        user = config.inlists["user_inlists"].targets
        assert user == {"star1_controls": [tmp_path / "my_inlist1"],
                        "binary_controls": [tmp_path / "my_inlist_project"],
                        "binary_job": [tmp_path / "my_inlist_project"]}
        assert config.output.star1 == ("history", "final_model")
        assert config.output.star2 == ()
        assert config.output.binary_history is True
        assert config.submission.backend == "slurm_array"
        assert config.submission.mem_per_cpu == "4G"

    def test_minimal_defaults(self, write_ini, tmp_path):
        config = totest.read_config(write_ini(MINIMAL))
        assert config.grid.run_directory == tmp_path
        assert config.scenario.is_set is False
        assert config.output == totest.OutputSection()
        assert config.submission == totest.SubmissionSection()
        assert config.extras == {}

    def test_missing_file(self, tmp_path):
        with pytest.raises(totest.ConfigError, match="not found"):
            totest.read_config(tmp_path / "nope.ini")

    def test_unknown_section(self, write_ini):
        with pytest.raises(totest.ConfigError,
                           match=r"\[mesa_extras\]: unknown section"):
            totest.read_config(write_ini(MINIMAL + "[mesa_extras]\na = 1\n"))

    def test_missing_grid_section(self, write_ini):
        with pytest.raises(totest.ConfigError, match=r"\[grid\]: missing"):
            totest.read_config(write_ini("[output]\nhistory_interval = 1\n"))

    def test_absolute_and_relative_paths(self, write_ini, tmp_path):
        config = totest.read_config(write_ini(
            "[grid]\ntype = HMS-HMS\npoints = /data/grid.csv\n"
            "[extras]\nstar_run = src/run.f\n"))
        assert str(config.grid.points) == "/data/grid.csv"
        assert config.extras == {"star_run": tmp_path / "src" / "run.f"}

    @pytest.mark.parametrize("key", ["makefile_binary", "makefile_star",
                                     "binary_run", "star_run"])
    def test_extras_takes_mesa_work_files(self, write_ini, key):
        config = totest.read_config(write_ini(
            MINIMAL + "[extras]\n{} = x\n".format(key)))
        assert list(config.extras) == [key]

    @pytest.mark.parametrize("key", ["user_star1_extras", "mesa_binary_extras",
                                     "posydon_star_binary_extras",
                                     "run_star_extras", "my_extras"])
    def test_extras_refuses_sources(self, write_ini, key):
        """Each Fortran source has one place: [user_inlists]."""
        with pytest.raises(totest.ConfigError, match=r"\[extras\] {}: not "
                           r"an \[extras\] key: \[extras\] takes MESA's "
                           r"work files \(makefile_binary, makefile_star, "
                           r"binary_run, star_run\); Fortran sources go in "
                           r"\[user_inlists\]".format(key)):
            totest.read_config(write_ini(
                MINIMAL + "[extras]\n{} = x\n".format(key)))

    def test_inline_comments(self, write_ini):
        config = totest.read_config(write_ini(
            "[grid]\ntype = HMS-HMS   ; the grid type\npoints = g.csv\n"))
        assert config.grid.type == "HMS-HMS"

    def test_interpolation(self, write_ini, tmp_path):
        config = totest.read_config(write_ini(
            MINIMAL + "[user_inlists]\nstar1_job = ${grid:points}.inlist\n"))
        assert config.inlists["user_inlists"].targets["star1_job"] == \
            [tmp_path / "grid.csv.inlist"]


class TestGridSection:
    @pytest.mark.parametrize("text, message", [
        ("[grid]\npoints = g.csv\n", r"\[grid\] type: missing"),
        ("[grid]\ntype = HMS-HMS\n", r"\[grid\] points: missing"),
        ("[grid]\ntype = MS-MS\npoints = g.csv\n",
         r"\[grid\] type: unknown grid type 'MS-MS'"),
        (MINIMAL + "mode = sample\n", r"\[grid\] mode: unknown value"),
        (MINIMAL + "colour = red\n", r"\[grid\] colour: unknown key"),
        (MINIMAL + "keep_photos = maybe\n",
         r"\[grid\] keep_photos: expected true or false"),
        (MINIMAL + "mode = dynamic\n",
         r"\[grid\] psycris_inifile: needed for mode = dynamic"),
    ])
    def test_errors(self, write_ini, text, message):
        with pytest.raises(totest.ConfigError, match=message):
            totest.read_config(write_ini(text))

    def test_optional_keys(self, write_ini, tmp_path):
        config = totest.read_config(write_ini(
            MINIMAL + "mode = dynamic\npsycris_inifile = p.ini\n"
            "zams = z.data\nkeep_profiles = yes\nkeep_photos = off\n"))
        grid = config.grid
        assert grid.mode == "dynamic"
        assert grid.psycris_inifile == tmp_path / "p.ini"
        assert grid.zams_star1 == grid.zams_star2 == tmp_path / "z.data"
        assert (grid.keep_profiles, grid.keep_photos) == (True, False)

    def test_zams_per_star(self, write_ini, tmp_path):
        config = totest.read_config(write_ini(
            MINIMAL + "zams_star2 = z2.data\n"))
        assert config.grid.zams_star1 is None
        assert config.grid.zams_star2 == tmp_path / "z2.data"


class TestScenarioSection:
    def test_all_keys(self, write_ini, tmp_path):
        config = totest.read_config(write_ini(
            MINIMAL + "[scenario]\nrelease = v2.0\nrerun = hot\n"
            "mesa = r11701\ninlist_repo = repo\n"))
        scenario = config.scenario
        assert (scenario.release, scenario.rerun, scenario.mesa) == \
            ("v2.0", "hot", "r11701")
        assert scenario.inlist_repo == tmp_path / "repo"
        # what a legacy scenario uses to check the repository out
        assert (scenario.branch, scenario.commit) == (None, None)

    def test_release_alone(self, write_ini):
        """No ref or commit: the clone is read as it is."""
        config = totest.read_config(write_ini(
            MINIMAL + "[scenario]\nrelease = v2.0\n"))
        assert config.scenario.is_set
        assert config.scenario.inlist_repo is None

    @pytest.mark.parametrize("text, message", [
        ("source = path\n", r"\[scenario\] source: no longer used: the "
         "setup reads inlist_repo as it is, and the release names the "
         "version. Remove this key, and update the clone yourself"),
        ("release = v2.0\nref = v2_sources\n",
         r"\[scenario\] ref: no longer used"),
        ("release = v2.0\ncommit = 8482d7d\n",
         r"\[scenario\] commit: no longer used"),
        ("branch = main\n", r"\[scenario\] branch: no longer used"),
        ("rerun = dedt\n", r"\[scenario\] rerun: a rerun needs a release"),
        ("inlist_repo = repo\n", r"\[scenario\] release: missing"),
        ("mesa = r11701\n", r"\[scenario\] release: missing"),
        ("system = HMS-HMS\n", r"\[scenario\] system: unknown key"),
    ])
    def test_errors(self, write_ini, text, message):
        with pytest.raises(totest.ConfigError, match=message):
            totest.read_config(write_ini(MINIMAL + "[scenario]\n" + text))


class TestInlistSections:
    def test_posydon_inlists_is_refused(self, write_ini):
        with pytest.raises(totest.ConfigError, match=r"\[posydon_inlists\]:"
                           " no longer used: POSYDON's inlists come from "
                           r"\[scenario\]"):
            totest.read_config(write_ini(
                MINIMAL + "[posydon_inlists]\nstar1 = x\n"))

    @pytest.mark.parametrize("key", ["run_star_extras", "run_binary_extras",
                                     "run_star1_extras", "run_star2_extras"])
    def test_sources_only_in_user_inlists(self, write_ini, tmp_path, key):
        config = totest.read_config(write_ini(
            MINIMAL + "[user_inlists]\n{} = mine.f\n".format(key)))
        assert config.inlists["user_inlists"].whole_files == {
            key: tmp_path / "mine.f"}
        with pytest.raises(totest.ConfigError, match=r"\[mesa_inlists\] "
                           "{}: Fortran sources go in ".format(key) +
                           r"\[user_inlists\]; MESA's own come from "
                           r"\$MESA_DIR"):
            totest.read_config(write_ini(
                MINIMAL + "[mesa_inlists]\n{} = mine.f\n".format(key)))

    @pytest.mark.parametrize("layer", ["mesa_inlists",
                                       "user_inlists"])
    def test_every_layer_is_read(self, write_ini, tmp_path, layer):
        config = totest.read_config(write_ini(
            MINIMAL + "[{}]\nbinary_job = a, b\n".format(layer)))
        assert config.inlists[layer].targets == {
            "binary_job": [tmp_path / "a", tmp_path / "b"]}

    def test_shorthands(self, write_ini, tmp_path):
        config = totest.read_config(write_ini(
            MINIMAL + "[user_inlists]\nstar2 = i2\n"))
        assert config.inlists["user_inlists"].targets == {
            "star2_controls": [tmp_path / "i2"],
            "star2_job": [tmp_path / "i2"]}

    def test_formation_steps(self, write_ini, tmp_path):
        config = totest.read_config(write_ini(
            MINIMAL + "[user_inlists]\nstar1_formation = s1, s2\n"
            "star1_formation_step1 = tweak\n"))
        targets = config.inlists["user_inlists"].targets
        assert targets == {
            "star1_formation_step0_controls": [tmp_path / "s1"],
            "star1_formation_step0_job": [tmp_path / "s1"],
            "star1_formation_step1_controls": [tmp_path / "s2",
                                               tmp_path / "tweak"],
            "star1_formation_step1_job": [tmp_path / "s2",
                                          tmp_path / "tweak"]}

    def test_whole_files(self, write_ini, tmp_path):
        config = totest.read_config(write_ini(
            MINIMAL + "[user_inlists]\nprofile_columns = p.list\n"))
        assert config.inlists["user_inlists"].whole_files == {
            "profile_columns": tmp_path / "p.list"}

    @pytest.mark.parametrize("text, message", [
        ("star1_controls = ,\n", "no file given"),
        ("profile_columns = a.list, b.list\n", "expected one file"),
        ("star3 = x\n", "unknown key"),
    ])
    def test_errors(self, write_ini, text, message):
        with pytest.raises(totest.ConfigError, match=message):
            totest.read_config(write_ini(MINIMAL + "[user_inlists]\n" + text))


class TestOutputSection:
    def test_values(self, write_ini):
        config = totest.read_config(write_ini(
            MINIMAL + "[output]\nhistory_interval = 5\n"
            "star1 = history, final_profile\nbinary = none\n"))
        assert config.output.history_interval == 5
        assert config.output.star1 == ("history", "final_profile")
        assert config.output.star2 == ("history", "final_model")
        assert config.output.binary_history is False

    @pytest.mark.parametrize("text, message", [
        ("history_interval = often\n", "expected an integer"),
        ("star2 = everything\n", r"\[output\] star2: unknown value"),
        ("binary = final_model\n", "expected history or none"),
        ("star3 = none\n", "unknown key"),
    ])
    def test_errors(self, write_ini, text, message):
        with pytest.raises(totest.ConfigError, match=message):
            totest.read_config(write_ini(MINIMAL + "[output]\n" + text))


class TestSubmissionSection:
    def test_values(self, write_ini):
        config = totest.read_config(write_ini(
            MINIMAL + "[submission]\nbackend = slurm_mpi\nnodes = 2\n"
            "job_array = false\n"
            "mpi_tasks = 8\ncpus_per_task = 4\naccount = test\n"
            "email = a@b.c\nwork_dir = /scratch\nnewgroup = posydon\n"))
        submission = config.submission
        assert (submission.backend, submission.nodes, submission.mpi_tasks,
                submission.cpus_per_task) == ("slurm_mpi", 2, 8, 4)
        assert submission.job_array is False
        assert submission.work_dir == "/scratch"

    @pytest.mark.parametrize("text, message", [
        ("backend = pbs\n", "unknown value 'pbs'"),
        ("nodes = two\n", "expected an integer"),
        ("queue = long\n", "unknown key"),
    ])
    def test_errors(self, write_ini, text, message):
        with pytest.raises(totest.ConfigError, match=message):
            totest.read_config(write_ini(MINIMAL + "[submission]\n" + text))
