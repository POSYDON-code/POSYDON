"""Unit tests of posydon/grids/setup/legacy.py
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

# import other needed code for the tests, which is not already imported in the
# module you like to test
import doctest
from inspect import isroutine
from pathlib import Path

import pytest

# import the module which will be tested
import posydon.grids.setup.legacy as totest
from posydon.grids.setup.config import ConfigError
from posydon.grids.setup.layers import section_for

# shaped like the CO-HeMS_yggdrasil.ini files of the v2 grids
LEGACY_INI = """
[slurm]
number_of_nodes=1
number_of_mpi_tasks=1
number_of_cpus_per_task=4
job_array=True
user=tester
partition=public-cpu
account=test
walltime='2-00:00:00'
work_dir=/scratch
email=tester@posydon.invalid
newgroup=posydon

[mesa_inlists]
posydon_github_root=/opt/POSYDON
scenario = ['posydon', 'fix_iniHe_pre_rerun5-0e736a853c3c2e522b76299a54aea8f5eea15d08', 'CO-He_star']
single_star_grid = False
star1_formation_job_mesa_defaults = ${posydon_github_root}/defaults/star_job.defaults
star1_formation_controls_mesa_defaults = ${posydon_github_root}/defaults/controls.defaults
binary_controls_mesa_defaults = ${posydon_github_root}/defaults/binary_controls.defaults
binary_job_mesa_defaults = ${posydon_github_root}/defaults/binary_job.defaults
star1_job_mesa_defaults = ${posydon_github_root}/defaults/star_job.defaults
star1_controls_mesa_defaults = ${posydon_github_root}/defaults/controls.defaults
star2_job_mesa_defaults = ${posydon_github_root}/defaults/star_job.defaults
star2_controls_mesa_defaults = ${posydon_github_root}/defaults/controls.defaults
star_history_columns = ${posydon_github_root}/columns/history_columns.list
binary_history_columns = ${posydon_github_root}/columns/binary_history_columns.list
profile_columns = ${posydon_github_root}/columns/profile_columns.list
history_interval = 1
binary_history = True
history_star1 = True
final_profile_star1 = False
final_model_star1 = True
history_star2 = False
final_profile_star2 = False
final_model_star2 = False

[mesa_extras]
makefile_binary = ${MESA_DIR}/binary/work/make/makefile
mesa_binary_extras = ${MESA_DIR}/binary/work/src/run_binary_extras.f
binary_run = ${MESA_DIR}/binary/work/src/binary_run.f

[run_parameters]
grid = grid.csv
"""


@pytest.fixture
def write_ini(tmp_path):
    def write(text=LEGACY_INI, replace=()):
        for old, new in replace:
            assert old in text
            text = text.replace(old, new)
        path = tmp_path / "grid.ini"
        path.write_text(text)
        return path
    return write


def read(path, **kwargs):
    kwargs.setdefault("mesa_dir", "/mesa")
    return totest.read_legacy_config(path, **kwargs)


class TestElements:
    def test_dir(self):
        elements = {'COLUMN_KEYS', 'ConfigError', 'GridConfig',
                    'INLIST_KEY', 'KNOWN_KEYS', 'LEGACY_GRID_TYPES',
                    'LEGACY_LAYERS', 'OUTPUT_KEYS', 'Path',
                    'SubmissionSection', '__authors__', '__builtins__',
                    '__cached__', '__doc__', '__file__', '__loader__',
                    '__name__', '__package__', '__spec__', '_as_list',
                    '_read_inlist_keys', '_read_with_mesa_dir', 'os',
                    'parse_inifile', 're', 'read_legacy_config',
                    'rewrite_grid_params', 'split_scenario_ref',
                    'SUPPORT_COMMIT', 'map_grid_params_path'}
        assert set(dir(totest)) == elements

    def test_instance_LEGACY_GRID_TYPES(self):
        assert totest.LEGACY_GRID_TYPES[("CO-H_star", False)] == "CO-HMS_RLO"
        assert totest.LEGACY_GRID_TYPES[("MS-MS", False)] == "HMS-HMS"
        assert totest.LEGACY_GRID_TYPES[("CO-He_star", True)] == \
            "single_HeMS"

    def test_instance_LEGACY_LAYERS(self):
        assert totest.LEGACY_LAYERS == {"mesa_defaults": "mesa_inlists",
                                        "posydon_defaults": "posydon_inlists",
                                        "special": "posydon_inlists",
                                        "user": "user_inlists"}

    @pytest.mark.parametrize("name", ['split_scenario_ref',
                                      'read_legacy_config'])
    def test_instance_functions(self, name):
        assert isroutine(getattr(totest, name))

    def test_doctests(self):
        result = doctest.testmod(totest)
        assert result.failed == 0
        assert result.attempted > 0


class TestInlistKey:
    @pytest.mark.parametrize("key, groups", [
        ("binary_controls_mesa_defaults",
         ("binary", None, "controls", "mesa_defaults")),
        ("star2_job_user", ("star2", None, "job", "user")),
        ("star1_formation_controls_posydon_defaults",
         ("star1", "formation_", "controls", "posydon_defaults")),
        ("star1_controls_special", ("star1", None, "controls", "special")),
    ])
    def test_matches(self, key, groups):
        match = totest.INLIST_KEY.match(key)
        assert match.group("who", "formation", "part", "suffix") == groups

    @pytest.mark.parametrize("key", ["star_history_columns", "zams_filename",
                                     "star3_job_user", "binary_controls"])
    def test_no_match(self, key):
        assert totest.INLIST_KEY.match(key) is None


class TestSplitScenarioRef:
    @pytest.mark.parametrize("ref, expected", [
        ("main-18710e943edd926a5653c4cdb7d6e18e5bdb35a2",
         ("main", "18710e943edd926a5653c4cdb7d6e18e5bdb35a2")),
        ("CO-HeMS_v1_rerun-0e736a8", ("CO-HeMS_v1_rerun", "0e736a8")),
        ("a-b-c-1234567", ("a-b-c", "1234567")),
    ])
    def test_split(self, ref, expected):
        assert totest.split_scenario_ref(ref) == expected

    @pytest.mark.parametrize("ref", ["main", "-18710e9", "main-xyz1234",
                                     "main-12345"])
    def test_invalid(self, ref):
        with pytest.raises(ConfigError, match="expected 'branch-commit'"):
            totest.split_scenario_ref(ref)


class TestRewriteGridParams:
    GRID_PARAMS_INI = LEGACY_INI.replace(
        "${posydon_github_root}/defaults/", "${posydon_github_root}/"
        "grid_params/defaults/r11701/star/").replace(
        "history_interval = 1", "history_interval = 1\nzams_filename = "
        "${posydon_github_root}/grid_params/POSYDON-MESA-INLISTS/r11701/"
        "ZAMS_models/zams.data")

    def test_missing_paths_are_rewritten(self, write_ini, tmp_path):
        calls = []

        def support():
            calls.append(1)
            return tmp_path / "cache"
        config = read(write_ini(self.GRID_PARAMS_INI), support=support)
        assert calls == [1]  # extracted once, lazily
        assert config.inlists["mesa_inlists"].targets["binary_job"] == [
            tmp_path / "cache/r11701/mesa/default/defaults/star/"
            "binary_job.defaults"]
        assert config.legacy["mesa_inlists"]["zams_filename"] == \
            str(tmp_path / "cache/r11701/ZAMS_models/zams.data")
        assert config.grid.zams_star1 == \
            tmp_path / "cache/r11701/ZAMS_models/zams.data"
        assert any("9 grid_params/ paths are read from POSYDON-MESA-"
                   "INLISTS at 8482d7d" in note for note in config.notes)

    def test_existing_paths_and_others_stay(self, write_ini, tmp_path):
        root = tmp_path / "POSYDON"
        defaults = root / "grid_params/defaults/r11701/star"
        defaults.mkdir(parents=True)
        for name in ("star_job", "controls", "binary_controls",
                     "binary_job"):
            (defaults / (name + ".defaults")).write_text("")
        text = self.GRID_PARAMS_INI.replace("/opt/POSYDON", str(root))
        config = read(write_ini(text), support=lambda: tmp_path / "cache")
        assert config.inlists["mesa_inlists"].targets["binary_job"] == [
            defaults / "binary_job.defaults"]
        assert config.extras["binary_run"] == \
            Path("/mesa/binary/work/src/binary_run.f")
        assert any("1 grid_params/ paths" in n for n in config.notes)

    def test_lists_and_unknown_paths(self):
        sections = [{"a": ["/p/grid_params/pipeline_quest.ini",
                           "/p/grid_params/unknown/file"],
                     "b": None}]
        count = totest.rewrite_grid_params(sections, lambda: Path("/c"))
        assert count == 1
        assert sections[0]["a"] == ["/c/templates/pipeline_quest.ini",
                                    "/p/grid_params/unknown/file"]

    def test_no_support(self, write_ini):
        config = read(write_ini(self.GRID_PARAMS_INI))
        assert config.legacy["mesa_inlists"]["zams_filename"].startswith(
            "/opt/POSYDON/grid_params/")


class TestReadLegacyConfig:
    def test_scenario_and_grid_type(self, write_ini):
        config = read(write_ini())
        assert config.source_format == "legacy"
        assert config.grid.type == "CO-HeMS"
        scenario = config.scenario
        assert (scenario.source, scenario.branch, scenario.system_type) == \
            ("posydon", "fix_iniHe_pre_rerun5", "CO-He_star")
        assert scenario.commit == "0e736a853c3c2e522b76299a54aea8f5eea15d08"
        assert config.notes == [config.notes[0]]
        assert "legacy ini format" in config.notes[0]

    def test_layers_and_targets(self, write_ini):
        config = read(write_ini())
        mesa = config.inlists["mesa_inlists"].targets
        assert mesa["binary_controls"] == \
            [Path("/opt/POSYDON/defaults/binary_controls.defaults")]
        assert mesa["star1_formation_every_step_controls"] == \
            [Path("/opt/POSYDON/defaults/controls.defaults")]
        assert set(mesa) == {"binary_controls", "binary_job",
                             "star1_controls", "star1_job", "star2_controls",
                             "star2_job",
                             "star1_formation_every_step_controls",
                             "star1_formation_every_step_job"}
        for target in mesa:
            section_for(target)  # every target is one layers.py knows
        assert config.inlists["user_inlists"].targets == {}
        assert config.inlists["posydon_inlists"].whole_files == {
            "history_columns": Path("/opt/POSYDON/columns/"
                                    "history_columns.list"),
            "binary_history_columns": Path("/opt/POSYDON/columns/"
                                           "binary_history_columns.list"),
            "profile_columns": Path("/opt/POSYDON/columns/"
                                    "profile_columns.list")}

    def test_output(self, write_ini):
        output = read(write_ini()).output
        assert output.history_interval == 1
        assert output.binary_history is True
        assert output.star1 == ("history", "final_model")
        assert output.star2 == ()

    def test_grid_and_command_line(self, write_ini, tmp_path):
        config = read(write_ini(), grid_type="dynamic", run_directory="run")
        assert config.grid.points == Path("grid.csv")
        assert config.grid.mode == "dynamic"
        assert config.grid.run_directory == Path("run")
        assert config.grid.zams_star1 is None
        assert (config.grid.keep_profiles, config.grid.keep_photos) == \
            (False, False)

    def test_run_parameters_options(self, write_ini):
        config = read(write_ini(replace=[(
            "grid = grid.csv", "grid = grid.csv\npsycris_inifile = p.ini\n"
            "keep_profiles = True\nkeep_photos = True")]))
        assert config.grid.psycris_inifile == Path("p.ini")
        assert (config.grid.keep_profiles, config.grid.keep_photos) == \
            (True, True)

    def test_extras_and_mesa_dir(self, write_ini):
        config = read(write_ini(), mesa_dir="/opt/mesa")
        assert config.extras == {
            "makefile_binary": Path("/opt/mesa/binary/work/make/makefile"),
            "mesa_binary_extras": Path("/opt/mesa/binary/work/src/"
                                       "run_binary_extras.f"),
            "binary_run": Path("/opt/mesa/binary/work/src/binary_run.f")}

    def test_submission(self, write_ini):
        submission = read(write_ini()).submission
        assert submission.backend == "slurm_array"
        assert (submission.partition, submission.account,
                submission.walltime) == ("public-cpu", "test", "2-00:00:00")
        assert (submission.nodes, submission.mpi_tasks,
                submission.cpus_per_task) == (1, 1, 4)
        assert (submission.work_dir, submission.newgroup) == \
            ("/scratch", "posydon")

    def test_submission_types(self, write_ini):
        assert read(write_ini(), submission_type="shell").submission\
            .backend == "shell"
        no_array = write_ini(replace=[("job_array=True", "job_array=False")])
        shell = read(no_array, submission_type="shell").submission
        assert (shell.backend, shell.job_array) == ("shell", False)
        assert read(no_array, submission_type="slurm").submission\
            .backend == "slurm_mpi"

    def test_raw_sections_are_kept(self, write_ini):
        legacy = read(write_ini()).legacy
        assert set(legacy) == {"run_parameters", "slurm", "mesa_inlists",
                               "mesa_extras"}
        keys = list(legacy["mesa_inlists"])
        assert keys.index("star1_formation_job_mesa_defaults") < \
            keys.index("binary_controls_mesa_defaults")


class TestGridTypes:
    @pytest.mark.parametrize("system_type, single, grid_type", [
        ("HMS-HMS", "False", "HMS-HMS"),
        ("CO-H_star", "False", "CO-HMS_RLO"),
        ("HeMS-HMS", "False", "HeMS-HMS"),
        ("HMS-HMS", "True", "single_HMS"),
        ("CO-He_star", "True", "single_HeMS"),
    ])
    def test_grid_type(self, write_ini, system_type, single, grid_type):
        config = read(write_ini(replace=[
            ("'CO-He_star']", "'{}']".format(system_type)),
            ("single_star_grid = False",
             "single_star_grid = {}".format(single))]))
        assert config.grid.type == grid_type

    def test_MS_MS(self, write_ini):
        config = read(write_ini(replace=[("'CO-He_star']", "'MS-MS']")]))
        assert config.grid.type == "HMS-HMS"
        assert "system type 'MS-MS' does not exist; read as HMS-HMS" in \
            config.notes

    def test_unknown_system_type(self, write_ini):
        config = read(write_ini(replace=[("'CO-He_star']", "'XX-YY']")]))
        assert config.grid.type is None
        assert any("unknown system type 'XX-YY'" in note
                   for note in config.notes)

    def test_no_scenario(self, write_ini):
        config = read(write_ini(replace=[(
            "scenario = ['posydon', 'fix_iniHe_pre_rerun5-0e736a853c3c2e522b"
            "76299a54aea8f5eea15d08', 'CO-He_star']\n", "")]))
        assert config.grid.type is None
        assert config.scenario.is_set is False
        assert any("no scenario" in note for note in config.notes)


class TestInlistKeys:
    def test_explicit_layers_and_lists(self, write_ini):
        config = read(write_ini(replace=[("single_star_grid = False", (
            "single_star_grid = False\n"
            "binary_controls_user = /u/inlist_project\n"
            "star1_controls_posydon_defaults = /p/inlist1\n"
            "star1_formation_controls_user = /u/step1, /u/step2\n"
            "star1_controls_special = /p/flag"))]))
        user = config.inlists["user_inlists"].targets
        assert user["binary_controls"] == [Path("/u/inlist_project")]
        assert user["star1_formation_step0_controls"] == [Path("/u/step1")]
        assert user["star1_formation_step1_controls"] == [Path("/u/step2")]
        posydon = config.inlists["posydon_inlists"].targets
        assert posydon["star1_controls"] == [Path("/p/inlist1"),
                                             Path("/p/flag")]

    def test_single_star_uses_formation_steps(self, write_ini):
        config = read(write_ini(replace=[
            ("single_star_grid = False", "single_star_grid = True\n"
             "star1_controls_user = /u/step1, /u/step2"),
            ("'CO-He_star']", "'CO-He_star']")]))
        targets = config.inlists["user_inlists"].targets
        assert set(targets) == {"star1_formation_step0_controls",
                                "star1_formation_step1_controls"}
        assert "star1_formation_every_step_controls" in \
            config.inlists["mesa_inlists"].targets
        assert "star2_controls" in config.inlists["mesa_inlists"].targets

    def test_none_values_are_skipped(self, write_ini):
        config = read(write_ini(replace=[(
            "single_star_grid = False",
            "single_star_grid = False\nstar2_job_user = None")]))
        assert config.inlists["user_inlists"].targets == {}

    def test_zams(self, write_ini):
        config = read(write_ini(replace=[(
            "single_star_grid = False",
            "single_star_grid = False\nzams_filename_2 = /z/two.data")]))
        assert config.grid.zams_star1 is None
        assert config.grid.zams_star2 == Path("/z/two.data")

    def test_zams_filename_wins_over_per_star(self, write_ini):
        """As today: zams_filename sets both stars."""
        config = read(write_ini(replace=[(
            "single_star_grid = False",
            "single_star_grid = False\nzams_filename_1 = /z/one.data\n"
            "zams_filename = /z/both.data")]))
        assert config.grid.zams_star1 == config.grid.zams_star2 == \
            Path("/z/both.data")

    def test_missing_column_key(self, write_ini):
        config = read(write_ini(replace=[(
            "profile_columns = ${posydon_github_root}/columns/"
            "profile_columns.list\n", "")]))
        assert "profile_columns" not in \
            config.inlists["posydon_inlists"].whole_files

    def test_unused_keys_are_noted(self, write_ini):
        config = read(write_ini(replace=[(
            "single_star_grid = False",
            "single_star_grid = False\nmystery_key = 1")]))
        assert "[mesa_inlists] keys not used: mystery_key" in config.notes


class TestErrors:
    def test_missing_file(self, tmp_path):
        with pytest.raises(ConfigError, match="not found"):
            read(tmp_path / "nope.ini")

    def test_missing_section(self, write_ini):
        path = write_ini(replace=[("[mesa_extras]", "[extras]")])
        with pytest.raises(ConfigError,
                           match=r"\[mesa_extras\]: missing section"):
            read(path)

    def test_missing_grid(self, write_ini):
        with pytest.raises(ConfigError,
                           match=r"\[run_parameters\] grid: missing"):
            read(write_ini(replace=[("grid = grid.csv", "")]))

    def test_bad_scenario(self, write_ini):
        with pytest.raises(ConfigError, match="expected 'branch-commit'"):
            read(write_ini(replace=[(
                "fix_iniHe_pre_rerun5-"
                "0e736a853c3c2e522b76299a54aea8f5eea15d08",
                "fix_iniHe_pre_rerun5")]))


class TestMesaDir:
    def test_missing_mesa_dir(self, write_ini, monkeypatch):
        monkeypatch.delenv("MESA_DIR", raising=False)
        with pytest.raises(ConfigError, match="MESA_DIR is not set"):
            totest.read_legacy_config(write_ini())

    def test_environment_is_used(self, write_ini, monkeypatch):
        monkeypatch.setenv("MESA_DIR", "/env/mesa")
        config = totest.read_legacy_config(write_ini())
        assert config.extras["binary_run"] == \
            Path("/env/mesa/binary/work/src/binary_run.f")

    def test_environment_is_restored(self, write_ini, monkeypatch):
        monkeypatch.setenv("MESA_DIR", "/env/mesa")
        read(write_ini(), mesa_dir="/other")
        assert totest.os.environ["MESA_DIR"] == "/env/mesa"
        monkeypatch.delenv("MESA_DIR")
        read(write_ini(), mesa_dir="/other")
        assert "MESA_DIR" not in totest.os.environ
