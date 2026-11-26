"""Unit tests of posydon/grids/setup/keycheck.py
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>",
]

# import other needed code for the tests, which is not already imported in the
# module you like to test
from inspect import isroutine
from pathlib import Path

import pytest

# import the module which will be tested
import posydon.grids.setup.keycheck as totest
from posydon.grids.setup.config import GridConfig
from posydon.grids.setup.defaults import fill_defaults
from posydon.grids.setup.manifest import MESA_DEFAULTS, Release
from posydon.grids.setup.resolve import ResolvedScenario
from posydon.grids.setup.scenario import ScenarioFiles

# the defaults of a small MESA installation, written like MESA's own
DEFAULTS = {
    "binary_controls": "      mdot_scheme = 'Kolb'\n",
    "binary_job": "      evolve_both_stars = .true.\n",
    "star_controls": """\
      !### max_model_number
      ! Negative means no maximum.
      max_model_number = -1
      x_ctrl(1:num_x_ctrls) = 0d0
      x_logical_ctrl(1:num_x_ctrls) = .false.
      !### xa_central_lower_limit_species
      ! xa_central_lower_limit_species(1) = ''
      !### num_accretion_species
      !### accretion_species_id
      ! E.g., if `min_overshoot_q = 0.1`, then only the outer 90%
""",
    "star_job": "      pgstar_flag = .false.\n      ! `s% dt` = 1\n",
}


def write(path, text=""):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


@pytest.fixture
def release(tmp_path):
    """A release whose MESA installation has the DEFAULTS above."""
    folder = tmp_path / "r99999" / "mesa" / "other" / "defaults"
    for name, rel in MESA_DEFAULTS.items():
        write(folder / rel, DEFAULTS[name])
    return Release(name="v9", folder=tmp_path / "r99999" / "releases" / "v9",
                   mesa_version="r99999", mesa="other", frozen=False)


def resolved_for(release):
    return ResolvedScenario(files=ScenarioFiles(), repo=Path("/repo"),
                            release=release)


def new_config(grid_type="HMS-HMS"):
    config = GridConfig()
    config.grid.type = grid_type
    return config


def filled(config, release):
    resolved = resolved_for(release)
    fill_defaults(config, resolved, "/mesa")
    return resolved


class TestElements:
    def test_dir(self):
        elements = {'CHECKED_LAYERS', 'DEFAULTS_NAMELISTS', 'Path', 'Pwarn',
                    'UnknownKeys', '_ASSIGNMENT', '_HEADING', '__authors__',
                    '__builtins__', '__cached__', '__doc__', '__file__',
                    '__loader__', '__name__', '__package__', '__spec__',
                    'check_keys', 'dataclass', 'defaults_keys',
                    'drop_wiring_keys', 'key_name', 'layer_targets',
                    'mesa_installation', 're', 'read_inlist',
                    'reference_files', 'section_for', 'unknown_keys'}
        # Pwarn's warnings of other tests may add a registry
        assert set(dir(totest)) - {'__warningregistry__'} == elements

    @pytest.mark.parametrize("name", ['key_name', 'defaults_keys',
                                      'reference_files', 'unknown_keys',
                                      'check_keys'])
    def test_instance_functions(self, name):
        assert isroutine(getattr(totest, name))

    def test_every_defaults_file_has_a_namelist(self):
        assert set(totest.DEFAULTS_NAMELISTS) == set(MESA_DEFAULTS)


class TestKeyName:
    @pytest.mark.parametrize("key, name", [
        ("max_model_number", "max_model_number"),
        ("Max_Model_Number", "max_model_number"),
        ("x_logical_ctrl(5)", "x_logical_ctrl"),
        ("accretion_species_id(1)", "accretion_species_id"),
        ("x_ctrl(1:num_x_ctrls)", "x_ctrl"),
        ("xa_central_lower_limit_species (1)",
         "xa_central_lower_limit_species")])
    def test_names(self, key, name):
        assert totest.key_name(key) == name


class TestDefaultsKeys:
    def test_set_commented_and_headings(self, tmp_path):
        path = write(tmp_path / "controls.defaults",
                     DEFAULTS["star_controls"])
        assert totest.defaults_keys(path) == {
            "max_model_number", "x_ctrl", "x_logical_ctrl",
            "xa_central_lower_limit_species", "num_accretion_species",
            "accretion_species_id"}

    def test_upper_case_and_prose(self, tmp_path):
        path = write(tmp_path / "star_job.defaults",
                     "  PGSTAR_FLAG = .false.  ! a = b\n"
                     "  ! s% x(k) = 1\n  ! not a key\n  !### two words\n")
        assert totest.defaults_keys(path) == {"pgstar_flag"}


class TestReferenceFiles:
    def test_the_release_installation(self, release):
        folder = release.folder.parents[1] / "mesa" / "other" / "defaults"
        references = totest.reference_files(new_config(),
                                            resolved_for(release))
        assert references == {
            "binary_controls": [folder / "binary/binary_controls.defaults"],
            "binary_job": [folder / "binary/binary_job.defaults"],
            "controls": [folder / "star/controls.defaults"],
            "star_job": [folder / "star/star_job.defaults"]}

    def test_without_a_release_the_mesa_layer_defaults(self, tmp_path):
        controls = write(tmp_path / "controls.defaults")
        config = new_config()
        layer = config.inlists["mesa_inlists"]
        layer.add("star1_controls", str(controls))
        layer.add("star2_controls", controls)  # the same file, listed once
        layer.add("star1_job", write(tmp_path / "my_star_job"))  # no
        layer.add("binary_job", tmp_path / "missing.defaults")  # no
        resolved = ResolvedScenario(files=ScenarioFiles(), repo=tmp_path)
        assert totest.reference_files(config, resolved) == \
            {"controls": [controls]}
        assert totest.reference_files(config) == {"controls": [controls]}


class TestUnknownKeys:
    def test_misspelled_and_misplaced(self, tmp_path, release):
        mine = write(tmp_path / "mine", """\
&controls
    max_modle_number = 100
    max_model_number = 100
    x_logical_ctrl(5) = .true.
    X_CTRL(2) = 1d0
    xa_central_lower_limit_species(1) = 'he4'
    accretion_species_id(1) = 'h1'
/
&star_job
    max_model_number = 100
    pgstar_flag = .true.
/
""")
        config = new_config()
        config.inlists["user_inlists"].add("star1_controls", mine)
        config.inlists["user_inlists"].add("star1_job", mine)
        resolved = filled(config, release)
        found = totest.unknown_keys(config, resolved)
        assert [(u.layer, u.path, u.namelist, u.keys) for u in found] == [
            ("user_inlists", mine, "controls", ["max_modle_number"]),
            ("user_inlists", mine, "star_job", ["max_model_number"])]
        folder = release.folder.parents[1] / "mesa" / "other" / "defaults"
        assert found[0].references == [folder / "star/controls.defaults"]

    def test_a_file_of_several_targets_is_listed_once(self, tmp_path,
                                                      release):
        mine = write(tmp_path / "mine", "&controls\n  a = 1\n  b = 2\n/\n")
        config = new_config("CO-HeMS")
        for target in ("star1_controls", "star2_controls",
                       "star1_formation_step0_controls",
                       "star1_formation_step1_controls"):
            config.inlists["user_inlists"].add(target, mine)
        resolved = filled(config, release)
        found = totest.unknown_keys(config, resolved)
        assert [(u.path, u.keys) for u in found] == [(mine, ["a", "b"])]

    def test_wiring_keys_and_good_files_pass(self, tmp_path, release):
        mine = write(tmp_path / "mine", """\
&binary_job
    read_extra_binary_job_inlist1 = .true.
    inlist_names(1) = 'x'
    evolve_both_stars = .false.
/
""")
        config = new_config()
        config.inlists["user_inlists"].add("binary_job", mine)
        assert totest.unknown_keys(config, filled(config, release)) == []

    def test_mesa_inlists_replacement_checked_against_the_release(
            self, tmp_path, release):
        """With a release, a [mesa_inlists] file that replaces MESA's
        defaults for a target is checked against the installation's."""
        newer = write(tmp_path / "newer" / "binary_controls.defaults",
                      "mdot_scheme = 'Kolb'\nnew_in_a_later_mesa = 1\n")
        config = new_config()
        config.inlists["mesa_inlists"].add("binary_controls", newer)
        found = totest.unknown_keys(config, filled(config, release))
        assert [(u.layer, u.path, u.keys) for u in found] == [
            ("mesa_inlists", newer, ["new_in_a_later_mesa"])]

    def test_without_a_release(self, tmp_path):
        """The mesa_inlists layer's .defaults files are the reference; a
        namelist without one, and missing files, are not checked."""
        controls = write(tmp_path / "controls.defaults", "alpha = 1\n")
        partial = write(tmp_path / "partial_controls", "beta = 2\n")
        mine = write(tmp_path / "mine", "&controls\n  Alpha = 2\n"
                     "  gamma = 3\n/\n&star_job\n  anything = 1\n/\n")
        config = new_config()
        mesa = config.inlists["mesa_inlists"]
        mesa.add("star1_controls", controls)
        mesa.add("star2_controls", partial)
        user = config.inlists["user_inlists"]
        user.add("star1_controls", mine)
        user.add("star1_job", mine)  # no star_job reference: not checked
        user.add("star2_controls", tmp_path / "missing")
        found = totest.unknown_keys(config)
        assert [(u.layer, u.path, u.keys) for u in found] == [
            ("mesa_inlists", partial, ["beta"]),
            ("user_inlists", mine, ["gamma"])]

    def test_release_files_are_not_checked(self, tmp_path, release):
        theirs = write(tmp_path / "release_file", "&controls\n  zz = 1\n/\n")
        config = new_config()
        resolved = filled(config, release)
        resolved.files.inlists.add("star1_controls", theirs)
        assert totest.unknown_keys(config, resolved) == []


class TestCheckKeys:
    def test_one_warning_lists_every_file_and_namelist(self, tmp_path,
                                                       release):
        mine = write(tmp_path / "mine", "&controls\n  max_modle_number = 1\n"
                     "/\n&star_job\n  a = 1\n  b = 2\n/\n")
        config = new_config()
        config.inlists["user_inlists"].add("star1_controls", mine)
        config.inlists["user_inlists"].add("star1_job", mine)
        resolved = filled(config, release)
        with pytest.warns(Warning) as record:
            found = totest.check_keys(config, resolved)
        assert len(found) == 2
        assert len(record) == 1
        message = str(record[0].message)
        folder = release.folder.parents[1] / "mesa" / "other" / "defaults"
        assert "MESA stops at a key it does not know, when the run " \
            "starts: [user_inlists] {0}: max_modle_number is not a key of " \
            "&controls in this MESA version ({1}); [user_inlists] {0}: " \
            "a, b are not keys of &star_job in this MESA version " \
            "({2})".format(mine, folder / "star/controls.defaults",
                           folder / "star/star_job.defaults") in message

    def test_silent_when_all_is_known(self, release, recwarn):
        config = new_config()
        assert totest.check_keys(config, filled(config, release)) == []
        assert len(recwarn) == 0
