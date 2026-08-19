"""Unit tests for posydon.CLI.grids.inlists."""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

import os
from pathlib import Path

import pytest

from posydon.CLI.grids import inlists
from posydon.utils.posydonwarning import OverwriteWarning

totest = inlists


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


class TestInlistSection:
    """Tests for the InlistSection abstraction."""

    def test_to_fortran_float_d_notation(self):
        section = totest.InlistSection("controls", {"m1": 10.5})
        assert (
            section.to_fortran()
            == "&controls\n\tm1 = 1.0500000000d+01\n/ ! end of controls namelist\n"
        )

    def test_to_fortran_float_int(self):
        section = totest.InlistSection("controls", {"m1": 10.0})
        assert (
            section.to_fortran() == "&controls\n\tm1 = 10\n/ ! end of controls namelist\n"
        )

    def test_to_fortran_bool(self):
        section = totest.InlistSection("controls", {"on": True, "off": False})
        assert (
            section.to_fortran()
            == "&controls\n\ton = .true.\n\toff = .false.\n/ ! end of controls namelist\n"
        )

    def test_to_fortran_str_and_int_passthrough(self):
        section = totest.InlistSection("controls", {"s": "foo", "i": 5})
        fortran = section.to_fortran()
        assert "s = foo" in fortran
        assert "i = 5" in fortran

    def test_format_value(self):
        assert totest.InlistSection._format_value(True) == ".true."
        assert totest.InlistSection._format_value(False) == ".false."
        assert totest.InlistSection._format_value(10.5) == "1.0500000000d+01"
        assert totest.InlistSection._format_value(10.0) == "10"
        assert totest.InlistSection._format_value("abc") == "abc"
        assert totest.InlistSection._format_value(7) == 7

    def test_merge(self):
        a = totest.InlistSection("controls", {"x": "1"})
        b = totest.InlistSection("controls", {"y": "2"})
        assert a.merge(b) is a
        assert a.parameters == {"x": "1", "y": "2"}
        c = totest.InlistSection("controls", {"x": "9"})
        a.merge(c)
        assert a.parameters == {"x": "9", "y": "2"}

    def test_merge_different_names_raises(self):
        with pytest.raises(
            ValueError, match="Cannot merge sections with different names"
        ):
            totest.InlistSection("controls", {}).merge(
                totest.InlistSection("star_job", {})
            )

    def test_parse_value_strips(self):
        assert totest.InlistSection._parse_value("  foo  ") == "foo"

    def test_from_string(self):
        section = totest.InlistSection.from_string(
            "controls", "&controls\n  m1 = 1\n/ ! end\n"
        )
        assert section.name == "controls"
        assert section.parameters == {"m1": "1"}

    def test_from_string_bare_line_skipped(self):
        section = totest.InlistSection.from_string(
            "controls", "&controls\nfoo\n  m1 = 1\n/ ! end\n"
        )
        assert section.parameters == {"m1": "1"}

    def test_from_file(self, tmp_path):
        f = tmp_path / "c.inlist"
        f.write_text("&controls\n  m1 = 2\n/ ! end\n")
        section = totest.InlistSection.from_file("controls", str(f))
        assert section.parameters == {"m1": "2"}

    def test_repr(self):
        assert repr(totest.InlistSection("controls", {})) == (
            "InlistSection(name='controls', parameters=0)"
        )
        assert repr(totest.InlistSection("controls", {"a": "1"})) == (
            "InlistSection(name='controls', parameters=1: [a])"
        )
        section = totest.InlistSection(
            "controls", {"a": "1", "b": "2", "c": "3", "d": "4"}
        )
        assert "parameters=4" in repr(section)
        assert "..." in repr(section)


class TestInlist:
    """Tests for the Inlist abstraction."""

    def test_attribute_access(self):
        inlist = totest.Inlist("inlist")
        inlist.add_section(totest.InlistSection("controls", {"m1": "1"}))
        assert inlist.controls.parameters == {"m1": "1"}
        assert inlist.sections["controls"].name == "controls"
        with pytest.raises(AttributeError):
            inlist.star_job

    def test_setattr_adds_section(self):
        inlist = totest.Inlist("inlist")
        inlist.my_section = totest.InlistSection("my_section", {"a": "1"})
        assert inlist.sections["my_section"].parameters == {"a": "1"}
        assert inlist.my_section.parameters == {"a": "1"}

    def test_setattr_plain(self):
        inlist = totest.Inlist("inlist")
        inlist.extra = "something"
        assert inlist.extra == "something"

    def test_set_name(self):
        inlist = totest.Inlist("a")
        inlist.name = "b"
        assert inlist.name == "b"

    def test_add_section_merges_existing(self):
        inlist = totest.Inlist("inlist")
        inlist.add_section(totest.InlistSection("controls", {"x": "1"}))
        inlist.add_section(totest.InlistSection("controls", {"y": "2"}))
        assert inlist.sections["controls"].parameters == {"x": "1", "y": "2"}

    def test_merge_returns_new(self):
        a = totest.Inlist("inlist")
        a.add_section(totest.InlistSection("controls", {"x": "1"}))
        b = totest.Inlist("inlist")
        b.add_section(totest.InlistSection("controls", {"y": "2"}))
        b.add_section(totest.InlistSection("star_job", {"s": "0"}))
        merged = a.merge(b)
        assert isinstance(merged, totest.Inlist)
        assert merged is not a
        assert merged.controls.parameters == {"x": "1", "y": "2"}
        assert a.controls.parameters == {"x": "1"}
        assert merged.star_job.parameters == {"s": "0"}

    def test_to_file_ordering(self, tmp_path):
        inlist = totest.Inlist("inlist")
        inlist.add_section(totest.InlistSection("star_job", {"s": "0"}))
        inlist.add_section(totest.InlistSection("controls", {"m": "1"}))
        inlist.add_section(totest.InlistSection("binary_controls", {"b": "1"}))
        inlist.add_section(totest.InlistSection("extra_section", {"e": "1"}))
        out = tmp_path / "inlist"
        inlist.to_file(str(out))
        content = out.read_text()
        assert content.index("&controls") < content.index("&star_job")
        assert content.index("&star_job") < content.index("&binary_controls")
        assert content.index("&binary_controls") < content.index("&extra_section")

    def test_to_file_empty_sections(self, tmp_path):
        inlist = totest.Inlist("empty")
        out = tmp_path / "empty_inlist"
        inlist.to_file(str(out))
        assert out.read_text() == ""

    def test_from_string_multi_section(self):
        content = (
            "&controls\n  m1 = 1\n/ ! end\n&star_job\n  load = .true.\n/ ! end\n"
        )
        inlist = totest.Inlist.from_string(content, name="inlist")
        assert inlist.name == "inlist"
        assert inlist.controls.parameters == {"m1": "1"}
        assert inlist.star_job.parameters == {"load": ".true."}

    def test_from_string_sectionless_with_section(self):
        inlist = totest.Inlist.from_string(
            "  m1 = 1\n  m2 = ''\n", name="d", section="controls"
        )
        assert inlist.controls.parameters == {"m1": "1"}
        assert "m2" not in inlist.controls.parameters

    def test_from_string_sectionless_no_section(self):
        inlist = totest.Inlist.from_string("  m1 = 1\n")
        assert inlist.sections == {}

    def test_from_string_param_before_header(self):
        content = "  foo = 1\n&controls\n  m1 = 2\n/ ! end\n"
        inlist = totest.Inlist.from_string(content, name="inlist")
        assert inlist.controls.parameters == {"m1": "2"}

    def test_from_string_skips_empty_and_comment_lines(self):
        content = (
            "&controls\n\n  m1 = 1\n! comment\n  m2 = '.'\n  m3 = ''\nfoo\n"
            "/ ! end\n&star_job\n  s = ''\n  load = .true.\n/ ! end\n"
        )
        inlist = totest.Inlist.from_string(content, name="inlist")
        assert inlist.controls.parameters == {"m1": "1"}
        assert inlist.star_job.parameters == {"load": ".true."}

    def test_from_file(self, tmp_path):
        f = tmp_path / "x.inlist"
        f.write_text("&controls\n  m1 = 1\n/ ! end\n")
        inlist = totest.Inlist.from_file(str(f))
        assert inlist.name == "x.inlist"
        assert inlist.controls.parameters == {"m1": "1"}

    def test_from_file_with_name(self, tmp_path):
        f = tmp_path / "x.inlist"
        f.write_text("&controls\n  m1 = 1\n/ ! end\n")
        inlist = totest.Inlist.from_file(str(f), name="custom")
        assert inlist.name == "custom"

    def test_repr(self):
        assert repr(totest.Inlist("empty")) == "Inlist(name='empty', sections=0)"
        inlist = totest.Inlist("inlist")
        inlist.add_section(totest.InlistSection("controls", {}))
        assert repr(inlist) == "Inlist(name='inlist', sections=1: [controls])"


class TestMESAInlists:
    """Tests for the MESAInlists class."""

    def _write_mesa_dir(self, tmp_path):
        """Write a synthetic MESA defaults directory."""
        mesa = tmp_path / "mesa"
        star = mesa / "star"
        binary = mesa / "binary"
        star.mkdir(parents=True)
        binary.mkdir(parents=True)
        (star / "controls.defaults").write_text(
            "  read_extra_controls_inlist1 = .true.\n"
            "  extra_controls_inlist1_name = 'inlist1'\n"
            "  num_x_ctrls = 3\n"
            "  initial_mass = 10.0\n"
            "  empty_val = ''\n"
        )
        (star / "star_job.defaults").write_text(
            "  load_saved_model = .false.\n"
            "  read_extra_star_job_inlist1 = .true.\n"
        )
        (binary / "binary_job.defaults").write_text(
            "  evolve_both_stars = .true.\n"
            "  read_extra_binary_job_inlist1 = .true.\n"
        )
        (binary / "binary_controls.defaults").write_text(
            "  mdot_scheme = 'Kolb'\n"
            "  read_extra_controls_inlist1 = .true.\n"
        )
        return str(mesa)

    def test_constructor_builds_base_inlists(self, tmp_path):
        path = self._write_mesa_dir(tmp_path)
        mi = totest.MESAInlists(path)
        assert mi.path == path
        assert mi.base_star_inlist.controls.parameters == {
            "1": "3",
            "initial_mass": "10.0",
        }
        assert mi.base_star_inlist.star_job.parameters == {
            "load_saved_model": ".false."
        }
        assert mi.base_binary_inlist.binary_controls.parameters == {
            "mdot_scheme": "'Kolb'"
        }
        assert mi.base_binary_inlist.binary_job.parameters == {
            "evolve_both_stars": ".true."
        }

    def test_clean_parameters(self):
        mi = totest.MESAInlists.__new__(totest.MESAInlists)
        assert mi._clean_parameters(
            {"a": "1", "read_extra_x": "2", "inlist_y": "3"}
        ) == {"a": "1"}
        assert mi._clean_parameters({"num_x_ctrls": "3", "b": "2"}) == {
            "1": "3",
            "b": "2",
        }

    def test_repr(self, tmp_path):
        path = self._write_mesa_dir(tmp_path)
        mi = totest.MESAInlists(path)
        assert repr(mi).startswith("MESAInlists(path='{0}'".format(path))


class TestInlistManager:
    """Tests for the InlistManager class."""

    def _inlist(self, name, section, params):
        inlist = totest.Inlist(name)
        inlist.add_section(totest.InlistSection(section, params))
        return inlist

    def test_append_and_getitem(self):
        mgr = totest.InlistManager()
        a = self._inlist("a", "controls", {"m1": "1"})
        b = self._inlist("b", "star_job", {"s": "1"})
        mgr.append_binary_inlist(a)
        mgr.append_binary_star1_inlist(b)
        mgr.append_binary_star2_inlist(a)
        mgr.append_star1_inlist(b)
        mgr.append_star2_inlist(a)
        assert mgr.binary_inlists == [a]
        assert mgr.binary_star1_inlists == [b]
        assert mgr.binary_star2_inlists == [a]
        assert mgr.star1_inlists == [b]
        assert mgr.star2_inlists == [a]
        assert mgr["binary_inlists"] == [a]
        assert mgr.keys() == [
            "binary_inlists",
            "binary_star1_inlists",
            "binary_star2_inlists",
            "star1_inlists",
            "star2_inlists",
        ]
        with pytest.raises(KeyError):
            mgr["nope"]

    def test_repr(self):
        mgr = totest.InlistManager()
        assert "binary_inlists=0" in repr(mgr)
        assert "star2_inlists=0)" in repr(mgr)

    def test_write_inlists_single_and_multi(self, tmp_path):
        mgr = totest.InlistManager()
        mgr.append_binary_inlist(self._inlist("p", "binary_controls", {"x": "1"}))
        mgr.append_binary_star1_inlist(self._inlist("i1", "controls", {"m": "1"}))
        mgr.append_binary_star2_inlist(self._inlist("i2a", "controls", {"m": "1"}))
        mgr.append_binary_star2_inlist(self._inlist("i2b", "controls", {"m": "2"}))
        mgr.append_star1_inlist(self._inlist("s1", "controls", {"s": "1"}))
        mgr.append_star2_inlist(self._inlist("s2a", "controls", {"s": "1"}))
        mgr.append_star2_inlist(self._inlist("s2b", "controls", {"s": "2"}))
        out = tmp_path / "out"
        mgr.write_inlists(str(out))
        assert (out / "binary" / "inlist_project").read_text().startswith(
            "&binary_controls"
        )
        assert (out / "binary" / "inlist1").read_text().startswith("&controls")
        assert (out / "binary" / "inlist2_step0").exists()
        assert (out / "binary" / "inlist2_step1").exists()
        assert (out / "star1" / "inlist_step0").read_text().startswith("&controls")
        assert (out / "star2" / "inlist_step0").exists()
        assert (out / "star2" / "inlist_step1").exists()

    def test_write_inlists_empty_noop(self, tmp_path):
        mgr = totest.InlistManager()
        out = tmp_path / "out"
        mgr.write_inlists(str(out))
        assert not (out / "binary").exists()


class TestMergeInlistLayers:
    """Tests for merge_inlist_layers()."""

    def test_single_star_grid_returns_empty(self, capsys):
        merged = totest.merge_inlist_layers({"single_star_grid": True}, [], "/work")
        assert merged.final_binary_controls == {}
        assert merged.final_binary_job == {}
        assert merged.final_star1_binary_controls == {}
        assert merged.final_star1_binary_job == {}
        assert merged.final_star2_binary_controls == {}
        assert merged.final_star2_binary_job == {}
        assert merged.grid_params_binary_controls == []
        assert merged.grid_params_binary_job == []
        assert merged.grid_params_star1_binary_controls == []
        assert merged.grid_params_star1_binary_job == []
        assert merged.grid_params_star2_binary_controls == []
        assert merged.grid_params_star2_binary_job == []
        assert capsys.readouterr().out == ""

    def test_all_sections_merged(self, tmp_path, capsys):
        d = tmp_path / "inlists"
        d.mkdir()
        bc = d / "binary_controls.inlist"
        bc.write_text(
            "&binary_controls\n  mdot_scheme = 'Kolb'\n  read_extra_controls_inlist1 = .true.\n/ ! end\n"
        )
        bj = d / "binary_job.inlist"
        bj.write_text(
            "&binary_job\n  evolve_both_stars = .true.\n  extra_inlist_name = 'x'\n/ ! end\n"
        )
        s1c = d / "star1_controls.inlist"
        s1c.write_text(
            "&controls\n  m1 = 10.0\n  num_x_ctrls = 3\n  read_extra_controls_inlist1 = .true.\n/ ! end\n"
        )
        s1j = d / "star1_job.inlist"
        s1j.write_text(
            "&star_job\n  load_saved_model = .false.\n  read_extra_star_job_inlist1 = .true.\n/ ! end\n"
        )
        s2c = d / "star2_controls.inlist"
        s2c.write_text(
            "&controls\n  m2 = 8.0\n  num_x_ctrls = 5\n  extra_controls_inlist1_name = 'y'\n/ ! end\n"
        )
        s2j = d / "star2_job.inlist"
        s2j.write_text(
            "&star_job\n  save_model_when_terminate = .false.\n  read_extra_star_job_inlist1 = .true.\n/ ! end\n"
        )

        mesa_inlists = {
            "single_star_grid": False,
            "binary_controls_user": str(bc),
            "binary_job_user": str(bj),
            "star1_controls_user": str(s1c),
            "star1_job_user": str(s1j),
            "star2_controls_user": str(s2c),
            "star2_job_user": str(s2j),
            "zams_filename_1": None,
            "zams_filename_2": None,
            "final_profile_star1": False,
            "final_profile_star2": False,
            "final_model_star1": False,
            "final_model_star2": False,
            "history_star1": True,
            "history_star2": True,
            "history_interval": 100,
            "binary_history": True,
        }
        grid_parameters = [
            "mdot_scheme",
            "evolve_both_stars",
            "m1",
            "m2",
            "load_saved_model",
            "save_model_when_terminate",
        ]
        merged = totest.merge_inlist_layers(mesa_inlists, grid_parameters, str(tmp_path))

        assert merged.final_binary_controls == {"mdot_scheme": "'Kolb'"}
        assert merged.final_binary_job == {
            "evolve_both_stars": ".true.",
            "inlist_names(1)": "'{0}'".format(
                os.path.join(tmp_path, "binary", "inlist1")
            ),
            "inlist_names(2)": "'{0}'".format(
                os.path.join(tmp_path, "binary", "inlist2")
            ),
        }
        assert merged.final_star1_binary_controls["m1"] == "10.0"
        assert merged.final_star1_binary_controls["1"] == "3"
        assert merged.final_star1_binary_job["load_saved_model"] == ".false."
        assert merged.final_star2_binary_controls["m2"] == "8.0"
        assert merged.final_star2_binary_controls["1"] == "5"
        assert merged.final_star2_binary_job["save_model_when_terminate"] == ".false."

        assert merged.grid_params_binary_controls == ["mdot_scheme"]
        assert merged.grid_params_binary_job == ["evolve_both_stars"]
        assert merged.grid_params_star1_binary_controls == ["m1"]
        assert merged.grid_params_star1_binary_job == ["load_saved_model"]
        assert merged.grid_params_star2_binary_controls == ["m2"]
        assert merged.grid_params_star2_binary_job == ["save_model_when_terminate"]

        assert merged.final_star1_binary_controls["read_extra_controls_inlist1"] == ".true."
        assert (
            merged.final_star1_binary_controls["extra_controls_inlist1_name"]
            == "'inlist_grid_star1_binary_controls'"
        )
        assert merged.final_star2_binary_controls["read_extra_controls_inlist1"] == ".true."
        assert (
            merged.final_star2_binary_controls["extra_controls_inlist1_name"]
            == "'inlist_grid_star2_binary_controls'"
        )
        assert merged.final_star1_binary_job["read_extra_star_job_inlist1"] == ".true."
        assert (
            merged.final_star1_binary_job["extra_star_job_inlist1_name"]
            == "'inlist_grid_star1_binary_job'"
        )
        assert merged.final_star2_binary_job["read_extra_star_job_inlist1"] == ".true."
        assert (
            merged.final_star2_binary_job["extra_star_job_inlist1_name"]
            == "'inlist_grid_star2_binary_job'"
        )

        captured = capsys.readouterr()
        assert "Grid parameters that effect binary_controls: mdot_scheme" in captured.out
        assert "Grid parameters that effect binary_job: evolve_both_stars" in captured.out
        assert "Grid parameters that effect star1_binary_controls: m1" in captured.out
        assert "Grid parameters that effect star1_binary_job: load_saved_model" in captured.out
        assert "Grid parameters that effect star2_binary_controls: m2" in captured.out
        assert "Grid parameters that effect star2_binary_job: save_model_when_terminate" in captured.out


class TestBuildStarFormation:
    """Tests for build_star_formation()."""

    def test_star1_returns_none_when_zams(self, tmp_path, inlist_files):
        mesa_inlists = _binary_mesa_inlists(inlist_files)
        mesa_inlists["zams_filename_1"] = "/some/path/zams1.data"
        result = totest.build_star_formation("star1", mesa_inlists, str(tmp_path), {})
        assert result is None

    def test_star2_returns_none_when_zams(self, tmp_path, inlist_files):
        mesa_inlists = _binary_mesa_inlists(inlist_files)
        mesa_inlists["zams_filename_2"] = "/some/path/zams2.data"
        result = totest.build_star_formation("star2", mesa_inlists, str(tmp_path), {})
        assert result is None

    def test_star2_returns_none_when_no_formation(self, tmp_path, inlist_files):
        mesa_inlists = _binary_mesa_inlists(inlist_files)
        result = totest.build_star_formation("star2", mesa_inlists, str(tmp_path), {})
        assert result is None

    def test_star2_formation_steps(self, tmp_path):
        step0 = tmp_path / "formation2" / "step0"
        step0.parent.mkdir(parents=True)
        step0.write_text(
            "&controls\n  m2 = 8.0\n  num_x_ctrls = 5\n  read_extra_controls_inlist1 = .true.\n"
            "/ ! end\n&star_job\n  load_saved_model = .false.\n  extra_inlist_name = 'x'\n/ ! end\n"
        )
        os.makedirs(os.path.join(tmp_path, "star2"), exist_ok=True)
        mesa_inlists = {
            "single_star_grid": False,
            "zams_filename_1": None,
            "zams_filename_2": None,
            "star2_formation_controls_user": str(step0),
            "star2_formation_job_user": str(step0),
        }
        counterpart = {}
        result = totest.build_star_formation(
            "star2", mesa_inlists, str(tmp_path), counterpart
        )
        step_file = os.path.join(tmp_path, "star2", "inlist_step0")
        assert result == " {0}".format(step_file)
        content = Path(step_file).read_text()
        assert "\tm2 = 8.0" in content
        assert "\t1 = 5" in content
        assert "\tload_saved_model = .false." in content
        assert "/ ! end of star2_controls namelist" in content
        assert "/ ! end of star2_job namelist" in content
        assert counterpart["saved_model_name"] == "'initial_star2_step0.mod'"
        assert "\tsave_model_when_terminate = .true." in content

    def test_star1_single_star_zams_sets_controls(self, tmp_path):
        special = tmp_path / "special.inlist"
        special.write_text(
            "&controls\n  initial_mass = 10.0\n/ ! end\n&star_job\n/ ! end\n"
        )
        os.makedirs(os.path.join(tmp_path, "star1"), exist_ok=True)
        mesa_inlists = {
            "single_star_grid": True,
            "zams_filename_1": "/path/zams1.data",
            "zams_filename_2": None,
            "star1_controls_special": str(special),
            "star1_job_special": str(special),
        }
        result = totest.build_star_formation("star1", mesa_inlists, str(tmp_path), {})
        step_file = os.path.join(tmp_path, "star1", "inlist_step0")
        assert result == " {0}".format(step_file)
        content = Path(step_file).read_text()
        assert "\tzams_filename = '/path/zams1.data'" in content

    def test_star1_single_star_pops_zams_filename(self, tmp_path):
        special = tmp_path / "special.inlist"
        special.write_text(
            "&controls\n  zams_filename = 'orig'\n  zams_filename_1 = 'x'\n/ ! end\n&star_job\n/ ! end\n"
        )
        os.makedirs(os.path.join(tmp_path, "star1"), exist_ok=True)
        mesa_inlists = {
            "single_star_grid": True,
            "zams_filename_1": None,
            "zams_filename_2": None,
            "star1_controls_special": str(special),
            "star1_job_special": str(special),
        }
        result = totest.build_star_formation("star1", mesa_inlists, str(tmp_path), {})
        step_file = os.path.join(tmp_path, "star1", "inlist_step0")
        assert result == " {0}".format(step_file)
        content = Path(step_file).read_text()
        assert "\tzams_filename_1 = 'x'" in content
        assert "\tzams_filename = 'orig'" not in content

    def test_step_file_overwrite_pwarn(self, tmp_path):
        special = tmp_path / "special.inlist"
        special.write_text(
            "&controls\n  initial_mass = 10.0\n/ ! end\n&star_job\n/ ! end\n"
        )
        os.makedirs(os.path.join(tmp_path, "star1"), exist_ok=True)
        step_file = os.path.join(tmp_path, "star1", "inlist_step0")
        Path(step_file).write_text("old")
        mesa_inlists = {
            "single_star_grid": True,
            "zams_filename_1": None,
            "zams_filename_2": None,
            "star1_controls_special": str(special),
            "star1_job_special": str(special),
        }
        with pytest.warns(OverwriteWarning):
            result = totest.build_star_formation(
                "star1", mesa_inlists, str(tmp_path), {}
            )
        assert result == " {0}".format(step_file)
        content = Path(step_file).read_text()
        assert "\tinitial_mass = 10.0" in content


class TestApplyOutputControls:
    """Tests for apply_output_controls()."""

    def test_all_flags(self):
        bc = {}
        s1c = {}
        s2c = {}
        s1j = {}
        s2j = {}
        mesa_inlists = {
            "final_profile_star1": True,
            "final_profile_star2": True,
            "final_model_star1": True,
            "final_model_star2": True,
            "history_star1": False,
            "history_star2": False,
            "history_interval": 50,
            "binary_history": False,
            "zams_filename_1": "/z1",
            "zams_filename_2": "/z2",
        }
        totest.apply_output_controls(bc, s1c, s2c, s1j, s2j, mesa_inlists)
        assert s1j["write_profile_when_terminate"] == ".true."
        assert s1j["filename_for_profile_when_terminate"] == "'final_profile_star1.data'"
        assert s2j["write_profile_when_terminate"] == ".true."
        assert s2j["filename_for_profile_when_terminate"] == "'final_profile_star2.data'"
        assert s1j["save_model_when_terminate"] == ".true."
        assert s1j["save_model_filename"] == "'final_star1.mod'"
        assert s2j["save_model_when_terminate"] == ".true."
        assert s2j["save_model_filename"] == "'final_star2.mod'"
        assert s1c["do_history_file"] == ".false."
        assert s2c["do_history_file"] == ".false."
        assert bc["history_interval"] == "-1"
        assert s1c["history_interval"] == 50
        assert s2c["history_interval"] == 50
        assert s1c["zams_filename"] == "'/z1'"
        assert s2c["zams_filename"] == "'/z2'"


class TestWriteBinaryInlists:
    """Tests for write_binary_inlists()."""

    def test_overwrite_pwarn(self, tmp_path):
        os.makedirs(os.path.join(tmp_path, "binary"), exist_ok=True)
        for name in ["inlist_project", "inlist1", "inlist2"]:
            Path(tmp_path, "binary", name).write_text("old")
        with pytest.warns(OverwriteWarning):
            result = totest.write_binary_inlists(
                str(tmp_path),
                {"a": "1"},
                {"b": "2"},
                {"c": "3"},
                {"d": "4"},
                {"e": "5"},
                {"f": "6"},
            )
        assert result == (
            os.path.join(tmp_path, "binary", "inlist_project"),
            os.path.join(tmp_path, "binary", "inlist1"),
            os.path.join(tmp_path, "binary", "inlist2"),
        )


class TestConstructStaticInlistBackCompat:
    """Tests for the back-compat branches of construct_static_inlist()."""

    def test_zams_filename_backcompat(self, tmp_path, inlist_files):
        for sub in ["binary", "star1", "star2"]:
            os.makedirs(os.path.join(tmp_path, sub), exist_ok=True)
        mesa_inlists = _binary_mesa_inlists(inlist_files)
        del mesa_inlists["zams_filename_1"]
        del mesa_inlists["zams_filename_2"]
        del mesa_inlists["single_star_grid"]
        mesa_inlists["zams_filename"] = "/backcompat/zams.data"
        totest.construct_static_inlist(mesa_inlists, [], working_directory=str(tmp_path))
        assert mesa_inlists["zams_filename_1"] == "/backcompat/zams.data"
        assert mesa_inlists["zams_filename_2"] == "/backcompat/zams.data"
        inlist1 = Path(tmp_path, "binary", "inlist1").read_text()
        assert "\tzams_filename = '/backcompat/zams.data'" in inlist1

    def test_missing_keys_defaulted(self, tmp_path, inlist_files):
        for sub in ["binary", "star1", "star2"]:
            os.makedirs(os.path.join(tmp_path, sub), exist_ok=True)
        mesa_inlists = _binary_mesa_inlists(inlist_files)
        del mesa_inlists["zams_filename_1"]
        del mesa_inlists["zams_filename_2"]
        del mesa_inlists["single_star_grid"]
        totest.construct_static_inlist(mesa_inlists, [], working_directory=str(tmp_path))
        assert mesa_inlists["zams_filename_1"] is None
        assert mesa_inlists["zams_filename_2"] is None
        assert mesa_inlists["single_star_grid"] is False
