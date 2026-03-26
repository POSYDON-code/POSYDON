"""Unit tests of posydon/popsyn/synthetic_population.py
"""

__authors__ = [
    "Elizabeth Teng <elizabethteng@u.northwestern.edu>"
]

# import the module which will be tested
import posydon.popsyn.synthetic_population as totest
from posydon.utils.constants import Zsun

# aliases
np = totest.np
pd = totest.pd

import warnings
from inspect import isclass, isroutine

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import approx, fixture, raises, warns

warnings.simplefilter("always")
import os
import shutil

from posydon.unit_tests._helper_functions_for_tests.population import (
    make_ini,
    make_test_pop,
)


# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['DFInterface','History','Oneline',
                    'Population','PopulationIO','PopulationRunner',
                    'Rates','TransientPopulation',
                    '__authors__','__builtins__', '__cached__', '__doc__',
                    '__file__','__loader__', '__name__', '__package__', '__spec__',
                    'np', 'pd', 'tqdm', 'os', 'shutil','plt',
                    'Zsun', 'binarypop_kwargs_from_ini','plot_pop','SimulationProperties',
                    'calculate_model_weights','saved_ini_parameters',
                    'convert_metallicity_to_string','Pwarn','cosmology','const',
                    'get_shell_comoving_volume', 'get_comoving_distance_from_redshift',
                    'get_cosmic_time_from_redshift', 'redshift_from_cosmic_time_interpolator',
                    'DEFAULT_SFH_MODEL', 'get_redshift_bin_edges',
                    'get_redshift_bin_centers', 'SFR_per_met_at_z',
                    'BinaryPopulation', 'HISTORY_MIN_ITEMSIZE','ONELINE_MIN_ITEMSIZE'
                   ]
        totest_elements = set(dir(totest))
        missing_in_test = set(elements) - totest_elements
        assert len(missing_in_test) == 0, "There are missing objects in "\
                                          +f"{totest.__name__}: "\
                                          +f"{missing_in_test}. Please "\
                                          +"check, whether they have been "\
                                          +"removed on purpose and update "\
                                          +"this unit test."
        new_in_test = totest_elements - set(elements)
        assert len(new_in_test) == 0, "There are new objects in "\
                                      +f"{totest.__name__}: {new_in_test}. "\
                                      +"Please check, whether they have been "\
                                      +"added on purpose and update this "\
                                      +"unit test."

class TestPopulationRunner:

    def test_init(self):
        # missing argument
        with raises(TypeError,match="missing 1 required positional argument: 'path_to_ini'"):
            totest.PopulationRunner()
        # bad input
        with raises(ValueError, match="You did not provide a valid path_to_ini!"):
            totest.PopulationRunner("test")

    def test_evolve(self,tmp_path,monkeypatch):
        # mock dependencies
        class DummyPop:
            def __init__(self, **kwargs):
                self.kwargs = kwargs
                self.comm = None
                self.metallicity = kwargs["metallicity"]
            def evolve(self,**kwargs):
                self.evolved = True
            def combine_saved_files(self, *args):
                self.combined = True
        def dummy_kwargs(path):
            return {
                "metallicity": 0.1,
                "temp_directory": "tmp_dir",
                "verbose": False}
        def dummy_kwargs_list(path):
            return {
                "metallicity": [0.1,1.],
                "temp_directory": "tmp_dir",
                "verbose": False}
        def dummy_merge(pop,overwrite):
            pop.merged = True

        # Mock out functions
        monkeypatch.setattr(totest, "binarypop_kwargs_from_ini", dummy_kwargs)
        monkeypatch.setattr(totest, "BinaryPopulation", DummyPop)
        monkeypatch.setattr(totest, "convert_metallicity_to_string", lambda x: "0.1")
        run = totest.PopulationRunner(make_ini(tmp_path))
        # overwrite=False, directory doesn't exist
        monkeypatch.setattr(os.path, "exists", lambda path: False)
        run.merge_parallel_runs = dummy_merge
        run.evolve()
        assert run.binary_populations[0].evolved is True
        assert run.binary_populations[0].merged is True
        # overwrite=False, directory exists
        monkeypatch.setattr(os.path, "exists", lambda path: True)
        monkeypatch.setattr(totest, "binarypop_kwargs_from_ini", dummy_kwargs_list)
        run = totest.PopulationRunner(make_ini(tmp_path), verbose=True)
        with raises(FileExistsError, match="tmp_dir"):
            run.evolve(overwrite=False)
        # overwrite=True, directory exists
        removed = {}
        monkeypatch.setattr(shutil, "rmtree", lambda path: removed.setdefault("called", path))
        run.merge_parallel_runs = dummy_merge
        run.evolve(overwrite=True)
        assert removed["called"] == "0.1_Zsun_tmp_dir"
        assert run.binary_populations[0].evolved is True
        assert run.binary_populations[0].merged is True

    def test_merge_parallel_runs(self, tmp_path, monkeypatch, capsys):
        class DummyPop:
            def __init__(self, metallicity, temp_directory,**kwargs):
                self.metallicity = metallicity
                self.kwargs = {"temp_directory": temp_directory}
                self.combine_args = None
                self.combined = False

            def combine_saved_files(self, out_path, files):
                self.combine_args = (out_path, files)
                self.combined = True

        def dummy_kwargs(path):
            return {
                "metallicity": 0.1,
                "temp_directory": "tmp_dir",
                "verbose": False}

        monkeypatch.setattr(totest, "binarypop_kwargs_from_ini", dummy_kwargs)
        monkeypatch.setattr(totest, "BinaryPopulation", DummyPop)
        monkeypatch.setattr(totest, "convert_metallicity_to_string",
                            lambda x: str(os.path.join(tmp_path, "0.1")))

        # 1) File exists case: should raise FileExistsError
        pop = DummyPop(metallicity=0.1, temp_directory=str(tmp_path))
        output_file = os.path.join(tmp_path,"0.1_Zsun_population.h5")
        with open(output_file, "w") as f:
            f.write("test")
        run = totest.PopulationRunner(make_ini(tmp_path.parent))
        run.verbose = False
        with raises(FileExistsError, match="Files were not merged"):
            run.merge_parallel_runs(pop)

        # 2) Normal merge case
        file1 = os.path.join(tmp_path,"file1.tmp")
        file2 = os.path.join(tmp_path,"file2.tmp")
        output_file = os.path.join(tmp_path,"0.1_Zsun_population.h5")
        with open(file1, "w") as f:
            f.write("test")
        with open(file2, "w") as f:
            f.write("test")
        pop = DummyPop(metallicity=0.1, temp_directory=str(tmp_path))
        run = totest.PopulationRunner(make_ini(tmp_path.parent))
        run.verbose = True
        monkeypatch.setattr(totest, "convert_metallicity_to_string", lambda x: "0.1")
        run.merge_parallel_runs(pop)
        assert pop.combined is True
        out_path, files = pop.combine_args
        assert out_path == "0.1_Zsun_population.h5"
        # Filter out output file if somehow included (defensive)
        filtered_files = [f for f in files if os.path.basename(f) != out_path]
        assert set(os.path.basename(f) for f in filtered_files) == {"file1.tmp", "file2.tmp"}
        captured = capsys.readouterr()
        assert "Merging" in captured.out
        assert "Files merged!" in captured.out
        assert f"Removing files in {tmp_path}" in captured.out
        # cleanup
        for f in [file1, file2, output_file]:
            if os.path.exists(f):
                os.remove(f)
        assert len(os.listdir(tmp_path)) == 0

        run.verbose = False
        run.merge_parallel_runs(pop)
        assert not os.path.exists(pop.kwargs["temp_directory"])

class TestDFInterface:

    def test_head_tail_select(self, tmp_path):
        # Setup test HDF5 file
        data = pd.DataFrame({
            "index": np.repeat(np.arange(5), 2),
            "time": np.random.rand(10),
            "value": np.random.rand(10)
        })
        hdf_path = os.path.join(tmp_path,"test.h5")
        data.to_hdf(hdf_path, key="history", format="table", index=False)

        dfi = totest.DFInterface()
        dfi.filename = str(hdf_path)
        dfi.chunksize = 3

        head = dfi.head("history", n=3)
        tail = dfi.tail("history", n=2)
        subset = dfi.select("history", columns=["time"])

        assert len(head) == 3
        assert len(tail) == 2
        assert "time" in subset.columns
        assert subset.shape[1] == 1

    def test_repr_methods(self, tmp_path):
        df = pd.DataFrame({"index": range(10), "x": np.random.rand(10)})
        path = os.path.join(tmp_path, "test_repr.h5")
        df.to_hdf(path, key="history", format="table", index=False)

        dfi = totest.DFInterface()
        dfi.filename = str(path)

        s = dfi.get_repr("history")
        html = dfi.get_html_repr("history")

        assert isinstance(s, str)
        assert "x" in s
        assert isinstance(html, str)
        assert "<table" in html

class TestHistory:

    def test_init(self, tmp_path):
        with raises(FileNotFoundError, match="does not exist!"):
            totest.History("nonexistent_file.h5")

        df = pd.DataFrame({
            "binary_index": np.repeat(np.arange(3), 2),
            "a": np.random.rand(6),
        })
        file_path = os.path.join(tmp_path, "test_history.h5")
        df.set_index("binary_index", inplace=True)
        df.to_hdf(file_path, key="history", format="table")

        hist = totest.History(str(file_path), verbose=True, chunksize=2)

        assert hist.filename == str(file_path)
        assert hist.lengths is not None
        assert hist.columns == ["a"]
        assert hist.number_of_systems == 3
        assert isinstance(hist.indices, np.ndarray)

    def test_getitem_and_len(self, tmp_path):
        df = pd.DataFrame({
            "binary_index": np.repeat(np.arange(3), 2),
            "a": np.random.rand(6),
        })
        file_path = os.path.join(tmp_path, "test_getitem.h5")
        df.set_index("binary_index", inplace=True)
        df.to_hdf(file_path, key="history", format="table")

        lengths_df = pd.DataFrame({'lengths': [2, 2, 2]}, index=[0, 1, 2])
        lengths_df.to_hdf(file_path, key="history_lengths")

        hist = totest.History(str(file_path))

        # adding history_lengths
        assert hist.lengths.equals(lengths_df)

        # __len__
        assert len(hist) == 6

        # __getitem__ with "none" slice indices
        out_startnone = hist[:3]
        out_stopnone = hist[2:]
        assert not out_startnone.empty
        assert all(i in out_startnone.index for i in range(3))
        assert not out_stopnone.empty
        assert all(i>=2 for i in out_stopnone.index)

        # __getitem__ with int
        out = hist[0]
        assert isinstance(out, pd.DataFrame)

        # __getitem__ with list of int
        out = hist[[0, 1]]
        assert isinstance(out, pd.DataFrame)
        out_none = hist[[]]
        assert isinstance(out_none,pd.DataFrame)
        assert out_none.empty

        # __getitem__ with numpy array of int
        out = hist[np.array([0, 1])]
        assert isinstance(out, pd.DataFrame)
        out_none = hist[np.array([], dtype=int)]
        assert isinstance(out_none,pd.DataFrame)
        assert out_none.empty

        # __getitem__ with bool array
        full_data = pd.read_hdf(file_path, key="history")
        mask = full_data["a"] > -1
        empty_mask = np.array([],dtype=bool)
        out = hist[mask.to_numpy()]
        out_none = hist[empty_mask]
        assert not out.empty
        assert (out["a"] > -1).all()
        assert isinstance(out_none,pd.DataFrame)
        assert out_none.empty

        # __getitem__ with str column
        out = hist["a"]
        assert "a" in out.columns

        # __getitem__ with list of str columns
        out = hist[["a"]]
        assert "a" in out.columns

        # Invalid column
        with raises(ValueError, match="is not a valid column name"):
            hist["bad_column"]

        # Invalid list of column names
        with raises(ValueError, match="Not all columns in"):
            hist[["a", "bad"]]

        # Invalid type
        with raises(ValueError, match="Invalid key type"):
            hist[None]

    def test_slice(self, tmp_path):
        df = pd.DataFrame({
            "index": np.repeat(np.arange(5), 2),
            "val": np.random.rand(10)
        })
        path = os.path.join(tmp_path, "test_slice.h5")
        df.to_hdf(path, key="history", format="table", index=False)
        hist = totest.History(str(path))
        sliced = hist[1:3]
        assert isinstance(sliced, pd.DataFrame)

    def test_head_tail_repr(self, tmp_path):
        df = pd.DataFrame({
            "index": np.repeat(np.arange(5), 2),
            "val": np.random.rand(10)
        })
        path = os.path.join(tmp_path, "test_repr2.h5")
        df.to_hdf(path, key="history", format="table", index=False)
        hist = totest.History(str(path))

        head = hist.head(n=3)
        tail = hist.tail(n=2)
        rep = repr(hist)
        html = hist._repr_html_()

        assert isinstance(head, pd.DataFrame)
        assert isinstance(tail, pd.DataFrame)
        assert isinstance(rep, str)
        assert isinstance(html, str)
        assert "val" in rep
        assert "<table" in html

class TestOneline:

    def test_init(self, tmp_path):
        with raises(FileNotFoundError, match="does not exist!"):
            totest.Oneline("nonexistent_file.h5")

        df = pd.DataFrame({
            "index": np.arange(5),
            "a": np.random.rand(5),
            "b": np.random.rand(5),
        })
        file_path = os.path.join(tmp_path, "test_oneline.h5")
        df.set_index("index", inplace=True)
        df.to_hdf(file_path, key="oneline", format="table")

        one = totest.Oneline(str(file_path), verbose=True, chunksize=2)

        assert one.filename == str(file_path)
        assert one.columns == ["a", "b"]
        assert one.number_of_systems == 5
        assert isinstance(one.indices, np.ndarray)

    def test_getitem_variants(self, tmp_path):
        df = pd.DataFrame({
            "index": np.arange(6),
            "a": np.random.rand(6),
            "b": np.random.rand(6),
        })
        file_path = os.path.join(tmp_path, "test_getitem_oneline.h5")
        df.set_index("index", inplace=True)
        df.to_hdf(file_path, key="oneline", format="table")

        one = totest.Oneline(str(file_path))

        # slice with start/stop
        out = one[1:4]
        assert isinstance(out, pd.DataFrame)
        assert all(i in out.index for i in range(1, 4))

        # start=None
        out = one[:3]
        assert isinstance(out, pd.DataFrame)
        assert all(i in out.index for i in range(3))

        # stop=None
        out = one[3:]
        assert isinstance(out, pd.DataFrame)
        assert all(i >= 3 for i in out.index)

        # int
        out = one[2]
        assert isinstance(out, pd.DataFrame)
        assert 2 in out.index

        # list of int
        out = one[[1, 3, 5]]
        assert set(out.index) == {1, 3, 5}

        # numpy array of int
        out = one[np.array([0, 2, 4])]
        assert set(out.index) == {0, 2, 4}

        # numpy array of bool
        mask = np.array([True, False, True, False, True, False])
        out = one[mask]
        assert set(out.index) == {0, 2, 4}

        # pandas DataFrame mask of bool
        mask_df = pd.DataFrame({"mask": [True, False, True, False, True, False]})
        out = one[mask_df]
        assert set(out.index) == {0, 2, 4}

        # str column
        out = one["a"]
        assert list(out.columns) == ["a"]

        # list of str columns
        out = one[["a", "b"]]
        assert list(out.columns) == ["a", "b"]

        # invalid column
        with raises(ValueError, match="is not a valid column"):
            one["bad"]

        # invalid list of str columns
        with raises(ValueError, match="Not all columns in"):
            one[["a", "bad"]]

        # invalid type
        with raises(ValueError, match="Invalid key type"):
            one[None]

    def test_invalid_float_list(self, tmp_path):
        df = pd.DataFrame({
            "index": np.arange(3),
            "val": np.random.rand(3)
        })
        file_path = os.path.join(tmp_path, "test_invalid_float_list.h5")
        df.set_index("index", inplace=True)
        df.to_hdf(file_path, key="oneline", format="table")

        one = totest.Oneline(str(file_path))

        with raises(ValueError, match="elements in list are not integers"):
            one[[1.1, 2.2]]

class TestPopulationIO:

    @fixture
    def popio(self):
        p = totest.PopulationIO()
        p.mass_per_metallicity = pd.DataFrame({"Z": [0.02], "mass": [1.0]})
        p.ini_params = {"param1": 42, "param2": "abc"}
        return p

    def test_load_metadata(self,monkeypatch,popio):
        # bad input
        with raises(ValueError,match='does not contain .h5'):
            popio._load_metadata("not_pop.txt")
        # examples
        called={}
        monkeypatch.setattr(popio, "_load_ini_params", lambda f: called.setdefault("ini", f))
        monkeypatch.setattr(popio, "_load_mass_per_metallicity", lambda f: called.setdefault("mass", f))
        popio._load_metadata("file.h5")
        assert called == {"ini": "file.h5", "mass": "file.h5"}

    def test_save_mass_per_metallicity(self, tmp_path, popio):
        filename = os.path.join(tmp_path, "mass.h5")
        popio._save_mass_per_metallicity(filename)

        with pd.HDFStore(filename, "r") as store:
            df = store["mass_per_metallicity"]

        pd.testing.assert_frame_equal(df, popio.mass_per_metallicity)

    def test_save_ini_params(self,tmp_path,popio,monkeypatch):
        filename = os.path.join(tmp_path, "ini_out.h5")
        monkeypatch.setattr("posydon.popsyn.synthetic_population.saved_ini_parameters", ["param1", "param3"])

        popio._save_ini_params(filename)

        with pd.HDFStore(filename, "r") as store:
            df = store["ini_parameters"]

        assert list(df.columns) == ["param1"]
        assert df["param1"][0] == 42
    def test_load_ini_params(self,tmp_path,popio,monkeypatch):
        filename = os.path.join(tmp_path, "ini.h5")
        monkeypatch.setattr("posydon.popsyn.synthetic_population.saved_ini_parameters", ["param1", "param2", "param3"])

        df = pd.DataFrame({"param1": [1], "param2": ["x"]})
        with pd.HDFStore(filename, "w") as store:
            store.put("ini_parameters", df)

        popio._load_ini_params(filename)

        assert popio.ini_params["param1"] == 1
        assert popio.ini_params["param2"] == "x"
        assert "param3" not in popio.ini_params

class TestPopulation:

    def test_population_init(self, tmp_path, monkeypatch):

        # bad input
        with raises(ValueError, match="does not contain .h5"):
            totest.Population("hello.txt")

        # missing /history
        filename = os.path.join(tmp_path, "pop_missing.h5")
        with pd.HDFStore(filename, "w") as store:
            store.put("ini_parameters", pd.DataFrame({"Parameter": [], "Value": []}), format="table")
        with raises(ValueError, match="does not contain a history table"):
            totest.Population(str(filename))

        # /history exists, /oneline missing
        history_df = pd.DataFrame({"binary_index": [0], "event": [0], "time": [0.0]})
        with pd.HDFStore(filename, "a") as store:
            store.put("history", history_df, format="table")
        with raises(ValueError, match="does not contain an oneline table"):
            totest.Population(str(filename))

        # /history and /oneline exist, no ini_parameters
        oneline_df = pd.DataFrame({
            "binary_index": [0],
            "S1_mass_i": [1],
            "S2_mass_i": [1],
            "state_i": ["initially_single_star"],
            "metallicity": [0.02]
            })
        with pd.HDFStore(filename, "a") as store:
            store.put("oneline", oneline_df, format="table")
            store.put("mass_per_metallicity", pd.DataFrame({"simulated_mass": [0]}, index=[0.02]), format="table")
        with raises(ValueError, match='does not contain an ini_parameters table'):
            totest.Population(str(filename))

        # /history and /oneline exist, yes ini_parameters, no mass_per_metallicity
        filename_no_mass = os.path.join(tmp_path, "pop_no_mass.h5")
        with pd.HDFStore(filename_no_mass, "w") as store:
            store.put("history", history_df, format="table")
            store.put("oneline", oneline_df, format="table")
            store.put("ini_parameters", pd.DataFrame({"Parameter": ["metallicity"], "Value": [0.02]}), format="table")
        with raises(ValueError, match='does not contain a mass_per_metallicity table'):
            totest.Population(str(filename_no_mass))

        # metallicity specified
        monkeypatch.setattr(
            "posydon.popsyn.synthetic_population.binarypop_kwargs_from_ini",
            lambda ini_file: {"dummy_param": 1},
        )

        pop_with_metallicity = totest.Population(
            str(filename_no_mass), metallicity=0.02, ini_file=str(tmp_path / "dummy.ini")
        )
        assert pop_with_metallicity.mass_per_metallicity is not None
        assert pop_with_metallicity.solar_metallicities[0] == 0.02
        assert pop_with_metallicity.metallicities[0] == 0.02 * Zsun

        # everything exists
        pop = make_test_pop(tmp_path, filename="full_pop.h5")
        assert pop.number_of_systems > 0
        assert isinstance(pop.history, totest.History)
        assert isinstance(pop.oneline, totest.Oneline)

        pop_with_metallicity = totest.Population(
            str(pop.filename), metallicity=0.02, ini_file=str(tmp_path / "dummy.ini")
        )
        assert pop_with_metallicity.mass_per_metallicity is not None
        assert pop_with_metallicity.solar_metallicities[0] == 0.02
        assert pop_with_metallicity.metallicities[0] == 0.02 * Zsun

    def test_export_selection(self, tmp_path, monkeypatch):
        pop = make_test_pop(tmp_path)
        export_file = tmp_path / "exp.h5"

        # bad input
        with raises(ValueError, match='does not contain .h5'):
            pop.export_selection([0], 'hello.txt')

        with raises(ValueError, match="Both overwrite and append cannot be True!"):
            pop.export_selection([0], str(export_file), append=True, overwrite=True)

        dummy_file = tmp_path / "exists.h5"
        pd.DataFrame({"a": [1]}).to_hdf(dummy_file, "dummy", format="table")
        with raises(FileExistsError, match='Set overwrite or append to True'):
            pop.export_selection([0], str(dummy_file), overwrite=False, append=False)

        # overwrite
        out_file = tmp_path / "out.h5"
        pop.export_selection([0], str(out_file), overwrite=True, history_chunksize=1)

        # append
        pop.export_selection([0], str(out_file), append=True, history_chunksize=1)

        # write export
        pop.export_selection(
            [0], os.path.join(tmp_path, 'new.h5'), append=False, overwrite=False, history_chunksize=1
        )

        # test case: oneline missing metallicity
        class DummyOnelineNoMetal:
            columns = ["S1_mass_i", "S2_mass_i", "state_i"]
            number_of_systems = 1
            def __getitem__(self, cols):
                return pd.DataFrame({
                    "S1_mass_i": [1], "S2_mass_i": [1], "state_i": ["initial"]
                }, index=[0])
            def __len__(self): return 1

        pop.oneline = DummyOnelineNoMetal()
        pop.export_selection([0], str(tmp_path / "out2.h5"), overwrite=True)

        # mass_per_metallicity updated
        df = pd.read_hdf(out_file, "mass_per_metallicity")
        assert "number_of_systems" in df.columns

        # multiple metallicities error
        class DummyNoMet:
            columns = ["foo"]
            number_of_systems = 1
            def __getitem__(self, idx): return pd.DataFrame({"foo": [1]})
            def __len__(self): return 1

        pop.oneline = DummyNoMet()
        pop.metallicities = [0.02, 0.01]

        with raises(ValueError, match="multiple metallicities"):
            pop.export_selection([0], str(tmp_path / "multi_met.h5"), overwrite=True)

    def test_calculate_formation_channels(self, tmp_path):
        pop = make_test_pop(tmp_path)

        class DummyOneline:
            columns = ["interp_class_HMS_HMS", "mt_history_HMS_HMS"]
            number_of_systems = 4

            def select(self, start=None, stop=None, columns=None):
                data = [
                    {"interp_class_HMS_HMS": "initial_MT", "mt_history_HMS_HMS": "Stable contact phase"},
                    {"interp_class_HMS_HMS": "stable_MT", "mt_history_HMS_HMS": None},
                    {"interp_class_HMS_HMS": "stable_reverse_MT", "mt_history_HMS_HMS": None},
                    {"interp_class_HMS_HMS": "no_MT", "mt_history_HMS_HMS": None},
                ]
                selected = data[start:stop]
                while len(selected) < (stop - start):
                    selected.append(data[-1])
                df = pd.DataFrame(selected)
                if columns is not None:
                    df = df[columns]
                return df

        pop.oneline = DummyOneline()
        pop.chunksize = 2

        pop.calculate_formation_channels(mt_history=True)
        assert hasattr(pop, "formation_channels")
        assert all(col in pop.formation_channels.columns for col in ["channel", "channel_debug"])
        assert any("contact" in str(c) for c in pop.formation_channels["channel"])

        pop.calculate_formation_channels(mt_history=False)
        assert hasattr(pop, "formation_channels")
        assert "channel" in pop.formation_channels.columns

        pop.calculate_formation_channels(mt_history=True)
        with pd.HDFStore(pop.filename, "r") as store:
            assert "/formation_channels" in store.keys()

    def test_create_transient_population(self):
        # missing argument
        # bad input
        # examples
        pass

class TestTransientPopulation:

    def test_select(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_calculate_model_weights(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_calculate_cosmic_weights(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_efficiency(self):
        # missing argument
        # bad input
        # examples
        pass

class TestRates:

    def test_select_rate_slice(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_calculate_intrinsic_rate_density(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_calculate_observable_population(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_observable_population(self):
        # missing argument
        # bad input
        # examples
        pass
    def test_edges_metallicity_bins(self):
        # TODO: needs Rates object setup
        pass
