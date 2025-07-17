"""Unit tests of posydon/popsyn/transient_select_funcs.py
"""

__authors__ = [
    "Elizabeth Teng <elizabethteng@u.northwestern.edu>"
]

# import the module which will be tested
import posydon.popsyn.transient_select_funcs as totest
# aliases
np = totest.np
pd = totest.pd
PATH_TO_POSYDON_DATA = totest.PATH_TO_POSYDON_DATA

# import other needed code for the tests, which is not already imported in the
# module you like to test
from pytest import fixture, raises, warns, approx
from inspect import isroutine, isclass
import posydon.popsyn.selection_effects as selection_effects
from posydon.utils.posydonwarning import ReplaceValueWarning
import warnings
warnings.simplefilter("always")

# define test classes collecting several test functions
class TestElements:
    # check for objects, which should be an element of the tested module
    def test_dir(self):
        elements = ['PATH_TO_PDET_GRID', 'GRB_selection', 'chi_eff', 'm_chirp', \
                    'mass_ratio', 'BBH_selection_function','DCO_detectability', \
                    '__builtins__', '__cached__', '__doc__', \
                    '__file__','__loader__', '__name__', '__package__', '__spec__', \
                    'np', 'pd', 'PATH_TO_POSYDON_DATA', \
                    'os', 'tqdm', 'warnings', 'Pwarn','selection_effects']
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

    def test_instance_GRB_selection(self):
        assert isroutine(totest.GRB_selection)

    def test_instance_chi_eff(self):
        assert isroutine(totest.chi_eff)

    def test_instance_m_chirp(self):
        assert isroutine(totest.m_chirp)

    def test_instance_mass_ratio(self):
        assert isroutine(totest.mass_ratio)

    def test_instance_BBH_selection_function(self):
        assert isroutine(totest.BBH_selection_function)

    def test_instance_DCO_detectability(self):
        assert isroutine(totest.DCO_detectability)
        
class TestFunctions:
    
    @fixture
    def history_chunk(self):
        return pd.DataFrame({
            'binary_index': [10, 10, 10],
            'S1_state': ['MS', 'HG', 'BH'],
            'S2_state': ['MS', 'HG', 'BH'],
            'step_names': ['step_RLO', 'step_RLO', 'step_SN'],
            'orbital_period': [1.0, 1.1, 1.2],
            'eccentricity': [0.1, 0.15, 0.2],
            'S1_spin': [0.3, 0.35, 0.4],
            'S2_spin': [0.5, 0.55, 0.6],
            'S1_mass': [10.0, 9.5, 9.0],
            'S2_mass': [8.0, 7.5, 7.0],
            'time': [1.0e6, 2.0e6, 3.0e6]})

    @fixture
    def oneline_chunk(self):
        return pd.DataFrame({
            'metallicity': [0.02],
            'S1_m_disk_radiated': [0.5],
            'S2_m_disk_radiated': [0.0],
        }, index=[10])
    
    @fixture
    def formation_channels_chunk(self):
        return pd.DataFrame({
            'channel': ['foo_CC1','bar_CC2']
        }, index=[10,11])
    
    @fixture
    def history_BBH(self):
        return pd.DataFrame({
            'event': ['MID', 'END', 'END'],
            'time': [1e6, 5e6, 6e6],
            'S1_state': ['BH', 'BH', 'BH'],
            'S2_state': ['BH', 'BH', 'BH'],
            'step_names': ['step_SN', 'step_SN', 'step_SN'],
            'state': ['detached', 'detached', 'detached'],
            'S1_mass': [30, 35, 40],
            'S2_mass': [20, 25, 30],
            'S1_spin': [0.5, 0.6, 0.7],
            'S2_spin': [0.4, 0.3, 0.2],
            'orbital_period': [0.5, 0.6, 0.7],
            'eccentricity': [0.1, 0.2, 0.3],
        }, index=[0,1,2])

    @fixture
    def oneline_BBH(self):
        return pd.DataFrame({
            'metallicity': [0.01, 0.02, 0.03],
            'S1_spin_orbit_tilt_second_SN': [0.1, 0.2, 0.3],
            'S2_spin_orbit_tilt_second_SN': [0.4, 0.5, 0.6],
        }, index=[0,1,2])

    @fixture
    def formation_channels_BBH(self):
        return pd.DataFrame({
            'channel': ['foo', 'bar', 'baz'],
        }, index=[0,1,2])

    @fixture
    def array(self):
        return np.array([1.0,2.0,3.0])
    
    @fixture
    def nan_array(self):
        return np.array([np.nan,np.nan,np.nan])
    
    @fixture
    def wrong_array(self):
        return np.array(['1.0','2.0','3.0'])
    
    @fixture
    def transient_pop_chunk(self):
        return pd.DataFrame({
            'S1_mass': [30, 35],
            'S2_mass': [25, 30],
            'S1_spin': [0.1, 0.2],
            'S2_spin': [0.1, 0.2],
            'S1_spin_orbit_tilt_at_merger': [0.5, 0.6],
            'S2_spin_orbit_tilt_at_merger': [0.4, 0.5],
            'q': [0.83, 0.86],
            'chi_eff': [0.1, 0.2]})

    @fixture
    def z_events_chunk(self):
        return pd.DataFrame({
            'event_1': [0.1, np.nan],
            'event_2': [0.2, 0.3]})

    @fixture
    def z_events_chunk_with_nan(self):
        return pd.DataFrame({
            'event_1': [1.0, np.nan],   
            'event_2': [np.nan, np.nan] 
        }, index=[0,1])

    @fixture
    def z_weights_chunk(self):
        return pd.DataFrame({
            'event_1': [1.0, 1.0],
            'event_2': [1.0, 1.0]
        }, index=[0, 1])
    

    def test_GRB_selection(self,history_chunk,oneline_chunk,
                           formation_channels_chunk):
        # missing argument
        with raises(TypeError,match="missing 2 required positional arguments"):
            totest.GRB_selection()
        # bad input
        with raises(TypeError,match='string indices must be integers'):
            totest.GRB_selection("1.1", "1.2")
        with raises(AttributeError,match="'float' object has no attribute 'index'"):
            totest.GRB_selection(1.1, 1.2)
        with raises(ValueError,match='S1_S2 must be either S1 or S2'):
            totest.GRB_selection(history_chunk, oneline_chunk.copy(),
                                 S1_S2='test')
        # example with S1
        df = totest.GRB_selection(history_chunk, oneline_chunk.copy(),
                           formation_channels_chunk, S1_S2='S1')
        assert not df.empty
        assert df.index[0] == 10
        assert 'S1_mass_preSN' in df.columns
        assert 'S1_mass_postSN' in df.columns
        assert df['time'].iloc[0] == 3.0  # 3 Myr = 3e6 years * 1e-6
        assert df['channel'].iloc[0] == 'foo_CC1'
        # example with no formation channels
        df = totest.GRB_selection(history_chunk, oneline_chunk.copy(),
                           formation_channels_chunk=None, S1_S2='S1')
        assert not df.empty
        assert 'channel' not in df.columns 
        # example with S2
        chunk = oneline_chunk.copy()
        chunk['S1_m_disk_radiated'] = [0.0]
        chunk['S2_m_disk_radiated'] = [0.5]
        df = totest.GRB_selection(history_chunk, chunk,
                           formation_channels_chunk, S1_S2='S2')
        assert not df.empty
        assert 'S2_mass_postSN' in df.columns
        assert 'metallicity' in df.columns
        assert 'channel' in df.columns      
        # example with no disk radiation
        chunk = oneline_chunk.copy()
        chunk['S1_m_disk_radiated'] = [0.0]
        df = totest.GRB_selection(history_chunk, chunk, formation_channels_chunk=None,S1_S2='S1')
        assert df.empty
        
    def test_chi_eff(self,array,nan_array,wrong_array):
        # missing argument
        with raises(TypeError,match="missing 6 required positional arguments"):
            totest.chi_eff()
        # bad input
        with raises(TypeError,match="ufunc 'cos' not supported for the input types"):
            totest.chi_eff(array,array,array,array,array,wrong_array)
        # undefined values
        with warns(ReplaceValueWarning,match="a_1 contains undefined values"):
            totest.chi_eff(array,array,nan_array.copy(),array,array,array)
        with warns(ReplaceValueWarning,match="a_2 contains undefined values"):
            totest.chi_eff(array,array,array,nan_array.copy(),array,array)
        with warns(ReplaceValueWarning,match="tilt_1 contains undefined values"):
            totest.chi_eff(array,array,array,array,nan_array.copy(),array)
        with warns(ReplaceValueWarning,match="tilt_2 contains undefined values"):
            totest.chi_eff(array,array,array,array,array,nan_array.copy())
        # example
        assert totest.chi_eff(array,array,array,
                              array,array,array)[0] == 0.5403023058681398
        
    def test_m_chirp(self):
        # missing argument
        with raises(TypeError,match="missing 1 required positional argument: 'm_2'"):
            totest.m_chirp(3.)
        # bad input
        with raises(TypeError,match="can't multiply sequence by non-int of type 'str'"): 
            totest.m_chirp("3.","2.")
        # examples
        tests = [(4.,2.,2.433457367572823),
                 (40.,10.,16.65106414803746)]
        for (m1,m2,mc) in tests:
            assert totest.m_chirp(m1,m2) == mc
            
    def test_mass_ratio(self):
        # missing argument
        with raises(TypeError,match="missing 1 required positional argument: 'm_2'"):
            totest.mass_ratio(3.)
        # bad input
        with raises(TypeError,match="unsupported operand type"): 
            totest.mass_ratio("3.","2.")
        # examples
        tests = [(5.,1.,0.2),
                 (1.,5.,0.2),
                 (4.,2.,0.5)]
        for (m1,m2,q) in tests:
            assert totest.mass_ratio(np.array([m1]),
                                     np.array([m2])) == q
            
    def test_BBH_selection_function(self, history_BBH, oneline_BBH,
                                   formation_channels_BBH):
        # missing argument
        with raises(TypeError,match="missing 2 required positional arguments"):
            totest.BBH_selection_function()
        # bad input
        with raises(AttributeError,match="'float' object has no attribute 'index'"):
            totest.BBH_selection_function(1.1, 1.2)
        # example without formation channels 
        df = totest.BBH_selection_function(history_BBH, oneline_BBH)
        assert not df.empty
        assert all(col in df.columns for col in [
            'time', 't_inspiral', 'metallicity', 'S1_state', 'S2_state',
            'S1_mass', 'S2_mass', 'S1_spin', 'S2_spin',
            'S1_spin_orbit_tilt_at_merger', 'S2_spin_orbit_tilt_at_merger',
            'orbital_period', 'chirp_mass', 'mass_ratio', 'chi_eff', 'eccentricity'
        ])
        assert (df.index == oneline_BBH.index).all()
        assert df['t_inspiral'].iloc[1] == 0.0
        # example with formation channels 
        df = totest.BBH_selection_function(history_BBH, oneline_BBH, formation_channels_BBH)
        assert 'channel' in df.columns
        assert (df['channel'] == formation_channels_BBH['channel']).all()
    
    def test_DCO_detectability(self,
                               transient_pop_chunk,
                               z_events_chunk, 
                               z_events_chunk_with_nan,
                               z_weights_chunk,
                               monkeypatch):
        class FakeKNNmodel:
            def __init__(self, grid_path, sensitivity_key):
                pass
            def predict_pdet(self, df):
                # Return a fixed probability (e.g., 0.5) for each row in df
                return np.full(len(df), 0.5)

        monkeypatch.setattr('posydon.popsyn.selection_effects.KNNmodel',
                            FakeKNNmodel)

        # missing argument
        with raises(TypeError,match="missing 4 required positional arguments"):
            totest.DCO_detectability()
        # bad input
        with raises(ValueError,match='Unknown sensitivity sens_example'):
            totest.DCO_detectability("sens_example",
                                     transient_pop_chunk,
                                     z_events_chunk,
                                     z_weights_chunk)
        # example: basic functionality    
        out = totest.DCO_detectability('O3actual_H1L1V1', transient_pop_chunk, 
                                z_events_chunk, z_weights_chunk)
        assert isinstance(out, pd.DataFrame)
        assert out.shape == z_weights_chunk.shape
        assert (out.values <= 1.0).all()
        # example: missing q
        transient = transient_pop_chunk.drop(columns=['q'])
        out = totest.DCO_detectability('O3actual_H1L1V1', transient, 
                                z_events_chunk, z_weights_chunk)
        assert not out.empty
        # example: missing chi_eff
        transient = transient_pop_chunk.drop(columns=['chi_eff'])
        out = totest.DCO_detectability('O3actual_H1L1V1', transient, 
                                z_events_chunk, z_weights_chunk)
        assert not out.empty
        assert (out.values <= 1.0).all()
        # example: masking for nans in z_events_chunk
        out = totest.DCO_detectability('O3actual_H1L1V1',
                                       transient_pop_chunk,
                                       z_events_chunk_with_nan,
                                       z_weights_chunk)
        assert (out['event_2'] == 0.0).all()
