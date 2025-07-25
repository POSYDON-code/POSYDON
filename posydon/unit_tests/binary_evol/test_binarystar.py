"""Unit tests of posydon/binary_evol/binarystar.py
"""

__authors__ = [
    "Assistant <assistant@cline.bot>"
]

# import the module which will be tested
import posydon.binary_evol.binarystar as totest
from posydon.binary_evol.binarystar import BinaryStar, BINARYPROPERTIES, MAXIMUM_STEP_TIME
from posydon.binary_evol.singlestar import SingleStar, STARPROPERTIES
from posydon.binary_evol.simulationproperties import SimulationProperties

# import other needed code for the tests
import numpy as np
import pandas as pd
import pytest
from unittest.mock import Mock, patch, MagicMock
import copy
import signal
from posydon.utils.posydonerror import FlowError

class TestBinaryStarInit:
    """Test BinaryStar initialization."""
    
    def test_init_default(self):
        """Test BinaryStar initialization with default parameters."""
        binary = BinaryStar()
        
        # Check that binary has the correct attributes
        assert hasattr(binary, 'star_1')
        assert hasattr(binary, 'star_2')
        assert isinstance(binary.star_1, SingleStar)
        assert isinstance(binary.star_2, SingleStar)
        assert hasattr(binary, 'properties')
        assert isinstance(binary.properties, SimulationProperties)
        
        # Check binary properties are initialized
        for prop in BINARYPROPERTIES:
            assert hasattr(binary, prop)
            assert hasattr(binary, prop + '_history')
        
        # Check special default values
        assert np.array_equal(binary.V_sys, np.array([0.0, 0.0, 0.0]))
        assert binary.mass_transfer_case == 'None'
        assert binary.inspiral_time is None
        assert binary.index is None
    
    def test_init_with_stars(self):
        """Test BinaryStar initialization with provided stars."""
        star1 = SingleStar(mass=10.0, state='H-rich_Core_H_burning')
        star2 = SingleStar(mass=8.0, state='H-rich_Core_H_burning')
        
        binary = BinaryStar(star_1=star1, star_2=star2, index=42)
        
        assert binary.star_1 is star1
        assert binary.star_2 is star2
        assert binary.index == 42
    
    def test_init_with_binary_kwargs(self):
        """Test BinaryStar initialization with binary parameters."""
        binary = BinaryStar(
            time=1e6,
            separation=100.0,
            orbital_period=10.0,
            eccentricity=0.1,
            state='detached',
            event='ZAMS'
        )
        
        assert binary.time == 1e6
        assert binary.separation == 100.0
        assert binary.orbital_period == 10.0
        assert binary.eccentricity == 0.1
        assert binary.state == 'detached'
        assert binary.event == 'ZAMS'
        
        # Check history is initialized
        assert binary.time_history == [1e6]
        assert binary.separation_history == [100.0]
    
    def test_init_orbital_calculations(self):
        """Test automatic orbital parameter calculations."""
        star1 = SingleStar(mass=10.0)
        star2 = SingleStar(mass=8.0)
        
        # Test period calculation from separation
        binary1 = BinaryStar(star_1=star1, star_2=star2, separation=100.0)
        assert binary1.orbital_period is not None
        assert binary1.orbital_period > 0
        
        # Test separation calculation from period
        binary2 = BinaryStar(star_1=star1, star_2=star2, orbital_period=10.0)
        assert binary2.separation is not None
        assert binary2.separation > 0
    
    def test_init_grid_attributes(self):
        """Test initialization of grid-specific attributes."""
        binary = BinaryStar()
        
        grid_types = ['HMS_HMS', 'CO_HMS_RLO', 'CO_HeMS', 'CO_HeMS_RLO']
        for grid_type in grid_types:
            assert hasattr(binary, f'interp_class_{grid_type}')
            assert hasattr(binary, f'mt_history_{grid_type}')
            assert hasattr(binary, f'culmulative_mt_case_{grid_type}')
            assert getattr(binary, f'interp_class_{grid_type}') is None
    
    def test_init_with_custom_properties(self):
        """Test initialization with custom SimulationProperties."""
        properties = SimulationProperties()
        binary = BinaryStar(properties=properties)
        
        assert binary.properties is properties

class TestBinaryStarEvolution:
    """Test BinaryStar evolution methods."""
    
    def test_evolve_basic(self):
        """Test basic evolution setup and flow."""
        binary = BinaryStar(time=0.0, state='detached', event='ZAMS')
        
        # Mock properties to avoid actual evolution
        mock_properties = Mock()
        mock_properties.max_simulation_time = 1e10
        mock_properties.max_n_steps_per_binary = 5
        mock_properties.end_events = ['END']
        mock_properties.end_states = ['merged']
        mock_properties.pre_evolve = Mock()
        mock_properties.post_evolve = Mock()
        binary.properties = mock_properties
        
        # Mock run_step to set END event after one step
        def mock_run_step():
            binary.event = 'END'
        
        with patch.object(binary, 'run_step', side_effect=mock_run_step):
            binary.evolve()
        
        mock_properties.pre_evolve.assert_called_once_with(binary)
        mock_properties.post_evolve.assert_called_once_with(binary)
    
    def test_evolve_time_check(self):
        """Test evolution time validation."""
        binary = BinaryStar(time=1e15)  # Very large time
        binary.properties.max_simulation_time = 1e10
        
        with pytest.raises(ValueError, match="birth time.*greater than"):
            binary.evolve()
    
    def test_evolve_max_steps(self):
        """Test evolution with maximum steps limit."""
        binary = BinaryStar(time=0.0, state='detached', event='ZAMS')
        
        mock_properties = Mock()
        mock_properties.max_simulation_time = 1e10
        mock_properties.max_n_steps_per_binary = 2
        mock_properties.end_events = ['END']
        mock_properties.end_states = ['merged']
        mock_properties.pre_evolve = Mock()
        mock_properties.post_evolve = Mock()
        binary.properties = mock_properties
        
        # Mock run_step to never set END event
        step_count = 0
        def mock_run_step():
            nonlocal step_count
            step_count += 1
            binary.event = 'continuing'  # Never ends
        
        with patch.object(binary, 'run_step', side_effect=mock_run_step):
            with pytest.raises(RuntimeError, match="Exceeded maximum number of steps"):
                binary.evolve()
    
    def test_run_step_undefined_state(self):
        """Test run_step with undefined state."""
        binary = BinaryStar()
        binary.star_1.state = 'undefined_state'
        binary.star_2.state = 'another_undefined_state'
        binary.state = 'undefined_binary_state'
        binary.event = 'undefined_event'
        
        # Mock UNDEFINED_STATES to include our test state
        undefined_state = ('undefined_state', 'another_undefined_state', 
                          'undefined_binary_state', 'undefined_event')
        
        with patch('posydon.binary_evol.binarystar.UNDEFINED_STATES', [undefined_state]):
            with pytest.raises(FlowError, match="known undefined state"):
                binary.run_step()
    
    def test_run_step_no_next_step(self):
        """Test run_step when no next step is defined."""
        binary = BinaryStar()
        binary.star_1.state = 'H-rich_Core_H_burning'
        binary.star_2.state = 'H-rich_Core_H_burning'
        binary.state = 'detached'
        binary.event = 'ZAMS'
        
        # Mock properties with empty flow
        mock_properties = Mock()
        mock_properties.flow = {}
        binary.properties = mock_properties
        
        with pytest.raises(ValueError, match="Undefined next step"):
            binary.run_step()
    
    def test_run_step_invalid_step_name(self):
        """Test run_step with invalid step name."""
        binary = BinaryStar()
        binary.star_1.state = 'H-rich_Core_H_burning'
        binary.star_2.state = 'H-rich_Core_H_burning'
        binary.state = 'detached'
        binary.event = 'ZAMS'
        
        total_state = ('H-rich_Core_H_burning', 'H-rich_Core_H_burning', 'detached', 'ZAMS')
        
        # Mock properties with invalid step name that returns None
        mock_properties = Mock()
        mock_properties.flow = {total_state: 'invalid_step_name'}
        mock_properties.pre_step = Mock()
        mock_properties.post_step = Mock()
        
        # Mock getattr to return None for the invalid step name
        def mock_getattr(obj, name, default=None):
            if name == 'invalid_step_name':
                return None
            return getattr(obj, name, default)
        
        binary.properties = mock_properties
        
        with patch('builtins.getattr', side_effect=mock_getattr):
            with pytest.raises(ValueError, match="does not correspond to a function"):
                binary.run_step()
    
    def test_run_step_success(self):
        """Test successful run_step execution."""
        binary = BinaryStar()
        binary.star_1.state = 'H-rich_Core_H_burning'
        binary.star_2.state = 'H-rich_Core_H_burning'
        binary.state = 'detached'
        binary.event = 'ZAMS'
        
        total_state = ('H-rich_Core_H_burning', 'H-rich_Core_H_burning', 'detached', 'ZAMS')
        
        # Mock properties and step function
        mock_step = Mock()
        mock_properties = Mock()
        mock_properties.flow = {total_state: 'test_step'}
        mock_properties.test_step = mock_step
        mock_properties.pre_step = Mock()
        mock_properties.post_step = Mock()
        binary.properties = mock_properties
        
        with patch.object(binary, 'append_state') as mock_append:
            binary.run_step()
        
        mock_properties.pre_step.assert_called_once_with(binary, 'test_step')
        mock_step.assert_called_once_with(binary)
        mock_append.assert_called_once()
        mock_properties.post_step.assert_called_once_with(binary, 'test_step')

class TestBinaryStarHistory:
    """Test BinaryStar history and state management."""
    
    def test_append_state(self):
        """Test appending state to history."""
        binary = BinaryStar(time=0.0, separation=100.0)
        
        # Change current values
        binary.time = 1e6
        binary.separation = 200.0
        
        binary.append_state()
        
        assert len(binary.time_history) == 2
        assert binary.time_history == [0.0, 1e6]
        assert len(binary.separation_history) == 2
        assert binary.separation_history == [100.0, 200.0]
    
    def test_append_state_redirect_event(self):
        """Test append_state with redirect event and history_verbose=False."""
        binary = BinaryStar(time=0.0, event='redirect_step')
        binary.history_verbose = False
        
        initial_length = len(binary.time_history)
        binary.append_state()
        
        # Should not append redirect events when history_verbose=False
        assert len(binary.time_history) == initial_length
    
    def test_append_state_redirect_event_verbose(self):
        """Test append_state with redirect event and history_verbose=True."""
        binary = BinaryStar(time=0.0, event='redirect_step')
        binary.history_verbose = True
        
        initial_length = len(binary.time_history)
        binary.append_state()
        
        # Should append redirect events when history_verbose=True
        assert len(binary.time_history) == initial_length + 1
    
    def test_switch_star(self):
        """Test switching stars."""
        star1 = SingleStar(mass=10.0, state='H-rich_Core_H_burning')
        star2 = SingleStar(mass=8.0, state='stripped_He_Core_He_burning')
        binary = BinaryStar(star_1=star1, star_2=star2)
        
        binary.switch_star()
        
        assert binary.star_1 is star2
        assert binary.star_2 is star1
    
    def test_restore(self):
        """Test restoring binary to previous state."""
        binary = BinaryStar(time=0.0, separation=100.0)
        
        # Add some history
        binary.time = 1e6
        binary.separation = 200.0
        binary.append_state()
        
        binary.time = 2e6
        binary.separation = 300.0
        binary.append_state()
        
        # Restore to index 1
        binary.restore(i=1)
        
        assert binary.time == 1e6
        assert binary.separation == 200.0
        assert len(binary.time_history) == 2
        assert len(binary.separation_history) == 2
    
    def test_restore_with_hooks(self):
        """Test restore with extra hook columns."""
        binary = BinaryStar(time=0.0, separation=100.0)
        
        # Add some history first
        binary.time = 1e6
        binary.separation = 200.0
        binary.append_state()
        
        # Mock properties with hooks
        mock_hook = Mock()
        mock_hook.extra_binary_col_names = ['custom_column']
        binary.properties.all_hooks_classes = [mock_hook]
        
        # Add custom column
        binary.custom_column = [1, 2, 3]
        
        binary.restore(i=1)
        
        assert binary.custom_column == [1, 2]
    
    def test_reset(self):
        """Test resetting binary to initial state."""
        binary = BinaryStar(time=0.0, separation=100.0)
        
        # Add some history
        binary.time = 1e6
        binary.separation = 200.0
        binary.append_state()
        
        binary.reset()
        
        assert binary.time == 0.0
        assert binary.separation == 100.0
        assert len(binary.time_history) == 1
        assert len(binary.separation_history) == 1
    
    def test_reset_with_new_properties(self):
        """Test reset with new SimulationProperties."""
        binary = BinaryStar()
        new_properties = SimulationProperties()
        
        with patch('posydon.binary_evol.binarystar.SimulationProperties') as mock_sim_props:
            mock_sim_props.return_value = new_properties
            binary.reset(properties=new_properties)
        
        # Note: The actual implementation has a bug - it should be SimulationProperties(new_properties)
        # but it passes the properties directly. We test the actual behavior.
        assert binary.properties == new_properties

class TestBinaryStarStates:
    """Test BinaryStar state management."""
    
    def test_update_star_states(self):
        """Test updating star states."""
        binary = BinaryStar()
        binary.star_1.state = 'some_state'  # Non-CO, non-massless state
        binary.star_2.state = 'massless_remnant'  # Should not be updated
        
        # Mock check_state_of_star
        with patch('posydon.binary_evol.binarystar.check_state_of_star') as mock_check:
            mock_check.return_value = 'H-rich_Core_H_burning'
            
            binary.update_star_states()
            
            # Only star 1 should be updated (star 2 is massless remnant)
            mock_check.assert_called_once_with(binary.star_1, star_CO=False)
            assert binary.star_1.state == 'H-rich_Core_H_burning'
    
    def test_update_star_states_massless_remnant(self):
        """Test update_star_states with massless remnant."""
        binary = BinaryStar()
        binary.star_1.state = 'massless_remnant'
        binary.star_2.state = 'massless_remnant'
        
        with patch('posydon.binary_evol.binarystar.check_state_of_star') as mock_check:
            binary.update_star_states()
            
            # Should not call check_state_of_star for massless remnants
            mock_check.assert_not_called()

class TestBinaryStarDataFrames:
    """Test BinaryStar DataFrame conversion methods."""
    
    def test_to_df_basic(self):
        """Test basic to_df functionality."""
        star1 = SingleStar(mass=10.0, state='H-rich_Core_H_burning')
        star2 = SingleStar(mass=8.0, state='H-rich_Core_H_burning')
        binary = BinaryStar(star_1=star1, star_2=star2, index=42, time=0.0)
        
        # Add some history
        binary.time = 1e6
        binary.append_state()
        
        df = binary.to_df()
        
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 2  # Two time steps
        assert 'time' in df.columns
        assert 'S1_mass' in df.columns
        assert 'S2_mass' in df.columns
        assert df.index.name == 'binary_index'
    
    def test_to_df_extra_columns(self):
        """Test to_df with extra columns."""
        binary = BinaryStar(index=42, time=0.0)
        binary.custom_param = [1.0, 2.0]
        
        extra_columns = {'custom_param': 'float64'}
        df = binary.to_df(extra_columns=extra_columns)
        
        assert 'custom_param' in df.columns
        assert len(df) == 2
    
    def test_to_df_ignore_columns(self):
        """Test to_df with ignored columns."""
        binary = BinaryStar(index=42, time=0.0)
        
        df = binary.to_df(ignore_columns=['time_history'])
        
        assert 'time' not in df.columns
    
    def test_to_df_only_select_columns(self):
        """Test to_df with only selected columns."""
        binary = BinaryStar(index=42, time=0.0, separation=100.0)
        
        df = binary.to_df(only_select_columns=['time'])
        
        assert 'time' in df.columns
        assert 'separation' not in df.columns
    
    def test_to_df_v_sys_expansion(self):
        """Test to_df V_sys column expansion."""
        binary = BinaryStar(index=42, V_sys=np.array([1.0, 2.0, 3.0]))
        
        df = binary.to_df()
        
        assert 'V_sys_x' in df.columns
        assert 'V_sys_y' in df.columns
        assert 'V_sys_z' in df.columns
        assert 'V_sys' not in df.columns
        assert df['V_sys_x'].iloc[0] == 1.0
        assert df['V_sys_y'].iloc[0] == 2.0
        assert df['V_sys_z'].iloc[0] == 3.0
    
    def test_to_df_unequal_lengths(self):
        """Test to_df with unequal history lengths (failed binary)."""
        binary = BinaryStar(index=42, time=0.0)
        
        # Simulate unequal history lengths
        binary.time_history = [0.0, 1e6]
        binary.separation_history = [100.0]  # Shorter history
        
        df = binary.to_df()
        
        # Should handle unequal lengths by padding with NaN
        assert len(df) == 2
        assert pd.isna(df['separation'].iloc[1])
    
    def test_from_df(self):
        """Test creating BinaryStar from DataFrame."""
        # Create a test DataFrame
        data = {
            'time': [0.0, 1e6],
            'separation': [100.0, 200.0],
            'S1_mass': [10.0, 9.8],
            'S2_mass': [8.0, 7.9],
            'S1_state': ['H-rich_Core_H_burning', 'H-rich_Core_He_burning'],
            'S2_state': ['H-rich_Core_H_burning', 'H-rich_Core_He_burning']
        }
        df = pd.DataFrame(data)
        df.index = [42, 42]
        df.index.name = 'binary_index'
        
        binary = BinaryStar.from_df(df)
        
        assert binary.index == 42
        assert len(binary.time_history) == 2
        assert binary.time == 1e6
        assert binary.star_1.mass == 9.8
        assert len(binary.star_1.mass_history) == 2
    
    def test_from_df_extra_columns(self):
        """Test from_df with extra columns."""
        data = {
            'time': [0.0, 1e6],
            'custom_param': [1.0, 2.0],
            'S1_mass': [10.0, 9.8],
            'S2_mass': [8.0, 7.9]
        }
        df = pd.DataFrame(data)
        df.index = [42, 42]
        df.index.name = 'binary_index'
        
        extra_columns = {'custom_param': 'float64'}
        binary = BinaryStar.from_df(df, extra_columns=extra_columns)
        
        assert hasattr(binary, 'custom_param')
        assert binary.custom_param == [1.0, 2.0]
    
    def test_to_oneline_df(self):
        """Test converting to oneline DataFrame."""
        binary = BinaryStar(index=42, time=0.0, separation=100.0)
        
        # Add some history
        binary.time = 1e6
        binary.separation = 200.0
        binary.append_state()
        
        df = binary.to_oneline_df()
        
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 1  # Single row
        assert 'time_i' in df.columns  # Initial time
        assert 'time_f' in df.columns  # Final time
        assert df['time_i'].iloc[0] == 0.0
        assert df['time_f'].iloc[0] == 1e6
    
    def test_to_oneline_df_no_history(self):
        """Test to_oneline_df with history=False."""
        binary = BinaryStar(index=42)
        
        df = binary.to_oneline_df(history=False)
        
        # Should only have FAILED/WARNING columns for failed binary check
        assert 'FAILED' in df.columns
        assert 'WARNING' in df.columns
        assert df['FAILED'].iloc[0] == 0
        assert df['WARNING'].iloc[0] == 0
    
    def test_to_oneline_df_with_traceback(self):
        """Test to_oneline_df with failed binary."""
        binary = BinaryStar(index=42)
        binary.traceback = "Some error"  # Simulate failure
        
        df = binary.to_oneline_df()
        
        assert 'FAILED' in df.columns
        assert df['FAILED'].iloc[0] == 1
    
    def test_to_oneline_df_with_warnings(self):
        """Test to_oneline_df with warning."""
        binary = BinaryStar(index=42)
        binary.warnings = "Some warning"  # Simulate warning
        
        df = binary.to_oneline_df()
        
        assert 'WARNING' in df.columns
        assert df['WARNING'].iloc[0] == 1
    
    def test_from_oneline_df(self):
        """Test creating BinaryStar from oneline DataFrame."""
        # Create test oneline DataFrame
        data = {
            'time_i': [0.0],
            'time_f': [1e6],
            'separation_i': [100.0],
            'separation_f': [200.0],
            'S1_mass_i': [10.0],
            'S1_mass_f': [9.8],
            'S2_mass_i': [8.0],
            'S2_mass_f': [7.9],
            'FAILED': [0],
            'WARNING': [0]
        }
        df = pd.DataFrame(data)
        df.index = [42]
        df.index.name = 'binary_index'
        
        binary = BinaryStar.from_oneline_df(df)
        
        assert binary.index == 42
        assert len(binary.time_history) == 1
        assert binary.time == 0.0  # Takes initial value
        assert binary.star_1.mass == 10.0
    
    def test_from_oneline_df_failed(self):
        """Test from_oneline_df with failed binary."""
        data = {
            'time_i': [0.0],
            'S1_mass_i': [10.0],
            'S2_mass_i': [8.0],
            'FAILED': [1],
            'WARNING': [0]
        }
        df = pd.DataFrame(data)
        df.index = [42]
        df.index.name = 'binary_index'
        
        binary = BinaryStar.from_oneline_df(df)
        
        assert binary.event == 'FAILED'

# TestBinaryStarFromRun class removed - from_run method is marked with pragma no cover

class TestBinaryStarUtilities:
    """Test BinaryStar utility methods."""
    
    def test_repr(self):
        """Test __repr__ method."""
        binary = BinaryStar(index=42, state='detached', event='ZAMS')
        
        repr_str = repr(binary)
        
        assert 'BinaryStar' in repr_str
        assert 'state: detached' in repr_str
        assert 'event: ZAMS' in repr_str
    
    def test_repr_with_error(self):
        """Test __repr__ with error message."""
        binary = BinaryStar(index=42)
        binary.error_message = "Test error"
        
        repr_str = repr(binary)
        
        assert 'BINARY FAILED: Test error' in repr_str
    
    def test_repr_with_warning(self):
        """Test __repr__ with warning message."""
        binary = BinaryStar(index=42)
        binary.warning_message = "Test warning"
        
        repr_str = repr(binary)
        
        assert 'WARNING FOUND: Test warning' in repr_str
    
    def test_str(self):
        """Test __str__ method."""
        star1 = SingleStar(mass=10.0, state='H-rich_Core_H_burning')
        star2 = SingleStar(mass=8.0, state='stripped_He_Core_He_burning')
        binary = BinaryStar(star_1=star1, star_2=star2, 
                           state='detached', event='ZAMS', orbital_period=10.0)
        
        str_repr = str(binary)
        
        assert 'BinaryStar(' in str_repr
        assert 'detached' in str_repr
        assert 'ZAMS' in str_repr
        assert 'p=10.00' in str_repr
        assert 'M=10.00' in str_repr
        assert 'M=8.00' in str_repr
    
    def test_str_with_nan_values(self):
        """Test __str__ with non-numeric values."""
        binary = BinaryStar(state='detached', event='ZAMS')
        binary.orbital_period = None
        binary.star_1.mass = "not_a_number"
        binary.star_2.mass = [1, 2, 3]  # List instead of number
        
        str_repr = str(binary)
        
        assert 'p=nan' in str_repr
        assert 'M=nan' in str_repr
    
    def test_initial_condition_message(self):
        """Test initial_condition_message method."""
        star1 = SingleStar(mass=10.0, state='H-rich_Core_H_burning')
        star1.natal_kick_array = [1, 2, 3, 4]
        star2 = SingleStar(mass=8.0, state='stripped_He_Core_He_burning')
        star2.natal_kick_array = [5, 6, 7, 8]
        
        binary = BinaryStar(star_1=star1, star_2=star2, 
                           orbital_period=10.0, eccentricity=0.1,
                           state='detached', event='ZAMS')
        
        message = binary.initial_condition_message()
        
        assert 'Failed Binary Initial Conditions' in message
        assert 'S1 mass: 10.0' in message
        assert 'S2 mass: 8.0' in message
        assert 'orbital period: 10.0' in message
        assert 'eccentricity: 0.1' in message
    
    def test_initial_condition_message_custom(self):
        """Test initial_condition_message with custom parameters."""
        binary = BinaryStar()
        custom_params = ["Custom message 1\n", "Custom message 2\n"]
        
        message = binary.initial_condition_message(ini_params=custom_params)
        
        assert message == "Custom message 1\nCustom message 2\n"

class TestModuleConstants:
    """Test module-level constants and functions."""
    
    def test_binaryproperties_list(self):
        """Test BINARYPROPERTIES constant."""
        assert isinstance(BINARYPROPERTIES, list)
        assert len(BINARYPROPERTIES) > 0
        
        # Check for essential properties
        essential_props = ['state', 'event', 'time', 'separation', 'orbital_period']
        for prop in essential_props:
            assert prop in BINARYPROPERTIES
    
    def test_maximum_step_time(self):
        """Test MAXIMUM_STEP_TIME constant."""
        assert isinstance(MAXIMUM_STEP_TIME, int)
        assert MAXIMUM_STEP_TIME > 0
    
    def test_signal_handler(self):
        """Test signal_handler function."""
        with pytest.raises(RuntimeError, match="Binary Step Exceeded Allotted Time"):
            totest.signal_handler(signal.SIGALRM, None)
    
    def test_star_attributes_dictionaries(self):
        """Test star attribute mapping dictionaries."""
        # Test STAR_ATTRIBUTES_FROM_BINARY_HISTORY
        assert isinstance(totest.STAR_ATTRIBUTES_FROM_BINARY_HISTORY, dict)
        assert 'mass' in totest.STAR_ATTRIBUTES_FROM_BINARY_HISTORY
        assert len(totest.STAR_ATTRIBUTES_FROM_BINARY_HISTORY['mass']) == 2
        
        # Test STAR_ATTRIBUTES_FROM_STAR_HISTORY
        assert isinstance(totest.STAR_ATTRIBUTES_FROM_STAR_HISTORY, dict)
        assert 'state' in totest.STAR_ATTRIBUTES_FROM_STAR_HISTORY
        
        # Test BINARY_ATTRIBUTES_FROM_HISTORY
        assert isinstance(totest.BINARY_ATTRIBUTES_FROM_HISTORY, dict)
        assert 'time' in totest.BINARY_ATTRIBUTES_FROM_HISTORY

class TestSignalHandling:
    """Test signal handling for step timeout."""
    
    def test_signal_alarm_setup(self):
        """Test that signal alarm is properly set up."""
        # The signal should be configured in the module
        current_handler = signal.signal(signal.SIGALRM, signal.SIG_DFL)
        signal.signal(signal.SIGALRM, current_handler)  # Restore
        
        # Verify our handler is set
        assert current_handler == totest.signal_handler

class TestIntegrationScenarios:
    """Test integrated scenarios and edge cases."""
    
    def test_full_binary_lifecycle(self):
        """Test a complete binary evolution scenario."""
        # Create initial binary
        star1 = SingleStar(mass=10.0, state='H-rich_Core_H_burning')
        star2 = SingleStar(mass=8.0, state='H-rich_Core_H_burning')
        binary = BinaryStar(star_1=star1, star_2=star2, 
                           time=0.0, separation=100.0, 
                           state='detached', event='ZAMS', index=1)
        
        # Simulate evolution steps
        binary.time = 1e6
        binary.state = 'RLO1'
        binary.event = 'oRLO1'
        binary.append_state()
        
        binary.time = 2e6
        binary.state = 'merged'
        binary.event = 'END'
        binary.append_state()
        
        # Test conversion to DataFrame and back
        df = binary.to_df()
        binary_restored = BinaryStar.from_df(df)
        
        assert binary_restored.index == binary.index
        assert len(binary_restored.time_history) == len(binary.time_history)
        assert binary_restored.time == binary.time
        
        # Test oneline conversion
        oneline_df = binary.to_oneline_df()
        binary_from_oneline = BinaryStar.from_oneline_df(oneline_df)
        
        assert binary_from_oneline.index == binary.index
    
    def test_error_handling_scenarios(self):
        """Test various error scenarios."""
        binary = BinaryStar()
        
        # Test DataFrame conversion with missing attributes
        with pytest.raises(AttributeError):
            binary.to_df(extra_columns={'nonexistent': 'float64'})
    
    def test_edge_case_values(self):
        """Test edge case values in binary properties."""
        # Test with None values
        binary = BinaryStar(time=None, separation=None)
        df = binary.to_df()
        assert pd.isna(df['time'].iloc[0])
        assert pd.isna(df['separation'].iloc[0])
        
        # Test with very large/small values
        binary = BinaryStar(time=1e50, separation=1e-10)
        df = binary.to_df()
        assert df['time'].iloc[0] == 1e50
        assert df['separation'].iloc[0] == 1e-10
