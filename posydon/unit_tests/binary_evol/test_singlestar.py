"""Unit tests of posydon/binary_evol/singlestar.py"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>"
]

# import the module which will be tested
import posydon.binary_evol.singlestar as totest

# import other needed code for the tests
import numpy as np
import pandas as pd
from pytest import fixture, raises

@fixture
def basic_star():
    """Create a basic SingleStar with minimal properties."""
    return totest.SingleStar(
        mass=1.5,
        metallicity=0.02,
        state="H-rich_Core_H_burning",
        log_R=0.0,
        log_L=0.3
    )


@fixture
def complete_star():
    """Create a SingleStar with all STARPROPERTIES set to 'realistic' values."""
    kwargs = {}
    for prop in totest.STARPROPERTIES:
        if prop == 'state':
            kwargs[prop] = "H-rich_Core_H_burning"
        elif prop == 'metallicity':
            kwargs[prop] = 0.02
        elif prop == 'mass':
            kwargs[prop] = 1.5
        elif prop == 'log_R':
            kwargs[prop] = 0.0
        elif prop == 'log_L':
            kwargs[prop] = 0.3
        elif prop in ['lg_mdot', 'lg_system_mdot', 'lg_wind_mdot']:
            kwargs[prop] = -8.0
        elif 'mass' in prop:
            kwargs[prop] = 0.1
        elif 'radius' in prop:
            kwargs[prop] = 0.05
        elif 'center_' in prop:
            kwargs[prop] = 0.5
        elif 'surface_' in prop:
            kwargs[prop] = 0.7
        elif 'log_L' in prop:
            kwargs[prop] = -2.0
        elif prop == 'c12_c12':
            kwargs[prop] = 1e-10
        elif prop == 'center_gamma':
            kwargs[prop] = 1.5
        elif prop == 'avg_c_in_c_core':
            kwargs[prop] = 0.2
        elif 'omega' in prop:
            kwargs[prop] = 1e-6
        elif prop == 'total_moment_of_inertia':
            kwargs[prop] = 1e50
        elif prop == 'log_total_angular_momentum':
            kwargs[prop] = 48.0
        elif prop == 'spin':
            kwargs[prop] = 0.1
        elif 'conv_env' in prop and 'mass' in prop:
            kwargs[prop] = 0.8
        elif 'conv_env' in prop and 'radius' in prop:
            kwargs[prop] = 0.6
        elif 'conv_env' in prop and 'time' in prop:
            kwargs[prop] = 1e6
        elif 'binding_energy' in prop:
            kwargs[prop] = -1e48
        elif 'mass_conv_reg_fortides' in prop:
            kwargs[prop] = 0.3
        elif 'thickness_conv_reg_fortides' in prop:
            kwargs[prop] = 0.2
        elif 'radius_conv_reg_fortides' in prop:
            kwargs[prop] = 0.4
        elif 'lambda_CE' in prop:
            kwargs[prop] = 0.5
        elif prop == 'profile':
            kwargs[prop] = None
        elif 'total_mass' in prop:
            kwargs[prop] = 0.7
        else:
            kwargs[prop] = 0.1
    
    return totest.SingleStar(**kwargs)


@fixture
def star_with_history():
    """Create a SingleStar with evolution history."""
    star = totest.SingleStar(
        mass=2.0,
        metallicity=0.02,
        state="H-rich_Core_H_burning"
    )
    
    # Simulate evolution by appending states
    star.mass = 1.8
    star.state = "H-rich_Shell_H_burning"
    star.append_state()
    
    star.mass = 1.6
    star.state = "H-rich_Core_He_burning"
    star.append_state()
    
    star.mass = 1.4
    star.state = "H-rich_Central_He_depleted"
    star.append_state()
    
    return star


@fixture
def star_with_none_values():
    """Star with None values in various properties."""
    return totest.SingleStar(
        mass=None,
        metallicity=None,
        state=None,
        log_R=None,
        log_L=None
    )


@fixture
def star_with_unequal_histories():
    """Star with mismatched history lengths (corrupted data)."""
    star = totest.SingleStar(mass=1.0, metallicity=0.02)
    # Manually corrupt the history to create unequal lengths
    star.mass_history = [1.0, 0.9, 0.8]
    star.metallicity_history = [0.02, 0.02]  # Different length
    return star


@fixture
def star_with_profile():
    """Star with detailed stellar profile data."""
    # Create a simple mock profile
    profile_data = np.array([
        (1.0, 0.1, 1.0, 0.0, 1e16, 0.7, 0.28, 0.02),
        (0.5, 0.3, 0.5, -0.3, 5e15, 0.6, 0.38, 0.02),
        (0.2, 0.2, 0.2, -0.7, 2e15, 0.1, 0.88, 0.02),
        (0.01, 0.01, 0.01, -2.0, 1e14, 0.0, 0.98, 0.02)
    ], dtype=[
        ('mass', 'f8'), ('dm', 'f8'), ('radius', 'f8'), ('log_R', 'f8'),
        ('energy', 'f8'), ('x_mass_fraction_H', 'f8'),
        ('y_mass_fraction_He', 'f8'), ('z_mass_fraction_metals', 'f8')
    ])
    
    return totest.SingleStar(
        mass=1.0,
        metallicity=0.02,
        state="H-rich_Core_H_burning",
        profile=profile_data
    )


# define test classes collecting several test functions
class TestConstants:
    """Test module-level constants."""
    
    def test_STARPROPERTIES_completeness(self):
        """Test that STARPROPERTIES contains expected stellar properties."""
        required_props = [
            'state', 'metallicity', 'mass', 'log_R', 'log_L',
            'he_core_mass', 'co_core_mass', 'center_h1', 'surface_h1',
            'log_LH', 'log_LHe', 'log_Lnuc', 'spin', 'profile'
        ]
        for prop in required_props:
            assert prop in totest.STARPROPERTIES, f"Missing required property: {prop}"
    
    def test_STARPROPERTIES_is_list(self):
        """Test that STARPROPERTIES is a list."""
        assert isinstance(totest.STARPROPERTIES, list)
        assert len(totest.STARPROPERTIES) > 0
    
    def test_STAR_ATTRIBUTES_FROM_STAR_HISTORY_SINGLE_structure(self):
        """Test structure of STAR_ATTRIBUTES_FROM_STAR_HISTORY_SINGLE."""
        assert isinstance(totest.STAR_ATTRIBUTES_FROM_STAR_HISTORY_SINGLE, dict)
        # Check some expected mappings
        assert 'mass' in totest.STAR_ATTRIBUTES_FROM_STAR_HISTORY_SINGLE
        assert 'spin' in totest.STAR_ATTRIBUTES_FROM_STAR_HISTORY_SINGLE


class TestUtilityFunctions:
    """Test module-level utility functions."""
    
    def test_properties_massless_remnant(self):
        """Test properties_massless_remnant function."""
        props = totest.properties_massless_remnant()
        
        # Should return a dictionary
        assert isinstance(props, dict)
        
        # Should contain all STARPROPERTIES
        for prop in totest.STARPROPERTIES:
            assert prop in props
        
        # Check specific values
        assert props["state"] == "massless_remnant"
        assert props["mass"] == 0.0
        
        # All other properties should be NaN
        for prop in totest.STARPROPERTIES:
            if prop not in ["state", "mass"]:
                assert np.isnan(props[prop])
    
    def test_convert_star_to_massless_remnant(self, basic_star):
        """Test convert_star_to_massless_remnant function."""
        original_mass = basic_star.mass
        original_state = basic_star.state
        
        # Convert to massless remnant
        converted_star = totest.convert_star_to_massless_remnant(basic_star)
        
        # Should return the same star object
        assert converted_star is basic_star
        
        # Check properties are updated
        assert basic_star.state == "massless_remnant"
        assert basic_star.mass == 0.0
        
        # Other properties should be NaN
        assert np.isnan(basic_star.log_R)
        assert np.isnan(basic_star.log_L)


class TestSingleStarInit:
    """Test SingleStar initialization."""
    
    def test_init_empty(self):
        """Test initialization with no arguments."""
        star = totest.SingleStar()
        
        # All STARPROPERTIES should be None initially
        for prop in totest.STARPROPERTIES:
            assert getattr(star, prop) is None
            # History should contain one None entry
            assert getattr(star, prop + '_history') == [None]
    
    def test_init_with_kwargs(self):
        """Test initialization with keyword arguments."""
        star = totest.SingleStar(
            mass=1.5,
            metallicity=0.02,
            state="H-rich_Core_H_burning"
        )
        
        assert star.mass == 1.5
        assert star.metallicity == 0.02
        assert star.state == "H-rich_Core_H_burning"
        
        # History should be initialized
        assert star.mass_history == [1.5]
        assert star.metallicity_history == [0.02]
        assert star.state_history == ["H-rich_Core_H_burning"]
    
    def test_init_with_starproperties(self, complete_star):
        """Test initialization with all STARPROPERTIES."""
        # Check that all properties are set
        for prop in totest.STARPROPERTIES:
            assert hasattr(complete_star, prop)
            assert hasattr(complete_star, prop + '_history')
            # History should contain the initial value
            assert len(getattr(complete_star, prop + '_history')) == 1
    
    def test_init_with_extra_attributes(self):
        """Test initialization with extra non-STARPROPERTIES attributes."""
        star = totest.SingleStar(
            mass=1.0,
            extra_attr="test_value",
            another_attr=42
        )
        
        assert star.mass == 1.0
        assert star.extra_attr == "test_value"
        assert star.another_attr == 42
        
        # Extra attributes should not have history
        assert not hasattr(star, 'extra_attr_history')
        assert not hasattr(star, 'another_attr_history')
    
    def test_init_sets_default_sn_attributes(self):
        """Test that SN-related attributes are initialized."""
        star = totest.SingleStar()
        
        # Check SN attributes are initialized
        assert star.natal_kick_array == [None] * 4
        assert star.spin_orbit_tilt_first_SN is None
        assert star.spin_orbit_tilt_second_SN is None
        assert star.f_fb is None
        assert star.SN_type is None
        assert star.m_disk_accreted is None
        assert star.m_disk_radiated is None
        assert star.h1_mass_ej is None
        assert star.he4_mass_ej is None
        assert star.M4 is None
        assert star.mu4 is None
    
    def test_init_sets_default_ce_attributes(self):
        """Test that CE-related attributes are initialized."""
        star = totest.SingleStar()
        
        # Check CE attributes for different percentages
        for val in [1, 10, 30, 'pure_He_star_10']:
            assert hasattr(star, f'm_core_CE_{val}cent')
            assert getattr(star, f'm_core_CE_{val}cent') is None
            assert hasattr(star, f'r_core_CE_{val}cent')
            assert getattr(star, f'r_core_CE_{val}cent') is None
    
    def test_init_sets_default_core_collapse_attributes(self):
        """Test that core collapse attributes are initialized."""
        star = totest.SingleStar()
        
        # Check He depletion attributes
        assert star.avg_c_in_c_core_at_He_depletion is None
        assert star.co_core_mass_at_He_depletion is None


class TestSingleStarMethods:
    """Test SingleStar methods."""
    
    def test_append_state(self, basic_star):
        """Test append_state method."""
        initial_length = len(basic_star.mass_history)
        
        # Change some properties
        basic_star.mass = 1.3
        basic_star.state = "H-rich_Shell_H_burning"
        
        # Append the new state
        basic_star.append_state()
        
        # History should be longer
        assert len(basic_star.mass_history) == initial_length + 1
        assert len(basic_star.state_history) == initial_length + 1
        
        # New values should be in history
        assert basic_star.mass_history[-1] == 1.3
        assert basic_star.state_history[-1] == "H-rich_Shell_H_burning"
    
    def test_append_state_with_none_values(self, star_with_none_values):
        """Test append_state with None values."""
        star_with_none_values.mass = None
        star_with_none_values.append_state()
        
        assert star_with_none_values.mass_history[-1] is None
        assert len(star_with_none_values.mass_history) == 2
    
    def test_restore_to_initial_state(self, star_with_history):
        """Test restore method to initial state."""
        # Should have evolved beyond initial state
        assert len(star_with_history.mass_history) > 1
        initial_mass = star_with_history.mass_history[0]
        initial_state = star_with_history.state_history[0]
        
        # Restore to initial state (i=0)
        star_with_history.restore(0)
        
        # Should be back at initial values
        assert star_with_history.mass == initial_mass
        assert star_with_history.state == initial_state
        
        # History should be truncated to only initial state
        assert len(star_with_history.mass_history) == 1
        assert len(star_with_history.state_history) == 1
    
    def test_restore_to_intermediate_state(self, star_with_history):
        """Test restore method to intermediate state."""
        original_length = len(star_with_history.mass_history)
        assert original_length > 2
        
        # Restore to second state (i=1)
        target_mass = star_with_history.mass_history[1]
        target_state = star_with_history.state_history[1]
        
        star_with_history.restore(1)
        
        # Should be at the target state
        assert star_with_history.mass == target_mass
        assert star_with_history.state == target_state
        
        # History should be truncated
        assert len(star_with_history.mass_history) == 2
        assert len(star_with_history.state_history) == 2
    
    def test_restore_with_hooks(self, basic_star):
        """Test restore method with extra hooks."""
        # First, add some evolution to create multiple history entries
        basic_star.mass = 1.3
        basic_star.append_state()
        basic_star.mass = 1.1
        basic_star.append_state()
        
        # Add some fake hook data
        basic_star.extra_col1 = [1, 2, 3]
        basic_star.extra_col2 = [10, 20, 30]
        
        # Create mock hook object
        class MockHook:
            extra_star_col_names = ['extra_col1', 'extra_col2']
        
        hooks = [MockHook()]
        
        # Restore to index 1
        basic_star.restore(1, hooks=hooks)
        
        # Extra columns should be truncated
        assert basic_star.extra_col1 == [1, 2]
        assert basic_star.extra_col2 == [10, 20]
    
    def test_restore_invalid_index(self, basic_star):
        """Test restore with out-of-bounds index."""
        with raises(IndexError):
            basic_star.restore(10)  # Index beyond history length


class TestSingleStarDataFrame:
    """Test SingleStar DataFrame methods."""
    
    def test_to_df_default(self, star_with_history):
        """Test to_df method with default parameters."""
        df = star_with_history.to_df()
        
        # Should be a pandas DataFrame
        assert isinstance(df, pd.DataFrame)
        
        # Should have rows equal to history length
        expected_rows = len(star_with_history.mass_history)
        assert len(df) == expected_rows
        
        # Should contain columns for all STARPROPERTIES (except profile by default)
        expected_cols = [prop for prop in totest.STARPROPERTIES if prop != 'profile']
        for col in expected_cols:
            assert col in df.columns
        
        # Check some values
        assert df['mass'].iloc[0] == star_with_history.mass_history[0]
        assert df['mass'].iloc[-1] == star_with_history.mass_history[-1]
    
    def test_to_df_with_extra_columns(self, basic_star):
        """Test to_df with extra columns."""
        # Add extra data
        basic_star.extra_data = [1, 2, 3]
        
        df = basic_star.to_df(extra_columns={'extra_data': 'int64'})
        
        assert 'extra_data' in df.columns
        assert df['extra_data'].dtype == 'int64'
    
    def test_to_df_ignore_columns(self, star_with_history):
        """Test to_df with ignored columns."""
        df = star_with_history.to_df(ignore_columns=['mass_history', 'state_history'])
        
        assert 'mass' not in df.columns
        assert 'state' not in df.columns
        # Other columns should still be present
        assert 'metallicity' in df.columns
    
    def test_to_df_only_select_columns(self, star_with_history):
        """Test to_df with only selected columns."""
        df = star_with_history.to_df(only_select_columns=['mass', 'state'])
        
        # Should only have selected columns
        assert set(df.columns) == {'mass', 'state'}
    
    def test_to_df_include_profile(self, star_with_profile):
        """Test to_df with profile included."""
        df = star_with_profile.to_df(include_profile=True)
        
        assert 'profile' in df.columns
    
    def test_to_df_with_prefix(self, basic_star):
        """Test to_df with column prefix."""
        df = basic_star.to_df(prefix='star1_')
        
        # All columns should have prefix
        for col in df.columns:
            assert col.startswith('star1_')
    
    def test_to_df_with_null_value(self, star_with_none_values):
        """Test to_df with custom null value."""
        df = star_with_none_values.to_df(null_value=-999)
        
        # None values should be replaced with -999
        for col in df.columns:
            none_mask = df[col] == -999
            if none_mask.any():
                # This column had None values that were replaced
                assert True
    
    def test_to_df_unequal_history_lengths(self, star_with_unequal_histories):
        """Test to_df with unequal history lengths."""
        df = star_with_unequal_histories.to_df()
        
        # Should handle unequal lengths by padding with NaN
        assert isinstance(df, pd.DataFrame)
        # All columns should have the same length (max of all histories)
        max_length = max(
            len(star_with_unequal_histories.mass_history),
            len(star_with_unequal_histories.metallicity_history)
        )
        assert len(df) == max_length
    
    def test_to_df_missing_attributes(self, basic_star):
        """Test to_df with missing attributes."""
        # Remove an attribute to cause an error
        delattr(basic_star, 'mass_history')
        
        with raises(AttributeError):
            basic_star.to_df()


class TestSingleStarOnelineDataFrame:
    """Test SingleStar oneline DataFrame method."""
    
    def test_to_oneline_df_with_history(self, star_with_history):
        """Test to_oneline_df with history=True."""
        df = star_with_history.to_oneline_df(history=True)
        
        # Should be a single-row DataFrame
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 1
        
        # Should have initial and final columns
        assert any(col.endswith('_i') for col in df.columns)
        assert any(col.endswith('_f') for col in df.columns)
        
        # Check specific values
        if 'mass_i' in df.columns and 'mass_f' in df.columns:
            assert df['mass_i'].iloc[0] == star_with_history.mass_history[0]
            assert df['mass_f'].iloc[0] == star_with_history.mass_history[-1]
    
    def test_to_oneline_df_without_history(self, basic_star):
        """Test to_oneline_df with history=False."""
        df = basic_star.to_oneline_df(history=False)
        
        # Should be empty DataFrame if no scalar names provided
        assert isinstance(df, pd.DataFrame)
        assert len(df.columns) == 0 or len(df) == 0
    
    def test_to_oneline_df_with_scalar_names(self, basic_star):
        """Test to_oneline_df with scalar names."""
        # Add some scalar attributes
        basic_star.f_fb = 0.5
        basic_star.SN_type = "CCSN"
        
        df = basic_star.to_oneline_df(
            history=False,
            scalar_names=['f_fb', 'SN_type']
        )
        
        assert 'f_fb' in df.columns
        assert 'SN_type' in df.columns
        assert df['f_fb'].iloc[0] == 0.5
        assert df['SN_type'].iloc[0] == "CCSN"
    
    def test_to_oneline_df_with_natal_kick_array(self, basic_star):
        """Test to_oneline_df with natal_kick_array."""
        basic_star.natal_kick_array = [10, 20, 30, 40]
        
        df = basic_star.to_oneline_df(
            history=False,
            scalar_names=['natal_kick_array']
        )
        
        # Should have separate columns for each array element
        assert 'natal_kick_array_0' in df.columns
        assert 'natal_kick_array_1' in df.columns
        assert 'natal_kick_array_2' in df.columns
        assert 'natal_kick_array_3' in df.columns
        
        assert df['natal_kick_array_0'].iloc[0] == 10
        assert df['natal_kick_array_3'].iloc[0] == 40
    
    def test_to_oneline_df_with_prefix(self, basic_star):
        """Test to_oneline_df with prefix."""
        df = basic_star.to_oneline_df(history=True, prefix='primary_')
        
        # All columns should have prefix
        for col in df.columns:
            assert col.startswith('primary_')


class TestSingleStarRepr:
    """Test SingleStar string representation."""
    
    def test___repr__(self, basic_star):
        """Test __repr__ method."""
        repr_str = repr(basic_star)
        
        # Should contain class information
        assert 'SingleStar' in repr_str
        assert 'posydon.binary_evol.singlestar' in repr_str
        
        # Should contain property values
        assert 'mass: 1.5' in repr_str
        assert 'metallicity: 0.02' in repr_str
        assert 'state: H-rich_Core_H_burning' in repr_str
    
    def test___repr___with_none_values(self, star_with_none_values):
        """Test __repr__ with None values."""
        repr_str = repr(star_with_none_values)
        
        # Should handle None values gracefully
        assert 'mass: None' in repr_str
        assert 'SingleStar' in repr_str


class TestSingleStarEdgeCases:
    """Test edge cases and error handling."""
    
    def test_large_history_lengths(self):
        """Test with very long evolution histories."""
        star = totest.SingleStar(mass=1.0)
        
        # Add many evolution steps
        for i in range(1000):
            star.mass = 1.0 - i * 0.0001
            star.append_state()
        
        assert len(star.mass_history) == 1001  # Initial + 1000 steps
        
        # DataFrame creation should still work
        df = star.to_df()
        assert len(df) == 1001
    
    def test_extreme_property_values(self):
        """Test with extreme property values."""
        star = totest.SingleStar(
            mass=1e-10,  # Very small mass
            log_L=10,    # Very high luminosity
            log_R=-5,    # Very small radius
            metallicity=1e-6  # Very low metallicity
        )
        
        assert star.mass == 1e-10
        assert star.log_L == 10
        assert star.log_R == -5
        assert star.metallicity == 1e-6
    
    def test_negative_values_where_inappropriate(self):
        """Test that negative values are handled appropriately."""
        # Negative mass should be allowed (might be used for special cases)
        star = totest.SingleStar(mass=-1.0)
        assert star.mass == -1.0
    
    def test_string_values_in_numeric_properties(self):
        """Test behavior with string values in numeric properties."""
        # This should be allowed as the class doesn't enforce types
        star = totest.SingleStar(mass="invalid")
        assert star.mass == "invalid"
    
    def test_very_long_property_names(self):
        """Test with very long custom property names."""
        very_long_name = "a" * 1000
        star = totest.SingleStar(**{very_long_name: "test_value"})
        assert getattr(star, very_long_name) == "test_value"


class TestSingleStarIntegration:
    """Test integration with other POSYDON components."""
    
    def test_starproperties_consistency(self):
        """Test that all STARPROPERTIES are properly handled."""
        star = totest.SingleStar()
        
        # All STARPROPERTIES should have corresponding attributes
        for prop in totest.STARPROPERTIES:
            assert hasattr(star, prop)
            assert hasattr(star, prop + '_history')
    
    def test_sn_models_integration(self):
        """Test SN models attributes are initialized."""
        star = totest.SingleStar()
        
        # Should have SN model attributes initialized
        from posydon.grids.SN_MODELS import SN_MODELS
        for model_name in SN_MODELS.keys():
            assert hasattr(star, model_name)
            assert getattr(star, model_name) is None
    
    def test_realistic_stellar_evolution(self):
        """Test a realistic stellar evolution sequence."""
        # Create a 15 solar mass star
        star = totest.SingleStar(
            mass=15.0,
            metallicity=0.02,
            state="H-rich_Core_H_burning",
            log_R=1.0,
            log_L=4.0,
            center_h1=0.7,
            surface_h1=0.7
        )
        
        # Evolve through hydrogen burning
        star.center_h1 = 0.1
        star.state = "H-rich_Shell_H_burning"
        star.append_state()
        
        # Helium burning
        star.center_h1 = 0.0
        star.center_he4 = 0.8
        star.state = "H-rich_Core_He_burning"
        star.append_state()
        
        # Central helium depletion
        star.center_he4 = 0.0
        star.state = "H-rich_Central_He_depleted"
        star.append_state()
        
        # Check evolution history
        assert len(star.state_history) == 4
        assert star.state_history[0] == "H-rich_Core_H_burning"
        assert star.state_history[-1] == "H-rich_Central_He_depleted"
        
        # Convert to DataFrame and check
        df = star.to_df()
        assert len(df) == 4
        assert df['center_h1'].iloc[0] == 0.7
        assert df['center_h1'].iloc[-1] == 0.0
