"""Unit tests of posydon/binary_evol/flow_chart.py
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>"
]

# import the module which will be tested
import posydon.binary_evol.flow_chart as totest

# import other needed code for the tests
import numpy as np
from pytest import fixture, raises, warns, approx
from inspect import isroutine, isfunction

# define test classes collecting several test functions
class TestConstants:
    """Test module-level constants and state lists."""
    
    def test_STAR_STATES_ALL_completeness(self):
        """Test that STAR_STATES_ALL contains expected stellar states."""
        required_states = [
            'WD', 'NS', 'BH', 'massless_remnant',
            'H-rich_Core_H_burning', 'H-rich_Core_He_burning',
            'stripped_He_Core_He_burning', 'stripped_He_Central_He_depleted'
        ]
        for state in required_states:
            assert state in totest.STAR_STATES_ALL, f"Missing required state: {state}"
    
    def test_STAR_STATES_ALL_is_list(self):
        """Test that STAR_STATES_ALL is a list."""
        assert isinstance(totest.STAR_STATES_ALL, list)
        assert len(totest.STAR_STATES_ALL) > 0
    
    def test_STAR_STATES_CO_subset(self):
        """Test that STAR_STATES_CO contains only compact objects."""
        expected_co_states = ['BH', 'NS', 'WD']
        assert totest.STAR_STATES_CO == expected_co_states
        
        # All CO states should be in STAR_STATES_ALL
        for state in totest.STAR_STATES_CO:
            assert state in totest.STAR_STATES_ALL
    
    def test_STAR_STATES_NOT_NORMALSTAR_consistency(self):
        """Test STAR_STATES_NOT_NORMALSTAR includes CO states and massless remnant."""
        assert 'massless_remnant' in totest.STAR_STATES_NOT_NORMALSTAR
        for co_state in totest.STAR_STATES_CO:
            assert co_state in totest.STAR_STATES_NOT_NORMALSTAR
    
    def test_STAR_STATES_NOT_CO_consistency(self):
        """Test STAR_STATES_NOT_CO excludes compact objects."""
        for co_state in totest.STAR_STATES_CO:
            assert co_state not in totest.STAR_STATES_NOT_CO
        
        # Should contain normal stellar states
        assert 'H-rich_Core_H_burning' in totest.STAR_STATES_NOT_CO
        assert 'massless_remnant' in totest.STAR_STATES_NOT_CO
    
    def test_STAR_STATES_NORMALSTAR_consistency(self):
        """Test STAR_STATES_NORMALSTAR excludes compact objects and massless remnant."""
        for excluded_state in totest.STAR_STATES_NOT_NORMALSTAR:
            assert excluded_state not in totest.STAR_STATES_NORMALSTAR
        
        # Should contain normal stellar evolution states
        assert 'H-rich_Core_H_burning' in totest.STAR_STATES_NORMALSTAR
        assert 'stripped_He_Core_He_burning' in totest.STAR_STATES_NORMALSTAR
    
    def test_STAR_STATES_H_RICH_consistency(self):
        """Test STAR_STATES_H_RICH contains only hydrogen-rich states."""
        for state in totest.STAR_STATES_H_RICH:
            assert 'H-rich' in state or 'accreted_He' in state
            assert 'stripped_He' not in state
        
        # Should not contain He-rich stripped states
        assert 'stripped_He_Core_He_burning' not in totest.STAR_STATES_H_RICH
    
    def test_STAR_STATES_HE_RICH_consistency(self):
        """Test STAR_STATES_HE_RICH contains only helium-rich states."""
        for state in totest.STAR_STATES_HE_RICH:
            assert ('stripped_He' in state or 
                    'accreted_He_Core_He_burning' in state or
                    state in ['accreted_He_non_burning', 'H-rich_non_burning'])  # Special cases
        
        # Should not contain most H-rich states
        assert 'H-rich_Core_H_burning' not in totest.STAR_STATES_HE_RICH
        assert 'accreted_He_Core_He_burning' in totest.STAR_STATES_HE_RICH
    
    def test_STAR_STATES_C_DEPLETION_consistency(self):
        """Test STAR_STATES_C_DEPLETION contains only carbon depletion states."""
        for state in totest.STAR_STATES_C_DEPLETION:
            assert 'C_depletion' in state
        
        expected_c_depletion = [
            'H-rich_Central_C_depletion',
            'stripped_He_Central_C_depletion'
        ]
        for state in expected_c_depletion:
            assert state in totest.STAR_STATES_C_DEPLETION
    
    def test_STAR_STATES_EVOLVABLE_consistency(self):
        """Test evolvable states exclude carbon depletion states."""
        # H-rich evolvable should not contain C depletion states
        for state in totest.STAR_STATES_H_RICH_EVOLVABLE:
            assert 'C_depletion' not in state
            assert state in totest.STAR_STATES_H_RICH
        
        # He-rich evolvable should not contain C depletion states (except special case)
        for state in totest.STAR_STATES_HE_RICH_EVOLVABLE:
            if state != 'H-rich_non_burning':  # Special case
                assert 'C_depletion' not in state
        
        # Should include H-rich_non_burning as special case
        assert 'H-rich_non_burning' in totest.STAR_STATES_HE_RICH_EVOLVABLE
    
    def test_STAR_STATES_CC_includes_core_collapse_states(self):
        """Test STAR_STATES_CC contains core collapse states."""
        expected_cc_states = [
            'H-rich_Central_C_depletion',
            'H-rich_Central_He_depleted',
            'stripped_He_Central_He_depleted',
            'stripped_He_Central_C_depletion'
        ]
        for state in expected_cc_states:
            assert state in totest.STAR_STATES_CC
    
    def test_BINARY_STATES_ALL_completeness(self):
        """Test BINARY_STATES_ALL contains expected binary states."""
        expected_states = [
            'initially_single_star', 'detached', 'RLO1', 'RLO2',
            'contact', 'disrupted', 'merged', 'initial_RLOF'
        ]
        assert totest.BINARY_STATES_ALL == expected_states
    
    def test_BINARY_STATES_CC_consistency(self):
        """Test BINARY_STATES_CC equals BINARY_STATES_ALL."""
        assert totest.BINARY_STATES_CC == totest.BINARY_STATES_ALL
    
    def test_BINARY_EVENTS_ALL_completeness(self):
        """Test BINARY_EVENTS_ALL contains expected events."""
        expected_events = [
            None, 'CC1', 'CC2', 'ZAMS', 'oRLO1', 'oRLO2',
            'oCE1', 'oCE2', 'oDoubleCE1', 'oDoubleCE2',
            'CO_contact', 'MaxTime_exceeded', 'maxtime',
            'oMerging1', 'oMerging2'
        ]
        for event in expected_events:
            assert event in totest.BINARY_EVENTS_ALL
    
    def test_BINARY_EVENTS_OF_SN_OR_AFTER_DETACHED_excludes_core_collapse(self):
        """Test that SN events exclude core collapse events."""
        excluded_events = ['CC1', 'CC2', 'MaxTime_exceeded', 'maxtime']
        for event in excluded_events:
            assert event not in totest.BINARY_EVENTS_OF_SN_OR_AFTER_DETACHED
        
        # Should contain other events
        assert 'oRLO1' in totest.BINARY_EVENTS_OF_SN_OR_AFTER_DETACHED
        assert None in totest.BINARY_EVENTS_OF_SN_OR_AFTER_DETACHED
    
    def test_UNDEFINED_STATES_structure(self):
        """Test UNDEFINED_STATES contains expected structure."""
        assert isinstance(totest.UNDEFINED_STATES, list)
        # Should contain tuples for CO states with specific binary states
        if len(totest.UNDEFINED_STATES) > 0:
            for state_tuple in totest.UNDEFINED_STATES:
                assert isinstance(state_tuple, tuple)
                assert len(state_tuple) == 4  # (s1, s2, binary_state, event)

class TestPOSYDONFlowChart:
    """Test the POSYDON_FLOW_CHART construction."""
    
    def test_POSYDON_FLOW_CHART_is_dict(self):
        """Test that POSYDON_FLOW_CHART is a dictionary."""
        assert isinstance(totest.POSYDON_FLOW_CHART, dict)
        assert len(totest.POSYDON_FLOW_CHART) > 0
    
    def test_POSYDON_FLOW_CHART_key_structure(self):
        """Test that flow chart keys have correct structure."""
        sample_keys = list(totest.POSYDON_FLOW_CHART.keys())[:10]
        for key in sample_keys:
            assert isinstance(key, tuple)
            assert len(key) == 4  # (s1_state, s2_state, binary_state, event)
            s1, s2, binary_state, event = key
            
            # States should be valid
            assert s1 in totest.STAR_STATES_ALL
            assert s2 in totest.STAR_STATES_ALL
            assert binary_state in totest.BINARY_STATES_ALL
            assert event in totest.BINARY_EVENTS_ALL
    
    def test_POSYDON_FLOW_CHART_value_structure(self):
        """Test that flow chart values are step names."""
        expected_steps = [
            'step_HMS_HMS', 'step_detached', 'step_CO_HMS_RLO',
            'step_CO_HeMS', 'step_CO_HeMS_RLO', 'step_merged',
            'step_CE', 'step_SN', 'step_dco', 'step_end',
            'step_initially_single', 'step_disrupted'
        ]
        
        all_values = set(totest.POSYDON_FLOW_CHART.values())
        for value in all_values:
            assert isinstance(value, str)
            assert value in expected_steps
    
    def test_ZAMS_flow_chart_entries(self):
        """Test ZAMS entries in flow chart."""
        # Check that HMS+HMS ZAMS entries exist
        zams_hms_hms_found = False
        for key, value in totest.POSYDON_FLOW_CHART.items():
            s1, s2, binary_state, event = key
            if (event == 'ZAMS' and s1 in totest.STAR_STATES_ZAMS and 
                s2 in totest.STAR_STATES_ZAMS and binary_state in totest.BINARY_STATES_ZAMS):
                zams_hms_hms_found = True
                assert value == 'step_HMS_HMS'
        
        assert zams_hms_hms_found, "No HMS+HMS ZAMS entries found in flow chart"
    
    def test_detached_flow_chart_entries(self):
        """Test detached binary entries in flow chart."""
        detached_step_detached_found = False
        for key, value in totest.POSYDON_FLOW_CHART.items():
            s1, s2, binary_state, event = key
            if (binary_state == 'detached' and event is None and
                s1 in totest.STAR_STATES_NORMALSTAR and s2 in totest.STAR_STATES_NORMALSTAR):
                detached_step_detached_found = True
                assert value == 'step_detached'
                break
        
        assert detached_step_detached_found, "No normal star detached entries found in flow chart"
    
    def test_CO_HMS_RLO_flow_chart_entries(self):
        """Test CO+HMS RLO entries in flow chart."""
        co_hms_rlo_found = False
        for key, value in totest.POSYDON_FLOW_CHART.items():
            s1, s2, binary_state, event = key
            if ((s1 in totest.STAR_STATES_CO and s2 in totest.STAR_STATES_H_RICH_EVOLVABLE and
                 binary_state == 'RLO2' and event == 'oRLO2') or
                (s2 in totest.STAR_STATES_CO and s1 in totest.STAR_STATES_H_RICH_EVOLVABLE and
                 binary_state == 'RLO1' and event == 'oRLO1')):
                co_hms_rlo_found = True
                # The implementation uses step_CO_HeMS_RLO for both H-rich and He-rich RLO
                assert value in ['step_CO_HMS_RLO', 'step_CO_HeMS_RLO']
        
        assert co_hms_rlo_found, "No CO+HMS RLO entries found in flow chart"
    
    def test_CO_HeMS_flow_chart_entries(self):
        """Test CO+HeMS entries in flow chart."""
        co_hems_found = False
        for key, value in totest.POSYDON_FLOW_CHART.items():
            s1, s2, binary_state, event = key
            if ((s1 in totest.STAR_STATES_CO and s2 in totest.STAR_STATES_HE_RICH_EVOLVABLE and
                 binary_state == 'detached' and event is None) or
                (s2 in totest.STAR_STATES_CO and s1 in totest.STAR_STATES_HE_RICH_EVOLVABLE and
                 binary_state == 'detached' and event is None)):
                co_hems_found = True
                assert value == 'step_CO_HeMS'
        
        assert co_hems_found, "No CO+HeMS entries found in flow chart"
    
    def test_HeMS_HeMS_RLO_merger_entries(self):
        """Test He-rich star + He-rich star RLO leads to merger."""
        he_he_rlo_found = False
        for key, value in totest.POSYDON_FLOW_CHART.items():
            s1, s2, binary_state, event = key
            if (s1 in totest.STAR_STATES_HE_RICH and s2 in totest.STAR_STATES_HE_RICH and
                binary_state in ['RLO1', 'RLO2'] and event in ['oRLO1', 'oRLO2']):
                he_he_rlo_found = True
                assert value == 'step_merged'
        
        assert he_he_rlo_found, "No He+He RLO merger entries found in flow chart"
    
    def test_CE_flow_chart_entries(self):
        """Test common envelope entries in flow chart."""
        ce_found = False
        for key, value in totest.POSYDON_FLOW_CHART.items():
            s1, s2, binary_state, event = key
            if event in ['oCE1', 'oCE2', 'oDoubleCE1', 'oDoubleCE2']:
                ce_found = True
                # CE events may go to step_CE, step_end, or step_disrupted depending on configuration
                assert value in ['step_CE', 'step_end', 'step_disrupted']
        
        assert ce_found, "No CE entries found in flow chart"
    
    def test_SN_flow_chart_entries(self):
        """Test supernova entries in flow chart."""
        sn_found = False
        for key, value in totest.POSYDON_FLOW_CHART.items():
            s1, s2, binary_state, event = key
            if event in ['CC1', 'CC2']:
                sn_found = True
                # Some SN events may go to step_end in certain configurations
                assert value in ['step_SN', 'step_end']
        
        assert sn_found, "No SN entries found in flow chart"
    
    def test_DCO_flow_chart_entries(self):
        """Test double compact object entries in flow chart."""
        dco_found = False
        for key, value in totest.POSYDON_FLOW_CHART.items():
            s1, s2, binary_state, event = key
            if (s1 in totest.STAR_STATES_CO and s2 in totest.STAR_STATES_CO and
                binary_state == 'detached' and event is None):
                dco_found = True
                assert value == 'step_dco'
        
        assert dco_found, "No DCO entries found in flow chart"
    
    def test_end_flow_chart_entries(self):
        """Test step_end entries in flow chart."""
        end_found = False
        for key, value in totest.POSYDON_FLOW_CHART.items():
            s1, s2, binary_state, event = key
            if value == 'step_end':
                end_found = True
                # Should include cases like initial_RLOF, maxtime, etc.
                assert (binary_state == 'initial_RLOF' or 
                        event in ['maxtime', 'MaxTime_exceeded', 'CO_contact'] or
                        (s1 == 'massless_remnant' and s2 == 'massless_remnant') or
                        (s1 in totest.STAR_STATES_CO and s2 == 'massless_remnant') or
                        (s1 == 'massless_remnant' and s2 in totest.STAR_STATES_CO) or
                        (s1 in totest.STAR_STATES_CO and s2 in totest.STAR_STATES_CO and binary_state == 'disrupted'))
        
        assert end_found, "No step_end entries found in flow chart"

class TestFlowChartFunction:
    """Test the flow_chart function."""
    
    def test_flow_chart_default_behavior(self):
        """Test flow_chart function with default parameters."""
        chart = totest.flow_chart()
        
        assert isinstance(chart, dict)
        assert len(chart) > 0
        
        # Should be a copy of POSYDON_FLOW_CHART
        assert chart == totest.POSYDON_FLOW_CHART
        
        # But not the same object
        assert chart is not totest.POSYDON_FLOW_CHART
    
    def test_flow_chart_with_custom_flow_chart(self):
        """Test flow_chart function with custom FLOW_CHART parameter."""
        custom_chart = {('NS', 'NS', 'detached', None): 'step_custom'}
        
        result_chart = totest.flow_chart(FLOW_CHART=custom_chart)
        
        assert result_chart == custom_chart
        # Note: The function may or may not return a copy depending on implementation
    
    def test_flow_chart_with_change_flow_chart(self):
        """Test flow_chart function with CHANGE_FLOW_CHART parameter."""
        change_chart = {('NS', 'NS', 'detached', None): 'step_custom'}
        
        result_chart = totest.flow_chart(CHANGE_FLOW_CHART=change_chart)
        
        # Should start with POSYDON default
        assert len(result_chart) >= len(totest.POSYDON_FLOW_CHART)
        
        # Should have the changed entry
        assert result_chart[('NS', 'NS', 'detached', None)] == 'step_custom'
        
        # Other entries should remain unchanged
        for key, value in totest.POSYDON_FLOW_CHART.items():
            if key != ('NS', 'NS', 'detached', None):
                assert result_chart[key] == value
    
    def test_flow_chart_with_both_parameters(self):
        """Test flow_chart function with both FLOW_CHART and CHANGE_FLOW_CHART."""
        custom_chart = {
            ('BH', 'BH', 'detached', None): 'step_original',
            ('NS', 'NS', 'detached', None): 'step_original'
        }
        change_chart = {('NS', 'NS', 'detached', None): 'step_changed'}
        
        result_chart = totest.flow_chart(FLOW_CHART=custom_chart, 
                                        CHANGE_FLOW_CHART=change_chart)
        
        # Should have both entries, with changed one overridden
        assert result_chart[('BH', 'BH', 'detached', None)] == 'step_original'
        assert result_chart[('NS', 'NS', 'detached', None)] == 'step_changed'
    
    def test_flow_chart_change_nonexistent_key(self):
        """Test flow_chart function changing a key that doesn't exist."""
        change_chart = {('custom_state', 'custom_state', 'custom_binary', 'custom_event'): 'step_custom'}
        
        result_chart = totest.flow_chart(CHANGE_FLOW_CHART=change_chart)
        
        # Should contain all original entries
        for key, value in totest.POSYDON_FLOW_CHART.items():
            assert result_chart[key] == value
        
        # The nonexistent key should not be added (function only changes existing keys)
        assert ('custom_state', 'custom_state', 'custom_binary', 'custom_event') not in result_chart

class TestInitialEccentricityFlowChart:
    """Test the initial_eccentricity_flow_chart function."""
    
    def test_initial_eccentricity_flow_chart_default(self):
        """Test initial_eccentricity_flow_chart with default parameters."""
        chart = totest.initial_eccentricity_flow_chart()
        
        assert isinstance(chart, dict)
        assert len(chart) > 0
        
        # Should modify ZAMS entries to go to step_detached
        zams_modified = False
        for key, value in chart.items():
            s1, s2, binary_state, event = key
            if event == 'ZAMS' and binary_state in totest.BINARY_STATES_ZAMS:
                zams_modified = True
                assert value == 'step_detached'
        
        assert zams_modified, "ZAMS entries not modified correctly"
    
    def test_initial_eccentricity_flow_chart_HMS_HMS_RLO_added(self):
        """Test that HMS+HMS RLO entries are added."""
        chart = totest.initial_eccentricity_flow_chart()
        
        # Should add HMS+HMS RLO entries
        hms_hms_rlo_found = False
        for key, value in chart.items():
            s1, s2, binary_state, event = key
            if (s1 in totest.STAR_STATES_H_RICH_EVOLVABLE and 
                s2 in totest.STAR_STATES_H_RICH_EVOLVABLE and
                binary_state in ['RLO1', 'RLO2'] and event in ['oRLO1', 'oRLO2']):
                hms_hms_rlo_found = True
                assert value == 'step_HMS_HMS_RLO'
        
        assert hms_hms_rlo_found, "HMS+HMS RLO entries not added"
    
    def test_initial_eccentricity_flow_chart_HeMS_HMS_RLO_added(self):
        """Test that HeMS+HMS RLO entries are added."""
        chart = totest.initial_eccentricity_flow_chart()
        
        # Should add HeMS+HMS RLO entries
        hems_hms_rlo_found = False
        for key, value in chart.items():
            s1, s2, binary_state, event = key
            if ((s1 in totest.STAR_STATES_HE_RICH_EVOLVABLE and 
                 s2 in totest.STAR_STATES_H_RICH_EVOLVABLE) or
                (s1 in totest.STAR_STATES_H_RICH_EVOLVABLE and 
                 s2 in totest.STAR_STATES_HE_RICH_EVOLVABLE)) and \
                binary_state in ['RLO1', 'RLO2'] and event in ['oRLO1', 'oRLO2']:
                hems_hms_rlo_found = True
                assert value == 'step_HMS_HMS_RLO'
        
        assert hems_hms_rlo_found, "HeMS+HMS RLO entries not added"
    
    def test_initial_eccentricity_flow_chart_with_custom_flow_chart(self):
        """Test initial_eccentricity_flow_chart with custom FLOW_CHART."""
        custom_chart = {('BH', 'BH', 'detached', None): 'step_custom'}
        
        result_chart = totest.initial_eccentricity_flow_chart(FLOW_CHART=custom_chart)
        
        # Should still have the custom entry
        assert result_chart[('BH', 'BH', 'detached', None)] == 'step_custom'
        
        # Should have added HMS+HMS RLO entries
        hms_hms_rlo_found = any(
            value == 'step_HMS_HMS_RLO' 
            for key, value in result_chart.items()
            if key[2] in ['RLO1', 'RLO2'] and key[3] in ['oRLO1', 'oRLO2']
        )
        assert hms_hms_rlo_found
    
    def test_initial_eccentricity_flow_chart_with_change_flow_chart(self):
        """Test initial_eccentricity_flow_chart with CHANGE_FLOW_CHART."""
        change_chart = {('NS', 'NS', 'detached', None): 'step_changed'}
        
        result_chart = totest.initial_eccentricity_flow_chart(CHANGE_FLOW_CHART=change_chart)
        
        # Should have the changed entry
        assert result_chart[('NS', 'NS', 'detached', None)] == 'step_changed'
        
        # Should still modify ZAMS entries
        zams_modified = any(
            value == 'step_detached' 
            for key, value in result_chart.items()
            if key[3] == 'ZAMS' and key[2] in totest.BINARY_STATES_ZAMS
        )
        assert zams_modified
    
    def test_initial_eccentricity_flow_chart_with_both_parameters(self):
        """Test initial_eccentricity_flow_chart with both parameters."""
        custom_chart = {('BH', 'BH', 'detached', None): 'step_original'}
        change_chart = {('BH', 'BH', 'detached', None): 'step_changed'}
        
        result_chart = totest.initial_eccentricity_flow_chart(
            FLOW_CHART=custom_chart, 
            CHANGE_FLOW_CHART=change_chart
        )
        
        # Should have the changed entry
        assert result_chart[('BH', 'BH', 'detached', None)] == 'step_changed'
        
        # Should still have modifications
        assert len(result_chart) > 1  # Should have added HMS+HMS RLO entries

class TestEdgeCasesAndValidation:
    """Test edge cases and validation."""
    
    def test_state_list_no_duplicates(self):
        """Test that state lists don't contain duplicates."""
        assert len(totest.STAR_STATES_ALL) == len(set(totest.STAR_STATES_ALL))
        assert len(totest.BINARY_STATES_ALL) == len(set(totest.BINARY_STATES_ALL))
        assert len(totest.BINARY_EVENTS_ALL) == len(set(totest.BINARY_EVENTS_ALL))
    
    def test_state_consistency_between_lists(self):
        """Test consistency between related state lists."""
        # All elements in derived lists should be in STAR_STATES_ALL
        all_derived_states = (totest.STAR_STATES_CO + 
                             totest.STAR_STATES_NOT_CO + 
                             totest.STAR_STATES_NORMALSTAR)
        
        for state in all_derived_states:
            assert state in totest.STAR_STATES_ALL
        
        # Check that partitions are complete and non-overlapping
        co_and_not_co = set(totest.STAR_STATES_CO + totest.STAR_STATES_NOT_CO)
        assert co_and_not_co == set(totest.STAR_STATES_ALL)
    
    def test_flow_chart_comprehensive_coverage(self):
        """Test that flow chart has reasonable coverage of state space."""
        # Count entries by step type
        step_counts = {}
        for value in totest.POSYDON_FLOW_CHART.values():
            step_counts[value] = step_counts.get(value, 0) + 1
        
        # Should have entries for major step types
        expected_major_steps = ['step_detached', 'step_end', 'step_SN', 'step_CE']
        for step in expected_major_steps:
            assert step in step_counts
            assert step_counts[step] > 0
        
        # step_end should be common (many termination conditions)
        assert step_counts['step_end'] > 50
    
    def test_empty_flow_chart_handling(self):
        """Test functions with empty flow charts."""
        empty_chart = {}
        
        # flow_chart should handle empty input
        result = totest.flow_chart(FLOW_CHART=empty_chart)
        assert result == {}
        
        # With changes to empty chart
        change_chart = {('NS', 'NS', 'detached', None): 'step_test'}
        result = totest.flow_chart(FLOW_CHART=empty_chart, CHANGE_FLOW_CHART=change_chart)
        assert result == {}  # Change should not add new keys
    
    def test_function_signatures(self):
        """Test that functions have expected signatures."""
        # flow_chart function should be callable
        assert isfunction(totest.flow_chart)
        assert isfunction(totest.initial_eccentricity_flow_chart)
        
        # Functions should accept the expected parameters
        import inspect
        
        flow_chart_sig = inspect.signature(totest.flow_chart)
        expected_params = ['FLOW_CHART', 'CHANGE_FLOW_CHART']
        for param in expected_params:
            assert param in flow_chart_sig.parameters
        
        initial_ecc_sig = inspect.signature(totest.initial_eccentricity_flow_chart)
        for param in expected_params:
            assert param in initial_ecc_sig.parameters

class TestIntegrationAndRealism:
    """Test integration aspects and realistic scenarios."""
    
    def test_realistic_binary_evolution_paths(self):
        """Test that realistic binary evolution paths exist in flow chart."""
        chart = totest.flow_chart()
        
        # Test HMS+HMS -> RLO -> CE path
        hms_hms_zams = ('H-rich_Core_H_burning', 'H-rich_Core_H_burning', 'detached', 'ZAMS')
        assert hms_hms_zams in chart
        assert chart[hms_hms_zams] == 'step_HMS_HMS'
        
        # Test CO+HMS detached path
        co_hms_detached = ('NS', 'H-rich_Core_H_burning', 'detached', None)
        assert co_hms_detached in chart
        assert chart[co_hms_detached] == 'step_detached'
        
        # Test DCO path
        dco_detached = ('NS', 'BH', 'detached', None)
        assert dco_detached in chart
        assert chart[dco_detached] == 'step_dco'
    
    def test_supernova_handling(self):
        """Test that supernova events are properly handled."""
        chart = totest.flow_chart()
        
        # Find SN entries
        sn_entries = [(k, v) for k, v in chart.items() if k[3] in ['CC1', 'CC2']]
        assert len(sn_entries) > 0
        
        # SN entries should go to step_SN or step_end depending on configuration
        for key, value in sn_entries:
            assert value in ['step_SN', 'step_end']
            
            # SN can happen from various states - this is just checking that events exist
            # The flow chart handles all possible combinations including edge cases
            s1, s2, binary_state, event = key
            # Basic validation that we have valid stellar states
            assert s1 in totest.STAR_STATES_ALL and s2 in totest.STAR_STATES_ALL
    
    def test_common_envelope_handling(self):
        """Test that common envelope events are properly handled."""
        chart = totest.flow_chart()
        
        # Find CE entries
        ce_events = ['oCE1', 'oCE2', 'oDoubleCE1', 'oDoubleCE2']
        ce_entries = [(k, v) for k, v in chart.items() if k[3] in ce_events]
        assert len(ce_entries) > 0
        
        # CE entries should go to step_CE, step_end, or step_disrupted depending on configuration
        for key, value in ce_entries:
            assert value in ['step_CE', 'step_end', 'step_disrupted']
    
    def test_initial_eccentricity_modifications_realistic(self):
        """Test that initial eccentricity modifications are realistic."""
        chart = totest.initial_eccentricity_flow_chart()
        
        # ZAMS should now go to detached
        zams_entries = [(k, v) for k, v in chart.items() if k[3] == 'ZAMS']
        for key, value in zams_entries:
            if key[2] in totest.BINARY_STATES_ZAMS:
                assert value == 'step_detached'
        
        # HMS+HMS RLO should exist
        hms_hms_rlo_entries = [
            (k, v) for k, v in chart.items() 
            if (k[0] in totest.STAR_STATES_H_RICH_EVOLVABLE and 
                k[1] in totest.STAR_STATES_H_RICH_EVOLVABLE and
                k[2] in ['RLO1', 'RLO2'] and k[3] in ['oRLO1', 'oRLO2'])
        ]
        assert len(hms_hms_rlo_entries) > 0
        
        for key, value in hms_hms_rlo_entries:
            assert value == 'step_HMS_HMS_RLO'
