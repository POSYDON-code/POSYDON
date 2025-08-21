#!/usr/bin/env python3
"""
Script to evolve a few binaries.
Used for validation of the branch.

Author: Max Briel
"""

import sys
import os
import argparse
import warnings
import signal
from posydon.binary_evol.binarystar import BinaryStar, SingleStar
from posydon.binary_evol.simulationproperties import SimulationProperties
from posydon.popsyn.io import simprop_kwargs_from_ini
from posydon.utils.common_functions import orbital_separation_from_period



target_rows = 12
line_length = 140
columns_to_show = ['step_names', 'state', 'event', 'S1_state', 'S1_mass', 'S2_state', 'S2_mass', 'orbital_period']

def load_inlist(verbose):

    sim_kwargs = simprop_kwargs_from_ini('script_data/population_params_default.ini', verbose=verbose)
    metallicity = {'metallicity':1, 'verbose':verbose}

    sim_kwargs['step_HMS_HMS'][1].update(metallicity)
    sim_kwargs['step_CO_HeMS'][1].update(metallicity)
    sim_kwargs['step_CO_HMS_RLO'][1].update(metallicity)
    sim_kwargs['step_CO_HeMS_RLO'][1].update(metallicity)
    sim_kwargs['step_detached'][1].update(metallicity)
    sim_kwargs['step_disrupted'][1].update(metallicity)
    sim_kwargs['step_merged'][1].update(metallicity)
    sim_kwargs['step_initially_single'][1].update(metallicity)

    sim_prop = SimulationProperties(**sim_kwargs)

    sim_prop.load_steps(verbose=verbose)
    return sim_prop

def write_binary_to_screen(binary):
    """Writes a binary DataFrame prettily to the screen
    
    Args:
        binary: BinaryStar object with evolved data
    """
    df = binary.to_df(**{'extra_columns':{'step_names':'str'}})
    
    # Filter to only existing columns
    available_columns = [col for col in columns_to_show if col in df.columns]
    df_filtered = df[available_columns]
    
    # Reset index to use a counter instead of NaN
    df_filtered = df_filtered.reset_index(drop=True)
    
    print("=" * line_length)
    
    # Print the DataFrame
    df_string = df_filtered.to_string(index=True, float_format='%.3f')
    print(df_string)
    
    # Add empty lines to reach exactly 10 rows of output
    current_rows = len(df_filtered) + 1 # add one for header
    
    if current_rows < target_rows:
        # Calculate the width of the output to print empty lines of the same width
        lines = df_string.split('\n')
        if len(lines) > 1:
            # Use the width of the data lines (skip header)
            empty_lines_needed = target_rows - current_rows
            for i in range(empty_lines_needed):
                print("")
    
    print("-" * line_length)


def print_failed_binary(binary,e,  max_error_lines=3):

    print("=" * line_length)
    print(f"🚨 Binary Evolution Failed!")
    print(f"Exception: {type(e).__name__}")
    print(f"Message: {e}")

    # Get the binary's current state and limit output
    try:
        df = binary.to_df(**{'extra_columns':{'step_names':'str'}})
        if len(df) > 0:
            # Select only the desired columns

            available_columns = [col for col in columns_to_show if col in df.columns]
            df_filtered = df[available_columns]
            
            # Reset index to use a counter instead of NaN
            df_filtered = df_filtered.reset_index(drop=True)
    
            # Limit to max_error_lines
            if len(df_filtered) > max_error_lines:
                df_filtered = df_filtered.tail(max_error_lines)
                print(f"\nShowing last {max_error_lines} evolution steps before failure:")
            else:
                print(f"\nEvolution steps before failure ({len(df_filtered)} steps):")
            
            df_string = df_filtered.to_string(index=True, float_format='%.3f')
            print(df_string)

            current_rows = len(df_filtered) + 1 + 5  # add one for header
            empty_lines_needed = target_rows - current_rows
            for i in range(empty_lines_needed):
                print("")
        else:
            print("\nNo evolution steps recorded before failure.")
    except Exception as inner_e:
        print(f"\nCould not retrieve binary state: {inner_e}")

    print("-" * line_length)

def evolve_binary(binary):
    
    # Capture warnings during evolution
    captured_warnings = []
    
    def warning_handler(message, category, filename, lineno, file=None, line=None):
        captured_warnings.append({
            'message': str(message),
            'category': category.__name__,
            'filename': filename,
            'lineno': lineno
        })
    
    # Set up warning capture
    old_showwarning = warnings.showwarning
    warnings.showwarning = warning_handler
    
    try:
        binary.evolve()
        # Display the evolution summary for successful evolution
        write_binary_to_screen(binary)
        
        # Show warnings if any were captured
        if captured_warnings:
            print(f"⚠️  {len(captured_warnings)} warning(s) raised during evolution:")
            for i, warning in enumerate(captured_warnings[:3], 1):  # Show max 3 warnings
                print(f"   {i}. {warning['category']}: {warning['message']}")
            if len(captured_warnings) > 3:
                print(f"   ... and {len(captured_warnings) - 3} more warning(s)")
            elif len(captured_warnings) <= 3:
                for i in range(4-len(captured_warnings)):
                    print("")
        else:
            print(f"No warning(s) raised during evolution\n\n\n\n")
        print("=" * line_length)

    except Exception as e:
        
        # turn off binary alarm in case of exception
        signal.alarm(0)
        
        print_failed_binary(binary, e)
        
        # Show warnings if any were captured before the exception
        if captured_warnings:
            print(f"\n⚠️  {len(captured_warnings)} warning(s) raised before failure:")
            for i, warning in enumerate(captured_warnings[:3], 1):  # Show max 3 warnings
                print(f"   {i}. {warning['category']}: {warning['message']}")
            if len(captured_warnings) > 3:
                print(f"   ... and {len(captured_warnings) - 3} more warning(s)")
        else:
            print(f"No warning(s) raised during evolution\n\n\n\n")
            
        print("=" * line_length)
    finally:
        # Always turn off binary alarm and restore warning handler
        signal.alarm(0)
        warnings.showwarning = old_showwarning
    

def evolve_binaries(verbose):
    """Evolves a few binaries to validate their output
    """
    sim_prop = load_inlist(verbose)

    ########################################
    # Failing binary in matching
    ########################################
    star_1 = SingleStar(**{'mass': 11.948472796094759, 'state': 'H-rich_Core_H_burning','metallicity':1,
                    'natal_kick_array': [231.97383621190582, 5.927334890264575, 1.5990566013567014, 6.137994236518587]})
    star_2 = SingleStar(**{'mass': 7.636958434479617, 'state': 'H-rich_Core_H_burning','metallicity':1,
                    'natal_kick_array': [None, None, None, None]})
    binary = BinaryStar(star_1, star_2, **{'time': 0.0,  'state': 'detached',  'event': 'ZAMS', 
                            'orbital_period':  190925.99636740884,'eccentricity': 0.0}, properties = sim_prop)
    evolve_binary(binary)
    ########################################
    # Failing binary in matching
    ########################################
    star_1 = SingleStar(**{'mass': 30.169861921689556, 'state': 'H-rich_Core_H_burning','metallicity':1,
                'natal_kick_array':  [77.96834852144123, 0.05021460132555987, 2.3146518208348152, 1.733054979982291]})
    star_2 = SingleStar(**{'mass': 10.972734402996027, 'state': 'H-rich_Core_H_burning','metallicity':1,
                    'natal_kick_array': [None, None, None, None]})
    binary = BinaryStar(star_1, star_2, **{'time': 0.0,  'state': 'detached',  'event': 'ZAMS', 
                            'orbital_period':  20479.71919353725,'eccentricity': 0.0}, properties = sim_prop)
    evolve_binary(binary)
    ########################################
    # flipped S1 and S2 ?
    ########################################
    star_1 = SingleStar(**{'mass': 9.474917413943635, 'state': 'H-rich_Core_H_burning','metallicity':1,
                'natal_kick_array':  [133.5713935237759, 4.398754864537542, 2.703102872841114, 1.4633904612711142]})
    star_2 = SingleStar(**{'mass': 9.311073918196263, 'state': 'H-rich_Core_H_burning','metallicity':1,
                    'natal_kick_array': [None, None, None, None]})
    binary = BinaryStar(star_1, star_2, **{'time': 0.0,  'state': 'detached',  'event': 'ZAMS', 
                            'orbital_period':  18.605997832086413,'eccentricity': 0.0}, properties = sim_prop)

    evolve_binary(binary)
    ########################################
    # flipped S1 and S2
    ########################################
    star_1 = SingleStar(**{'mass': 10.438541, 'state': 'H-rich_Core_H_burning','metallicity':1,
                'natal_kick_array':  [None, None, None, None]})
    star_2 = SingleStar(**{'mass': 1.400713, 'state': 'H-rich_Core_H_burning','metallicity':1,
                    'natal_kick_array': [0, 4.190728383757787, 1.1521129697118118, 5.015343794234789]})
    binary = BinaryStar(star_1, star_2, **{'time': 0.0,  'state': 'detached',  'event': 'ZAMS', 
                            'orbital_period':  9.824025,'eccentricity': 0.0}, properties = sim_prop)
    evolve_binary(binary)
    ########################################
    # flipped S1 and S2
    ########################################
    star_1= SingleStar(**{'mass': 9.845907 , 'state': 'H-rich_Core_H_burning','metallicity':1,
                'natal_kick_array':  [0, 4.190728383757787, 1.1521129697118118, 5.015343794234789]})
    star_2 = SingleStar(**{'mass': 9.611029, 'state': 'H-rich_Core_H_burning','metallicity':1,
                    'natal_kick_array': [None, None, None, None]})
    binary = BinaryStar(star_1, star_2, **{'time': 0.0,  'state': 'detached',  'event': 'ZAMS', 
                            'orbital_period':  3.820571,'eccentricity': 0.0}, properties = sim_prop)
    evolve_binary(binary)
    ########################################
    # Normal binary evolution
    ########################################
    star_1= SingleStar(**{'mass': 30.845907 , 'state': 'H-rich_Core_H_burning','metallicity':1,
                'natal_kick_array':  [0, 4.190728383757787, 1.1521129697118118, 5.015343794234789]})
    star_2 = SingleStar(**{'mass': 30.611029, 'state': 'H-rich_Core_H_burning','metallicity':1,
                    'natal_kick_array': [None, None, None, None]})
    binary = BinaryStar(star_1, star_2, **{'time': 0.0,  'state': 'detached',  'event': 'ZAMS', 
                            'orbital_period':  30.820571,'eccentricity': 0.0}, properties = sim_prop)
    evolve_binary(binary)
    ########################################
    # Normal binary
    ########################################
    star_1= SingleStar(**{'mass': 9.213534679594247 , 'state': 'H-rich_Core_H_burning','metallicity':1,
                'natal_kick_array':  [327.5906384501521, 1.7707176050073297, 1.573225822966838, 1.6757313876001914]})
    star_2 = SingleStar(**{'mass': 7.209878522799272, 'state': 'H-rich_Core_H_burning','metallicity':1,
                    'natal_kick_array': [None, None, None, None]})
    binary = BinaryStar(star_1, star_2, **{'time': 0.0,  'state': 'detached',  'event': 'ZAMS', 
                            'orbital_period':  63123.74544474666,'eccentricity': 0.0}, properties = sim_prop)
    evolve_binary(binary)
    ########################################
    # Normal binary
    ########################################
    star_1= SingleStar(**{'mass': 9.561158487732602 , 'state': 'H-rich_Core_H_burning','metallicity':1,
                'natal_kick_array':  [317.5423844462847, 2.9095984678057603, 1.754121288652108, 2.3693917842468784]})
    star_2 = SingleStar(**{'mass': 9.382732464319286, 'state': 'H-rich_Core_H_burning','metallicity':1,
                    'natal_kick_array': [None, None, None, None]})
    binary = BinaryStar(star_1, star_2, **{'time': 0.0,  'state': 'detached',  'event': 'ZAMS', 
                            'orbital_period':  27.77657038557851,'eccentricity': 0.0}, properties = sim_prop)
    evolve_binary(binary)
    ########################################
    #  Normal binary
    ########################################
    star1 = SingleStar(**{'mass': 7.552858,#29829485,
                      'state': 'H-rich_Core_H_burning',
                      'natal_kick_array': [40.91509926587841, 2.6295454150818256, 1.6718337470964977, 6.0408769315244895]})
    star2 = SingleStar(**{'mass': 6.742063, #481560266,
                        'state': 'H-rich_Core_H_burning',
                        'natal_kick_array': [None, None, None, None]})
    binary = BinaryStar(star1, star2, **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 
                        'orbital_period': 17.957531550841225, 'eccentricity': 0.0,},
                        properties=sim_prop)
    evolve_binary(binary)
    ########################################
    #  High BH spin options
    ########################################
    star_1 = SingleStar(**{'mass': 31.616785, 'state': 'H-rich_Core_H_burning',
                'natal_kick_array':  [10, 4.190728383757787, 1.1521129697118118, 5.015343794234789]})
    star_2 = SingleStar(**{'mass': 26.874267, 'state': 'H-rich_Core_H_burning',
                    'natal_kick_array': [None, None, None, None]})
    binary = BinaryStar(star_1, star_2, **{'time': 0.0,  'state': 'detached',  'event': 'ZAMS', 
                            'orbital_period':  501.99252706449792,'eccentricity': 0.0}, properties = sim_prop)
    evolve_binary(binary)
    ########################################
    #  Original a>1 spin error
    ########################################
    star_1 = SingleStar(**{'mass': 18.107506844123645, 'state': 'H-rich_Core_H_burning',
                'natal_kick_array':  [528.2970725443025, 4.190728383757787, 1.1521129697118118, 5.015343794234789]})
    star_2 = SingleStar(**{'mass': 15.641392951875442, 'state': 'H-rich_Core_H_burning',
                    'natal_kick_array': [None, None, None, None]})
    binary = BinaryStar(star_1, star_2, **{'time': 0.0,  'state': 'detached',  'event': 'ZAMS', 
                        'orbital_period':  151.99252706449792,'eccentricity': 0.0}, properties = sim_prop)
    ########################################
    # FIXED disrupted crash
    ########################################
    STAR1 = SingleStar(**{'mass': 52.967313,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':  [0.0, 0.0, 0.0, 0.0]})
    STAR2 = SingleStar(**{'mass': 36.306444,
                        'state': 'H-rich_Core_H_burning',
                        'natal_kick_array':  [0.0, 0.0, 0.0, 0.0]})
    BINARY = BinaryStar(STAR1, STAR2,  
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':12.877004, 'eccentricity': 0.0},
                        properties = sim_prop)
    evolve_binary(BINARY)
    ########################################
    # FIXED error with SN type
    ########################################
    STAR1 = SingleStar(**{'mass': 17.782576, 
                      'state': 'H-rich_Core_H_burning',
                      'natal_kick_array':  [0.0, 0.0, 0.0, 0.0]})
    STAR2 = SingleStar(**{'mass':3.273864,
                      'state': 'H-rich_Core_H_burning',
                      'natal_kick_array':  [0.0, 0.0, 0.0, 0.0]})
    BINARY = BinaryStar(STAR1, STAR2,  
                    **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':4513.150157, 'eccentricity': 0.0},
                    properties = sim_prop)
    evolve_binary(BINARY)
    ########################################
    # FIXED oRLO2 looping 
    ########################################
    STAR1 = SingleStar(**{'mass': 	170.638207, 
                      'state': 'H-rich_Core_H_burning',
                      'natal_kick_array': [4.921294, 4.31745, 1.777768, 3.509656]})
    STAR2 = SingleStar(**{'mass':37.917852,
                        'state': 'H-rich_Core_H_burning',
                        'natal_kick_array': [0.0, 0.0, 0.0, 0.0]})

    BINARY = BinaryStar(STAR1, STAR2,  
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':113.352736, 'eccentricity': 0.0},
                        properties = sim_prop)
    evolve_binary(BINARY)
    ########################################
    # Redirect to step_CO_HeMS (H-rich non-burning?)
    ########################################
    star_1 = SingleStar(**{'mass': 8.333579, 'state': 'H-rich_Core_H_burning',\
                                     'natal_kick_array': [17.125568, 4.101834, 0.917541, 3.961291]})
    star_2 = SingleStar(**{'mass' : 8.208376, 'state' : 'H-rich_Core_H_burning'})
    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period': 66.870417, 'eccentricity': 0.0},
                        properties = sim_prop)
    evolve_binary(binary)
    ########################################
    # FIXED oRLO2 looping 
    ########################################
    star_1 = SingleStar(**{'mass': 16.921378, 'state': 'H-rich_Core_H_burning',\
                                     'natal_kick_array': [268.837139, 5.773527, 2.568105, 2.519068]})
    star_2 = SingleStar(**{'mass' : 16.286318, 'state' : 'H-rich_Core_H_burning'})
    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period': 37.958768, 'eccentricity': 0.0},
                        properties = sim_prop)
    evolve_binary(binary)
    ########################################
    # FIXED? step_detached failure
    ########################################
    STAR1 = SingleStar(**{'mass': 19.787769, 'state': 'H-rich_Core_H_burning',
                                     'natal_kick_array': [24.464803, 0.666314, 1.954698, 5.598975]})
    STAR2 = SingleStar(**{'mass': 7.638741, 'state': 'H-rich_Core_H_burning'})
    BINARY = BinaryStar(STAR1, STAR2,
                    **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':3007.865561, 'eccentricity': 0.0},
                    properties = sim_prop)
    evolve_binary(BINARY)
    ########################################
    # Disrupted binary
    ########################################
    star_1 = SingleStar(**{'mass': 16.921378, 'state': 'H-rich_Core_H_burning',\
                                     'natal_kick_array': [268.837139, 5.773527, 2.568105, 2.519068]})
    star_2 = SingleStar(**{'mass' : 16.286318, 'state' : 'H-rich_Core_H_burning'})

    BINARY = BinaryStar(star_1, star_2,
                    **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':3007.865561, 'eccentricity': 0.0},
                    properties = sim_prop)
    evolve_binary(BINARY)
    ########################################
    # FIXED Detached binary failure (low mass)
    ########################################    
    STAR1 = SingleStar(**{'mass': 9,
                      'state': 'H-rich_Core_H_burning'})
    STAR2 = SingleStar(**{'mass':0.8,
                        'state': 'H-rich_Core_H_burning'})

    BINARY = BinaryStar(STAR1, STAR2,  
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':4513.150157, 'eccentricity': 0.0},
                        properties = sim_prop)
    evolve_binary(BINARY)
    ########################################
    # FIXED SN_TYPE = None crash
    ########################################    
    STAR1 = SingleStar(**{'mass': 17.782576, 
                      'state': 'H-rich_Core_H_burning'})
    STAR2 = SingleStar(**{'mass':3.273864,
                        'state': 'H-rich_Core_H_burning'})

    BINARY = BinaryStar(STAR1, STAR2,  
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':4513.150157, 'eccentricity': 0.0},
                        properties = sim_prop)
    evolve_binary(BINARY)
    ########################################
    # FIXED SN_TYPE errors
    ########################################  
    STAR1 = SingleStar(**{'mass': 6.782576, 
                      'state': 'H-rich_Core_H_burning'})
    STAR2 = SingleStar(**{'mass':3.273864,
                        'state': 'H-rich_Core_H_burning'})

    BINARY = BinaryStar(STAR1, STAR2,  
                    **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':4513.150157, 'eccentricity': 0.0},
                    properties = sim_prop)
    evolve_binary(BINARY)
    ########################################
    # FIXED SN_TYPE errors
    ########################################  
    STAR1 = SingleStar(**{'mass': 	40.638207, 
                      'state': 'H-rich_Core_H_burning',
                      'natal_kick_array': [30.921294, 4.31745, 1.777768, 3.509656]})
    STAR2 = SingleStar(**{'mass':37.917852,
                        'state': 'H-rich_Core_H_burning'})

    BINARY = BinaryStar(STAR1, STAR2,  
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':2113.352736, 'eccentricity': 0.0},
                        properties = sim_prop)
    evolve_binary(BINARY)
    ########################################
    # FIXED ECSN errors?
    ########################################
    STAR1 = SingleStar(**{'mass':  12.376778,
                      'state': 'H-rich_Core_H_burning',
                      'natal_kick_array': [80, 4.31745, 1.777768, 3.509656]})
    STAR2 = SingleStar(**{'mass': 9.711216,
                        'state': 'H-rich_Core_H_burning'})

    BINARY = BinaryStar(STAR1, STAR2,  
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':79.83702, 'eccentricity': 0.0},
                        properties = sim_prop)
    evolve_binary(BINARY)
    ########################################
    # Interpolator masses??
    ########################################
    STAR1 = SingleStar(**{'mass': 7.592921,
                      'state': 'H-rich_Core_H_burning'})
    STAR2 = SingleStar(**{'mass':5.038679 ,
                        'state': 'H-rich_Core_H_burning'})

    BINARY = BinaryStar(STAR1, STAR2,
                    **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':5.537807, 'eccentricity': 0.0},
                    properties = sim_prop)
    evolve_binary(BINARY)
    ########################################
    # Interpolator masses?
    ########################################
    star_1 = SingleStar(**{'mass': 38.741115, 'state': 'H-rich_Core_H_burning',\
                'natal_kick_array': [21.113771, 2.060135, 2.224789, 4.089729]})
    star_2 = SingleStar(**{'mass': 27.776178, 'state': 'H-rich_Core_H_burning',\
                    'natal_kick_array': [282.712103, 0.296252, 1.628433, 5.623812]})

    BINARY = BinaryStar(star_1, star_2,
                    **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period': 93.387072, 'eccentricity': 0.0},
                        properties = sim_prop)
    evolve_binary(BINARY)
    ########################################
    # FIXED NaN spin
    ########################################
    STAR1 = SingleStar(**{'mass': 	70.066924, 
                      'state': 'H-rich_Core_H_burning',
                      'metallicity':0.2})
    STAR2 = SingleStar(**{'mass':   34.183110,
                        'state': 'H-rich_Core_H_burning', 
                        'metallicity':0.2})
    BINARY = BinaryStar(STAR1, STAR2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':5.931492e+03, 
                        'separation': orbital_separation_from_period(5.931492e+03, STAR1.mass, STAR2.mass),
                        'eccentricity': 0.0},
                        properties = sim_prop)
    evolve_binary(BINARY)
    ########################################
    # FIXED NaN spin
    ########################################
    STAR1 = SingleStar(**{'mass': 	28.837286, 
                      'state': 'H-rich_Core_H_burning',})
    STAR2 = SingleStar(**{'mass':   6.874867,
                        'state': 'H-rich_Core_H_burning',})

    BINARY = BinaryStar(STAR1, STAR2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':35.609894, 
                        'separation': orbital_separation_from_period(35.609894, STAR1.mass, STAR2.mass),
                        'eccentricity': 0.0},
                        properties = sim_prop)
    evolve_binary(BINARY)
    ########################################
    # oRLO2 issue
    ########################################
    STAR1 = SingleStar(**{'mass':29.580210,
                      'state': 'H-rich_Core_H_burning'})
    STAR2 = SingleStar(**{'mass': 28.814626,
                        'state': 'H-rich_Core_H_burning'})

    BINARY = BinaryStar(STAR1, STAR2,  
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':40.437993, 'eccentricity': 0.0},
                        properties = sim_prop)
    evolve_binary(BINARY)
    
    ########################################
    # oRLO2 issue
    ########################################
    STAR1 = SingleStar(**{'mass':67.126795,
                      'state': 'H-rich_Core_H_burning'})
    STAR2 = SingleStar(**{'mass': 19.622908,
                        'state': 'H-rich_Core_H_burning'})

    BINARY = BinaryStar(STAR1, STAR2,  
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':1484.768582, 'eccentricity': 0.0},
                        properties = sim_prop)
    evolve_binary(BINARY)

    ########################################
    # oRLO2 issue
    ########################################
    STAR1 = SingleStar(**{'mass': 58.947503,
                        'state': 'H-rich_Core_H_burning'})
    STAR2 = SingleStar(**{'mass': 56.660506,
                        'state': 'H-rich_Core_H_burning'})

    BINARY = BinaryStar(STAR1, STAR2,  
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':2011.300659, 'eccentricity': 0.0},
                        properties = sim_prop)
    evolve_binary(BINARY)
    
    ########################################
    # oRLO2 issue
    ########################################
    STAR1 = SingleStar(**{'mass': 170.638207,
                        'state': 'H-rich_Core_H_burning',
                        'natal_kick_array':[47.979957374424956, 5.317304576107798, 2.7259013166068145, 4.700929589520818]})
    STAR2 = SingleStar(**{'mass': 37.917852,
                        'state': 'H-rich_Core_H_burning',}, )

    BINARY = BinaryStar(STAR1, STAR2,  
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':113.352736, 'eccentricity': 0.0},
                        properties = sim_prop)
    evolve_binary(BINARY)

    ########################################
    # oRLO2 issue
    ########################################
    STAR1 = SingleStar(**{'mass': 109.540207,
                        'state': 'H-rich_Core_H_burning',})
    STAR2 = SingleStar(**{'mass': 84.344530,
                        'state': 'H-rich_Core_H_burning',}, )

    BINARY = BinaryStar(STAR1, STAR2,  
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':5.651896, 'eccentricity': 0.0},
                        properties = sim_prop)
    evolve_binary(BINARY)
    
    ########################################
    # redirect
    ########################################
    STAR1 = SingleStar(**{'mass': 13.889634,
                      'state': 'H-rich_Core_H_burning'})
    STAR2 = SingleStar(**{'mass':0.490231,
                        'state': 'H-rich_Core_H_burning'})

    BINARY = BinaryStar(STAR1, STAR2,  
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':14513.150157, 'eccentricity': 0.0},
                        properties = sim_prop)
    evolve_binary(BINARY)
    ########################################
    # redirect
    ########################################
    STAR1 = SingleStar(**{'mass': 9,
                      'state': 'H-rich_Core_H_burning'})
    STAR2 = SingleStar(**{'mass':0.8,
                        'state': 'H-rich_Core_H_burning'})

    BINARY = BinaryStar(STAR1, STAR2,  
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':4513.150157, 'eccentricity': 0.0},
                        properties = sim_prop)
    evolve_binary(BINARY)
    ########################################
    # Max time
    ########################################
    star_1 = SingleStar(**{'mass': 103.07996766780799, 'state': 'H-rich_Core_H_burning',
                'natal_kick_array':  [0.0, 0.2965418610971261, 2.0789170290719117, 3.207488023705968]})
    star_2 = SingleStar(**{'mass': 83.66522615073987, 'state': 'H-rich_Core_H_burning',
                    'natal_kick_array': [None, None, None, None]})
    binary = BinaryStar(star_1, star_2, **{'time': 0.0,  'state': 'detached',  'event': 'ZAMS', 
                            'orbital_period':  1449.1101985875678,'eccentricity': 0.0}, properties = sim_prop)
    evolve_binary(binary)
    ########################################
    # Max time
    ########################################
    star_1 = SingleStar(**{'mass': 8.860934140643465, 'state': 'H-rich_Core_H_burning',
                'natal_kick_array':  [11.818027275431337, 2.812412688633058, 0.4998731824233789, 2.9272630485628643]})
    star_2 = SingleStar(**{'mass': 8.584716012668551, 'state': 'H-rich_Core_H_burning',
                    'natal_kick_array': [None, None, None, None]})
    binary = BinaryStar(star_1, star_2, **{'time': 0.0,  'state': 'detached',  'event': 'ZAMS', 
                            'orbital_period':  20.82030114750744,'eccentricity': 0.0}, properties = sim_prop)
    evolve_binary(binary)
    ########################################
    # PR421
    ########################################
    STAR1 = SingleStar(**{'mass': 24.035366,
                      'state': 'H-rich_Core_H_burning'})
    STAR2 = SingleStar(**{'mass':  23.187355,
                        'state': 'H-rich_Core_H_burning'})

    BINARY = BinaryStar(STAR1, STAR2,  
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':18.865029, 'eccentricity': 0.0},
                        properties = sim_prop)
    evolve_binary(BINARY)
    
    ########################################
    # CE class
    ########################################
    STAR1 = SingleStar(**{'mass':33.964274,
                      'state': 'H-rich_Core_H_burning'})
    STAR2 = SingleStar(**{'mass': 28.98149,
                        'state': 'H-rich_Core_H_burning'})

    BINARY = BinaryStar(STAR1, STAR2,  
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':82.370989, 'eccentricity': 0.0},
                        properties = sim_prop)
    evolve_binary(BINARY)
    
    ########################################
    # PR574 - stepCE fix
    ########################################
    STAR1 = SingleStar(**{'mass':29.580210,
                      'state': 'H-rich_Core_H_burning'})
    STAR2 = SingleStar(**{'mass': 28.814626*0.4,
                        'state': 'H-rich_Core_H_burning'})

    BINARY = BinaryStar(STAR1, STAR2,  
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':300.437993, 'eccentricity': 0.0},
                        properties = sim_prop)
    evolve_binary(BINARY)
    
    ########################################
    # e_ZAMS error
    ########################################
    star1 = SingleStar(**{'mass': 8.161885721822461,
                      'state': 'H-rich_Core_H_burning',
                      'natal_kick_array': [None, None, None, None]})
    star2 = SingleStar(**{'mass': 3.5907829421526154,
                        'state': 'H-rich_Core_H_burning',
                        'natal_kick_array': [None, None, None, None]})

    binary = BinaryStar(star1, star2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period': 36.873457164644144, 'eccentricity': 0.0},
                        properties = sim_prop)
    evolve_binary(binary)
    
    ########################################
    # e_ZAMS error
    ########################################
    star1 = SingleStar(**{'mass': 35.24755025317775,
                      'state': 'H-rich_Core_H_burning',
                      'natal_kick_array': [19.755993125895806, 0.37149222852233904, 1.6588846085306563,
                                          1.434617029858906]})
    star2 = SingleStar(**{'mass': 30.000450298072902,
                        'state': 'H-rich_Core_H_burning',
                        'natal_kick_array': [None, None, None, None]})

    binary = BinaryStar(star1, star2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period': 24060.02101364665, 'eccentricity': 0.8085077857996965},
                        properties = sim_prop)
    evolve_binary(binary)
    
    ########################################
    # e_ZAMS error
    ########################################
    star1 = SingleStar(**{'mass': 11.862930493162692,
                      'state': 'H-rich_Core_H_burning',
                      'natal_kick_array': [None, None, None, None]})
    star2 = SingleStar(**{'mass': 1.4739109294156703,
                        'state': 'H-rich_Core_H_burning',
                        'natal_kick_array': [None, None, None, None]})

    binary = BinaryStar(star1, star2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period': 4111.083887312003, 'eccentricity':0.0},
                        properties = sim_prop)
    evolve_binary(binary)
    
    ########################################
    # e_ZAMS error
    ########################################
    star1 = SingleStar(**{'mass': 8.527361341212108,
                      'state': 'H-rich_Core_H_burning',
                      'natal_kick_array': [None, None, None, None]})
    star2 = SingleStar(**{'mass': 0.7061748406821822,
                        'state': 'H-rich_Core_H_burning',
                        'natal_kick_array': [None, None, None, None]})

    binary = BinaryStar(star1, star2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period': 2521.1927287891444, 'eccentricity':0.0},
                        properties = sim_prop)    
    evolve_binary(binary)
    
    ########################################
    # e_ZAMS error
    ########################################
    star1 = SingleStar(**{'mass': 13.661942533447398   ,#29829485,
                      'state': 'H-rich_Core_H_burning',
                      'natal_kick_array': [None, None, None, None]})
    star2 = SingleStar(**{'mass': 4.466151109802313    , #481560266,
                        'state': 'H-rich_Core_H_burning',
                        'natal_kick_array': [None, None, None, None]})

    binary = BinaryStar(star1, star2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':3110.1346707516914, 'eccentricity':0.0},
                        properties = sim_prop)
    evolve_binary(binary)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Evolve binaries for validation.')
    parser.add_argument('--verbose', '-v', action='store_true', default=False,
                        help='Enable verbose output (default: False)')
    args = parser.parse_args()

    evolve_binaries(verbose=args.verbose)
