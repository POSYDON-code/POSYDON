from posydon.binary_evol.binarystar import BinaryStar, SingleStar
from posydon.utils.common_functions import orbital_separation_from_period


def get_test_binaries(sim_prop):

    test_binaries = []

    ########################################
    # Failing binary in matching
    ########################################
    star_1 = SingleStar(**{'mass': 11.948472796094759, 'state': 'H-rich_Core_H_burning','metallicity':1,
                    'natal_kick_array': [231.97383621190582, 5.927334890264575, 1.5990566013567014, 6.137994236518587]})
    star_2 = SingleStar(**{'mass': 7.636958434479617, 'state': 'H-rich_Core_H_burning','metallicity':1,
                    'natal_kick_array': [0.0, 0.0, 0.0, 0.0]})
    binary = BinaryStar(star_1, star_2, **{'time': 0.0,  'state': 'detached',  'event': 'ZAMS',
                            'orbital_period':  190925.99636740884,'eccentricity': 0.0}, properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # Failing binary in matching
    ########################################
    star_1 = SingleStar(**{'mass': 30.169861921689556, 'state': 'H-rich_Core_H_burning','metallicity':1,
                'natal_kick_array':  [77.96834852144123, 0.05021460132555987, 2.3146518208348152, 1.733054979982291]})
    star_2 = SingleStar(**{'mass': 10.972734402996027, 'state': 'H-rich_Core_H_burning','metallicity':1,
                    'natal_kick_array': [0.0, 0.0, 0.0, 0.0]})
    binary = BinaryStar(star_1, star_2, **{'time': 0.0,  'state': 'detached',  'event': 'ZAMS',
                            'orbital_period':  20479.71919353725,'eccentricity': 0.0}, properties = sim_prop)

    ########################################
    # flipped S1 and S2 ?
    ########################################
    star_1 = SingleStar(**{'mass': 9.474917413943635, 'state': 'H-rich_Core_H_burning','metallicity':1,
                'natal_kick_array':  [133.5713935237759, 4.398754864537542, 2.703102872841114, 1.4633904612711142]})
    star_2 = SingleStar(**{'mass': 9.311073918196263, 'state': 'H-rich_Core_H_burning','metallicity':1,
                    'natal_kick_array': [0.0, 0.0, 0.0, 0.0]})
    binary = BinaryStar(star_1, star_2, **{'time': 0.0,  'state': 'detached',  'event': 'ZAMS',
                            'orbital_period':  18.605997832086413,'eccentricity': 0.0}, properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # flipped S1 and S2
    ########################################
    star_1 = SingleStar(**{'mass': 10.438541, 'state': 'H-rich_Core_H_burning','metallicity':1,
                'natal_kick_array':  [0.0, 0.0, 0.0, 0.0]})
    star_2 = SingleStar(**{'mass': 1.400713, 'state': 'H-rich_Core_H_burning','metallicity':1,
                    'natal_kick_array': [0.0, 4.190728383757787, 1.1521129697118118, 5.015343794234789]})
    binary = BinaryStar(star_1, star_2, **{'time': 0.0,  'state': 'detached',  'event': 'ZAMS',
                            'orbital_period':  9.824025,'eccentricity': 0.0}, properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # flipped S1 and S2
    ########################################
    star_1= SingleStar(**{'mass': 9.845907 , 'state': 'H-rich_Core_H_burning','metallicity':1,
                'natal_kick_array':  [0.0, 4.190728383757787, 1.1521129697118118, 5.015343794234789]})
    star_2 = SingleStar(**{'mass': 9.611029, 'state': 'H-rich_Core_H_burning','metallicity':1,
                    'natal_kick_array': [0.0, 0.0, 0.0, 0.0]})
    binary = BinaryStar(star_1, star_2, **{'time': 0.0,  'state': 'detached',  'event': 'ZAMS',
                            'orbital_period':  3.820571,'eccentricity': 0.0}, properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # Normal binary evolution
    ########################################
    star_1= SingleStar(**{'mass': 30.845907 , 'state': 'H-rich_Core_H_burning','metallicity':1,
                'natal_kick_array':  [0.0, 4.190728383757787, 1.1521129697118118, 5.015343794234789]})
    star_2 = SingleStar(**{'mass': 30.611029, 'state': 'H-rich_Core_H_burning','metallicity':1,
                    'natal_kick_array': [0.0, 0.0, 0.0, 0.0]})
    binary = BinaryStar(star_1, star_2, **{'time': 0.0,  'state': 'detached',  'event': 'ZAMS',
                            'orbital_period':  30.820571,'eccentricity': 0.0}, properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # Normal binary
    ########################################
    star_1= SingleStar(**{'mass': 9.213534679594247 , 'state': 'H-rich_Core_H_burning','metallicity':1,
                'natal_kick_array':  [327.5906384501521, 1.7707176050073297, 1.573225822966838, 1.6757313876001914]})
    star_2 = SingleStar(**{'mass': 7.209878522799272, 'state': 'H-rich_Core_H_burning','metallicity':1,
                    'natal_kick_array': [0.0, 0.0, 0.0, 0.0]})
    binary = BinaryStar(star_1, star_2, **{'time': 0.0,  'state': 'detached',  'event': 'ZAMS',
                            'orbital_period':  63123.74544474666,'eccentricity': 0.0}, properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # Normal binary
    ########################################
    star_1= SingleStar(**{'mass': 9.561158487732602 , 'state': 'H-rich_Core_H_burning','metallicity':1,
                'natal_kick_array':  [317.5423844462847, 2.9095984678057603, 1.754121288652108, 2.3693917842468784]})
    star_2 = SingleStar(**{'mass': 9.382732464319286, 'state': 'H-rich_Core_H_burning','metallicity':1,
                    'natal_kick_array': [0.0, 0.0, 0.0, 0.0]})
    binary = BinaryStar(star_1, star_2, **{'time': 0.0,  'state': 'detached',  'event': 'ZAMS',
                            'orbital_period':  27.77657038557851,'eccentricity': 0.0}, properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    #  Normal binary
    ########################################
    star_1 = SingleStar(**{'mass': 7.552858,#29829485,
                      'state': 'H-rich_Core_H_burning',
                      'natal_kick_array': [40.91509926587841, 2.6295454150818256, 1.6718337470964977, 6.0408769315244895]})
    star_2 = SingleStar(**{'mass': 6.742063, #481560266,
                        'state': 'H-rich_Core_H_burning',
                        'natal_kick_array': [0.0, 0.0, 0.0, 0.0]})
    binary = BinaryStar(star_1, star_2, **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS',
                        'orbital_period': 17.957531550841225, 'eccentricity': 0.0,},
                        properties=sim_prop)
    test_binaries.append(binary)

    ########################################
    #  High BH spin options
    ########################################
    star_1 = SingleStar(**{'mass': 31.616785, 'state': 'H-rich_Core_H_burning',
                'natal_kick_array':  [10, 4.190728383757787, 1.1521129697118118, 5.015343794234789]})
    star_2 = SingleStar(**{'mass': 26.874267, 'state': 'H-rich_Core_H_burning',
                    'natal_kick_array': [0.0, 0.0, 0.0, 0.0]})
    binary = BinaryStar(star_1, star_2, **{'time': 0.0,  'state': 'detached',  'event': 'ZAMS',
                            'orbital_period':  501.99252706449792,'eccentricity': 0.0}, properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    #  Original a>1 spin error
    ########################################
    star_1 = SingleStar(**{'mass': 18.107506844123645, 'state': 'H-rich_Core_H_burning',
                'natal_kick_array':  [528.2970725443025, 4.190728383757787, 1.1521129697118118, 5.015343794234789]})
    star_2 = SingleStar(**{'mass': 15.641392951875442, 'state': 'H-rich_Core_H_burning',
                    'natal_kick_array': [0.0, 0.0, 0.0, 0.0]})
    binary = BinaryStar(star_1, star_2, **{'time': 0.0,  'state': 'detached',  'event': 'ZAMS',
                        'orbital_period':  151.99252706449792,'eccentricity': 0.0}, properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # FIXED disrupted crash
    ########################################
    star_1 = SingleStar(**{'mass': 52.967313,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':  [0.0, 0.0, 0.0, 0.0]})
    star_2 = SingleStar(**{'mass': 36.306444,
                        'state': 'H-rich_Core_H_burning',
                        'natal_kick_array':  [0.0, 0.0, 0.0, 0.0]})
    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':12.877004, 'eccentricity': 0.0},
                        properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # FIXED error with SN type
    ########################################
    star_1 = SingleStar(**{'mass': 17.782576,
                      'state': 'H-rich_Core_H_burning',
                      'natal_kick_array':  [0.0, 0.0, 0.0, 0.0]})
    star_2 = SingleStar(**{'mass':3.273864,
                      'state': 'H-rich_Core_H_burning',
                      'natal_kick_array':  [0.0, 0.0, 0.0, 0.0]})
    binary = BinaryStar(star_1, star_2,
                    **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':4513.150157, 'eccentricity': 0.0},
                    properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # FIXED oRLO2 looping
    ########################################
    star_1 = SingleStar(**{'mass': 	170.638207,
                      'state': 'H-rich_Core_H_burning',
                      'natal_kick_array': [4.921294, 4.31745, 1.777768, 3.509656]})
    star_2 = SingleStar(**{'mass':37.917852,
                        'state': 'H-rich_Core_H_burning',
                        'natal_kick_array': [0.0, 0.0, 0.0, 0.0]})

    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':113.352736, 'eccentricity': 0.0},
                        properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # Redirect to step_CO_HeMS (H-rich non-burning?)
    ########################################
    star_1 = SingleStar(**{'mass': 8.333579, 'state': 'H-rich_Core_H_burning',\
                                     'natal_kick_array': [17.125568, 4.101834, 0.917541, 3.961291]})
    star_2 = SingleStar(**{'mass' : 8.208376, 'state' : 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})
    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period': 66.870417, 'eccentricity': 0.0},
                        properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # FIXED oRLO2 looping
    ########################################
    star_1 = SingleStar(**{'mass': 16.921378, 'state': 'H-rich_Core_H_burning',\
                                     'natal_kick_array': [268.837139, 5.773527, 2.568105, 2.519068]})
    star_2 = SingleStar(**{'mass' : 16.286318, 'state' : 'H-rich_Core_H_burning',
                            'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})
    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period': 37.958768, 'eccentricity': 0.0},
                        properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # FIXED? step_detached failure
    ########################################
    star_1 = SingleStar(**{'mass': 19.787769, 'state': 'H-rich_Core_H_burning',
                                     'natal_kick_array': [24.464803, 0.666314, 1.954698, 5.598975]})
    star_2 = SingleStar(**{'mass': 7.638741, 'state': 'H-rich_Core_H_burning',
                           'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})
    binary = BinaryStar(star_1, star_2,
                    **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':3007.865561, 'eccentricity': 0.0},
                    properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # Disrupted binary
    ########################################
    star_1 = SingleStar(**{'mass': 16.921378, 'state': 'H-rich_Core_H_burning',\
                                     'natal_kick_array': [268.837139, 5.773527, 2.568105, 2.519068]})
    star_2 = SingleStar(**{'mass' : 16.286318, 'state' : 'H-rich_Core_H_burning',
                            'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})

    binary = BinaryStar(star_1, star_2,
                    **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':3007.865561, 'eccentricity': 0.0},
                    properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # FIXED Detached binary failure (low mass)
    ########################################
    star_1 = SingleStar(**{'mass': 9,
                      'state': 'H-rich_Core_H_burning',
                           'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})
    star_2 = SingleStar(**{'mass':0.8,
                        'state': 'H-rich_Core_H_burning',
                           'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})

    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':4513.150157, 'eccentricity': 0.0},
                        properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # FIXED SN_TYPE = None crash
    ########################################
    star_1 = SingleStar(**{'mass': 17.782576,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})
    star_2 = SingleStar(**{'mass':3.273864,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})

    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':4513.150157, 'eccentricity': 0.0},
                        properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # FIXED SN_TYPE errors
    ########################################
    star_1 = SingleStar(**{'mass': 6.782576,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})
    star_2 = SingleStar(**{'mass':3.273864,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})

    binary = BinaryStar(star_1, star_2,
                    **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':4513.150157, 'eccentricity': 0.0},
                    properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # FIXED SN_TYPE errors
    ########################################
    star_1 = SingleStar(**{'mass': 	40.638207,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array': [30.921294, 4.31745, 1.777768, 3.509656]})
    star_2 = SingleStar(**{'mass':37.917852,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})

    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':2113.352736, 'eccentricity': 0.0},
                        properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # FIXED ECSN errors?
    ########################################
    star_1 = SingleStar(**{'mass':  12.376778,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array': [80, 4.31745, 1.777768, 3.509656]})
    star_2 = SingleStar(**{'mass': 9.711216,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})


    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':79.83702, 'eccentricity': 0.0},
                        properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # Interpolator masses??
    ########################################
    star_1 = SingleStar(**{'mass': 7.592921,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})
    star_2 = SingleStar(**{'mass':5.038679 ,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})

    binary = BinaryStar(star_1, star_2,
                    **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':5.537807, 'eccentricity': 0.0},
                    properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # Interpolator masses?
    ########################################
    star_1 = SingleStar(**{'mass': 38.741115,
                           'state': 'H-rich_Core_H_burning',\
                           'natal_kick_array': [21.113771, 2.060135, 2.224789, 4.089729]})
    star_2 = SingleStar(**{'mass': 27.776178,
                           'state': 'H-rich_Core_H_burning',\
                           'natal_kick_array': [282.712103, 0.296252, 1.628433, 5.623812]})

    binary = BinaryStar(star_1, star_2,
                    **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period': 93.387072, 'eccentricity': 0.0},
                        properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # FIXED NaN spin
    ########################################
    star_1 = SingleStar(**{'mass': 	70.066924,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0],
                          'metallicity':1})
    star_2 = SingleStar(**{'mass':   34.183110,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0],
                          'metallicity':1})
    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':5.931492e+03,
                        'separation': orbital_separation_from_period(5.931492e+03, star_1.mass, star_2.mass),
                        'eccentricity': 0.0},
                        properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # FIXED NaN spin
    ########################################
    star_1 = SingleStar(**{'mass': 	28.837286,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})
    star_2 = SingleStar(**{'mass':   6.874867,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})

    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':35.609894,
                        'separation': orbital_separation_from_period(35.609894, star_1.mass, star_2.mass),
                        'eccentricity': 0.0},
                        properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # oRLO2 issue
    ########################################
    star_1 = SingleStar(**{'mass':29.580210,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})
    star_2 = SingleStar(**{'mass': 28.814626,
                          'state': 'H-rich_Core_H_burning',
                           'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})

    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':40.437993, 'eccentricity': 0.0},
                        properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # oRLO2 issue
    ########################################
    star_1 = SingleStar(**{'mass':67.126795,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})
    star_2 = SingleStar(**{'mass': 19.622908,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})

    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':1484.768582, 'eccentricity': 0.0},
                        properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # oRLO2 issue
    ########################################
    star_1 = SingleStar(**{'mass': 58.947503,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})
    star_2 = SingleStar(**{'mass': 56.660506,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})

    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':2011.300659, 'eccentricity': 0.0},
                        properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # oRLO2 issue
    ########################################
    star_1 = SingleStar(**{'mass': 170.638207,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[47.979957374424956, 5.317304576107798, 2.7259013166068145, 4.700929589520818]})
    star_2 = SingleStar(**{'mass': 37.917852,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})

    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':113.352736, 'eccentricity': 0.0},
                        properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # oRLO2 issue
    ########################################
    star_1 = SingleStar(**{'mass': 109.540207,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})
    star_2 = SingleStar(**{'mass': 84.344530,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})

    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':5.651896, 'eccentricity': 0.0},
                        properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # redirect
    ########################################
    star_1 = SingleStar(**{'mass': 13.889634,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})
    star_2 = SingleStar(**{'mass':0.490231,
                        'state': 'H-rich_Core_H_burning',
                           'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})

    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':14513.150157, 'eccentricity': 0.0},
                        properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # redirect
    ########################################
    star_1 = SingleStar(**{'mass': 9,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})
    star_2 = SingleStar(**{'mass':0.8,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})

    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':4513.150157, 'eccentricity': 0.0},
                        properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # Max time
    ########################################
    star_1 = SingleStar(**{'mass': 103.07996766780799,
                           'state': 'H-rich_Core_H_burning',
                           'natal_kick_array':  [0.0, 0.2965418610971261, 2.0789170290719117, 3.207488023705968]})
    star_2 = SingleStar(**{'mass': 83.66522615073987,
                           'state': 'H-rich_Core_H_burning',
                           'natal_kick_array': [0.0, 0.0, 0.0, 0.0]})
    binary = BinaryStar(star_1, star_2, **{'time': 0.0,  'state': 'detached',  'event': 'ZAMS',
                            'orbital_period':  1449.1101985875678,'eccentricity': 0.0}, properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # Max time
    ########################################
    star_1 = SingleStar(**{'mass': 8.860934140643465,
                           'state': 'H-rich_Core_H_burning',
                           'natal_kick_array':  [11.818027275431337, 2.812412688633058, 0.4998731824233789, 2.9272630485628643]})
    star_2 = SingleStar(**{'mass': 8.584716012668551,
                           'state': 'H-rich_Core_H_burning',
                           'natal_kick_array': [0.0, 0.0, 0.0, 0.0]})
    binary = BinaryStar(star_1, star_2, **{'time': 0.0,  'state': 'detached',  'event': 'ZAMS',
                            'orbital_period':  20.82030114750744,'eccentricity': 0.0}, properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # PR421
    ########################################
    star_1 = SingleStar(**{'mass': 24.035366,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})
    star_2 = SingleStar(**{'mass':  23.187355,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})

    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':18.865029, 'eccentricity': 0.0},
                        properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # CE class
    ########################################
    star_1 = SingleStar(**{'mass':33.964274,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})
    star_2 = SingleStar(**{'mass': 28.98149,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})

    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':82.370989, 'eccentricity': 0.0},
                        properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # PR574 - stepCE fix
    ########################################
    star_1 = SingleStar(**{'mass':29.580210,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})
    star_2 = SingleStar(**{'mass': 28.814626*0.4,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array':[0.0, 0.0, 0.0, 0.0]})

    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':300.437993, 'eccentricity': 0.0},
                        properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # e_ZAMS error
    ########################################
    star_1 = SingleStar(**{'mass': 8.161885721822461,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]})
    star_2 = SingleStar(**{'mass': 3.5907829421526154,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]})

    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period': 36.873457164644144, 'eccentricity': 0.0},
                        properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # e_ZAMS error
    ########################################
    star_1 = SingleStar(**{'mass': 35.24755025317775,
                      'state': 'H-rich_Core_H_burning',
                      'natal_kick_array': [19.755993125895806, 0.37149222852233904, 1.6588846085306563,
                                          1.434617029858906]})
    star_2 = SingleStar(**{'mass': 30.000450298072902,
                        'state': 'H-rich_Core_H_burning',
                        'natal_kick_array': [0.0, 0.0, 0.0, 0.0]})

    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period': 24060.02101364665, 'eccentricity': 0.8085077857996965},
                        properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # e_ZAMS error
    ########################################
    star_1 = SingleStar(**{'mass': 11.862930493162692,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]})
    star_2 = SingleStar(**{'mass': 1.4739109294156703,
                        'state': 'H-rich_Core_H_burning',
                        'natal_kick_array': [0.0, 0.0, 0.0, 0.0]})

    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period': 4111.083887312003, 'eccentricity':0.0},
                        properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # e_ZAMS error
    ########################################
    star_1 = SingleStar(**{'mass': 8.527361341212108,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]})
    star_2 = SingleStar(**{'mass': 0.7061748406821822,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]})

    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period': 2521.1927287891444, 'eccentricity':0.0},
                        properties = sim_prop)
    test_binaries.append(binary)

    ########################################
    # e_ZAMS error
    ########################################
    star_1 = SingleStar(**{'mass': 13.661942533447398   ,#29829485,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]})
    star_2 = SingleStar(**{'mass': 4.466151109802313    , #481560266,
                          'state': 'H-rich_Core_H_burning',
                          'natal_kick_array': [0.0, 0.0, 0.0, 0.0]})

    binary = BinaryStar(star_1, star_2,
                        **{'time': 0.0, 'state': 'detached', 'event': 'ZAMS', 'orbital_period':3110.1346707516914, 'eccentricity':0.0},
                        properties = sim_prop)
    test_binaries.append(binary)

    return test_binaries
