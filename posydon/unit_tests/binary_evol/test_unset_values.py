"""Unit tests of how unset properties of BinaryStar and SingleStar are
represented: numerical properties are NaN instead of None (see issue #236).
"""

__authors__ = [
    "Max Briel <max.briel@unige.ch>"
]

import numpy as np
import pandas as pd
from pytest import fixture, mark

from posydon.binary_evol.binarystar import (
    BINARYPROPERTIES,
    BinaryStar,
    unset_binary_value,
)
from posydon.binary_evol.singlestar import (
    STARPROPERTIES,
    SingleStar,
    unset_star_value,
)
from posydon.binary_evol.SN.step_SN import StepSN
from posydon.grids.SN_MODELS import SN_MODELS
from posydon.popsyn.io import BINARYPROPERTIES_DTYPES, STARPROPERTIES_DTYPES
from posydon.utils.common_functions import (
    MT_CASE_NO_RLO,
    get_binary_state_and_event_and_mt_case,
    infer_mass_transfer_case,
)

NATAL_KICK_KEYS = ['natal_kick_velocity', 'natal_kick_azimuthal_angle',
                   'natal_kick_polar_angle', 'natal_kick_mean_anomaly']


# properties which may still be None: categorical (string) and object
# properties, the binary index, and the SN_MODELS keys (a dict of
# core-collapse results per SN model)
ALLOWED_NONE = {'state', 'event', 'SN_type', 'profile', 'interp1d',
                'properties', 'index', *SN_MODELS.keys()}
ALLOWED_NONE_PREFIXES = ('interp_class_', 'mt_history_',
                         'cumulative_mt_case_')


def is_nan(value):
    return isinstance(value, float) and np.isnan(value)


class TestUnsetValues:
    def test_unset_binary_value(self):
        for key in BINARYPROPERTIES:
            if BINARYPROPERTIES_DTYPES.get(key) == 'float64':
                assert is_nan(unset_binary_value(key))
            else:
                assert unset_binary_value(key) is None
        # the flow chart relies on an unset event being None
        assert unset_binary_value('event') is None
        assert unset_binary_value('not_a_property') is None

    def test_unset_star_value(self):
        for key in STARPROPERTIES:
            if STARPROPERTIES_DTYPES.get(key) in ('float64', 'int'):
                assert is_nan(unset_star_value(key))
        assert unset_star_value('state') == ''
        assert unset_star_value('profile') is None
        assert unset_star_value('not_a_property') is None


class TestDefaults:
    def test_binary_numerical_defaults_are_nan(self):
        binary = BinaryStar()
        for key in BINARYPROPERTIES:
            if key == 'eccentricity':
                assert getattr(binary, key) == 0.0
            elif BINARYPROPERTIES_DTYPES.get(key) == 'float64':
                assert is_nan(getattr(binary, key)), key
                assert is_nan(getattr(binary, key + '_history')[0]), key
        assert binary.state is None
        assert binary.event is None
        assert is_nan(binary.inspiral_time)

    def test_binary_true_anomaly_defaults(self):
        # these are the names step_SN reads and writes
        binary = BinaryStar()
        assert is_nan(binary.true_anomaly_first_SN)
        assert is_nan(binary.true_anomaly_second_SN)
        assert not hasattr(binary, 'true_anomaly_SN1')
        assert not hasattr(binary, 'true_anomaly_SN2')

    def test_binary_numerical_defaults_support_arithmetic(self):
        # None would raise a TypeError here
        binary = BinaryStar()
        assert np.isnan(binary.time + 1.0)
        assert np.isnan(10**binary.lg_mtransfer_rate)

    def test_star_natal_kick_defaults_are_nan(self):
        star = SingleStar()
        for key in NATAL_KICK_KEYS:
            assert is_nan(getattr(star, key)), key

    def test_star_extra_numerical_defaults_are_nan(self):
        star = SingleStar()
        for key in ['spin_orbit_tilt_first_SN', 'spin_orbit_tilt_second_SN',
                    'f_fb', 'm_disk_accreted', 'm_disk_radiated',
                    'h1_mass_ej', 'he4_mass_ej', 'M4', 'mu4', 'Xi', 'sc',
                    'avg_c_in_c_core_at_He_depletion',
                    'co_core_mass_at_He_depletion']:
            assert is_nan(getattr(star, key)), key
        for quantity in ['m_core_CE', 'r_core_CE']:
            for val in [1, 10, 30, 'pure_He_star_10']:
                assert is_nan(getattr(star, f'{quantity}_{val}cent'))
        # categorical and object properties stay None
        assert star.SN_type is None
        assert star.interp1d is None

    def test_no_numerical_none_in_default_objects(self):
        binary = BinaryStar()
        for obj in [binary, binary.star_1, binary.star_2]:
            for key, value in vars(obj).items():
                if value is None:
                    assert key in ALLOWED_NONE or key.startswith(
                        ALLOWED_NONE_PREFIXES), key

    def test_given_values_are_kept(self):
        binary = BinaryStar(time=1.0, separation=2.0)
        assert binary.time == 1.0
        assert binary.separation == 2.0
        star = SingleStar(natal_kick_velocity=3.0)
        assert star.natal_kick_velocity == 3.0


class TestMassTransferCaseWithUnsetValues:
    @mark.parametrize("unset", [None, np.nan])
    def test_infer_mass_transfer_case(self, unset):
        state = 'H-rich_Core_H_burning'
        assert infer_mass_transfer_case(unset, -5.0, state) == MT_CASE_NO_RLO
        assert infer_mass_transfer_case(1.0, unset, state) == MT_CASE_NO_RLO

    @mark.parametrize("unset", [None, np.nan])
    def test_default_binary_is_detached(self, unset):
        star_1 = SingleStar(state='H-rich_Core_H_burning', mass=10.0)
        star_2 = SingleStar(state='H-rich_Core_H_burning', mass=8.0)
        binary = BinaryStar(star_1=star_1, star_2=star_2,
                            rl_relative_overflow_1=unset,
                            rl_relative_overflow_2=unset,
                            lg_mtransfer_rate=unset)
        assert get_binary_state_and_event_and_mt_case(binary) ==\
               ['detached', None, 'None']


class TestNatalKickWithUnsetValues:
    @fixture
    def SN_step(self):
        return StepSN(verbose=False)

    def kicked_binary(self, SN_step, kick):
        star_1 = SingleStar(state='stripped_He_Core_C_depleted', mass=5.0)
        # the collapse already happened, only the orbital kick is missing
        star_1.state, star_1.mass = 'NS', 1.4
        star_1.SN_type, star_1.f_fb = 'CCSN', 0.0
        for key in NATAL_KICK_KEYS:
            setattr(star_1, key, kick)
        star_2 = SingleStar(state='H-rich_Core_H_burning', mass=10.0)
        binary = BinaryStar(star_1=star_1, star_2=star_2, time=1e7,
                            state='detached', event='CC1', separation=100.0,
                            eccentricity=0.0)
        SN_step.RNG = np.random.default_rng(42)
        SN_step.orbital_kick(binary)
        return binary

    def test_nan_kick_is_drawn_like_none(self, SN_step):
        # a NaN kick (e.g., an empty cell in a kick file or a value that
        # went through a DataFrame) must be drawn, not used as a kick
        with_none = self.kicked_binary(SN_step, None)
        with_nan = self.kicked_binary(SN_step, np.nan)
        for key in NATAL_KICK_KEYS:
            value = getattr(with_nan.star_1, key)
            assert pd.notna(value), key
            assert value == getattr(with_none.star_1, key), key
        assert with_nan.state == with_none.state


class TestStepSNUnsetValues:
    def test_reset_other_star_properties(self):
        star = SingleStar(lg_mdot=-5.0, lg_system_mdot=-6.0)
        StepSN(verbose=False)._reset_other_star_properties(star)
        assert is_nan(star.lg_mdot)
        assert is_nan(star.lg_system_mdot)
