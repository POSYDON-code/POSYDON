"""Unit tests of posydon/utils/ini_defaults.py
"""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

import os

from pytest import fixture, raises

# import the module which will be tested
import posydon.utils.ini_defaults as totest
from posydon.binary_evol.CE.step_CEE import StepCEE
from posydon.binary_evol.DT.double_CO import DoubleCO
from posydon.binary_evol.DT.step_detached import (
    DETACHED_RUNTIME_KWARGS,
    detached_defaults_from_ini,
    detached_step,
)
from posydon.binary_evol.DT.step_disrupted import DisruptedStep
from posydon.binary_evol.DT.step_initially_single import InitiallySingleStep
from posydon.binary_evol.DT.step_merged import MergedStep
from posydon.binary_evol.MESA.step_mesa import (
    CO_HeMS_RLO_step,
    CO_HeMS_step,
    CO_HMS_RLO_step,
    HMS_HMS_RLO_step,
    MesaGridStep,
    MS_MS_step,
)
from posydon.binary_evol.SN.step_SN import StepSN
from posydon.binary_evol.track_match import MATCHER_ONLY_KWARGS, TrackMatcher
from posydon.popsyn.defaults import default_kwargs
from posydon.popsyn.io import DEFAULT_INI_PATH


class TestElements:
    def test_dir(self):
        elements = {'DEFAULT_INI_PATH', 'IMPORT_KEYS', 'section_defaults',
                    '_read_default_sections', '__authors__', '__builtins__',
                    '__cached__', '__doc__', '__file__', '__loader__',
                    '__name__', '__package__', '__spec__', 'ast', 'copy',
                    'os', 'ConfigParser', 'lru_cache'}
        assert set(dir(totest)) == elements

    def test_default_ini_path(self):
        assert os.path.isfile(totest.DEFAULT_INI_PATH)
        # the population io uses the same file
        assert DEFAULT_INI_PATH == totest.DEFAULT_INI_PATH


class TestSectionDefaults:
    @fixture
    def clear_cache(self):
        totest._read_default_sections.cache_clear()
        yield
        totest._read_default_sections.cache_clear()

    def test_values_are_evaluated(self):
        defaults = totest.section_defaults('step_SN')
        assert defaults['mechanism'] == 'Fryer+12-delayed'
        assert defaults['kick'] is True
        assert defaults['mean_kick_ECSN'] is None
        assert isinstance(defaults['sigma_kick_ECSN'], float)

    def test_import_keys_left_out(self):
        defaults = totest.section_defaults('step_SN')
        assert 'import' not in defaults
        assert 'absolute_import' not in defaults

    def test_exclude(self):
        defaults = totest.section_defaults('step_SN', exclude=['kick'])
        assert 'kick' not in defaults
        assert 'mechanism' in defaults

    def test_returns_copies(self):
        defaults = totest.section_defaults('BinaryPopulation_options')
        defaults['metallicities'].append(42.)
        defaults['number_of_binaries'] = -1
        fresh = totest.section_defaults('BinaryPopulation_options')
        assert 42. not in fresh['metallicities']
        assert fresh['number_of_binaries'] != -1

    def test_unknown_section(self):
        with raises(KeyError, match="not in the default inifile"):
            totest.section_defaults('step_does_not_exist')

    def test_missing_default_ini(self, monkeypatch, clear_cache, tmp_path):
        monkeypatch.setattr(totest, 'DEFAULT_INI_PATH',
                            str(tmp_path / 'missing.ini'))
        with raises(FileNotFoundError, match="Default inifile not found"):
            totest.section_defaults('step_SN')


class TestClassDefaultsFollowIni:
    """The default inifile is the single source of the class defaults."""

    @staticmethod
    def check(cls, section, runtime_keys, exclude=()):
        ini = totest.section_defaults(section, exclude=exclude)
        for key, value in ini.items():
            assert key in cls.DEFAULT_KWARGS, \
                f"{cls.__name__} does not accept [{section}] {key}"
            assert repr(cls.DEFAULT_KWARGS[key]) == repr(value), \
                f"{cls.__name__}.{key} differs from [{section}]"
        # everything else is set at run time
        extra = set(cls.DEFAULT_KWARGS) - set(ini)
        assert extra <= set(runtime_keys), \
            f"{cls.__name__} has defaults not in [{section}]: " \
            f"{extra - set(runtime_keys)}"

    def test_mesa_steps(self):
        runtime = MesaGridStep.RUNTIME_KWARGS
        for cls, section in [(MesaGridStep, 'step_HMS_HMS'),
                             (MS_MS_step, 'step_HMS_HMS'),
                             (CO_HMS_RLO_step, 'step_CO_HMS_RLO'),
                             (CO_HeMS_RLO_step, 'step_CO_HeMS_RLO'),
                             (CO_HeMS_step, 'step_CO_HeMS'),
                             # no own section, uses the one of step_HMS_HMS
                             (HMS_HMS_RLO_step, 'step_HMS_HMS')]:
            self.check(cls, section, runtime)

    def test_detached_steps(self):
        runtime = DETACHED_RUNTIME_KWARGS
        self.check(detached_step, 'step_detached', runtime,
                   MATCHER_ONLY_KWARGS)
        # subclasses also accept the parameters of detached_step
        inherited = set(runtime) | set(detached_step.DEFAULT_KWARGS)
        for cls, section in [(DisruptedStep, 'step_disrupted'),
                             (MergedStep, 'step_merged'),
                             (InitiallySingleStep, 'step_initially_single'),
                             (DoubleCO, 'step_dco')]:
            self.check(cls, section, inherited, MATCHER_ONLY_KWARGS)

    def test_detached_defaults_from_ini(self):
        defaults = detached_defaults_from_ini('step_merged')
        assert set(DETACHED_RUNTIME_KWARGS) <= set(defaults)
        assert not set(defaults) & MATCHER_ONLY_KWARGS
        defaults = detached_defaults_from_ini('step_merged', base={'a': 1})
        assert defaults['a'] == 1

    def test_CE_step(self):
        self.check(StepCEE, 'step_CE', {'metallicity', 'track_matcher'},
                   MATCHER_ONLY_KWARGS)

    def test_SN_step(self):
        self.check(StepSN, 'step_SN', {'RNG'})

    def test_population_defaults(self):
        assert default_kwargs == \
            totest.section_defaults('BinaryPopulation_options')

    def test_matcher_only_kwargs(self):
        assert MATCHER_ONLY_KWARGS == \
            set(TrackMatcher.DEFAULT_KWARGS) - {'metallicity', 'verbose'}
