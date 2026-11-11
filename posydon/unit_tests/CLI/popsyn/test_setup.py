"""Unit tests of posydon/utils/CLI/popsyn/setup.py"""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

import os
from unittest.mock import MagicMock, patch

import pytest

# import the module which will be tested
import posydon.CLI.popsyn.setup as totest
from posydon.binary_evol.SN.step_SN import StepSN
from posydon.popsyn.io import write_resolved_inifile
from posydon.utils.posydonwarning import StepWarning


@pytest.fixture(autouse=True)
def mock_write_resolved():
    """Do not write resolved ini files from mocked (non-existing) ini files."""
    with patch('posydon.CLI.popsyn.setup.write_resolved_inifile') as mock:
        yield mock


class TestSNModelValidation:
    """Test class for SN model validation functions."""

    @patch('posydon.CLI.popsyn.setup.simprop_kwargs_from_ini')
    @patch('posydon.CLI.popsyn.setup.get_SN_MODEL_NAME')
    def test_check_SN_MODEL_validity_with_interp_false(self, mock_get_sn_name, mock_simprop):
        """Test that function returns True when use_interp_values is False."""
        # Setup mock to return use_interp_values=False
        mock_simprop.return_value = {
            'step_SN': (StepSN, {'use_interp_values': False})
        }

        result = totest.check_SN_MODEL_validity("test.ini")
        assert result is True
        # get_SN_MODEL_NAME should not be called when use_interp_values is False
        mock_get_sn_name.assert_not_called()

    @patch('posydon.CLI.popsyn.setup.simprop_kwargs_from_ini')
    @patch('posydon.CLI.popsyn.setup.get_SN_MODEL_NAME')
    def test_check_SN_MODEL_validity_with_valid_model(self, mock_get_sn_name, mock_simprop):
        """Test that function returns True when model is valid."""
        # Setup mocks
        mock_simprop.return_value = {
            'step_SN': (StepSN, {'use_interp_values': True, 'model': 'valid_model'})
        }
        mock_get_sn_name.return_value = "VALID_MODEL"

        result = totest.check_SN_MODEL_validity("test.ini")
        assert result is True
        mock_get_sn_name.assert_called_once()

    @patch('posydon.CLI.popsyn.setup.simprop_kwargs_from_ini')
    @patch('posydon.CLI.popsyn.setup.get_SN_MODEL_NAME')
    def test_check_SN_MODEL_validity_with_invalid_model(self, mock_get_sn_name, mock_simprop):
        """Test that function returns False when model is invalid."""
        # Setup mocks
        mock_simprop.return_value = {
            'step_SN': (StepSN, {'use_interp_values': True, 'model': 'invalid_model'})
        }
        mock_get_sn_name.return_value = None

        result = totest.check_SN_MODEL_validity("test.ini", verbose_on_fail=False)
        assert result is False

    @patch('posydon.CLI.popsyn.setup.simprop_kwargs_from_ini')
    @patch('posydon.CLI.popsyn.setup.get_SN_MODEL_NAME')
    def test_check_SN_MODEL_validity_with_verbose_on_fail(self, mock_get_sn_name, mock_simprop):
        """Test that verbose mode is called on failure."""
        # Setup mocks
        mock_simprop.return_value = {
            'step_SN': (StepSN, {'use_interp_values': True, 'model': 'invalid_model'})
        }
        mock_get_sn_name.return_value = None

        result = totest.check_SN_MODEL_validity("test.ini", verbose_on_fail=True)
        assert result is False
        # Should be called twice: once normal, once with verbose=True
        assert mock_get_sn_name.call_count == 2

    @patch('posydon.CLI.popsyn.setup.simprop_kwargs_from_ini')
    @patch('posydon.CLI.popsyn.setup.get_SN_MODEL_NAME')
    def test_check_SN_MODEL_validity_custom_class(self, mock_get_sn_name,
                                                  mock_simprop):
        """A class that is not a StepSN skips the SN model check."""
        class CustomSN:
            pass

        mock_simprop.return_value = {'step_SN': (CustomSN, {'kick': 1.0})}
        assert totest.check_SN_MODEL_validity("test.ini") is True
        mock_get_sn_name.assert_not_called()

    @patch('posydon.CLI.popsyn.setup.simprop_kwargs_from_ini')
    @patch('posydon.CLI.popsyn.setup.get_SN_MODEL_NAME')
    def test_check_SN_MODEL_validity_subclass_uses_defaults(
            self, mock_get_sn_name, mock_simprop):
        """A StepSN subclass is checked with its defaults plus the section."""
        class MySN(StepSN):
            pass

        mock_simprop.return_value = {'step_SN': (MySN, {'engine': 'N20'})}
        mock_get_sn_name.return_value = "VALID_MODEL"
        assert totest.check_SN_MODEL_validity("test.ini") is True
        checked_model = mock_get_sn_name.call_args[0][0]
        assert checked_model['engine'] == 'N20'
        for key in ('use_interp_values', 'mechanism'):
            assert checked_model[key] == StepSN.DEFAULT_KWARGS[key]

    def test_validate_ini_file_with_custom_SN_step(self, tmp_path):
        """A [step_SN] section with another class passes the setup checks."""
        (tmp_path / "my_SN.py").write_text(
            "class MySNStep:\n"
            "    def __init__(self, **kwargs):\n"
            "        pass\n")
        ini_file = tmp_path / "population_params.ini"
        ini_file.write_text(
            "[step_SN]\n"
            "  import = None\n"
            f"  absolute_import = ['{tmp_path / 'my_SN.py'}', 'MySNStep']\n"
            "  kick_scale = 100.0\n")
        totest.validate_ini_file(str(ini_file))


class TestIniFileValidation:
    """Test class for INI file validation."""

    def test_validate_ini_file_not_found(self):
        """Test that FileNotFoundError is raised when file doesn't exist."""
        with pytest.raises(FileNotFoundError, match="File .* not found"):
            totest.validate_ini_file("/nonexistent/file.ini")

    @patch('posydon.CLI.popsyn.setup.check_flow_steps')
    @patch('posydon.CLI.popsyn.setup.check_SN_MODEL_validity')
    @patch('os.path.exists')
    def test_validate_ini_file_invalid_sn_model(self, mock_exists, mock_check_sn,
                                                mock_check_flow):
        """Test that ValueError is raised when SN model is invalid."""
        mock_exists.return_value = True
        mock_check_sn.return_value = False

        with pytest.raises(ValueError, match="The step_SN MODEL is not valid"):
            totest.validate_ini_file("test.ini")

    @patch('posydon.CLI.popsyn.setup.check_flow_steps')
    @patch('posydon.CLI.popsyn.setup.check_SN_MODEL_validity')
    @patch('os.path.exists')
    def test_validate_ini_file_success(self, mock_exists, mock_check_sn,
                                       mock_check_flow):
        """Test that validation passes with valid file and model."""
        mock_exists.return_value = True
        mock_check_sn.return_value = True

        # Should not raise any exception
        totest.validate_ini_file("test.ini")
        mock_check_flow.assert_called_once_with("test.ini")


DEFAULT_INI = os.path.join(os.path.dirname(totest.__file__), '..', '..',
                           'popsyn', 'population_params_default.ini')


def _flow(*step_names):
    """Return a flow-chart function whose flow uses the given steps."""
    def flow_function():
        return {(i, None, None, None): name
                for i, name in enumerate(step_names)}
    return flow_function


class TestFlowStepCheck:
    """Test class for the check of the flow against the configured steps."""

    @pytest.fixture(autouse=True)
    def no_sections_for_mocked_ini(self):
        """The mocked 'test.ini' does not exist, so it has no sections."""
        real_step_sections = totest._step_sections
        with patch('posydon.CLI.popsyn.setup._step_sections',
                   side_effect=lambda ini: (set() if ini == "test.ini"
                                            else real_step_sections(ini))):
            yield

    @patch('posydon.CLI.popsyn.setup.simprop_kwargs_from_ini')
    def test_no_flow_section(self, mock_simprop):
        mock_simprop.return_value = {'step_end': (None, {})}
        with pytest.raises(ValueError, match=r"no \[flow\] section"):
            totest.check_flow_steps("test.ini")

    @patch('posydon.CLI.popsyn.setup.simprop_kwargs_from_ini')
    def test_missing_step_raises(self, mock_simprop):
        mock_simprop.return_value = {
            'flow': (_flow('step_end', 'step_HMS_HMS_RLO'), {}),
            'step_end': (None, {}),
        }
        with pytest.raises(ValueError, match="step_HMS_HMS_RLO"):
            totest.check_flow_steps("test.ini")

    @patch('posydon.CLI.popsyn.setup.simprop_kwargs_from_ini')
    def test_unused_step_warns(self, mock_simprop):
        mock_simprop.return_value = {
            'flow': (_flow('step_end'), {}),
            'step_end': (None, {}),
            'step_ce': (None, {}),
            'HMS_HMS_path': '/some/path',
        }
        with pytest.warns(StepWarning, match="step_ce"):
            totest.check_flow_steps("test.ini")

    @patch('posydon.CLI.popsyn.setup.simprop_kwargs_from_ini')
    def test_matching_steps_pass(self, mock_simprop, recwarn):
        mock_simprop.return_value = {
            'flow': (_flow('step_end', 'step_SN'), {}),
            'step_end': (None, {}),
            'step_SN': (None, {}),
        }
        totest.check_flow_steps("test.ini")
        assert not [w for w in recwarn if "never uses" in str(w.message)]

    @patch('posydon.CLI.popsyn.setup.simprop_kwargs_from_ini')
    def test_unused_default_step_does_not_warn(self, mock_simprop, recwarn):
        # default sections are always stacked in, so a flow that skips one of
        # them (here step_dco) is fine
        mock_simprop.return_value = {
            'flow': (_flow('step_end'), {}),
            'step_end': (None, {}),
            'step_dco': (None, {}),
        }
        totest.check_flow_steps("test.ini")
        assert not [w for w in recwarn if "never uses" in str(w.message)]

    def test_misspelled_section_without_import_warns(self, tmp_path):
        # a typo in a section that only sets parameters, e.g. step_SNN for
        # step_SN, has no import and must still be found
        ini_file = tmp_path / "typo.ini"
        ini_file.write_text("[step_SNN]\n"
                            "  mechanism = 'Fryer+12-delayed'\n")
        with pytest.warns(StepWarning, match="step_SNN"):
            totest.check_flow_steps(str(ini_file))

    def test_default_ini_matches_default_flow(self, recwarn):
        totest.check_flow_steps(DEFAULT_INI)
        assert not [w for w in recwarn if "never uses" in str(w.message)]

    def test_eccentric_flow_needs_HMS_HMS_RLO_step(self, tmp_path):
        override = tmp_path / "eccentric_flow.ini"
        override.write_text(
            "[flow]\n"
            "  import = ['posydon.binary_evol.flow_chart', "
            "'initial_eccentricity_flow_chart']\n"
            "  absolute_import = None\n")
        with pytest.raises(ValueError, match="step_HMS_HMS_RLO"):
            totest.check_flow_steps([DEFAULT_INI, str(override)])


class TestSetupPopsynFunction:
    """Test class for setup_popsyn_function."""

    @pytest.fixture
    def mock_args(self):
        """Create mock command-line arguments."""
        args = MagicMock()
        args.ini_file = "test.ini"
        args.job_array = 10
        args.email = "test@example.com"
        args.partition = "normal"
        args.walltime = "24:00:00"
        args.merge_walltime = "12:00:00"
        args.mem_per_cpu = "4G"
        args.account = "test_account"
        return args

    @patch('posydon.CLI.popsyn.setup.validate_ini_file')
    @patch('posydon.CLI.popsyn.setup.binarypop_kwargs_from_ini')
    def test_setup_popsyn_function_too_few_binaries(self, mock_binarypop, mock_validate, mock_args):
        """Test that ValueError is raised when number of binaries is too small."""
        mock_binarypop.return_value = {
            'metallicities': [1.0],
            'number_of_binaries': 5,  # Less than job_array (10)
            'secondary_mass_scheme': 'q=1'
        }

        with pytest.raises(ValueError, match="number of binaries is less than the job array"):
            totest.setup_popsyn_function(mock_args)

    @patch('posydon.CLI.popsyn.setup.validate_ini_file')
    @patch('posydon.CLI.popsyn.setup.binarypop_kwargs_from_ini')
    @patch('posydon.CLI.popsyn.setup.create_python_scripts')
    @patch('posydon.CLI.popsyn.setup.create_slurm_scripts')
    @patch('posydon.CLI.popsyn.setup.create_bash_submit_script')
    @patch('os.makedirs')
    def test_setup_popsyn_function_success(
        self, mock_makedirs, mock_bash_submit_script,
        mock_slurm_scripts, mock_python_scripts, mock_binarypop,
        mock_validate, mock_args
    ):
        """Test successful setup of population synthesis."""
        # Setup mocks
        metallicities = [0.01, 1.0]
        mock_binarypop.return_value = {
            'metallicities': metallicities,
            'number_of_binaries': 1000,
            'secondary_mass_scheme': 'q=1'
        }

        # Call function
        totest.setup_popsyn_function(mock_args)

        # Verify calls
        mock_validate.assert_called_once_with(mock_args.ini_file)
        resolved_ini = os.path.abspath("test_resolved.ini")
        mock_python_scripts.assert_called_once_with(resolved_ini)
        mock_binarypop.assert_called_once_with(resolved_ini,
                                               use_defaults=False)
        assert mock_slurm_scripts.call_count == len(metallicities)
        mock_bash_submit_script.assert_called_once()

    @patch('posydon.CLI.popsyn.setup.validate_ini_file')
    @patch('posydon.CLI.popsyn.setup.binarypop_kwargs_from_ini')
    @patch('posydon.CLI.popsyn.setup.create_python_scripts')
    @patch('posydon.CLI.popsyn.setup.create_slurm_scripts')
    @patch('posydon.CLI.popsyn.setup.create_bash_submit_script')
    @patch('os.makedirs')
    def test_setup_popsyn_function_creates_log_directories(
        self, mock_makedirs, mock_bash_submit_script,
        mock_slurm_scripts, mock_python_scripts, mock_binarypop,
        mock_validate, mock_args
    ):
        """Test that log directories are created for each metallicity."""
        metallicities = [0.1, 1.0]
        mock_binarypop.return_value = {
            'metallicities': metallicities,
            'number_of_binaries': 1000,
            'secondary_mass_scheme': 'q=1'
        }

        totest.setup_popsyn_function(mock_args)

        # Check that makedirs was called for each metallicity
        assert mock_makedirs.call_count == len(metallicities)
        # Verify the directory names
        calls = [call[0][0] for call in mock_makedirs.call_args_list]
        assert any(['1e+00_logs' in call for call in calls])
        assert any(['1e-01_logs' in call for call in calls])


class TestFlatMassRatioWarning:
    """Test class for flat_mass_ratio mass ratio warning logic."""

    @pytest.fixture
    def mock_args(self):
        """Create mock command-line arguments."""
        args = MagicMock()
        args.ini_file = "test.ini"
        args.job_array = 10
        return args

    @patch('posydon.CLI.popsyn.setup.validate_ini_file')
    @patch('posydon.CLI.popsyn.setup.binarypop_kwargs_from_ini')
    @patch('posydon.CLI.popsyn.setup.create_python_scripts')
    @patch('posydon.CLI.popsyn.setup.create_slurm_scripts')
    @patch('posydon.CLI.popsyn.setup.create_bash_submit_script')
    @patch('posydon.CLI.popsyn.setup.Pwarn')
    @patch('os.makedirs')
    def test_flat_mass_ratio_low_q_issues_warning(
        self, mock_makedirs, mock_pwarn, mock_bash, mock_slurm,
        mock_python, mock_binarypop, mock_validate, mock_args
    ):
        """Test warning is issued when q_min < 0.05 with flat_mass_ratio."""
        mock_binarypop.return_value = {
            'metallicities': [1.0],
            'number_of_binaries': 1000,
            'secondary_mass_scheme': 'flat_mass_ratio',
            'secondary_mass_min': 0.1,
            'primary_mass_max': 150.0,
        }

        totest.setup_popsyn_function(mock_args)

        mock_pwarn.assert_called_once()
        call_args = mock_pwarn.call_args[0]
        assert 'InappropriateValueWarning' in call_args

    @patch('posydon.CLI.popsyn.setup.validate_ini_file')
    @patch('posydon.CLI.popsyn.setup.binarypop_kwargs_from_ini')
    @patch('posydon.CLI.popsyn.setup.create_python_scripts')
    @patch('posydon.CLI.popsyn.setup.create_slurm_scripts')
    @patch('posydon.CLI.popsyn.setup.create_bash_submit_script')
    @patch('posydon.CLI.popsyn.setup.Pwarn')
    @patch('os.makedirs')
    def test_flat_mass_ratio_acceptable_q_no_warning(
        self, mock_makedirs, mock_pwarn, mock_bash, mock_slurm,
        mock_python, mock_binarypop, mock_validate, mock_args
    ):
        """Test no warning is issued when q_min >= 0.05 with flat_mass_ratio."""
        mock_binarypop.return_value = {
            'metallicities': [1.0],
            'number_of_binaries': 1000,
            'secondary_mass_scheme': 'flat_mass_ratio',
            'secondary_mass_min': 5.0,
            'primary_mass_max': 80.0,
        }

        totest.setup_popsyn_function(mock_args)

        mock_pwarn.assert_not_called()

    @patch('posydon.CLI.popsyn.setup.validate_ini_file')
    @patch('posydon.CLI.popsyn.setup.binarypop_kwargs_from_ini')
    @patch('posydon.CLI.popsyn.setup.create_python_scripts')
    @patch('posydon.CLI.popsyn.setup.create_slurm_scripts')
    @patch('posydon.CLI.popsyn.setup.create_bash_submit_script')
    @patch('posydon.CLI.popsyn.setup.Pwarn')
    @patch('os.makedirs')
    def test_non_flat_mass_ratio_scheme_no_warning(
        self, mock_makedirs, mock_pwarn, mock_bash, mock_slurm,
        mock_python, mock_binarypop, mock_validate, mock_args
    ):
        """Test no warning is issued when secondary_mass_scheme is not flat_mass_ratio."""
        mock_binarypop.return_value = {
            'metallicities': [1.0],
            'number_of_binaries': 1000,
            'secondary_mass_scheme': 'q=1',
        }

        totest.setup_popsyn_function(mock_args)

        mock_pwarn.assert_not_called()


class TestIntegration:
    """Integration tests for the setup module."""

    @patch('posydon.CLI.popsyn.setup.binarypop_kwargs_from_ini')
    @patch('posydon.CLI.popsyn.setup.simprop_kwargs_from_ini')
    @patch('posydon.CLI.popsyn.setup.get_SN_MODEL_NAME')
    @patch('posydon.CLI.popsyn.setup.create_python_scripts')
    @patch('posydon.CLI.popsyn.setup.create_slurm_scripts')
    @patch('posydon.CLI.popsyn.setup.create_bash_submit_script')
    def test_full_setup_workflow(
        self, mock_bash, mock_slurm, mock_python,
        mock_get_sn, mock_simprop, mock_binarypop, tmp_path
    ):
        """Test full workflow of setting up a population synthesis run."""
        # Create temporary INI file
        ini_file = tmp_path / "test.ini"
        ini_file.write_text("[test]")

        # Setup mocks
        mock_simprop.return_value = {
            'flow': (lambda: {('a', 'b', 'c', 'd'): 'step_SN'}, {}),
            'step_SN': (StepSN, {'use_interp_values': False})
        }
        mock_binarypop.return_value = {
            'metallicities': [1.0],
            'number_of_binaries': 100,
            'secondary_mass_scheme': 'q=1'
        }

        # Create mock args
        args = MagicMock()
        args.ini_file = str(ini_file)
        args.job_array = 10
        args.email = None
        args.partition = None
        args.walltime = "24:00:00"
        args.merge_walltime = "12:00:00"
        args.mem_per_cpu = "4G"
        args.account = None

        # Change to temp directory
        original_dir = os.getcwd()
        os.chdir(tmp_path)

        try:
            # This should complete without errors
            totest.setup_popsyn_function(args)

            # Verify the workflow
            mock_python.assert_called_once()
            mock_slurm.assert_called_once()
            mock_bash.assert_called_once()
        finally:
            os.chdir(original_dir)


class TestResolvedIni:
    """The setup writes all parameters of the run to one resolved ini file."""

    @patch('posydon.CLI.popsyn.setup.create_slurm_scripts')
    @patch('posydon.CLI.popsyn.setup.create_bash_submit_script')
    def test_setup_writes_resolved_ini(self, mock_bash, mock_slurm,
                                       tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        monkeypatch.setattr(totest, 'write_resolved_inifile',
                            write_resolved_inifile)
        user_ini = tmp_path / "my_run.ini"
        user_ini.write_text("[BinaryPopulation_options]\n"
                            "  number_of_binaries = 100\n"
                            "[step_CE]\n"
                            "  common_envelope_efficiency = 0.5\n")
        args = MagicMock()
        args.ini_file = str(user_ini)
        args.job_array = 10

        totest.setup_popsyn_function(args)

        resolved = tmp_path / "my_run_resolved.ini"
        assert resolved.is_file()
        # the resolved file alone gives the same parameters as the user ini
        # file stacked on the defaults
        from posydon.popsyn.io import simprop_kwargs_from_ini
        assert totest.binarypop_kwargs_from_ini(str(resolved),
                                                use_defaults=False) \
            == totest.binarypop_kwargs_from_ini(str(user_ini))
        sim = simprop_kwargs_from_ini(str(resolved), use_defaults=False)
        assert sim['step_CE'][1]['common_envelope_efficiency'] == 0.5
        assert sim['step_CE'] == simprop_kwargs_from_ini(str(user_ini))['step_CE']
        # the job scripts read the resolved file without the defaults
        run_script = (tmp_path / "run_metallicity.py").read_text()
        assert f"'{resolved}', use_defaults=False" in run_script
