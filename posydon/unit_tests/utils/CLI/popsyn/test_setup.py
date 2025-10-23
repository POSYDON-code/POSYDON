"""Unit tests of posydon/utils/CLI/popsyn/setup.py"""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

import os
import tempfile
from unittest.mock import MagicMock, mock_open, patch

import pytest

# import the module which will be tested
import posydon.utils.CLI.popsyn.setup as totest


class TestSNModelValidation:
    """Test class for SN model validation functions."""

    @patch('posydon.utils.CLI.popsyn.setup.simprop_kwargs_from_ini')
    @patch('posydon.utils.CLI.popsyn.setup.get_SN_MODEL_NAME')
    def test_check_SN_MODEL_validity_with_interp_false(self, mock_get_sn_name, mock_simprop):
        """Test that function returns True when use_interp_values is False."""
        # Setup mock to return use_interp_values=False
        mock_simprop.return_value = {
            'step_SN': (None, {'use_interp_values': False})
        }

        result = totest.check_SN_MODEL_validity("test.ini")
        assert result is True
        # get_SN_MODEL_NAME should not be called when use_interp_values is False
        mock_get_sn_name.assert_not_called()

    @patch('posydon.utils.CLI.popsyn.setup.simprop_kwargs_from_ini')
    @patch('posydon.utils.CLI.popsyn.setup.get_SN_MODEL_NAME')
    def test_check_SN_MODEL_validity_with_valid_model(self, mock_get_sn_name, mock_simprop):
        """Test that function returns True when model is valid."""
        # Setup mocks
        mock_simprop.return_value = {
            'step_SN': (None, {'use_interp_values': True, 'model': 'valid_model'})
        }
        mock_get_sn_name.return_value = "VALID_MODEL"

        result = totest.check_SN_MODEL_validity("test.ini")
        assert result is True
        mock_get_sn_name.assert_called_once()

    @patch('posydon.utils.CLI.popsyn.setup.simprop_kwargs_from_ini')
    @patch('posydon.utils.CLI.popsyn.setup.get_SN_MODEL_NAME')
    def test_check_SN_MODEL_validity_with_invalid_model(self, mock_get_sn_name, mock_simprop):
        """Test that function returns False when model is invalid."""
        # Setup mocks
        mock_simprop.return_value = {
            'step_SN': (None, {'use_interp_values': True, 'model': 'invalid_model'})
        }
        mock_get_sn_name.return_value = None

        result = totest.check_SN_MODEL_validity("test.ini", verbose_on_fail=False)
        assert result is False

    @patch('posydon.utils.CLI.popsyn.setup.simprop_kwargs_from_ini')
    @patch('posydon.utils.CLI.popsyn.setup.get_SN_MODEL_NAME')
    def test_check_SN_MODEL_validity_with_verbose_on_fail(self, mock_get_sn_name, mock_simprop):
        """Test that verbose mode is called on failure."""
        # Setup mocks
        mock_simprop.return_value = {
            'step_SN': (None, {'use_interp_values': True, 'model': 'invalid_model'})
        }
        mock_get_sn_name.return_value = None

        result = totest.check_SN_MODEL_validity("test.ini", verbose_on_fail=True)
        assert result is False
        # Should be called twice: once normal, once with verbose=True
        assert mock_get_sn_name.call_count == 2


class TestIniFileValidation:
    """Test class for INI file validation."""

    def test_validate_ini_file_not_found(self):
        """Test that FileNotFoundError is raised when file doesn't exist."""
        with pytest.raises(FileNotFoundError, match="File .* not found"):
            totest.validate_ini_file("/nonexistent/file.ini")

    @patch('posydon.utils.CLI.popsyn.setup.check_SN_MODEL_validity')
    @patch('os.path.exists')
    def test_validate_ini_file_invalid_sn_model(self, mock_exists, mock_check_sn):
        """Test that ValueError is raised when SN model is invalid."""
        mock_exists.return_value = True
        mock_check_sn.return_value = False

        with pytest.raises(ValueError, match="The step_SN MODEL is not valid"):
            totest.validate_ini_file("test.ini")

    @patch('posydon.utils.CLI.popsyn.setup.check_SN_MODEL_validity')
    @patch('os.path.exists')
    def test_validate_ini_file_success(self, mock_exists, mock_check_sn):
        """Test that validation passes with valid file and model."""
        mock_exists.return_value = True
        mock_check_sn.return_value = True

        # Should not raise any exception
        totest.validate_ini_file("test.ini")


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

    @patch('posydon.utils.CLI.popsyn.setup.validate_ini_file')
    @patch('posydon.utils.CLI.popsyn.setup.binarypop_kwargs_from_ini')
    def test_setup_popsyn_function_too_few_binaries(self, mock_binarypop, mock_validate, mock_args):
        """Test that ValueError is raised when number of binaries is too small."""
        mock_binarypop.return_value = {
            'metallicity': [1.0],
            'number_of_binaries': 5  # Less than job_array (10)
        }

        with pytest.raises(ValueError, match="number of binaries is less than the job array"):
            totest.setup_popsyn_function(mock_args)

    @patch('posydon.utils.CLI.popsyn.setup.validate_ini_file')
    @patch('posydon.utils.CLI.popsyn.setup.binarypop_kwargs_from_ini')
    @patch('posydon.utils.CLI.popsyn.setup.create_python_scripts')
    @patch('posydon.utils.CLI.popsyn.setup.create_slurm_scripts')
    @patch('posydon.utils.CLI.popsyn.setup.create_bash_submit_script')
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
            'metallicity': metallicities,
            'number_of_binaries': 1000
        }

        # Call function
        totest.setup_popsyn_function(mock_args)

        # Verify calls
        mock_validate.assert_called_once_with(mock_args.ini_file)
        mock_python_scripts.assert_called_once_with(mock_args.ini_file)
        assert mock_slurm_scripts.call_count == len(metallicities)
        mock_bash_submit_script.assert_called_once()

    @patch('posydon.utils.CLI.popsyn.setup.validate_ini_file')
    @patch('posydon.utils.CLI.popsyn.setup.binarypop_kwargs_from_ini')
    @patch('posydon.utils.CLI.popsyn.setup.create_python_scripts')
    @patch('posydon.utils.CLI.popsyn.setup.create_slurm_scripts')
    @patch('posydon.utils.CLI.popsyn.setup.create_bash_submit_script')
    @patch('os.makedirs')
    def test_setup_popsyn_function_creates_log_directories(
        self, mock_makedirs, mock_bash_submit_script,
        mock_slurm_scripts, mock_python_scripts, mock_binarypop,
        mock_validate, mock_args
    ):
        """Test that log directories are created for each metallicity."""
        metallicities = [0.1, 1.0]
        mock_binarypop.return_value = {
            'metallicity': metallicities,
            'number_of_binaries': 1000
        }

        totest.setup_popsyn_function(mock_args)

        # Check that makedirs was called for each metallicity
        assert mock_makedirs.call_count == len(metallicities)
        # Verify the directory names
        calls = [call[0][0] for call in mock_makedirs.call_args_list]
        assert any(['1e+00_logs' in call for call in calls])
        assert any(['1e-01_logs' in call for call in calls])


class TestIntegration:
    """Integration tests for the setup module."""

    @patch('posydon.utils.CLI.popsyn.setup.binarypop_kwargs_from_ini')
    @patch('posydon.utils.CLI.popsyn.setup.simprop_kwargs_from_ini')
    @patch('posydon.utils.CLI.popsyn.setup.get_SN_MODEL_NAME')
    @patch('posydon.utils.CLI.popsyn.setup.create_python_scripts')
    @patch('posydon.utils.CLI.popsyn.setup.create_slurm_scripts')
    @patch('posydon.utils.CLI.popsyn.setup.create_bash_submit_script')
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
            'step_SN': (None, {'use_interp_values': False})
        }
        mock_binarypop.return_value = {
            'metallicity': [1.0],
            'number_of_binaries': 100
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
