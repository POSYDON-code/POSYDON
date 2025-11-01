"""Unit tests of posydon/utils/CLI/popsyn/check.py"""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

import os
from unittest.mock import MagicMock, patch

import pytest

# import the module which will be tested
import posydon.CLI.popsyn.check as totest
from posydon.utils.common_functions import convert_metallicity_to_string


class TestGetIniFile:
    """Test class for get_ini_file function."""

    @pytest.fixture
    def mock_args(self, tmp_path):
        """Create mock arguments with a temporary run folder."""
        args = MagicMock()
        args.run_folder = str(tmp_path)
        return args

    def test_get_ini_file_no_files(self, mock_args, tmp_path):
        """Test that FileNotFoundError is raised when no INI file exists."""
        with pytest.raises(FileNotFoundError, match="No INI file found"):
            totest.get_ini_file(mock_args)

    def test_get_ini_file_single_file(self, mock_args, tmp_path, capsys):
        """Test successful retrieval of single INI file."""
        # Create a single INI file
        ini_file = tmp_path / "test.ini"
        ini_file.write_text("[test]")

        result = totest.get_ini_file(mock_args)

        assert result == str(ini_file)
        captured = capsys.readouterr()
        assert "Using INI file" in captured.out

    def test_get_ini_file_multiple_files_valid_selection(
        self, mock_args, tmp_path, capsys, monkeypatch
    ):
        """Test selection when multiple INI files exist."""
        # Create multiple INI files
        ini1 = tmp_path / "test1.ini"
        ini2 = tmp_path / "test2.ini"
        ini1.write_text("[test1]")
        ini2.write_text("[test2]")

        # Mock user input to select index 1
        monkeypatch.setattr('builtins.input', lambda _: "1")

        result = totest.get_ini_file(mock_args)

        assert result == str(ini2)

    def test_get_ini_file_multiple_files_invalid_selection(
        self, mock_args, tmp_path, capsys, monkeypatch
    ):
        """Test that invalid selection defaults to first file."""
        # Create multiple INI files
        ini1 = tmp_path / "test1.ini"
        ini2 = tmp_path / "test2.ini"
        ini1.write_text("[test1]")
        ini2.write_text("[test2]")

        # Mock user input with invalid index
        monkeypatch.setattr('builtins.input', lambda _: "99")

        result = totest.get_ini_file(mock_args)

        # Should default to first file
        assert result == str(ini1)
        captured = capsys.readouterr()
        assert "Invalid index" in captured.out

    def test_get_ini_file_multiple_files_non_numeric_input(
        self, mock_args, tmp_path, capsys, monkeypatch
    ):
        """Test that non-numeric input (ValueError) defaults to first file."""
        # Create multiple INI files
        ini1 = tmp_path / "test1.ini"
        ini2 = tmp_path / "test2.ini"
        ini1.write_text("[test1]")
        ini2.write_text("[test2]")

        # Mock user input with non-numeric value (triggers ValueError on line 72)
        monkeypatch.setattr('builtins.input', lambda _: "abc")

        result = totest.get_ini_file(mock_args)

        # Should default to first file
        assert result == str(ini1)
        captured = capsys.readouterr()
        assert "Invalid input" in captured.out




class TestValidateRunFolder:
    """Test class for validate_run_folder function."""

    def test_validate_run_folder_not_exists(self):
        """Test that FileNotFoundError is raised for non-existent folder."""
        with pytest.raises(FileNotFoundError, match="does not exist"):
            totest.validate_run_folder("/nonexistent/folder")

    def test_validate_run_folder_empty(self, tmp_path):
        """Test that ValueError is raised for empty folder."""
        empty_folder = tmp_path / "empty"
        empty_folder.mkdir()

        with pytest.raises(ValueError, match="is empty"):
            totest.validate_run_folder(str(empty_folder))

    def test_validate_run_folder_success(self, tmp_path):
        """Test successful validation of valid folder."""
        # Create folder with a file
        test_file = tmp_path / "test.txt"
        test_file.write_text("test")

        # Should not raise any exception
        totest.validate_run_folder(str(tmp_path))


class TestGetBinaryParams:
    """Test class for get_binary_params function."""

    @patch('posydon.CLI.popsyn.check.binarypop_kwargs_from_ini')
    def test_get_binary_params(self, mock_binarypop):
        """Test that binary parameters are extracted correctly."""
        mock_binarypop.return_value = {
            'metallicity': [0.01, 1.0, 2.0],
            'number_of_binaries': 1000,
            'other_param': 'value'
        }

        num_met, num_bin, metallicities, params = totest.get_binary_params("test.ini")

        assert num_met == 3
        assert num_bin == 1000
        assert metallicities == [0.01, 1.0, 2.0]
        assert params == mock_binarypop.return_value


class TestGetRunConfiguration:
    """Test class for get_run_configuration function."""

    @pytest.fixture
    def mock_args(self, tmp_path):
        """Create mock arguments."""
        args = MagicMock()
        args.run_folder = str(tmp_path)
        return args

    @patch('posydon.CLI.popsyn.check.validate_run_folder')
    @patch('posydon.CLI.popsyn.check.get_ini_file')
    @patch('posydon.CLI.popsyn.check.get_binary_params')
    def test_get_run_configuration(
        self, mock_get_params, mock_get_ini, mock_validate, mock_args
    ):
        """Test successful retrieval of run configuration."""
        # Setup mocks
        mock_get_ini.return_value = "test.ini"
        mock_get_params.return_value = (2, 1000, [1.0, 2.0], {'test': 'params'})

        result = totest.get_run_configuration(mock_args)

        assert len(result) == 5
        assert result[0] == "test.ini"
        assert result[1] == 2
        assert result[2] == 1000
        assert result[3] == [1.0, 2.0]
        assert result[4] == {'test': 'params'}


class TestCheckPopulationFiles:
    """Test class for check_population_files function."""

    def test_check_population_files_all_exist(self, tmp_path, capsys):
        """Test when all population files exist."""
        metallicities = [0.01, 1.0]

        # Create population files
        for met in metallicities:
            str_met = convert_metallicity_to_string(met)
            pop_file = tmp_path / f"{str_met}_Zsun_population.h5"
            pop_file.write_text("dummy")

        all_exist, status = totest.check_population_files(str(tmp_path), metallicities)

        assert all_exist is True
        assert all(status.values())

    def test_check_population_files_some_missing(self, tmp_path, capsys):
        """Test when some population files are missing."""
        metallicities = [0.01, 1.0]

        # Create only one file
        str_met = convert_metallicity_to_string(metallicities[0])
        pop_file = tmp_path / f"{str_met}_Zsun_population.h5"
        pop_file.write_text("dummy")

        all_exist, status = totest.check_population_files(str(tmp_path), metallicities)

        assert all_exist is False
        assert status[metallicities[0]] is True
        assert status[metallicities[1]] is False


class TestCheckBinaryCounts:
    """Test class for check_binary_counts function."""

    @patch('posydon.CLI.popsyn.check.Population')
    def test_check_binary_counts_all_match(self, mock_population, tmp_path, capsys):
        """Test when all binary counts match expected."""
        metallicities = [0.01, 1.0]
        expected_count = 1000

        # Create mock population files
        for met in metallicities:
            str_met = convert_metallicity_to_string(met)
            pop_file = tmp_path / f"{str_met}_Zsun_population.h5"
            pop_file.write_text("dummy")

        # Mock Population to return expected count
        mock_pop_instance = MagicMock()
        mock_pop_instance.number_of_systems = expected_count
        mock_population.return_value = mock_pop_instance

        all_match, counts = totest.check_binary_counts(str(tmp_path), metallicities, expected_count)

        assert all_match is True
        assert all(count == expected_count for count in counts.values())

    @patch('posydon.CLI.popsyn.check.Population')
    def test_check_binary_counts_mismatch(self, mock_population, tmp_path, capsys):
        """Test when binary counts don't match expected."""
        metallicities = [0.01, 1.0]
        expected_count = 1000

        # Create mock population files
        for met in metallicities:
            str_met = convert_metallicity_to_string(met)
            pop_file = tmp_path / f"{str_met}_Zsun_population.h5"
            pop_file.write_text("dummy")

        # Mock Population to return different counts
        def mock_pop_constructor(file):
            mock_pop = MagicMock()
            # First file has correct count, second has wrong count
            if "7e-02" in file or "1e-02" in file:
                mock_pop.number_of_systems = expected_count
            else:
                mock_pop.number_of_systems = 500  # Wrong count
            return mock_pop

        mock_population.side_effect = mock_pop_constructor

        all_match, counts = totest.check_binary_counts(str(tmp_path), metallicities, expected_count)

        assert all_match is False
        captured = capsys.readouterr()
        assert "MISMATCH" in captured.out

    @patch('posydon.CLI.popsyn.check.Population')
    def test_check_binary_counts_error_loading(self, mock_population, tmp_path, capsys):
        """Test when there's an error loading a population file."""
        metallicities = [0.01]
        expected_count = 1000

        # Create mock population file
        str_met = convert_metallicity_to_string(metallicities[0])
        pop_file = tmp_path / f"{str_met}_Zsun_population.h5"
        pop_file.write_text("dummy")

        # Mock Population to raise an exception
        mock_population.side_effect = Exception("File corrupted")

        all_match, counts = totest.check_binary_counts(str(tmp_path), metallicities, expected_count)

        assert all_match is False
        captured = capsys.readouterr()
        assert "ERROR" in captured.out
        assert "File corrupted" in captured.out


class TestCheckRunStatus:
    """Test class for check_run_status function."""

    @patch('posydon.CLI.popsyn.check.check_binary_counts')
    @patch('posydon.CLI.popsyn.check.check_population_files')
    def test_check_run_status_all_pass(self, mock_check_files, mock_check_counts):
        """Test when all files exist and counts match."""
        metallicities = [1.0]

        mock_check_files.return_value = (True, {1.0: True})
        mock_check_counts.return_value = (True, {1.0: 1000})

        files_exist, counts_match, file_status = totest.check_run_status(
            "/test/folder", metallicities, 1000
        )

        assert files_exist is True
        assert counts_match is True
        assert file_status == {1.0: True}

    @patch('posydon.CLI.popsyn.check.check_binary_counts')
    @patch('posydon.CLI.popsyn.check.check_population_files')
    def test_check_run_status_files_missing(self, mock_check_files, mock_check_counts):
        """Test when some files are missing."""
        metallicities = [1.0]

        mock_check_files.return_value = (False, {1.0: False})
        # check_binary_counts should not be called when files don't exist

        files_exist, counts_match, file_status = totest.check_run_status(
            "/test/folder", metallicities, 1000
        )

        assert files_exist is False
        assert counts_match is False
        mock_check_counts.assert_not_called()

    @patch('posydon.CLI.popsyn.check.check_binary_counts')
    @patch('posydon.CLI.popsyn.check.check_population_files')
    def test_check_run_status_counts_mismatch(self, mock_check_files, mock_check_counts):
        """Test when files exist but counts don't match."""
        metallicities = [1.0]

        mock_check_files.return_value = (True, {1.0: True})
        mock_check_counts.return_value = (False, {1.0: 500})  # Wrong count

        files_exist, counts_match, file_status = totest.check_run_status(
            "/test/folder", metallicities, 1000
        )

        assert files_exist is True
        assert counts_match is False
        assert file_status == {1.0: True}


class TestGetExpectedBatchCount:
    """Test class for get_expected_batch_count function."""

    def test_get_expected_batch_count_file_not_found(self, tmp_path):
        """Test when SLURM script doesn't exist."""
        result = totest.get_expected_batch_count(str(tmp_path), "1e+00")
        assert result is None

    def test_get_expected_batch_count_success(self, tmp_path):
        """Test successful extraction of batch count."""
        # Create SLURM script
        metallicity = 1.0
        str_met = convert_metallicity_to_string(metallicity)
        slurm_file = tmp_path / f"{str_met}_Zsun_slurm_array.slurm"
        slurm_content = """#!/bin/bash
#SBATCH --array=0-9
#SBATCH --job-name=test
"""
        slurm_file.write_text(slurm_content)

        result = totest.get_expected_batch_count(str(tmp_path), str_met)
        assert result == 10  # 0-9 inclusive = 10 elements

    def test_get_expected_batch_count_no_array_line(self, tmp_path):
        """Test when SLURM script exists but has no array line (line 329)."""
        # Create SLURM script without #SBATCH --array= line
        metallicity = 1.0
        str_met = convert_metallicity_to_string(metallicity)
        slurm_file = tmp_path / f"{str_met}_Zsun_slurm_array.slurm"
        slurm_content = """#!/bin/bash
#SBATCH --job-name=test
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4G
"""
        slurm_file.write_text(slurm_content)

        result = totest.get_expected_batch_count(str(tmp_path), str_met)
        assert result is None  # Should return None when no array line found

    def test_get_expected_batch_count_array_without_dash(self, tmp_path):
        """Test when array line doesn't contain dash (line 326->323 branch)."""
        # Create SLURM script with array format without dash
        metallicity = 1.0
        str_met = convert_metallicity_to_string(metallicity)
        slurm_file = tmp_path / f"{str_met}_Zsun_slurm_array.slurm"
        slurm_content = """#!/bin/bash
#SBATCH --array=5
#SBATCH --job-name=test
"""
        slurm_file.write_text(slurm_content)

        result = totest.get_expected_batch_count(str(tmp_path), str_met)
        assert result is None  # Should return None when array doesn't have dash
class TestFindMissingBatchIndices:
    """Test class for find_missing_batch_indices function."""

    def test_find_missing_batch_indices(self, tmp_path):
        """Test finding missing batch indices."""
        batch_folder = tmp_path / "batches"
        batch_folder.mkdir()

        # Create some batch files (missing indices 1, 3)
        for idx in [0, 2, 4]:
            batch_file = batch_folder / f"evolution.combined.{idx}.h5"
            batch_file.write_text("dummy")

        expected_count = 5
        missing = totest.find_missing_batch_indices(str(batch_folder), expected_count)

        assert missing == {1, 3}


class TestPrintBatchStatus:
    """Test class for print_batch_status function."""

    def test_print_batch_status_unknown_expected(self, capsys):
        """Test when expected count is unknown."""
        status = totest.print_batch_status("1e+00", None, 5, set())
        assert status == "unknown_expected_count"

    def test_print_batch_status_incomplete(self, capsys):
        """Test when batches are incomplete."""
        status = totest.print_batch_status("1e+00", 10, 8, {2, 5})
        assert status == "incomplete"
        captured = capsys.readouterr()
        assert "MISSING" in captured.out

    def test_print_batch_status_incomplete_many_missing(self, capsys):
        """Test when batches are incomplete with >10 missing indices (lines 390-391)."""
        # Create more than 10 missing indices
        missing_indices = set(range(5, 20))  # 15 missing indices
        status = totest.print_batch_status("1e+00", 20, 5, missing_indices)
        assert status == "incomplete"
        captured = capsys.readouterr()
        assert "MISSING" in captured.out
        # Should show "and X more" message
        assert "and" in captured.out
        assert "more" in captured.out

    def test_print_batch_status_complete(self, capsys):
        """Test when batches are complete."""
        status = totest.print_batch_status("1e+00", 10, 10, set())
        assert status == "complete"

    def test_print_batch_status_extra_files(self, capsys):
        """Test when there are extra files."""
        status = totest.print_batch_status("1e+00", 10, 12, set())
        assert status == "extra_files"


class TestSelectJobId:
    """Test class for select_job_id function."""

    def test_select_job_id_no_logs(self, tmp_path, capsys):
        """Test when no log files exist."""
        result = totest.select_job_id(str(tmp_path), "1e+00")
        assert result is None

    def test_select_job_id_single_job(self, tmp_path):
        """Test with a single job ID."""
        # Create log directory and files
        log_dir = tmp_path / "1e+00_logs"
        log_dir.mkdir()

        (log_dir / "popsyn_12345_0.out").write_text("log")
        (log_dir / "popsyn_12345_1.out").write_text("log")

        result = totest.select_job_id(str(tmp_path), "1e+00")
        assert result == 12345

    def test_select_job_id_multiple_jobs(self, tmp_path, monkeypatch):
        """Test with multiple job IDs and valid selection (line 439->exit branch)."""
        # Create log directory and files
        log_dir = tmp_path / "1e+00_logs"
        log_dir.mkdir()

        (log_dir / "popsyn_12345_0.out").write_text("log")
        (log_dir / "popsyn_67890_0.out").write_text("log")

        # Mock user selection - valid input that returns immediately
        monkeypatch.setattr('builtins.input', lambda _: "1")

        result = totest.select_job_id(str(tmp_path), "1e+00")

        # Should return the second job ID (index 1)
        assert result == 67890

    def test_select_job_id_multiple_jobs_first_index(self, tmp_path, monkeypatch):
        """Test with multiple job IDs selecting first index (also tests 439->exit)."""
        # Create log directory and files
        log_dir = tmp_path / "1e+00_logs"
        log_dir.mkdir()

        (log_dir / "popsyn_12345_0.out").write_text("log")
        (log_dir / "popsyn_67890_0.out").write_text("log")

        # Mock user selection - select index 0
        monkeypatch.setattr('builtins.input', lambda _: "0")

        result = totest.select_job_id(str(tmp_path), "1e+00")

        # Should return the first job ID (index 0)
        assert result == 12345


class TestReadBatchLogFile:
    """Test class for read_batch_log_file function."""

    def test_read_batch_log_file_not_found(self, capsys):
        """Test when log file doesn't exist."""
        totest.read_batch_log_file("/nonexistent/file.out", 0, "1e+00", 12345)
        captured = capsys.readouterr()
        assert "Log file not found" in captured.out

    def test_read_batch_log_file_empty(self, tmp_path, capsys):
        """Test with empty log file."""
        log_file = tmp_path / "test.out"
        log_file.write_text("")

        totest.read_batch_log_file(str(log_file), 0, "1e+00", 12345)
        captured = capsys.readouterr()
        assert "empty file" in captured.out

    def test_read_batch_log_file_time_limit(self, tmp_path, capsys):
        """Test detecting time limit exceeded with exactly 3 lines."""
        log_file = tmp_path / "test.out"
        log_file.write_text("Line 1\nLine 2\nDUE TO TIME LIMIT\n")

        totest.read_batch_log_file(str(log_file), 5, "1e+00", 12345)
        captured = capsys.readouterr()
        assert "Wall time exceeded" in captured.out

    def test_read_batch_log_file_time_limit_many_lines(self, tmp_path, capsys):
        """Test detecting time limit exceeded with many lines (>3 lines, lines 471-476)."""
        log_file = tmp_path / "test.out"
        # Create a log file with many lines where the last 3 contain DUE TO TIME LIMIT
        content = "\n".join([f"Log line {i}" for i in range(20)])
        content += "\nSome error message\nMore context\nDUE TO TIME LIMIT\n"
        log_file.write_text(content)

        totest.read_batch_log_file(str(log_file), 5, "1e+00", 12345)
        captured = capsys.readouterr()
        # Should detect time limit in the last 3 lines
        assert "Wall time exceeded" in captured.out
        # Should NOT print individual lines when time limit is detected
        assert "1e+00_logs/popsyn_12345_5.out:" not in captured.out

    def test_read_batch_log_file_normal(self, tmp_path, capsys):
        """Test normal log file reading without time limit."""
        log_file = tmp_path / "test.out"
        log_file.write_text("Line 1\nLine 2\nLine 3\n")

        totest.read_batch_log_file(str(log_file), 7, "1e+00", 12345)
        captured = capsys.readouterr()
        # Should print the log file path and last 3 lines
        assert "1e+00_logs/popsyn_12345_7.out:" in captured.out
        assert "Line 1" in captured.out or "Line 2" in captured.out or "Line 3" in captured.out

    def test_read_batch_log_file_few_lines(self, tmp_path, capsys):
        """Test log file with fewer than 3 lines (lines 478-487)."""
        log_file = tmp_path / "test.out"
        log_file.write_text("Only one line\nAnother line\n")

        totest.read_batch_log_file(str(log_file), 3, "1e+00", 12345)
        captured = capsys.readouterr()
        assert "File contains 2 line(s):" in captured.out
        assert "Only one line" in captured.out
        assert "1e+00_logs/popsyn_12345_3.out:" in captured.out

    def test_read_batch_log_file_few_lines_time_limit(self, tmp_path, capsys):
        """Test log file with <3 lines containing time limit (lines 481-482)."""
        log_file = tmp_path / "test.out"
        log_file.write_text("DUE TO TIME LIMIT\n")

        totest.read_batch_log_file(str(log_file), 3, "1e+00", 12345)
        captured = capsys.readouterr()
        assert "File contains 1 line(s):" in captured.out
        assert "Wall time exceeded" in captured.out


class TestAnalyzeMissingBatchLogs:
    """Test class for analyze_missing_batch_logs function."""

    def test_analyze_missing_batch_logs_no_missing(self, capsys):
        """Test when there are no missing indices."""
        # Should return early without doing anything
        totest.analyze_missing_batch_logs("/test/folder", "1e+00", set())
        captured = capsys.readouterr()
        # Should not print anything about reading logs
        assert "READING THE LOGS" not in captured.out

    @patch('posydon.CLI.popsyn.check.select_job_id')
    @patch('posydon.CLI.popsyn.check.read_batch_log_file')
    def test_analyze_missing_batch_logs_with_job_id(
        self, mock_read_log, mock_select_job, tmp_path, capsys
    ):
        """Test analyzing missing batch logs with valid job ID."""
        missing_indices = {1, 3, 5}
        mock_select_job.return_value = 12345

        totest.analyze_missing_batch_logs(str(tmp_path), "1e+00", missing_indices)

        # Should call select_job_id
        mock_select_job.assert_called_once_with(str(tmp_path), "1e+00")

        # Should call read_batch_log_file for each missing index
        assert mock_read_log.call_count == len(missing_indices)

        captured = capsys.readouterr()
        assert "READING THE LOGS" in captured.out

    @patch('posydon.CLI.popsyn.check.select_job_id')
    def test_analyze_missing_batch_logs_no_job_id(
        self, mock_select_job, tmp_path, capsys
    ):
        """Test when no job ID is found."""
        missing_indices = {1, 2, 3}
        mock_select_job.return_value = None  # No job ID found

        totest.analyze_missing_batch_logs(str(tmp_path), "1e+00", missing_indices)

        # Should return early after select_job_id returns None
        captured = capsys.readouterr()
        # Should not try to read logs
        assert "READING THE LOGS" not in captured.out


class TestCheckBatch:
    """Test class for check_batch function."""

    def test_check_batch_folder_missing(self, tmp_path, capsys):
        """Test when batch folder doesn't exist."""
        result = totest.check_batch(str(tmp_path), 1.0, "batches")

        assert result['status'] == "folder_missing"
        assert result['found_count'] == 0
        assert result['metallicity'] == 1.0

    @patch('posydon.CLI.popsyn.check.get_expected_batch_count')
    @patch('posydon.CLI.popsyn.check.analyze_missing_batch_logs')
    def test_check_batch_incomplete(
        self, mock_analyze, mock_get_expected, tmp_path, capsys
    ):
        """Test when batch is incomplete."""
        # Create batch folder with some files
        metallicity = 1.0
        str_met = convert_metallicity_to_string(metallicity)
        batch_folder = tmp_path / f"{str_met}_Zsun_batches"
        batch_folder.mkdir()

        for idx in [0, 1, 2]:
            (batch_folder / f"evolution.combined.{idx}.h5").write_text("dummy")

        mock_get_expected.return_value = 5

        result = totest.check_batch(str(tmp_path), metallicity, "batches")

        assert result['status'] == "incomplete"
        assert result['found_count'] == 3
        assert result['expected_count'] == 5
        assert result['missing_indices'] == {3, 4}

    @patch('posydon.CLI.popsyn.check.get_expected_batch_count')
    @patch('posydon.CLI.popsyn.check.analyze_missing_batch_logs')
    def test_check_batch_complete(
        self, mock_analyze, mock_get_expected, tmp_path, capsys
    ):
        """Test when batch is complete (line 581->585 branch: condition is False)."""
        # Create batch folder with all expected files
        metallicity = 1.0
        str_met = convert_metallicity_to_string(metallicity)
        batch_folder = tmp_path / f"{str_met}_Zsun_batches"
        batch_folder.mkdir()

        # Create all 5 expected batch files
        for idx in [0, 1, 2, 3, 4]:
            (batch_folder / f"evolution.combined.{idx}.h5").write_text("dummy")

        mock_get_expected.return_value = 5

        result = totest.check_batch(str(tmp_path), metallicity, "batches")

        assert result['status'] == "complete"
        assert result['found_count'] == 5
        assert result['expected_count'] == 5
        assert result['missing_indices'] is None  # None because status is complete

    @patch('posydon.CLI.popsyn.check.get_expected_batch_count')
    @patch('posydon.CLI.popsyn.check.analyze_missing_batch_logs')
    def test_check_batch_expected_count_none(
        self, mock_analyze, mock_get_expected, tmp_path, capsys
    ):
        """Test when expected_batch_count is None (line 581->585 branch)."""
        # Create batch folder with some files
        metallicity = 1.0
        str_met = convert_metallicity_to_string(metallicity)
        batch_folder = tmp_path / f"{str_met}_Zsun_batches"
        batch_folder.mkdir()

        for idx in [0, 1, 2]:
            (batch_folder / f"evolution.combined.{idx}.h5").write_text("dummy")

        mock_get_expected.return_value = None  # Expected count is None

        result = totest.check_batch(str(tmp_path), metallicity, "batches")

        assert result['status'] == "unknown_expected_count"
        assert result['found_count'] == 3
        assert result['expected_count'] is None
        assert result['missing_indices'] is None  # None because expected count unknown


class TestGetBatchesStatus:
    """Test class for get_batches_status function."""

    @patch('posydon.CLI.popsyn.check.check_batch')
    def test_get_batches_status_no_missing_files(self, mock_check_batch):
        """Test when there are no missing files."""
        missing_files = {}
        synpop_params = {'temp_directory': 'batches'}

        result = totest.get_batches_status("/test/folder", missing_files, synpop_params)

        # Should not call check_batch if no missing files
        mock_check_batch.assert_not_called()
        assert result == {}

    @patch('posydon.CLI.popsyn.check.check_batch')
    def test_get_batches_status_with_missing_files(
        self, mock_check_batch, capsys
    ):
        """Test with multiple missing files."""
        missing_files = {0.01: False, 1.0: False}
        synpop_params = {'temp_directory': 'batches'}

        # Mock check_batch to return different statuses
        def mock_check_side_effect(run_folder, met, batch_dir):
            return {
                'metallicity': met,
                'status': 'incomplete' if met == 0.01 else 'complete',
                'missing_indices': {1, 2} if met == 0.01 else None
            }

        mock_check_batch.side_effect = mock_check_side_effect

        result = totest.get_batches_status("/test/folder", missing_files, synpop_params)

        # Should call check_batch for each missing file
        assert mock_check_batch.call_count == 2
        assert 0.01 in result
        assert 1.0 in result
        assert result[0.01]['status'] == 'incomplete'
        assert result[1.0]['status'] == 'complete'

    @patch('posydon.CLI.popsyn.check.check_batch')
    def test_get_batches_status_uses_temp_directory(
        self, mock_check_batch, capsys
    ):
        """Test that temp_directory from synpop_params is used."""
        missing_files = {1.0: False}
        synpop_params = {'temp_directory': 'custom_batches'}

        mock_check_batch.return_value = {
            'metallicity': 1.0,
            'status': 'complete'
        }

        totest.get_batches_status("/test/folder", missing_files, synpop_params)

        # Verify check_batch was called with the custom temp_directory
        mock_check_batch.assert_called_once_with(
            "/test/folder",
            1.0,
            'custom_batches'
        )


class TestGetUserConfirmation:
    """Test class for get_user_confirmation function."""

    def test_get_user_confirmation_yes(self, monkeypatch):
        """Test user confirms with 'yes'."""
        monkeypatch.setattr('builtins.input', lambda _: "yes")
        result = totest.get_user_confirmation("Confirm?")
        assert result is True

    def test_get_user_confirmation_y(self, monkeypatch):
        """Test user confirms with 'y'."""
        monkeypatch.setattr('builtins.input', lambda _: "y")
        result = totest.get_user_confirmation("Confirm?")
        assert result is True

    def test_get_user_confirmation_no(self, monkeypatch):
        """Test user declines with 'no'."""
        monkeypatch.setattr('builtins.input', lambda _: "no")
        result = totest.get_user_confirmation("Confirm?")
        assert result is False

    def test_get_user_confirmation_invalid(self, monkeypatch, capsys):
        """Test invalid input defaults to no."""
        monkeypatch.setattr('builtins.input', lambda _: "invalid")
        result = totest.get_user_confirmation("Confirm?")
        assert result is False
        captured = capsys.readouterr()
        assert "Unrecognized input" in captured.out

    def test_get_user_confirmation_custom_valid_yes(self, monkeypatch):
        """Test with custom valid_yes parameter (line 651->653 branch: valid_yes is not None)."""
        monkeypatch.setattr('builtins.input', lambda _: "yep")
        result = totest.get_user_confirmation("Confirm?", valid_yes=['yep', 'sure'])
        assert result is True

    def test_get_user_confirmation_custom_valid_no(self, monkeypatch):
        """Test with custom valid_no parameter (line 653-> branch: valid_no is not None)."""
        monkeypatch.setattr('builtins.input', lambda _: "nope")
        result = totest.get_user_confirmation("Confirm?", valid_no=['nope', 'never'])
        assert result is False

    def test_get_user_confirmation_both_custom(self, monkeypatch):
        """Test with both custom valid_yes and valid_no parameters."""
        monkeypatch.setattr('builtins.input', lambda _: "affirmative")
        result = totest.get_user_confirmation(
            "Confirm?",
            valid_yes=['affirmative', 'ok'],
            valid_no=['negative', 'cancel']
        )
        assert result is True


class TestSubmitSlurmJob:
    """Test class for submit_slurm_job function."""

    @patch('subprocess.run')
    def test_submit_slurm_job_success(self, mock_run, capsys):
        """Test successful SLURM job submission."""
        mock_result = MagicMock()
        mock_result.stdout = "Submitted batch job 12345\\n"
        mock_run.return_value = mock_result

        result = totest.submit_slurm_job("test.slurm")

        assert result is True
        captured = capsys.readouterr()
        assert "Job submitted" in captured.out

    @patch('subprocess.run')
    def test_submit_slurm_job_failure(self, mock_run, capsys):
        """Test failed SLURM job submission."""
        import subprocess
        mock_run.side_effect = subprocess.CalledProcessError(1, 'sbatch', stderr="Error")

        result = totest.submit_slurm_job("test.slurm")

        assert result is False
        captured = capsys.readouterr()
        assert "Failed to submit" in captured.out

    @patch('subprocess.run')
    def test_submit_slurm_job_not_found(self, mock_run, capsys):
        """Test when sbatch command is not found."""
        mock_run.side_effect = FileNotFoundError()

        result = totest.submit_slurm_job("test.slurm")

        assert result is False
        captured = capsys.readouterr()
        assert "sbatch command not found" in captured.out


class TestHandleBatchesComplete:
    """Test class for handle_batches_complete function."""

    @pytest.fixture
    def mock_args(self, tmp_path):
        """Create mock arguments."""
        args = MagicMock()
        args.run_folder = str(tmp_path)
        return args

    def test_handle_batches_complete_has_incomplete(self, mock_args):
        """Test when there are incomplete batches."""
        missing_files = {1.0: False}
        batch_status = {
            1.0: {'status': 'incomplete'}
        }

        result = totest.handle_batches_complete(mock_args, missing_files, batch_status)
        assert result is False

    @patch('posydon.CLI.popsyn.check.get_user_confirmation')
    def test_handle_batches_complete_user_declines(
        self, mock_confirm, mock_args, capsys
    ):
        """Test when user declines to resubmit merge jobs."""
        missing_files = {1.0: False}
        batch_status = {1.0: {'status': 'complete'}}
        mock_confirm.return_value = False

        result = totest.handle_batches_complete(mock_args, missing_files, batch_status)

        assert result is True
        captured = capsys.readouterr()
        assert "Merge jobs not resubmitted" in captured.out

    @patch('posydon.CLI.popsyn.check.get_user_confirmation')
    @patch('posydon.CLI.popsyn.check.submit_slurm_job')
    def test_handle_batches_complete_user_confirms(
        self, mock_submit, mock_confirm, mock_args, tmp_path, capsys
    ):
        """Test when user confirms to resubmit merge jobs."""
        # Create merge script
        merge_script = tmp_path / "7_Zsun_merge_popsyn.slurm"
        merge_script.write_text("#!/bin/bash\\n")

        missing_files = {1.0: False}
        batch_status = {1.0: {'status': 'complete'}}
        mock_confirm.return_value = True
        mock_submit.return_value = True

        result = totest.handle_batches_complete(mock_args, missing_files, batch_status)

        assert result is True
        mock_submit.assert_called_once()


class TestCheckPopsynFunction:
    """Test class for check_popsyn_function main function."""

    @pytest.fixture
    def mock_args(self, tmp_path):
        """Create mock arguments."""
        args = MagicMock()
        args.run_folder = str(tmp_path)
        args.job_array = 10
        return args

    @patch('posydon.CLI.popsyn.check.get_run_configuration')
    @patch('posydon.CLI.popsyn.check.check_run_status')
    def test_check_popsyn_function_all_passed(
        self, mock_check_status, mock_get_config, mock_args, capsys
    ):
        """Test when all checks pass."""
        # Setup mocks
        mock_get_config.return_value = (
            "test.ini", 1, 1000, [1.0], {'test': 'params'}
        )
        mock_check_status.return_value = (True, True, {1.0: True})

        result = totest.check_popsyn_function(mock_args)

        assert result == 0
        captured = capsys.readouterr()
        assert "All checks passed successfully" in captured.out

    @patch('posydon.CLI.popsyn.check.get_run_configuration')
    @patch('posydon.CLI.popsyn.check.check_run_status')
    @patch('posydon.CLI.popsyn.check.get_batches_status')
    @patch('posydon.CLI.popsyn.check.handle_batches_complete')
    def test_check_popsyn_function_needs_merge(
        self, mock_handle, mock_get_batches, mock_check_status,
        mock_get_config, mock_args
    ):
        """Test when only merge jobs need to be resubmitted."""
        # Setup mocks
        mock_get_config.return_value = (
            "test.ini", 1, 1000, [1.0], {'test': 'params'}
        )
        mock_check_status.return_value = (False, False, {1.0: False})
        mock_get_batches.return_value = {1.0: {'status': 'complete'}}
        mock_handle.return_value = True  # Merge handled

        result = totest.check_popsyn_function(mock_args)

        assert result == 2

    @patch('posydon.CLI.popsyn.check.get_run_configuration')
    @patch('posydon.CLI.popsyn.check.check_run_status')
    @patch('posydon.CLI.popsyn.check.get_batches_status')
    @patch('posydon.CLI.popsyn.check.handle_batches_complete')
    def test_check_popsyn_function_batch_folder_missing(
        self, mock_handle, mock_get_batches, mock_check_status,
        mock_get_config, mock_args, capsys
    ):
        """Test when batch folder is missing (return 1)."""
        # Setup mocks
        mock_get_config.return_value = (
            "test.ini", 1, 1000, [1.0], {'test': 'params'}
        )
        mock_check_status.return_value = (False, False, {1.0: False})
        mock_get_batches.return_value = {1.0: {'status': 'folder_missing'}}
        mock_handle.return_value = False  # Batches incomplete

        result = totest.check_popsyn_function(mock_args)

        assert result == 1
        captured = capsys.readouterr()
        assert "One or more batch folders are missing" in captured.out
        assert "Cannot generate rescue scripts" in captured.out

    @patch('posydon.CLI.popsyn.check.get_run_configuration')
    @patch('posydon.CLI.popsyn.check.check_run_status')
    @patch('posydon.CLI.popsyn.check.get_batches_status')
    @patch('posydon.CLI.popsyn.check.handle_batches_complete')
    @patch('posydon.CLI.popsyn.check.get_user_confirmation')
    def test_check_popsyn_function_incomplete_user_declines_rescue(
        self, mock_confirm, mock_handle, mock_get_batches,
        mock_check_status, mock_get_config, mock_args, capsys
    ):
        """Test when batches incomplete but user declines to create rescue scripts."""
        # Setup mocks
        mock_get_config.return_value = (
            "test.ini", 1, 1000, [1.0], {'test': 'params'}
        )
        mock_check_status.return_value = (False, False, {1.0: False})
        mock_get_batches.return_value = {
            1.0: {
                'status': 'incomplete',
                'missing_indices': {1, 2, 3}
            }
        }
        mock_handle.return_value = False  # Batches incomplete
        mock_confirm.return_value = False  # User declines

        result = totest.check_popsyn_function(mock_args)

        assert result == 2
        captured = capsys.readouterr()
        assert "Rescue scripts not created" in captured.out

    @patch('posydon.CLI.popsyn.check.get_run_configuration')
    @patch('posydon.CLI.popsyn.check.check_run_status')
    @patch('posydon.CLI.popsyn.check.get_batches_status')
    @patch('posydon.CLI.popsyn.check.handle_batches_complete')
    @patch('posydon.CLI.popsyn.check.get_user_confirmation')
    @patch('posydon.CLI.popsyn.check.create_batch_rescue_script')
    @patch('posydon.CLI.popsyn.check.create_bash_submit_rescue_script')
    @patch('os.system')
    def test_check_popsyn_function_rescue_created_and_submitted(
        self, mock_os_system, mock_create_bash, mock_create_rescue,
        mock_confirm, mock_handle, mock_get_batches,
        mock_check_status, mock_get_config, mock_args, capsys
    ):
        """Test when rescue scripts are created and submitted."""
        # Setup mocks
        mock_get_config.return_value = (
            "test.ini", 2, 1000, [0.01, 1.0], {'test': 'params'}
        )
        mock_check_status.return_value = (False, False, {0.01: False, 1.0: False})
        mock_get_batches.return_value = {
            0.01: {
                'status': 'incomplete',
                'missing_indices': {1, 2}
            },
            1.0: {
                'status': 'incomplete',
                'missing_indices': {3, 4}
            }
        }
        mock_handle.return_value = False  # Batches incomplete
        # User confirms both times (create rescue and submit)
        mock_confirm.side_effect = [True, True]
        mock_create_rescue.side_effect = ['rescue1.slurm', 'rescue2.slurm']
        mock_create_bash.return_value = '/test/folder/resubmit_slurm.sh'

        result = totest.check_popsyn_function(mock_args)

        assert result == 0
        captured = capsys.readouterr()
        assert "GENERATING A RESCUE SCRIPT" in captured.out
        assert "Rescue scripts submitted" in captured.out

        # Verify rescue scripts were created for both metallicities
        assert mock_create_rescue.call_count == 2
        mock_os_system.assert_called_once_with('sh /test/folder/resubmit_slurm.sh')

    @patch('posydon.CLI.popsyn.check.get_run_configuration')
    @patch('posydon.CLI.popsyn.check.check_run_status')
    @patch('posydon.CLI.popsyn.check.get_batches_status')
    @patch('posydon.CLI.popsyn.check.handle_batches_complete')
    @patch('posydon.CLI.popsyn.check.get_user_confirmation')
    @patch('posydon.CLI.popsyn.check.create_batch_rescue_script')
    @patch('posydon.CLI.popsyn.check.create_bash_submit_rescue_script')
    def test_check_popsyn_function_rescue_created_not_submitted(
        self, mock_create_bash, mock_create_rescue,
        mock_confirm, mock_handle, mock_get_batches,
        mock_check_status, mock_get_config, mock_args, capsys
    ):
        """Test when rescue scripts are created but user declines to submit."""
        # Setup mocks
        mock_get_config.return_value = (
            "test.ini", 1, 1000, [1.0], {'test': 'params'}
        )
        mock_check_status.return_value = (False, False, {1.0: False})
        mock_get_batches.return_value = {
            1.0: {
                'status': 'incomplete',
                'missing_indices': {1, 2, 3}
            }
        }
        mock_handle.return_value = False  # Batches incomplete
        # User confirms creation but declines submission
        mock_confirm.side_effect = [True, False]
        mock_create_rescue.return_value = 'rescue.slurm'
        mock_create_bash.return_value = '/test/folder/resubmit_slurm.sh'

        result = totest.check_popsyn_function(mock_args)

        assert result == 2
        captured = capsys.readouterr()
        assert "Please submit the rescue scripts" in captured.out
        assert "sh resubmit_slurm.sh" in captured.out
