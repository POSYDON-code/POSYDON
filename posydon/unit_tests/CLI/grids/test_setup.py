"""Unit tests of posydon/CLI/grids/setup.py"""

__authors__ = [
    "GitHub Copilot <copilot@github.com>"
]

import os
import subprocess
from unittest.mock import MagicMock, call, patch

import pytest

# import the module which will be tested
import posydon.CLI.grids.setup as totest


class TestSetupMESADefaults:
    """Test class for setup_MESA_defaults function."""

    @pytest.fixture
    def mock_mesa_dir(self, tmp_path):
        """Create a mock MESA directory structure for testing."""
        mesa_dir = tmp_path / "mesa"
        mesa_dir.mkdir()

        # Create binary work directory structure
        binary_make = mesa_dir / "binary" / "work" / "make"
        binary_make.mkdir(parents=True)
        (binary_make / "makefile").write_text("# Binary makefile")

        binary_src = mesa_dir / "binary" / "work" / "src"
        binary_src.mkdir(parents=True)
        (binary_src / "run_binary.f").write_text("! Binary run file")
        (binary_src / "run_binary_extras.f").write_text("! Binary extras")
        (binary_src / "run_star_extras.f").write_text("! Star binary extras")

        # Create star work directory structure
        star_make = mesa_dir / "star" / "work" / "make"
        star_make.mkdir(parents=True)
        (star_make / "makefile").write_text("# Star makefile")

        star_src = mesa_dir / "star" / "work" / "src"
        star_src.mkdir(parents=True)
        (star_src / "run.f").write_text("! Star run file")
        (star_src / "run_star_extras.f").write_text("! Star extras")

        # Create defaults directory with column files
        binary_defaults = mesa_dir / "binary" / "defaults"
        binary_defaults.mkdir(parents=True)
        (binary_defaults / "history_columns.list").write_text("# Binary history columns")

        star_defaults = mesa_dir / "star" / "defaults"
        star_defaults.mkdir(parents=True)
        (star_defaults / "history_columns.list").write_text("# Star history columns")
        (star_defaults / "profile_columns.list").write_text("# Profile columns")

        return str(mesa_dir)

    @pytest.fixture
    def mock_version_path(self, tmp_path):
        """Create a mock version path."""
        version_path = tmp_path / "r11701"
        version_path.mkdir()

        mesa_defaults = version_path / "MESA_defaults"
        mesa_defaults.mkdir()

        return str(version_path)

    def test_setup_MESA_defaults_success(self, mock_mesa_dir, mock_version_path):
        """Test successful setup of MESA defaults."""
        with patch.dict(os.environ, {'MESA_DIR': mock_mesa_dir}):
            inlists, extras, columns = totest.setup_MESA_defaults(mock_version_path)

            # Check that all three dictionaries are returned
            assert isinstance(inlists, dict)
            assert isinstance(extras, dict)
            assert isinstance(columns, dict)

            # Check extras dictionary has all expected keys
            expected_extras_keys = [
                'makefile_binary', 'makefile_star', 'star_run', 'binary_run',
                'binary_extras', 'star_binary_extras', 'star1_extras', 'star2_extras'
            ]
            for key in expected_extras_keys:
                assert key in extras, f"Missing key: {key}"
                assert extras[key], f"Empty path for key: {key}"

            # Check columns dictionary has all expected keys
            expected_column_keys = [
                'star_history_columns', 'binary_history_columns', 'profile_columns'
            ]
            for key in expected_column_keys:
                assert key in columns, f"Missing key: {key}"
                assert columns[key], f"Empty path for key: {key}"

    def test_setup_MESA_defaults_paths_correctness(self, mock_mesa_dir, mock_version_path):
        """Test that returned paths are correctly formatted."""
        with patch.dict(os.environ, {'MESA_DIR': mock_mesa_dir}):
            inlists, extras, columns = totest.setup_MESA_defaults(mock_version_path)

            # Check makefile paths
            assert extras['makefile_binary'].endswith('binary/work/make/makefile')
            assert extras['makefile_star'].endswith('star/work/make/makefile')

            # Check run file paths
            assert extras['star_run'].endswith('star/work/src/run.f')
            assert extras['binary_run'].endswith('binary/work/src/run_binary.f')

            # Check extras file paths
            assert extras['binary_extras'].endswith('binary/work/src/run_binary_extras.f')
            assert extras['star_binary_extras'].endswith('binary/work/src/run_star_extras.f')
            assert extras['star1_extras'].endswith('star/work/src/run_star_extras.f')
            assert extras['star2_extras'].endswith('star/work/src/run_star_extras.f')

            # Check that star1_extras and star2_extras point to the same file
            assert extras['star1_extras'] == extras['star2_extras']

            # Check column file paths
            assert columns['star_history_columns'].endswith('star/defaults/history_columns.list')
            assert columns['binary_history_columns'].endswith('binary/defaults/history_columns.list')
            assert columns['profile_columns'].endswith('star/defaults/profile_columns.list')

    def test_setup_MESA_defaults_all_paths_contain_mesa_dir(self, mock_mesa_dir, mock_version_path):
        """Test that all returned paths start with MESA_DIR."""
        with patch.dict(os.environ, {'MESA_DIR': mock_mesa_dir}):
            inlists, extras, columns = totest.setup_MESA_defaults(mock_version_path)

            # Check all extras paths
            for key, path in extras.items():
                assert path.startswith(mock_mesa_dir), f"{key} path doesn't start with MESA_DIR"

            # Check all column paths
            for key, path in columns.items():
                assert path.startswith(mock_mesa_dir), f"{key} path doesn't start with MESA_DIR"

    def test_setup_MESA_defaults_missing_mesa_dir_env(self, mock_version_path):
        """Test that function fails when MESA_DIR environment variable is not set."""
        with patch.dict(os.environ, {}, clear=True):
            with pytest.raises(KeyError):
                totest.setup_MESA_defaults(mock_version_path)

    def test_setup_MESA_defaults_missing_extras_file(self, mock_mesa_dir, mock_version_path):
        """Test that function raises ValueError when an extras file is missing."""
        # Remove one of the required files
        binary_run_path = os.path.join(mock_mesa_dir, "binary", "work", "src", "run_binary.f")
        os.remove(binary_run_path)

        with patch.dict(os.environ, {'MESA_DIR': mock_mesa_dir}):
            with pytest.raises(ValueError, match="does not exist"):
                totest.setup_MESA_defaults(mock_version_path)

    def test_setup_MESA_defaults_missing_column_file(self, mock_mesa_dir, mock_version_path):
        """Test that function raises ValueError when a column file is missing."""
        # Remove one of the required column files
        profile_columns_path = os.path.join(mock_mesa_dir, "star", "defaults", "profile_columns.list")
        os.remove(profile_columns_path)

        with patch.dict(os.environ, {'MESA_DIR': mock_mesa_dir}):
            with pytest.raises(ValueError, match="does not exist"):
                totest.setup_MESA_defaults(mock_version_path)

    def test_setup_MESA_defaults_empty_inlists_dict(self, mock_mesa_dir, mock_version_path):
        """Test that inlists dictionary is empty (as per current implementation)."""
        with patch.dict(os.environ, {'MESA_DIR': mock_mesa_dir}):
            inlists, extras, columns = totest.setup_MESA_defaults(mock_version_path)

            # Current implementation returns empty inlists dict
            assert inlists == {}

    def test_setup_MESA_defaults_file_existence_check(self, mock_mesa_dir, mock_version_path, capsys):
        """Test that files are checked for existence and errors are printed."""
        # Create a scenario with a missing file
        binary_extras_path = os.path.join(mock_mesa_dir, "binary", "work", "src", "run_binary_extras.f")
        os.remove(binary_extras_path)

        with patch.dict(os.environ, {'MESA_DIR': mock_mesa_dir}):
            with pytest.raises(ValueError):
                totest.setup_MESA_defaults(mock_version_path)

            # Check that error message was printed
            captured = capsys.readouterr()
            assert "does not exist" in captured.out

    def test_setup_MESA_defaults_mesa_path_helper(self, mock_mesa_dir, mock_version_path):
        """Test that the internal mesa_path helper function works correctly."""
        with patch.dict(os.environ, {'MESA_DIR': mock_mesa_dir}):
            inlists, extras, columns = totest.setup_MESA_defaults(mock_version_path)

            # Verify that paths are constructed using proper path joining
            # This indirectly tests the mesa_path helper function
            expected_binary_run = os.path.join(mock_mesa_dir, 'binary', 'work', 'src', 'run_binary.f')
            assert extras['binary_run'] == expected_binary_run

    def test_setup_MESA_defaults_with_special_characters_in_path(self, tmp_path):
        """Test that function handles paths with special characters."""
        # Create a MESA directory with spaces in the path
        mesa_dir = tmp_path / "mesa test dir"
        mesa_dir.mkdir()

        # Create minimal required structure
        for module in ['binary', 'star']:
            work_dirs = ['make', 'src']
            for work_dir in work_dirs:
                (mesa_dir / module / "work" / work_dir).mkdir(parents=True)
            (mesa_dir / module / "defaults").mkdir(parents=True)

        # Create all required files
        (mesa_dir / "binary" / "work" / "make" / "makefile").touch()
        (mesa_dir / "star" / "work" / "make" / "makefile").touch()
        (mesa_dir / "star" / "work" / "src" / "run.f").touch()
        (mesa_dir / "binary" / "work" / "src" / "run_binary.f").touch()
        (mesa_dir / "binary" / "work" / "src" / "run_binary_extras.f").touch()
        (mesa_dir / "binary" / "work" / "src" / "run_star_extras.f").touch()
        (mesa_dir / "star" / "work" / "src" / "run_star_extras.f").touch()
        (mesa_dir / "binary" / "defaults" / "history_columns.list").touch()
        (mesa_dir / "star" / "defaults" / "history_columns.list").touch()
        (mesa_dir / "star" / "defaults" / "profile_columns.list").touch()

        version_path = tmp_path / "version dir"
        version_path.mkdir()

        with patch.dict(os.environ, {'MESA_DIR': str(mesa_dir)}):
            inlists, extras, columns = totest.setup_MESA_defaults(str(version_path))

            # Should succeed without errors
            assert len(extras) > 0
            assert len(columns) > 0


class TestCheckFileExist:
    """Test class for check_file_exist function."""

    def test_check_file_exist_with_existing_file(self, tmp_path):
        """Test that function returns True for existing file."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("test content")

        result = totest.check_file_exist(str(test_file))
        assert result is True

    def test_check_file_exist_with_missing_file_raise_error(self, tmp_path):
        """Test that function raises ValueError for missing file when raise_error=True."""
        missing_file = tmp_path / "missing.txt"

        with pytest.raises(ValueError, match="does not exist"):
            totest.check_file_exist(str(missing_file), raise_error=True)

    def test_check_file_exist_with_missing_file_no_raise(self, tmp_path):
        """Test that function returns False for missing file when raise_error=False."""
        missing_file = tmp_path / "missing.txt"

        result = totest.check_file_exist(str(missing_file), raise_error=False)
        assert result is False

    def test_check_file_exist_prints_error_message(self, tmp_path, capsys):
        """Test that error message is printed when file doesn't exist."""
        missing_file = tmp_path / "missing.txt"

        with pytest.raises(ValueError):
            totest.check_file_exist(str(missing_file))

        captured = capsys.readouterr()
        assert "does not exist" in captured.out
        assert str(missing_file) in captured.out

    def test_check_file_exist_with_directory(self, tmp_path):
        """Test that function works with directories (os.path.exists returns True for dirs)."""
        test_dir = tmp_path / "test_dir"
        test_dir.mkdir()

        result = totest.check_file_exist(str(test_dir))
        assert result is True


class TestSetupInlistRepository:
    """Test class for setup_inlist_repository function."""

    @pytest.fixture
    def mock_subprocess_success(self):
        """Mock subprocess.run to return successful results."""
        mock_result = MagicMock()
        mock_result.stdout = "Success"
        mock_result.returncode = 0
        return mock_result

    def test_setup_inlist_repository_new_directory(self, tmp_path, mock_subprocess_success, capsys):
        """Test setup when inlist repository directory doesn't exist."""
        inlist_repo = tmp_path / "new_inlist_repo"
        mesa_version = "r11701"

        # Create the version folder after the repo would be cloned
        version_path = inlist_repo / mesa_version

        with patch('subprocess.run', return_value=mock_subprocess_success) as mock_run:
            # Create the directory structure as git clone would
            def side_effect(*args, **kwargs):
                if 'clone' in args[0]:
                    inlist_repo.mkdir(parents=True, exist_ok=True)
                    version_path.mkdir(parents=True, exist_ok=True)
                return mock_subprocess_success

            mock_run.side_effect = side_effect

            result = totest.setup_inlist_repository(str(inlist_repo), mesa_version)

            # Check that directory was created
            assert inlist_repo.exists()

            # Check that git clone was called
            assert any('clone' in str(call_args) for call_args in mock_run.call_args_list)

            # Check that git pull was called
            assert any('pull' in str(call_args) for call_args in mock_run.call_args_list)

            # Check return value
            assert result == str(version_path)

            # Check printed messages
            captured = capsys.readouterr()
            assert "We are setting up your inlist repository now" in captured.out
            assert "Creating inlist repository" in captured.out

    def test_setup_inlist_repository_existing_empty_directory(self, tmp_path, mock_subprocess_success, capsys):
        """Test setup when directory exists but is empty."""
        inlist_repo = tmp_path / "empty_repo"
        inlist_repo.mkdir()
        mesa_version = "r11701"
        version_path = inlist_repo / mesa_version

        with patch('subprocess.run', return_value=mock_subprocess_success) as mock_run:
            # Create version folder during clone
            def side_effect(*args, **kwargs):
                if 'clone' in args[0]:
                    version_path.mkdir(parents=True, exist_ok=True)
                return mock_subprocess_success

            mock_run.side_effect = side_effect

            result = totest.setup_inlist_repository(str(inlist_repo), mesa_version)

            # Check that git clone was called
            clone_calls = [c for c in mock_run.call_args_list if 'clone' in str(c)]
            assert len(clone_calls) == 1

            # Verify the clone command
            clone_call = clone_calls[0]
            command_list = clone_call[0][0]
            assert 'git' in command_list
            assert 'clone' in command_list
            # Check that the URL is in the command list
            assert any('POSYDON-MESA-INLISTS.git' in arg for arg in command_list)

            assert result == str(version_path)

    def test_setup_inlist_repository_existing_nonempty_directory(self, tmp_path, mock_subprocess_success, capsys):
        """Test setup when directory exists and contains files (repo already cloned)."""
        inlist_repo = tmp_path / "existing_repo"
        inlist_repo.mkdir()

        # Create some files to simulate existing repo
        (inlist_repo / "README.md").write_text("# POSYDON MESA INLISTS")

        mesa_version = "r11701"
        version_path = inlist_repo / mesa_version
        version_path.mkdir(parents=True)

        with patch('subprocess.run', return_value=mock_subprocess_success) as mock_run:
            result = totest.setup_inlist_repository(str(inlist_repo), mesa_version)

            # Check that git clone was NOT called
            clone_calls = [c for c in mock_run.call_args_list if 'clone' in str(c)]
            assert len(clone_calls) == 0

            # Check that git pull was still called
            pull_calls = [c for c in mock_run.call_args_list if 'pull' in str(c)]
            assert len(pull_calls) == 1

            # Check printed messages
            captured = capsys.readouterr()
            assert "Files found in inlist repository" in captured.out
            assert "assuming POSYDON inlist repository is already cloned" in captured.out

            assert result == str(version_path)

    def test_setup_inlist_repository_git_pull_with_correct_cwd(self, tmp_path, mock_subprocess_success):
        """Test that git pull is called with correct working directory."""
        inlist_repo = tmp_path / "repo"
        inlist_repo.mkdir()
        (inlist_repo / "README.md").write_text("test")

        mesa_version = "r11701"
        version_path = inlist_repo / mesa_version
        version_path.mkdir(parents=True)

        with patch('subprocess.run', return_value=mock_subprocess_success) as mock_run:
            totest.setup_inlist_repository(str(inlist_repo), mesa_version)

            # Find the git pull call
            pull_calls = [c for c in mock_run.call_args_list if 'pull' in str(c)]
            assert len(pull_calls) == 1

            # Check that cwd parameter was set correctly
            pull_call = pull_calls[0]
            assert pull_call[1].get('cwd') == str(inlist_repo)

    def test_setup_inlist_repository_missing_version_folder(self, tmp_path, mock_subprocess_success):
        """Test that ValueError is raised when MESA version folder doesn't exist."""
        inlist_repo = tmp_path / "repo"
        inlist_repo.mkdir()
        (inlist_repo / "README.md").write_text("test")

        mesa_version = "r99999"  # Non-existent version

        with patch('subprocess.run', return_value=mock_subprocess_success):
            with pytest.raises(ValueError, match="does not exist in the inlist repository"):
                totest.setup_inlist_repository(str(inlist_repo), mesa_version)

    def test_setup_inlist_repository_git_clone_failure(self, tmp_path):
        """Test handling of git clone failure."""
        inlist_repo = tmp_path / "repo"
        mesa_version = "r11701"

        # Mock subprocess to raise CalledProcessError for git clone
        with patch('subprocess.run') as mock_run:
            mock_run.side_effect = subprocess.CalledProcessError(1, 'git clone')

            with pytest.raises(subprocess.CalledProcessError):
                totest.setup_inlist_repository(str(inlist_repo), mesa_version)

    def test_setup_inlist_repository_git_pull_failure(self, tmp_path):
        """Test handling of git pull failure."""
        inlist_repo = tmp_path / "repo"
        inlist_repo.mkdir()
        (inlist_repo / "README.md").write_text("test")

        mesa_version = "r11701"
        version_path = inlist_repo / mesa_version
        version_path.mkdir(parents=True)

        # Mock subprocess to fail on git pull
        with patch('subprocess.run') as mock_run:
            mock_run.side_effect = subprocess.CalledProcessError(1, 'git pull')

            with pytest.raises(subprocess.CalledProcessError):
                totest.setup_inlist_repository(str(inlist_repo), mesa_version)

    def test_setup_inlist_repository_prints_git_output(self, tmp_path, capsys):
        """Test that git command output is printed."""
        inlist_repo = tmp_path / "repo"
        inlist_repo.mkdir()
        (inlist_repo / "README.md").write_text("test")

        mesa_version = "r11701"
        version_path = inlist_repo / mesa_version
        version_path.mkdir(parents=True)

        mock_result = MagicMock()
        mock_result.stdout = "Already up to date."
        mock_result.returncode = 0

        with patch('subprocess.run', return_value=mock_result):
            totest.setup_inlist_repository(str(inlist_repo), mesa_version)

            captured = capsys.readouterr()
            assert "Already up to date." in captured.out
            assert "Updating inlist repository" in captured.out

    def test_setup_inlist_repository_subprocess_parameters(self, tmp_path, mock_subprocess_success):
        """Test that subprocess.run is called with correct parameters."""
        inlist_repo = tmp_path / "repo"
        mesa_version = "r11701"
        version_path = inlist_repo / mesa_version

        with patch('subprocess.run', return_value=mock_subprocess_success) as mock_run:
            # Create directories during clone simulation
            def side_effect(*args, **kwargs):
                if 'clone' in args[0]:
                    inlist_repo.mkdir(parents=True, exist_ok=True)
                    version_path.mkdir(parents=True, exist_ok=True)
                return mock_subprocess_success

            mock_run.side_effect = side_effect

            totest.setup_inlist_repository(str(inlist_repo), mesa_version)

            # Check all subprocess calls have correct parameters
            for call_obj in mock_run.call_args_list:
                # Check that capture_output, text, and check parameters are set
                assert call_obj[1].get('capture_output') is True
                assert call_obj[1].get('text') is True
                assert call_obj[1].get('check') is True

    def test_setup_inlist_repository_correct_git_url(self, tmp_path, mock_subprocess_success):
        """Test that correct POSYDON MESA INLISTS URL is used."""
        inlist_repo = tmp_path / "repo"
        mesa_version = "r11701"
        version_path = inlist_repo / mesa_version

        with patch('subprocess.run', return_value=mock_subprocess_success) as mock_run:
            # Create directories during clone simulation
            def side_effect(*args, **kwargs):
                if 'clone' in args[0]:
                    inlist_repo.mkdir(parents=True, exist_ok=True)
                    version_path.mkdir(parents=True, exist_ok=True)
                return mock_subprocess_success

            mock_run.side_effect = side_effect

            totest.setup_inlist_repository(str(inlist_repo), mesa_version)

            # Find the clone call and verify URL
            clone_calls = [c for c in mock_run.call_args_list if 'clone' in str(c)]
            if clone_calls:
                clone_call = clone_calls[0]
                command = clone_call[0][0]
                assert 'https://github.com/POSYDON-code/POSYDON-MESA-INLISTS.git' in command

    def test_setup_inlist_repository_return_value_format(self, tmp_path, mock_subprocess_success):
        """Test that return value is correctly formatted path."""
        inlist_repo = tmp_path / "repo"
        inlist_repo.mkdir()
        (inlist_repo / "README.md").write_text("test")

        mesa_version = "r11701"
        version_path = inlist_repo / mesa_version
        version_path.mkdir(parents=True)

        with patch('subprocess.run', return_value=mock_subprocess_success):
            result = totest.setup_inlist_repository(str(inlist_repo), mesa_version)

            # Check that result is a string path
            assert isinstance(result, str)

            # Check that it ends with the MESA version
            assert result.endswith(mesa_version)

            # Check that it contains the inlist_repo path
            assert str(inlist_repo) in result

            # Check that the path exists
            assert os.path.exists(result)
