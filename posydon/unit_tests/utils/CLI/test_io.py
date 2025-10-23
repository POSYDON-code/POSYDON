"""Unit tests of posydon/utils/CLI/io.py"""

__authors__ = [
    "Max Briel <max.briel@gmail.com>",
]

import os
import textwrap
from unittest.mock import MagicMock, patch

import pytest

# import the module which will be tested
import posydon.utils.CLI.io as totest
from posydon.utils.common_functions import convert_metallicity_to_string


class TestOutputUtilityFunctions:
    """Test class for output utility functions."""

    def test_print_error(self, capsys):
        """Test print_error outputs red colored text."""
        test_message = "This is an error"
        totest.print_error(test_message)
        captured = capsys.readouterr()
        assert test_message in captured.out
        assert totest.COLOR_RED in captured.out
        assert totest.COLOR_RESET in captured.out

    def test_print_success(self, capsys):
        """Test print_success outputs green colored text."""
        test_message = "This is a success"
        totest.print_success(test_message)
        captured = capsys.readouterr()
        assert test_message in captured.out
        assert totest.COLOR_GREEN in captured.out
        assert totest.COLOR_RESET in captured.out

    def test_print_separator_line(self, capsys):
        """Test print_separator_line outputs dashes."""
        totest.print_separator_line()
        captured = capsys.readouterr()
        assert "-" * 80 in captured.out

    def test_print_separator_section(self, capsys):
        """Test print_separator_section outputs hashes."""
        totest.print_separator_section()
        captured = capsys.readouterr()
        assert "#" * 80 in captured.out

    def test_clear_previous_lines(self, capsys):
        """Test clear_previous_lines outputs ANSI codes."""
        num_lines = 3
        totest.clear_previous_lines(num_lines)
        captured = capsys.readouterr()
        # Check that ANSI escape codes are printed
        assert '\033[1A' in captured.out or '\033[K' in captured.out


class TestScriptCreation:
    """Test class for script creation functions."""

    def test_create_merge_script_text(self):
        """Test that create_merge_script_text returns valid Python code."""
        ini_file = "/path/to/test.ini"
        text = totest.create_merge_script_text(ini_file)

        # Check that the text contains essential imports and code
        assert "from posydon.popsyn.binarypopulation import BinaryPopulation" in text
        assert "binarypop_kwargs_from_ini" in text
        assert ini_file in text
        assert "synpop.combine_saved_files" in text
        assert "argparse" in text

    def test_create_run_script_text(self):
        """Test that create_run_script_text returns valid Python code."""
        ini_file = "/path/to/test.ini"
        text = totest.create_run_script_text(ini_file)

        # Check that the text contains essential imports and code
        assert "from posydon.popsyn.binarypopulation import BinaryPopulation" in text
        assert "from posydon.binary_evol.simulationproperties import SimulationProperties" in text
        assert ini_file in text
        assert "synpop.evolve()" in text

    def test_create_run_script(self, tmp_path):
        """Test that create_run_script creates a file."""
        # Create a temporary ini file
        ini_file = tmp_path / "test.ini"
        ini_file.write_text("[test]")

        # Change to temp directory to create the script there
        original_dir = os.getcwd()
        os.chdir(tmp_path)

        try:
            totest.create_run_script(str(ini_file))

            # Check that the file was created
            assert os.path.exists("run_metallicity.py")

            # Check file content
            with open("run_metallicity.py", "r") as f:
                content = f.read()
                assert str(ini_file) in content
        finally:
            os.chdir(original_dir)

    def test_create_merge_script(self, tmp_path):
        """Test that create_merge_script creates a file."""
        # Create a temporary ini file
        ini_file = tmp_path / "test.ini"
        ini_file.write_text("[test]")

        # Change to temp directory
        original_dir = os.getcwd()
        os.chdir(tmp_path)

        try:
            totest.create_merge_script(str(ini_file))

            # Check that the file was created
            assert os.path.exists("merge_metallicity.py")

            # Check file content
            with open("merge_metallicity.py", "r") as f:
                content = f.read()
                assert str(ini_file) in content
        finally:
            os.chdir(original_dir)

    def test_create_python_scripts(self, tmp_path, capsys):
        """Test that create_python_scripts creates both run and merge scripts."""
        # Create a temporary ini file
        ini_file = tmp_path / "test.ini"
        ini_file.write_text("[test]")

        # Change to temp directory
        original_dir = os.getcwd()
        os.chdir(tmp_path)

        try:
            totest.create_python_scripts(str(ini_file))

            # Check that both files were created
            assert os.path.exists("run_metallicity.py")
            assert os.path.exists("merge_metallicity.py")

            # Check console output
            captured = capsys.readouterr()
            assert "Created run script and merge script" in captured.out
        finally:
            os.chdir(original_dir)


class TestSlurmScriptCreation:
    """Test class for SLURM script creation functions."""

    def test_create_slurm_array(self, tmp_path):
        """Test that create_slurm_array creates a SLURM array script."""
        original_dir = os.getcwd()
        os.chdir(tmp_path)

        try:
            metallicity = 1.0
            job_array_length = 10
            partition = "normal"
            email = "test@example.com"
            walltime = "24:00:00"
            account = "test_account"
            mem_per_cpu = "4G"
            path_to_posydon = "/path/to/posydon"
            path_to_posydon_data = "/path/to/data"

            totest.create_slurm_array(
                metallicity, job_array_length, partition, email,
                walltime, account, mem_per_cpu,
                path_to_posydon, path_to_posydon_data
            )

            # Check that the file was created
            str_met = convert_metallicity_to_string(metallicity)
            filename = f"{str_met}_Zsun_slurm_array.slurm"
            assert os.path.exists(filename)

            # Check file content
            with open(filename, "r") as f:
                content = f.read()
                assert "#!/bin/bash" in content
                assert f"#SBATCH --array=0-{job_array_length-1}" in content
                assert f"#SBATCH --partition={partition}" in content
                assert f"#SBATCH --mail-user={email}" in content
                assert f"#SBATCH --time={walltime}" in content
                assert f"#SBATCH --mem-per-cpu={mem_per_cpu}" in content
                assert f"#SBATCH --account={account}" in content
                assert f"export PATH_TO_POSYDON={path_to_posydon}" in content
                assert "srun python ./run_metallicity.py" in content
        finally:
            os.chdir(original_dir)

    def test_create_slurm_array_minimal(self, tmp_path):
        """Test create_slurm_array with minimal options (no email, partition, account)."""
        original_dir = os.getcwd()
        os.chdir(tmp_path)

        try:
            metallicity = 1.0
            job_array_length = 5

            totest.create_slurm_array(
                metallicity, job_array_length,
                partition=None, email=None,
                walltime="12:00:00", account=None,
                mem_per_cpu="2G",
                path_to_posydon="/posydon",
                path_to_posydon_data="/data"
            )
            str_met = convert_metallicity_to_string(metallicity)
            filename = f"{str_met}_Zsun_slurm_array.slurm"
            assert os.path.exists(filename)

            with open(filename, "r") as f:
                content = f.read()
                assert "#SBATCH --partition" not in content
                assert "#SBATCH --mail-user" not in content
                assert "#SBATCH --account" not in content
        finally:
            os.chdir(original_dir)

    def test_create_slurm_merge(self, tmp_path):
        """Test that create_slurm_merge creates a merge SLURM script."""
        original_dir = os.getcwd()
        os.chdir(tmp_path)

        try:
            metallicity = 1.0
            partition = "normal"
            email = "test@example.com"
            merge_walltime = "12:00:00"
            account = "test_account"
            mem_per_cpu = "8G"
            path_to_posydon = "/path/to/posydon"
            path_to_posydon_data = "/path/to/data"

            totest.create_slurm_merge(
                metallicity, partition, email,
                merge_walltime, account, mem_per_cpu,
                path_to_posydon, path_to_posydon_data
            )

            # Check that the file was created
            str_met = convert_metallicity_to_string(metallicity)
            filename = f"{str_met}_Zsun_merge_popsyn.slurm"
            assert os.path.exists(filename)

            # Check file content
            with open(filename, "r") as f:
                content = f.read()
                assert "#!/bin/bash" in content
                assert f"#SBATCH --job-name={str_met}_Zsun_merge" in content
                assert f"#SBATCH --time={merge_walltime}" in content
                assert "srun python ./merge_metallicity.py" in content
        finally:
            os.chdir(original_dir)

    def test_create_slurm_merge_debug_partition(self, tmp_path):
        """Test that debug-cpu partition overrides walltime."""
        original_dir = os.getcwd()
        os.chdir(tmp_path)

        try:
            metallicity = 1.0
            totest.create_slurm_merge(
                metallicity=metallicity,
                partition="debug-cpu",
                email=None,
                merge_walltime="24:00:00",  # This should be overridden
                account=None,
                mem_per_cpu="4G",
                path_to_posydon="/posydon",
                path_to_posydon_data="/data"
            )

            str_met = convert_metallicity_to_string(metallicity)
            filename = f"{str_met}_Zsun_merge_popsyn.slurm"
            with open(filename, "r") as f:
                content = f.read()
                # Check that walltime was overridden to 14 minutes
                assert "#SBATCH --time=00:14:00" in content
        finally:
            os.chdir(original_dir)

    def test_create_slurm_rescue(self, tmp_path):
        """Test that create_slurm_rescue creates a rescue script."""
        original_dir = os.getcwd()
        os.chdir(tmp_path)

        try:
            metallicity = 1.0
            missing_indices = [1, 3, 5, 7]
            job_array_length = 10

            totest.create_slurm_rescue(
                metallicity=metallicity,
                missing_indices=missing_indices,
                job_array_length=job_array_length,
                partition="normal",
                email="test@example.com",
                walltime="20:00:00",
                account="test_account",
                mem_per_cpu="4G",
                path_to_posydon="/posydon",
                path_to_posydon_data="/data"
            )

            str_met = convert_metallicity_to_string(metallicity)
            filename = f"{str_met}_Zsun_rescue.slurm"
            assert os.path.exists(filename)

            with open(filename, "r") as f:
                content = f.read()
                assert "#SBATCH --array=1,3,5,7" in content
                assert f"#SBATCH --job-name={str_met}_popsyn_rescue" in content
                assert f"export SLURM_ARRAY_TASK_COUNT={job_array_length}" in content
                assert "export SLURM_ARRAY_TASK_MIN=0" in content
        finally:
            os.chdir(original_dir)

    def test_create_bash_submit_script(self, tmp_path):
        """Test that create_bash_submit_script creates a submission script."""
        original_dir = os.getcwd()
        os.chdir(tmp_path)

        try:
            metallicities = [0.01, 1.0, 2.0]
            filename = "test_submit.sh"

            totest.create_bash_submit_script(filename, metallicities)

            assert os.path.exists(filename)

            with open(filename, "r") as f:
                content = f.read()
                assert "#!/bin/bash" in content
                # Check that scripts for each metallicity are referenced
                assert "1e+00_Zsun_slurm_array.slurm" in content
                assert "1e-02_Zsun_slurm_array.slurm" in content
                assert "2e+00_Zsun_slurm_array.slurm" in content
                assert "sbatch --parsable" in content
                assert "--dependency=afterok" in content
        finally:
            os.chdir(original_dir)

    def test_create_bash_submit_rescue_script(self, tmp_path):
        """Test that create_bash_submit_rescue_script creates a rescue submission script."""
        original_dir = os.getcwd()
        os.chdir(tmp_path)

        try:
            rescue_scripts = [
                str(tmp_path / "1e+00_Zsun_rescue.slurm"),
                str(tmp_path / "1e-02_Zsun_rescue.slurm")
            ]

            # Create dummy rescue script files
            for script in rescue_scripts:
                with open(script, "w") as f:
                    f.write("#!/bin/bash\n")

            result = totest.create_bash_submit_rescue_script(str(tmp_path), rescue_scripts)

            resubmit_file = tmp_path / "resubmit_slurm.sh"
            assert os.path.exists(resubmit_file)
            assert result == str(resubmit_file)

            with open(resubmit_file, "r") as f:
                content = f.read()
                assert "#!/bin/bash" in content
                assert "1e+00_Zsun_rescue.slurm" in content
                assert "1e-02_Zsun_rescue.slurm" in content
                assert "1e+00_Zsun_merge_popsyn.slurm" in content
                assert "1e-02_Zsun_merge_popsyn.slurm" in content
        finally:
            os.chdir(original_dir)


class TestRescueScriptFunctions:
    """Test class for rescue script generation functions."""

    @pytest.fixture
    def mock_args(self):
        """Create mock arguments."""
        args = MagicMock()
        args.run_folder = "/test/run_folder"
        args.job_array = 10
        args.walltime = "24:00:00"
        args.mem_per_cpu = "4G"
        args.partition = "normal"
        args.account = "test_account"
        args.email = "test@example.com"
        return args

    @pytest.fixture
    def mock_batch_status(self):
        """Create mock batch status."""
        return {
            'metallicity': 1.0,
            'missing_indices': [1, 2, 5],
            'status': 'incomplete'
        }

    def test_create_batch_rescue_script_no_missing_indices(self, mock_args):
        """Test that create_batch_rescue_script raises error when no missing indices."""
        batch_status = {'metallicity': 1.0, 'missing_indices': []}

        with pytest.raises(ValueError, match="No missing indices found"):
            totest.create_batch_rescue_script(mock_args, batch_status)

    def test_create_batch_rescue_script(self, tmp_path, mock_args, mock_batch_status):
        """Test that create_batch_rescue_script generates a rescue script."""
        # Setup temporary directory structure
        run_folder = tmp_path / "run"
        run_folder.mkdir()
        mock_args.run_folder = str(run_folder)


        # Create a mock SLURM array script
        slurm_script = run_folder / "1e+00_Zsun_slurm_array.slurm"
        slurm_content = textwrap.dedent("""\
            #!/bin/bash
            #SBATCH --array=0-9
            #SBATCH --time=24:00:00
            #SBATCH --mem-per-cpu=4G
            #SBATCH --partition=normal
            #SBATCH --account=test_account
            #SBATCH --mail-user=test@example.com
            export PATH_TO_POSYDON=/path/to/posydon
            export PATH_TO_POSYDON_DATA=/path/to/data
            srun python ./run_metallicity.py 1.0
        """)
        slurm_script.write_text(slurm_content)

        original_dir = os.getcwd()
        os.chdir(run_folder)

        try:
            result = totest.create_batch_rescue_script(mock_args, mock_batch_status)

            # Check that rescue script was created
            rescue_script = run_folder / "1e+00_Zsun_rescue.slurm"
            assert rescue_script.exists()
            assert result == str(rescue_script)

            # Verify content
            content = rescue_script.read_text()
            assert "#SBATCH --array=1,2,5" in content
            assert "#SBATCH --job-name=1e+00_popsyn_rescue" in content
        finally:
            os.chdir(original_dir)
