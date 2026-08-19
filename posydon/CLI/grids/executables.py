"""Executables helpers for the posydon-grid setup CLI."""

import os
import shutil

from posydon.utils.posydonwarning import Pwarn


def _run_shell(command: str) -> None:
    """Run a shell command via os.system.

    Args:
        command: The shell command to execute.
    """
    os.system(command)


def check_file_exist(file_path: str, raise_error: bool = True) -> bool:
    """Check that a file exists, optionally raising a ValueError.

    Args:
        file_path: Path to the file to check.
        raise_error: If True, raise a ValueError when the file is missing.

    Returns:
        True if the file exists, False otherwise.
    """
    if os.path.exists(file_path):
        return True
    print("{0} does not exist".format(file_path))
    if raise_error:
        raise ValueError("{0} does not exist".format(file_path))
    return False


def make_executables(
    mesa_extras: dict, working_directory: str = os.getcwd()
) -> tuple[str, str, str]:
    """Pass mesa extra function and compile binary executable on the fly.

    Args:
        mesa_extras: Dictionary mapping extra/makefile names to their source
            file paths.
        working_directory: Directory under which the star and binary folders
            are created.

    Returns:
        Tuple of paths to the compiled binary, star1 and star2 executables.
    """
    star1_src_folder = os.path.join(working_directory, 'star1', 'src')
    if os.path.exists(star1_src_folder):  shutil.rmtree(star1_src_folder)
    os.makedirs(star1_src_folder)

    star1_make_folder = os.path.join(working_directory, 'star1', 'make')
    if os.path.exists(star1_make_folder):  shutil.rmtree(star1_make_folder)
    os.makedirs(star1_make_folder)

    star2_src_folder = os.path.join(working_directory, 'star2', 'src')
    if os.path.exists(star2_src_folder):  shutil.rmtree(star2_src_folder)
    os.makedirs(star2_src_folder)

    star2_make_folder = os.path.join(working_directory, 'star2', 'make')
    if os.path.exists(star2_make_folder):  shutil.rmtree(star2_make_folder)
    os.makedirs(star2_make_folder)

    binary_src_folder = os.path.join(working_directory, 'binary', 'src')
    if os.path.exists(binary_src_folder):  shutil.rmtree(binary_src_folder)
    os.makedirs(binary_src_folder)

    binary_make_folder = os.path.join(working_directory, 'binary', 'make')
    if os.path.exists(binary_make_folder):  shutil.rmtree(binary_make_folder)
    os.makedirs(binary_make_folder)

    if os.path.exists('mk'):
        Pwarn('Replace mk', "OverwriteWarning")
    with open('mk', "w") as f:
        for k, v in mesa_extras.items():
            if v is not None:
                if ('binary_extras' in k) or ('binary_run' in k):
                    shutil.copy(v, binary_src_folder)
                elif ('star_run' in k):
                    shutil.copy(v, star1_src_folder)
                    shutil.copy(v, star2_src_folder)
                elif ('star1_extras' in k):
                    shutil.copy(v, star1_src_folder)
                elif ('star2_extras' in k):
                    shutil.copy(v, star2_src_folder)
                elif 'makefile_binary' in k:
                    shutil.copy(v, os.path.join(binary_make_folder, k))
                    f.write('cd {0}\n'.format(binary_make_folder))
                    f.write('make -f {0}\n'.format(k))
                elif 'makefile_star' in k:
                    shutil.copy(v, os.path.join(star1_make_folder, k))
                    f.write('cd {0}\n'.format(star1_make_folder))
                    f.write('make -f {0}\n'.format(k))
                    shutil.copy(v, os.path.join(star2_make_folder, k))
                    f.write('cd {0}\n'.format(star2_make_folder))
                    f.write('make -f {0}\n'.format(k))
                elif 'mesa_dir' == k:
                    continue
                else:
                    shutil.copy(v, working_directory)

    _run_shell("chmod 755 mk")
    _run_shell('./mk')
    return os.path.join(working_directory,'binary','binary'), \
            os.path.join(working_directory,'star1','star'), \
            os.path.join(working_directory,'star2','star')
