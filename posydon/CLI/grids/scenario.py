"""Scenario-based inlist discovery for the posydon grid setup."""

import glob
import os
import subprocess

from posydon.utils.posydonwarning import Pwarn


def _run_git(args: list[str], cwd: str | None = None) -> subprocess.CompletedProcess:
    """Run a git command and return the resulting CompletedProcess.

    Args:
        args: The command and its arguments.
        cwd: Optional working directory for the command.

    Returns:
        The CompletedProcess produced by subprocess.run.
    """
    return subprocess.run(args, cwd=cwd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def setup_inlists_from_scenario(
    source: str, gitcommit: str, system_type: str, mesa_inlists: dict, mesa_extras: dict
) -> None:
    """Dynamically find the inlists the user wants from the supplied info.

    Args:
        source: Either 'posydon' or 'user'.
        gitcommit: The 'branch-githash' of the inlist repository to use.
        system_type: The MESA system type of the grid.
        mesa_inlists: The mesa_inlists section of the inifile, mutated in place.
        mesa_extras: The mesa_extras section of the inifile, mutated in place.

    Raises:
        ValueError: If the source is not understood or the user gitcommit has
            an invalid format.
    """
    where_am_i_now = os.getcwd()
    try:
        print("We are going to dynamically fetch the posydon inlists based on your scenario")
        if source == 'posydon':
            print("You have selected posydon as your source")
            print("checking if we have already cloned POSYDON-MESA-INLISTS for you")
            if not os.path.isdir('{0}/.posydon_mesa_inlists'.format(os.environ['HOME'])):
                print("We are clonining the repo for you")
                _run_git(['git', 'clone', 'https://github.com/POSYDON-code/POSYDON-MESA-INLISTS.git', '{0}/.posydon_mesa_inlists'.format(os.environ['HOME'])])
            else:
                Pwarn("git repository is already there, using that", "OverwriteWarning")

            inlists_dir = '{0}/.posydon_mesa_inlists'.format(os.environ['HOME'])
            branch = gitcommit.split('-')[0]
            githash = gitcommit.split('-')[1]

        elif source == 'user':
            print("You have selected user as your source "
                  "checking if we have already cloned USER-MESA-INLISTS for you "
                  "Validating the name of the git hash you want to use..."
                  "must be of format 'branch-githash'")

            if len(gitcommit.split('-')) != 2:
                raise ValueError("You have supplied an invalid user gitcommit format, must be of format 'branch-githash'")

            branch = gitcommit.split('-')[0]
            githash = gitcommit.split('-')[1]

            if not os.path.isdir('{0}/.user_mesa_inlists'.format(os.environ['HOME'])):
                print("We are clonining the repo for you")
                _run_git(['git', 'clone', 'https://github.com/POSYDON-code/USER-MESA-INLISTS.git', '{0}/.user_mesa_inlists'.format(os.environ['HOME'])])
            else:
                Pwarn("git repository is already there, using that", "OverwriteWarning")

            inlists_dir = '{0}/.user_mesa_inlists'.format(os.environ['HOME'])
            branch = gitcommit.split('-')[0]
            githash = gitcommit.split('-')[1]

        else:
            raise ValueError("supplied source is not valid/understood. Valid sources are user and posydon")

        os.chdir(inlists_dir)
        print("checking out branch: {0}".format(branch))

        _run_git(['git', 'checkout', '{0}'.format(branch)])

        print("For posterity we are pulling (specifically needed if you already have the repo clone)")
        _run_git(['git', 'pull'])

        print("checking out commit/tag: {0}".format(githash))

        _run_git(['git', 'checkout', '{0}'.format(githash)])

        if source == 'posydon':
            inlists_location_common = '{0}/{1}/{2}/'.format(inlists_dir, 'r11701', "default_common_inlists")
            print("Based on system_type {0} "
                  "We are populating the posydon inlists in the following directory: "
                  "{1}".format(system_type, inlists_location_common))

            inlist1 = os.path.join(inlists_location_common, 'binary', 'inlist1')
            if os.path.isfile(inlist1):
                with open(inlist1) as f:
                    for line in f.readlines():
                        if 'zams_filename' in line:
                            print("ZAMS_FILENAME detected, setting mesa_inlists['zams_filename'] for star 1")
                            zams_filename_1 = os.path.split(line.split("'")[1])[1]
                            zams_file_path = os.path.join(inlists_dir, 'r11701', "ZAMS_models", zams_filename_1)
                            if os.path.isfile(zams_file_path):
                                print("Verified locations of ZAMS data file, {0}".format(zams_file_path))
                                mesa_inlists['zams_filename_1'] = "{0}".format(zams_file_path)
                print("Running Single Grid: Setting mesa_star1_extras to {0}/binary/src/run_star_extras.f".format(inlists_location_common))

            inlist2 = os.path.join(inlists_location_common, 'binary', 'inlist2')
            if os.path.isfile(inlist2):
                with open(inlist2) as f:
                    for line in f.readlines():
                        if 'zams_filename' in line:
                            print("ZAMS_FILENAME detected, setting mesa_inlists['zams_filename'] for star 2")
                            zams_filename_2 = os.path.split(line.split("'")[1])[1]
                            zams_file_path = os.path.join(inlists_dir, 'r11701', "ZAMS_models", zams_filename_2)
                            if os.path.isfile(zams_file_path):
                                print("Verified locations of ZAMS data file, {0}".format(zams_file_path))
                                mesa_inlists['zams_filename_2'] = "{0}".format(zams_file_path)
                print("Running Single Grid: Setting mesa_star2_extras to {0}/binary/src/run_star_extras.f".format(inlists_location_common))

            print("Updating inifile values")
            mesa_inlists['star1_controls_posydon_defaults'] = '{0}/binary/inlist1'.format(inlists_location_common)
            mesa_inlists['star1_job_posydon_defaults'] = '{0}/binary/inlist1'.format(inlists_location_common)
            mesa_inlists['star2_controls_posydon_defaults'] = '{0}/binary/inlist2'.format(inlists_location_common)
            mesa_inlists['star2_job_posydon_defaults'] = '{0}/binary/inlist2'.format(inlists_location_common)
            mesa_inlists['binary_controls_posydon_defaults'] = '{0}/binary/inlist_project'.format(inlists_location_common)
            mesa_inlists['binary_job_posydon_defaults'] = '{0}/binary/inlist_project'.format(inlists_location_common)

            mesa_inlists['star_history_columns'] = '{0}/history_columns.list'.format(inlists_location_common)
            mesa_inlists['binary_history_columns'] = '{0}/binary_history_columns.list'.format(inlists_location_common)
            mesa_inlists['profile_columns'] = '{0}/profile_columns.list'.format(inlists_location_common)

            mesa_extras['posydon_binary_extras'] = '{0}/binary/src/run_binary_extras.f'.format(inlists_location_common)
            mesa_extras['posydon_star_binary_extras'] = '{0}/binary/src/run_star_extras.f'.format(inlists_location_common)

            mesa_extras["mesa_star1_extras"] = '{0}/binary/src/run_star_extras.f'.format(inlists_location_common)

            if system_type == 'HeMS-HMS':
                inlists_location = '{0}/{1}/{2}/'.format(inlists_dir, 'r11701', system_type)
                print("Based on system_type {0} "
                      "We are populating the user inlists in the following directory: "
                      "{1}".format(system_type, inlists_location))

                if os.path.isfile(os.path.join(inlists_location, "binary", "inlist1")):
                    mesa_inlists['star1_controls_user'] = '{0}'.format(os.path.join(inlists_location, "binary", "inlist1"))
                    mesa_inlists['star1_job_user'] = '{0}'.format(os.path.join(inlists_location, "binary", "inlist1"))

                if os.path.isfile(os.path.join(inlists_location, "binary", "inlist2")):
                    mesa_inlists['star2_controls_user'] = '{0}'.format(os.path.join(inlists_location, "binary", "inlist2"))
                    mesa_inlists['star2_job_user'] = '{0}'.format(os.path.join(inlists_location, "binary", "inlist2"))

                if os.path.isdir(os.path.join(inlists_location, "star1_formation")):
                    print("We are making star1 so we can unset the zams file we were going to use for star1")
                    mesa_inlists['zams_filename_1'] = None

                    star1_formation_scenario = sorted(glob.glob(os.path.join(inlists_location, "star1_formation", "*step*")))
                    print("These are the user we are using to make star1: {0}".format(star1_formation_scenario))
                    print("We are going to add a layer of posydon default common inlists to these user steps: {0}".format('{0}/binary/inlist1'.format(inlists_location_common)))
                    mesa_inlists['star1_formation_controls_posydon_defaults'] = []
                    mesa_inlists['star1_formation_job_posydon_defaults'] = []
                    for i in range(len(star1_formation_scenario)):
                        mesa_inlists['star1_formation_controls_posydon_defaults'].append('{0}/binary/inlist1'.format(inlists_location_common))
                        mesa_inlists['star1_formation_job_posydon_defaults'].append('{0}/binary/inlist1'.format(inlists_location_common))

                    mesa_inlists['star1_formation_controls_user'] = star1_formation_scenario
                    mesa_inlists['star1_formation_job_user'] = star1_formation_scenario

            elif system_type != "HMS-HMS" and not mesa_inlists['single_star_grid']:
                inlists_location = '{0}/{1}/{2}/'.format(inlists_dir, 'r11701', system_type)
                print("Based on system_type {0} "
                      "We are populating the user inlists in the following directory: "
                      "{1}".format(system_type, inlists_location))

                if os.path.isfile(os.path.join(inlists_location, "binary", "inlist_project")):
                    mesa_inlists['binary_controls_user'] = '{0}'.format(os.path.join(inlists_location, "binary", "inlist_project"))
                    mesa_inlists['binary_job_user'] = '{0}'.format(os.path.join(inlists_location, "binary", "inlist_project"))

                if os.path.isfile(os.path.join(inlists_location, "binary", "inlist1")):
                    mesa_inlists['star1_controls_user'] = '{0}'.format(os.path.join(inlists_location, "binary", "inlist1"))
                    mesa_inlists['star1_job_user'] = '{0}'.format(os.path.join(inlists_location, "binary", "inlist1"))

                if os.path.isfile(os.path.join(inlists_location, "binary", "inlist2")):
                    mesa_inlists['star2_controls_user'] = '{0}'.format(os.path.join(inlists_location, "binary", "inlist2"))
                    mesa_inlists['star2_job_user'] = '{0}'.format(os.path.join(inlists_location, "binary", "inlist2"))

                if os.path.isfile(os.path.join(inlists_location, "history_columns.list")):
                    mesa_inlists['star_history_columns'] = os.path.join(inlists_location, "history_columns.list")

                if os.path.isfile(os.path.join(inlists_location, "binary_history_columns.list")):
                    mesa_inlists['binary_history_columns'] = os.path.join(inlists_location, "binary_history_columns.list")

                if os.path.isfile(os.path.join(inlists_location, "profile_columns.list")):
                    mesa_inlists['profile_columns'] = os.path.join(inlists_location, "profile_columns.list")

                if os.path.isfile(os.path.join(inlists_location, "src", "run_binary_extras.f")):
                    mesa_extras['user_binary_extras'] = '{0}'.format(os.path.join(inlists_location, "src", "run_binary_extras.f"))

                if os.path.isfile(os.path.join(inlists_location, "src", "run_star_extras.f")):
                    mesa_extras['user_star_binary_extras'] = '{0}'.format(os.path.join(inlists_location, "src", "run_star_extras.f"))

                if os.path.isdir(os.path.join(inlists_location, "star1_formation")):
                    print("We are making star1 so we can unset the zams file we were going to use for star1")
                    mesa_inlists['zams_filename_1'] = None

                    star1_formation_scenario = sorted(glob.glob(os.path.join(inlists_location, "star1_formation", "*step*")))
                    print("These are the user we are using to make star1: {0}".format(star1_formation_scenario))
                    print("We are going to add a layer of posydon default common inlists to these user steps: {0}".format('{0}/binary/inlist1'.format(inlists_location_common)))
                    mesa_inlists['star1_formation_controls_posydon_defaults'] = []
                    mesa_inlists['star1_formation_job_posydon_defaults'] = []
                    for i in range(len(star1_formation_scenario)):
                        mesa_inlists['star1_formation_controls_posydon_defaults'].append('{0}/binary/inlist1'.format(inlists_location_common))
                        mesa_inlists['star1_formation_job_posydon_defaults'].append('{0}/binary/inlist1'.format(inlists_location_common))

                    mesa_inlists['star1_formation_controls_user'] = star1_formation_scenario
                    mesa_inlists['star1_formation_job_user'] = star1_formation_scenario

                if os.path.isdir(os.path.join(inlists_location, "star2_formation")):
                    print("We are making star2 so we can unset the zams file we were going to use for star2")
                    mesa_inlists['zams_filename_2'] = None

                    star2_formation_scenario = sorted(glob.glob(os.path.join(inlists_location, "star2_formation", "*step*")))
                    print("These are the user we are using to make star2: {0}".format(star2_formation_scenario))
                    print("We are going to add a layer of posydon default common inlists to these user steps: {0}".format('{0}/binary/inlist1'.format(inlists_location_common)))
                    mesa_inlists['star2_formation_controls_posydon_defaults'] = []
                    mesa_inlists['star2_formation_job_posydon_defaults'] = []
                    for i in range(len(star2_formation_scenario)):
                        mesa_inlists['star2_formation_controls_posydon_defaults'].append('{0}/binary/inlist1'.format(inlists_location_common))
                        mesa_inlists['star2_formation_job_posydon_defaults'].append('{0}/binary/inlist1'.format(inlists_location_common))

                    mesa_inlists['star2_formation_controls_user'] = star2_formation_scenario
                    mesa_inlists['star2_formation_job_user'] = star2_formation_scenario

            if system_type == "HMS-HMS" and mesa_inlists['single_star_grid']:
                print("You want a single star HMS grid, this means that we need to make a user inlist on the fly with a single line "
                      "x_logical_ctrl(1)=.true.")
                special_single_star_user_inlist = os.path.join(os.getcwd(), "special_single_star_user_inlist")
                if os.path.exists(special_single_star_user_inlist):
                    Pwarn('Replace ' + special_single_star_user_inlist, "OverwriteWarning")
                with open(special_single_star_user_inlist, 'wb') as f:
                    f.write(b'&controls\n\n')
                    f.write('\t{0} = {1}\n'.format("x_logical_ctrl(1)", ".true.").encode('utf-8'))

                    f.write(b'\n\n')

                    f.write(b"""
        / ! end of star1_controls namelist

        """)
                mesa_inlists['star1_controls_special'] = special_single_star_user_inlist
            elif system_type == "CO-He_star" and mesa_inlists['single_star_grid']:
                inlists_location = '{0}/{1}/{2}/'.format(inlists_dir, 'r11701', system_type)
                print("Based on system_type {0} "
                      "We are populating the user inlists in the following directory: "
                      "{1}".format(system_type, inlists_location))

                print("We are making star2 so we can unset the zams file we were going to use for star2")
                mesa_inlists['zams_filename_2'] = None

                single_star_scenario = sorted(glob.glob(os.path.join(inlists_location, "star1_formation", "*step*")))
                print("These are the user inlists used in the single star grid: {0}".format(single_star_scenario))
                mesa_inlists['star1_controls_user'] = single_star_scenario
                mesa_inlists['star1_job_user'] = single_star_scenario
                print("You want a single star He grid, "
                      "this means that we need to make the inlist that will be used to evolve the system "
                      "and make sure we layer on the line "
                      "x_logical_ctrl(1)=.true.")
                special_single_star_user_inlist = os.path.join(os.getcwd(), "special_single_star_user_inlist")
                if os.path.exists(special_single_star_user_inlist):
                    Pwarn('Replace ' + special_single_star_user_inlist, "OverwriteWarning")
                with open(special_single_star_user_inlist, 'wb') as f:
                    f.write(b'&controls\n\n')
                    f.write('\t{0} = {1}\n'.format("x_logical_ctrl(1)", ".true.").encode('utf-8'))

                    f.write(b'\n\n')
                    f.write(b"""
        / ! end of star1_controls namelist

        """)
                    f.write(b'&star_job\n\n')
                    f.write(b"""
        / / ! end of star_job namelist

        """)
                mesa_inlists['star1_controls_special'].append(special_single_star_user_inlist)
    finally:
        os.chdir(where_am_i_now)
