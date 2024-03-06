"""This class evolves a binary object according to a MESA grid."""


__authors__ = [
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Devina Misra <devina.misra@unige.ch>",
    "Emmanouil Zapartas <ezapartas@gmail.com>",
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
    "Jeffrey Andrews <jeffrey.andrews@northwestern.edu>",
    "Nam Tran <tranhn03@gmail.com>",
    "Zepei Xing <Zepei.Xing@unige.ch>",
    "Tassos Fragos <Anastasios.Fragkos@unige.ch>",
]


import os
import warnings
import numpy as np

from posydon.interpolation.interpolation import psyTrackInterp
from posydon.binary_evol.binarystar import BINARYPROPERTIES
from posydon.binary_evol.singlestar import STARPROPERTIES
from posydon.utils import common_functions as cf
from posydon.binary_evol.binarystar import BinaryStar
from posydon.interpolation.IF_interpolation import IFInterpolator
from posydon.utils.common_functions import (flip_stars,
                                            convert_metallicity_to_string,
                                            CO_radius, infer_star_state,
                                            set_binary_to_failed,)
from posydon.utils.data_download import data_download, PATH_TO_POSYDON_DATA
from posydon.grids.MODELS import MODELS


# left POSYDON, right MESA
POSYDON_TO_MESA = {
    'binary': {
        'state': True,
        'event': True,
        'time': 'age',
        'separation': 'binary_separation',
        'orbital_period': 'period_days',
        'eccentricity': True,
        'V_sys': True,
        'rl_relative_overflow_1': 'rl_relative_overflow_1',
        'rl_relative_overflow_2': 'rl_relative_overflow_2',
        'lg_mtransfer_rate': 'lg_mtransfer_rate',
        'trap_radius': 'trap_radius',
        'acc_radius': 'acc_radius',
        't_sync_rad_1': 't_sync_rad_1',
        't_sync_conv_1': 't_sync_conv_1',
        't_sync_rad_2': 't_sync_rad_2',
        't_sync_conv_2': 't_sync_conv_2',
        'mass_transfer_case': True,
        'nearest_neighbour_distance': True
    },
    'star': {
        'state': True,
        # TODO: 'Z', this should be included in the default columns of psygrid
        'metallicity': True,
        'mass': 'mass',                         # binary history
        'log_R': 'log_R',
        'log_L': 'log_L',
        'lg_mdot': 'lg_mstar_dot',              # binary history
        'lg_system_mdot': 'lg_system_mdot',     # binary history
        'lg_wind_mdot': 'lg_wind_mdot',         # binary history
        'he_core_mass': 'he_core_mass',
        'he_core_radius': 'he_core_radius',
        'c_core_mass': 'c_core_mass',
        'c_core_radius': 'c_core_radius',
        'o_core_mass': 'o_core_mass',
        'o_core_radius': 'o_core_radius',
        'co_core_mass': 'co_core_mass',
        'co_core_radius': 'co_core_radius',
        'center_h1': 'center_h1',
        'center_he4': 'center_he4',
        'center_c12': 'center_c12',
        'center_n14': 'center_n14',
        'center_o16': 'center_o16',
        'surface_h1': 'surface_h1',
        'surface_he4': 'surface_he4',
        'surface_c12': 'surface_c12',
        'surface_n14': 'surface_n14',
        'surface_o16': 'surface_o16',
        'log_LH': 'log_LH',
        'log_LHe': 'log_LHe',
        'log_LZ': 'log_LZ',
        'log_Lnuc': 'log_Lnuc',
        'c12_c12': 'c12_c12',
        'center_gamma': 'center_gamma',
        # 'avg_c_in_c_core', remove this for the moment because HMS-HMS did
        # not have this column
        'avg_c_in_c_core': None,
        'surf_avg_omega': 'surf_avg_omega',
        'surf_avg_omega_div_omega_crit': 'surf_avg_omega_div_omega_crit',
        'total_moment_of_inertia': 'total_moment_of_inertia',
        'log_total_angular_momentum': 'log_total_angular_momentum',
        'spin': 'spin_parameter',
        'conv_env_top_mass': 'conv_env_top_mass',
        'conv_env_bot_mass': 'conv_env_bot_mass',
        'conv_env_top_radius': 'conv_env_top_radius',
        'conv_env_bot_radius': 'conv_env_bot_radius',
        'conv_env_turnover_time_g': 'conv_env_turnover_time_g',
        'conv_env_turnover_time_l_b': 'conv_env_turnover_time_l_b',
        'conv_env_turnover_time_l_t': 'conv_env_turnover_time_l_t',
        'envelope_binding_energy': 'envelope_binding_energy',
        'mass_conv_reg_fortides': 'mass_conv_reg_fortides',
        'thickness_conv_reg_fortides': 'thickness_conv_reg_fortides',
        'radius_conv_reg_fortides': 'radius_conv_reg_fortides',
        'lambda_CE_1cent': 'lambda_CE_1cent',
        'lambda_CE_10cent': 'lambda_CE_10cent',
        'lambda_CE_30cent': 'lambda_CE_30cent',
        'lambda_CE_pure_He_star_10cent': 'lambda_CE_pure_He_star_10cent',
        'profile': True
    }
}


class MesaGridStep:
    """Superclass for steps using the POSYDON grids."""

    def __init__(
            self,
            metallicity,
            grid_name,
            path=PATH_TO_POSYDON_DATA,
            interpolation_path=None,
            interpolation_filename=None,
            interpolation_method="linear3c_kNN",
            save_initial_conditions=True,
            track_interpolation=False,
            stop_method='stop_at_max_time',     # "stop_at_end",
            stop_star="star_1",
            stop_var_name=None,
            stop_value=None,
            stop_interpolate=True,
            verbose=False):
        """Evolve a binary object given a MESA grid or interpolation object.

        Parameters
        ----------
        interpolation_method : 'nearest_neighbour', 'linear_knn_interpolation'
                               'linear_NN_interpolation'
             - nearest_neighbour requires a psygrid in the h5 format
             - linear_knn_interpolation requires a pkl file trained with
               interpolation class linear inteprolation, knn classification
             - linear_NN_interpolation requires a pkl file trained with
                 interpolation class linear inteprolation, neural network
                 classification
        track_interpolation: bool
             - True requires `nearest_neighbour` and will append the entire
               MESA history of properties.
        stop_method : stop_at_end or stop_at_max_time
             - stop_at_end will stop the binary at the end of the MESA track
             - stop_at_max_time binary will stop in the middle of a MESA track
               if doing so would exceed the binary's max_time
             - stop_at_condition will stop the binary when a five condition is
               met (see next parameters)
        stop_var_name: keys in STARPROPERTIES and
             BINARYPROPERTIES
             - key is one of the STARPROPERTIES and BINARYPROPERTIES
        stop_star: star_1 or star_2
             - if you choose stop_var_name from STARPROPERTIES you need to
               indicate which star you are referring too
        stop_value: float
             - when value is reached in the next MESA history model, the binary
               will stop the evolution
        stop_interpolate: True or False
             - True can only be chosen if stop_var_name='time', it means that
               the binary properties will be linearly interpolated between two
               timestamp to determine their values at the arbitrary time
               stop_value

        """
        # class variable
        self.path = path
        self.interpolation_method = interpolation_method
        self.save_initial_conditions = save_initial_conditions
        self.track_interpolation = track_interpolation
        self.stop_method = stop_method
        self.verbose = verbose

        if (self.track_interpolation
                and self.interpolation_method != 'nearest_neighbour'):
            raise ValueError('Track interpolation is currently supported only '
                             'by the nearest neighbour interpolation method!')

        # we load NN any time stop_at_max_time requested - regardless
        # of interp method
        if (self.stop_method == 'stop_at_max_time'
                or self.interpolation_method == 'nearest_neighbour'):
            self.load_psyTrackInterp(grid_name)

        grid_name = grid_name.replace('_%d', '')

        # Check interpolation method provided
        self.supported_interp_methods = ['linear_kNN', 'linear3c_kNN',
                                         '1NN_1NN']
        if self.interpolation_method in self.supported_interp_methods:

            # Set the interpolation path
            if interpolation_path is None:
                interpolation_path = (
                    self.path + os.path.split(grid_name)[0]
                    + '/interpolators/%s/' % self.interpolation_method)

            # Set the interpolation filename
            if interpolation_filename is None:
                interpolation_filename = (
                    interpolation_path
                    + os.path.split(grid_name)[1].replace('h5', 'pkl'))
            else:
                interpolation_filename = (interpolation_path
                                          + interpolation_filename)

            self.load_Interp(interpolation_filename)

            if (not (hasattr(self, '_psyTrackInterp')
                     or hasattr(self, '_Interp'))):
                raise ValueError(
                    'No interpolation methods specified in kwargs!\n'
                    'interpolation_method: {0}'.
                    format(self.interpolation_method))

        # this hook is used to avoid computing star and binary states in case
        # track_interpolation = False and stop = 'stop_at_condition' where
        # we drop the history
        self.flush_history = False
        self.flush_entries = None
        self.stop_star = stop_star
        self.stop_var_name = stop_var_name
        self.stop_value = stop_value
        self.stop_interpolate = stop_interpolate

    def load_psyTrackInterp(self, grid_name):
        """Load the interpolator that has been trained on the grid."""
        # Check if interpolation files exist
        filename = os.path.join(self.path,grid_name)
        if not (os.path.exists(filename.replace('%d','0')) or
                os.path.exists(filename.replace('_%d',''))):
            data_download()

        if self.verbose:
            print("loading psyTrackInterp: {}".format(filename))
        self._psyTrackInterp = psyTrackInterp(filename,
                                              interp_in_q=self.interp_in_q,
                                              verbose=self.verbose)
        self._psyTrackInterp.train()

    def load_Interp(self, filename):
        """Load the interpolator that has been trained on the grid."""
        if self.verbose:
            print("loading Interp: {}".format(filename))

        # Check if interpolation files exist
        if not os.path.exists(filename):
            data_download()

        # Load interpolator
        self._Interp = IFInterpolator()
        self._Interp.load(filename=filename)

    def close(self):
        """Close the inteprolator."""
        pass
        # if hasattr(self, '_Interp'):
        #     self._Interp.close()

    def get_final_MESA_step_time(self):
        """Infer the maximum MESA step time.

        Use the specified interpolation method to determine the maximum
        MESA simulation time about to be added to the binary. If this step
        exceeds the maximum age of the binary, a step(NN) is done.

        """
        if self.interpolation_method == 'nearest_neighbour':
            self.closest_binary, self.nearest_neighbour_distance, \
                self.termination_flags = self._psyTrackInterp.evaluate(
                    self.binary)
            if self.closest_binary.binary_history is None:
                return
            key = POSYDON_TO_MESA['binary']['time']
            max_MESA_sim_time = self.closest_binary.binary_history[key][-1]

        elif self.interpolation_method in self.supported_interp_methods:
            self.final_values, self.classes = self._Interp.evaluate(
                self.binary)

            max_MESA_sim_time = self.final_values[
                                        POSYDON_TO_MESA['binary']['time']]
        else:
            raise ValueError("unknown interpolation method: {}".
                             format(self.interpolation_method))
        return max_MESA_sim_time

    def __call__(self, binary):
        """Evolve a binary using the MESA step."""
        if not isinstance(binary, BinaryStar):
            raise ValueError("Must be an instance of BinaryStar")
        if not hasattr(self, 'step'):
            raise AttributeError("No step defined for {}".format(
                self.__name__))
        if self.flip_stars_before_step:
            flip_stars(binary)
        max_MESA_sim_time = self.get_final_MESA_step_time()

        if max_MESA_sim_time is None:
            if self.flip_stars_before_step:
                flip_stars(binary)
            binary.state = 'initial_RLOF'
            # binary.event = 'END'
            return

        binary_start_time = binary.time
        step_will_exceed_max_time = (binary.time+max_MESA_sim_time
                                     > binary.properties.max_simulation_time)
        if (step_will_exceed_max_time
                and self.stop_method == 'stop_at_max_time'):
            # self.step(binary, interp_method='nearest_neighbour')
            if self.interpolation_method != 'nearest_neighbour':
                self.closest_binary, self.nearest_neighbour_distance, \
                    self.termination_flags = self._psyTrackInterp.evaluate(
                                 self.binary)

            if self.track_interpolation:
                self.flush_history = False
            else:
                self.flush_history = True
            self.update_properties_NN(star_1_CO=self.star_1_CO,
                                      star_2_CO=self.star_2_CO,
                                      track_interpolation=True)
        else:
            self.step(binary, interp_method=self.interpolation_method)
        if (self.stop_method == 'stop_at_max_time'
                and binary.time >= binary.properties.max_simulation_time):
            
            # self.flush_history = True # needed???

            # stop_at_condition looks through the MESA output appended to the
            # binary object. The times have already been added to the binary's
            # start time, so the correct value to pass is max_simulation_time
            self.stop_at_condition(binary,
                                   star='star_1',
                                   property='time',
                                   value=binary.properties.max_simulation_time,
                                   interpolate=self.stop_interpolate,
                                   star_1_CO=self.star_1_CO,
                                   star_2_CO=self.star_2_CO)
        elif self.stop_method == 'stop_at_condition':
            self.stop_at_condition(binary,
                                   property=self.stop_var_name,
                                   value=self.stop_value,
                                   star=self.stop_star,
                                   interpolate=self.stop_interpolate,
                                   star_1_CO=self.star_1_CO,
                                   star_2_CO=self.star_2_CO)
        else:
            self.stop_at_end(binary,
                             property=self.stop_var_name,
                             value=self.stop_value,
                             star=self.stop_star,
                             interpolate=self.stop_interpolate,
                             star_1_CO=self.star_1_CO,
                             star_2_CO=self.star_2_CO)
        if self.flip_stars_before_step:
            flip_stars(binary)
        if binary.time > binary.properties.max_simulation_time:
            binary.event = 'MaxTime_exceeded'
        elif binary.time == binary.properties.max_simulation_time:
            binary.event = 'maxtime'
        return

    def step(self, binary, interp_method=None):
        """Take the step by calling the appropriate interpolator."""
        if interp_method is None:
            interp_method = self.interpolation_method

        if interp_method == 'nearest_neighbour':
            if self.track_interpolation:
                self.flush_history = False
            else:
                self.flush_history = True
            self.update_properties_NN(
                star_1_CO=self.star_1_CO, star_2_CO=self.star_2_CO,
                track_interpolation=self.track_interpolation)
        elif interp_method in self.supported_interp_methods:
            self.initial_final_interpolation(star_1_CO=self.star_1_CO,
                                             star_2_CO=self.star_2_CO)
        else:
            raise ValueError('Invalid interpolation method: {}'.
                             format(interp_method))

    def single_star(self, star_type):
        """Check if the type of the star is supported."""
        if star_type == 'single-HMS':
            pass
        elif star_type == 'single-HeMS':
            pass
        else:
            raise ValueError('Single star_type = %s unknown!' % star_type)

    # UPDATE STAR AND BINARY METHODS

    def update_properties_NN(self, star_1_CO=False, star_2_CO=False,
                             track_interpolation=False):
        """Update properites according to nearest neighbour interpolation.

        Parameters
        ----------
        star_1_CO : bool
            `True` if star_1 is a compact object.
        star_2_CO : bool
            `True` if star_2 is a compact object.
        track_interpolation : bool
            If True, uses track interpolation over initial-final interpolation

        """
        # simplify verbosity of code
        binary = self.binary
        stars = [binary.star_1, binary.star_2]
        stars_CO = [star_1_CO, star_2_CO]
        cb = self.closest_binary
        cb_bh = cb.binary_history
        cb_hs = [cb.history1, cb.history2]
        cb_fps = [cb.final_profile1, cb.final_profile2]

        # TOOD: I removed this which is now done in get_final_MESA_step_time
        # find the nearest_neighbour and the distance
        # self.closest_binary, self.nearest_neighbour_distance, \
        # self.termination_flags = self._psyTrackInterp.evaluate(self.binary)
        if (cb_bh['age'].size <= 1 or cb_bh['star_1_mass'].size <= 1):
            setattr(binary, "state", "initial_RLOF")
            # setattr(binary, "event", "END")
            return

        # check if the first interpolation gives 'initial_RLOF'
        interpolation_class = self.termination_flags[0]
        binary_state, binary_event, MT_case = (
            cf.get_binary_state_and_event_and_mt_case(
                binary, interpolation_class, verbose=self.verbose))
        setattr(binary, 'state', binary_state)
        setattr(binary, 'event', binary_event)
        setattr(binary, 'mass_transfer_case', MT_case)
        
        if binary.state == 'initial_RLOF':
            return
        
        if track_interpolation or self.save_initial_conditions:
            len_binary_hist = len(getattr(binary, "time_history"))
            length_binary_hist = len(cb_bh['age']) - 1
            length_star_hist = len(cb_hs[0]['center_he4']) - 1
            length_hist = min(length_binary_hist, length_star_hist)
            empy_h = [None] * length_hist
            MESA_history_bug_fix = False
            if length_binary_hist != length_star_hist:
                MESA_history_bug_fix = True
                warnings.warn(
                    'The MESA star_history and binary_history do not match '
                    'lenght %i != %i. This will cause errors, e.g. '
                    'get_binary_state_and_event_and_mt_case take '
                    'star.mdot_history - star.lg_wind_mdot, to '
                    'avoid the code to break we happened np.nan '
                    'to the missing values!' %
                    (length_binary_hist, length_star_hist))

        # update properties
        for key in BINARYPROPERTIES:
            key_h = key + "_history"
            if POSYDON_TO_MESA['binary'][key] is None:
                setattr(binary, key, None)
                if self.save_initial_conditions:
                    getattr(binary, key_h).append(empy_h[0])
                if track_interpolation:
                    getattr(binary, key_h).extend(empy_h)
            elif key == 'time':
                key_p = POSYDON_TO_MESA['binary'][key]
                current = getattr(self.binary, key)
                setattr(binary, key, current + cb_bh[key_p][-1])
                if self.save_initial_conditions:
                    history_of_attribute = current + cb_bh[key_p][:-1]
                    getattr(binary, key_h).append(history_of_attribute[0])
                if track_interpolation:
                    history_of_attribute = current + cb_bh[key_p][:-1]
                    getattr(binary, key_h).extend(history_of_attribute)
            elif key == 'nearest_neighbour_distance':
                NN_d = self.nearest_neighbour_distance
                setattr(binary, key, NN_d)
                if self.save_initial_conditions:
                    getattr(binary, key_h).append(NN_d)
                if track_interpolation:
                    getattr(binary, key_h).extend([NN_d]*length_hist)
            elif key in ['eccentricity', 'V_sys']:
                v_key = getattr(binary, key_h)[-1]
                setattr(binary, key, v_key)
                if self.save_initial_conditions:
                    getattr(binary, key_h).append(v_key)
                if track_interpolation:
                    getattr(binary, key_h).extend([v_key]*length_hist)
            elif key in ['state', 'event', 'mass_transfer_case']:
                continue
            else:
                key_p = POSYDON_TO_MESA['binary'][key]
                setattr(binary, key, cb_bh[key_p][-1])
                if self.save_initial_conditions:
                    history_of_attribute = cb_bh[key_p][:-1]
                    getattr(binary, key_h).append(history_of_attribute[0])
                if track_interpolation:
                    history_of_attribute = cb_bh[key_p][:-1]
                    getattr(binary, key_h).extend(history_of_attribute)
            if track_interpolation:
                if MESA_history_bug_fix:
                    real_len = max(length_binary_hist, length_star_hist)
                    missing_values = real_len + len_binary_hist - len(
                        getattr(binary, key_h))
                    if missing_values > 0:
                        getattr(binary, key_h).extend([np.nan]*missing_values)
                        # DEBUG
                        # print(key, missing_values)
                        # print('fixed', len(getattr(self.binary, key
                        #                            + "_history")))

        for k, star in enumerate(stars):
            for key in STARPROPERTIES:
                key_h = key + "_history"
                if not stars_CO[k]:
                    if POSYDON_TO_MESA['star'][key] is None:
                        setattr(star, key, None)
                        if self.save_initial_conditions:
                            getattr(star, key_h).append(empy_h[0])
                        if track_interpolation:
                            getattr(star, key_h).extend(empy_h)
                    elif key == 'mass':
                        key_p = 'star_%d_mass' % (k+1)
                        setattr(star, key, cb_bh[key_p][-1])
                        if self.save_initial_conditions:
                            history_of_attribute = cb_bh[key_p][:-1]
                            getattr(star, key_h).append(
                                history_of_attribute[0])
                        if track_interpolation:
                            history_of_attribute = cb_bh[key_p][:-1]
                            getattr(star, key_h).extend(history_of_attribute)
                    elif key == 'metallicity':
                        v_key = getattr(star, 'metallicity')
                        setattr(star, key, v_key)
                        if self.save_initial_conditions:
                            getattr(star, key_h).append(v_key)
                        if track_interpolation:
                            getattr(star, key_h).extend([v_key]*length_hist)
                    elif key == 'profile':
                        setattr(star, key, cb_fps[k])
                        if self.save_initial_conditions:
                            getattr(star, key_h).append(empy_h[0])
                        if track_interpolation:
                            getattr(star, key_h).extend(empy_h)
                    elif key in ['lg_mdot', 'lg_system_mdot', 'lg_wind_mdot']:
                        key_p = POSYDON_TO_MESA['star'][key]+'_%d' % (k+1)
                        setattr(star, key, cb_bh[key_p][-1])
                        if self.save_initial_conditions:
                            history_of_attribute = cb_bh[key_p][:-1]
                            getattr(star, key_h).append(
                                history_of_attribute[0])
                        if track_interpolation:
                            history_of_attribute = cb_bh[key_p][:-1]
                            getattr(star, key_h).extend(
                                history_of_attribute)
                    elif key == 'state':
                        continue
                    else:
                        key_p = POSYDON_TO_MESA['star'][key]
                        setattr(star, key, cb_hs[k][key_p][-1])
                        if self.save_initial_conditions:
                            history_of_attribute = cb_hs[k][key_p][:-1]
                            getattr(star, key_h).append(
                                history_of_attribute[0])
                        if track_interpolation:
                            history_of_attribute = cb_hs[k][key_p][:-1]
                            getattr(star, key_h).extend(history_of_attribute)
                else:   # star is a compact object
                    if key == 'mass':
                        key_p = 'star_%d_mass' % (k+1)
                        setattr(star, key, cb_bh[key_p][-1])
                        if self.save_initial_conditions:
                            history_of_attribute = cb_bh[key_p][:-1]
                            getattr(star, key_h).append(
                                history_of_attribute[0])
                        if track_interpolation:
                            history_of_attribute = cb_bh[key_p][:-1]
                            getattr(star, key_h).extend(history_of_attribute)
                    elif key == 'spin':
                        key_p = 'star_%d_mass' % (k+1)
                        mass_history = np.array(cb_bh[key_p])
                        spin_prev_step = getattr(star, 'spin')
                        spin = cf.spin_stable_mass_transfer(spin_prev_step, mass_history[0],
                                                            mass_history[-1])
                        setattr(star, key, spin)
                        if self.save_initial_conditions:
                            history_of_attribute =[
                                cf.spin_stable_mass_transfer(spin_prev_step ,mass_history[0],
                                                             mass_history[i])
                                for i in range(len(mass_history)-1)]
                            getattr(star, key_h).append(history_of_attribute[0])
                        if track_interpolation:
                            history_of_attribute =[
                                cf.spin_stable_mass_transfer(spin_prev_step, mass_history[0],
                                                             mass_history[i])
                                for i in range(len(mass_history)-1)]
                            getattr(star, key_h).extend(history_of_attribute)
                    elif key == 'log_R':
                        key_p = 'star_%d_mass' % (k+1)
                        mass = cb_bh[key_p][-1]
                        st = infer_star_state(star_mass=mass, star_CO=True)
                        setattr(star, key, np.log10(CO_radius(mass, st)))
                        if self.save_initial_conditions:
                            mass_history = np.array(cb_bh[key_p])
                            state_history = np.array([
                                infer_star_state(star_mass=x, star_CO=True)
                                for x in mass_history])
                            history_of_attribute = np.log10([
                                CO_radius(mass_history[i], state_history[i])
                                for i in range(len(mass_history)-1)])
                            getattr(star, key_h).append(
                                history_of_attribute[0])
                        if track_interpolation:
                            mass_history = np.array(cb_bh[key_p])
                            state_history = np.array([
                                infer_star_state(star_mass=x, star_CO=True)
                                for x in mass_history])
                            history_of_attribute = np.log10([
                                CO_radius(mass_history[i], state_history[i])
                                for i in range(len(mass_history)-1)])
                            getattr(star, key_h).extend(history_of_attribute)
                    elif key == 'metallicity':
                        v_key = getattr(star, 'metallicity')
                        setattr(star, key, v_key)
                        if self.save_initial_conditions:
                            getattr(star, key_h).append(v_key)
                        if track_interpolation:
                            getattr(star, key_h).extend([v_key]*length_hist)
                    elif key in ['lg_system_mdot', 'lg_wind_mdot']:
                        key_p = POSYDON_TO_MESA['star'][key]+'_%d' % (k+1)
                        setattr(star, key, cb_bh[key_p][-1])
                        if self.save_initial_conditions:
                            getattr(star, key_h).append(cb_bh[key_p][0])
                        if track_interpolation:
                            getattr(star, key_h).extend(cb_bh[key_p][:-1])
                    elif key in ['state', 'lg_mdot']:
                        continue
                    elif star.state == 'WD' and key in ['co_core_mass','he_core_mass','center_h1','center_he4','center_c12','center_n14','center_o16']:
                        continue
                    else:
                        setattr(star, key, None)
                        if self.save_initial_conditions:
                            getattr(star, key_h).append(empy_h[0])
                        if track_interpolation:
                            getattr(star, key_h).extend(empy_h)


                if track_interpolation:
                    if MESA_history_bug_fix:
                        real_len = max(length_binary_hist, length_star_hist)
                        missing_values_star = real_len + len_binary_hist - len(
                            getattr(star, key_h))
                        if missing_values_star > 0:
                            getattr(star, key_h).extend(
                                [np.nan]*missing_values)
                            # DEBUG
                            # print(key, missing_values_star_1)
                            # print('fixed', len(getattr(self.binary.star_1,
                            #                            key + "_history")))
        # convert these flags to default POSYDON star states
        setattr(stars[0], 'state',
                cf.check_state_of_star(stars[0], star_CO=stars_CO[0]))
        setattr(stars[1], 'state',
                cf.check_state_of_star(stars[1], star_CO=stars_CO[1]))
        interpolation_class = self.termination_flags[0]
        binary_state, binary_event, MT_case = (
            cf.get_binary_state_and_event_and_mt_case(
                binary, interpolation_class, verbose=self.verbose))
        setattr(binary, 'state', binary_state)
        setattr(binary, 'event', binary_event)
        setattr(binary, 'mass_transfer_case', MT_case)

        setattr(self.binary, f'interp_class_{self.grid_type}', interpolation_class)
        mt_history = self.termination_flags[2] # mass transfer history (TF12 plot label)
        setattr(self.binary, f'mt_history_{self.grid_type}', mt_history)

        if self.save_initial_conditions:
            # history N is how much to look back in the history
            # here N=1 as we only appended back the first entry
            state1_hist = cf.check_state_of_star(stars[0], i=len_binary_hist,
                                                 star_CO=stars_CO[0])
            getattr(stars[0], "state_history").append(state1_hist)
            state2_hist = cf.check_state_of_star(stars[1], i=len_binary_hist,
                                                 star_CO=stars_CO[1])
            getattr(stars[1], "state_history").append(state2_hist)
            binary_state, binary_event, MT_case = (
                cf.get_binary_state_and_event_and_mt_case(
                    binary, interpolation_class, i=len_binary_hist,
                    verbose=self.verbose))
            getattr(binary, "state_history").append(binary_state)
            getattr(binary, "event_history").append(None)
            getattr(binary, "mass_transfer_case_history").append(MT_case)

        if track_interpolation:
            if not self.flush_history:
                # history N is how much to look back in the history
                state1_hist = cf.check_state_of_star_history_array(
                    stars[0], N=length_hist, star_CO=stars_CO[0])
                getattr(stars[0], "state_history").extend(state1_hist)
                state2_hist = cf.check_state_of_star_history_array(
                    stars[1], N=length_hist, star_CO=stars_CO[1])
                getattr(stars[1], "state_history").extend(state2_hist)
                binary_state, binary_event, MT_case = (
                    cf.get_binary_state_and_event_and_mt_case_array(
                        binary, N=length_hist, verbose=self.verbose))
                getattr(binary, "state_history").extend(binary_state)
                getattr(binary, "event_history").extend(binary_event)
                getattr(binary, "mass_transfer_case_history").extend(MT_case)
                self.flush_entries = len_binary_hist   # this is needded!
                # this is to prevent the flushin of the initial value which is
                # appended twice
                if self.save_initial_conditions:
                    self.flush_entries += 1
            else:
                # the history is going to be flushed in self.stop
                # append None for a faster computation
                state1_hist = empy_h
                state2_hist = empy_h
                self.flush_entries = len_binary_hist
                # this is to prevent the flushin of the initial value which is
                # appended twice
                if self.save_initial_conditions:
                    self.flush_entries += 1
                binary_state = empy_h
                binary_event = empy_h
                MT_case = empy_h
                getattr(stars[0], "state_history").extend(state1_hist)
                getattr(stars[1], "state_history").extend(state2_hist)
                getattr(binary, "state_history").extend(binary_state)
                getattr(binary, "event_history").extend(binary_event)
                getattr(binary, "mass_transfer_case_history").extend(MT_case)

        if (star_2_CO or star_1_CO):
            # Updating Bondi-Hoyle accretion
            for k, star in enumerate(stars):
                if stars_CO[k]:
                    accretor = star
                    k_bh = k
                else:
                    donor = star
            key_bh = POSYDON_TO_MESA['star']['lg_mdot']+'_%d' % (k_bh+1)
            tmp_lg_mdot = np.log10(10**cb_bh[key_bh][-1] + cf.bondi_hoyle(
                binary, accretor, donor, idx=-1,
                wind_disk_criteria=True, scheme='Kudritzki+2000'))
            mdot_edd = cf.eddington_limit(binary, idx=-1)[0]
            if 10**tmp_lg_mdot > mdot_edd:
                tmp_lg_mdot = np.log10(mdot_edd)
            accretor.lg_mdot = tmp_lg_mdot
            if self.save_initial_conditions:
                mdot_history = np.array(cb_bh[key_bh])
                edd = cf.eddington_limit(binary, idx=len_binary_hist)[0]
                history_of_attribute = (np.log10(
                    10**cb_bh[key_bh][0] + cf.bondi_hoyle(
                        binary, accretor, donor, idx=len_binary_hist,
                        wind_disk_criteria=True, scheme='Kudritzki+2000')))
                if 10**history_of_attribute > edd:
                    history_of_attribute = np.log10(edd)
                accretor.lg_mdot_history.append(history_of_attribute)
            if track_interpolation:
                mdot_history = np.array(cb_bh[key_bh])
                # looping from range(-N,0) where 0 is excluded
                # note taht bondi_hoyle concatenates the current binary state
                # hence we loop one back range(-N-1,-1)
                tmp_h = [cf.bondi_hoyle(binary, accretor, donor, idx=i,
                                        wind_disk_criteria=True,
                                        scheme='Kudritzki+2000')
                         for i in range(-length_hist-1, -1)]
                tmp_edd = [cf.eddington_limit(binary, idx=i)[0]
                           for i in range(-length_hist-1, -1)]
                history_of_attribute = np.log10(10**cb_bh[key_bh][:-1] + tmp_h)
                history_of_attribute = [
                    y if x > y else x
                    for x, y in zip(history_of_attribute, tmp_edd)]
                accretor.lg_mdot_history.extend(history_of_attribute)

        # update post processed quanties
        key_post_processed = ['avg_c_in_c_core_at_He_depletion',
                              'co_core_mass_at_He_depletion',
                              'm_core_CE_1cent',
                              'm_core_CE_10cent',
                              'm_core_CE_30cent',
                              'm_core_CE_pure_He_star_10cent',
                              'r_core_CE_1cent',
                              'r_core_CE_10cent',
                              'r_core_CE_30cent',
                              'r_core_CE_pure_He_star_10cent']
        for k, star in enumerate(stars):
            for key in key_post_processed:
                setattr(star, key, cb.final_values['S%d_%s' % (k+1, key)])

        # update nearest neighbor core collapse quantites
        if interpolation_class != 'unstable_MT':
            for MODEL_NAME in MODELS.keys():
                for i, star in enumerate(stars):
                    if (not stars_CO[i] and
                        cb.final_values[f'S{i+1}_{MODEL_NAME}_CO_type'] != 'None'):
                        values = {}
                        for key in ['state', 'SN_type', 'f_fb', 'mass', 'spin',
                                    'm_disk_accreted', 'm_disk_radiated']:
                            if key == "state":
                                state = cb.final_values[f'S{i+1}_{MODEL_NAME}_CO_type']
                                values[key] = state
                            elif key == "SN_type":
                                values[key] = cb.final_values[f'S{i+1}_{MODEL_NAME}_{key}']
                            else:
                                values[key] = cb.final_values[f'S{i+1}_{MODEL_NAME}_{key}']
                        setattr(star, MODEL_NAME, values)
                    else:
                        setattr(star, key, None)

    def initial_final_interpolation(self, star_1_CO=False, star_2_CO=False):
        """Update the binary through initial-final interpolation."""
        # TODO: simplify verbosity of code
        binary = self.binary
        stars = [binary.star_1, binary.star_2]
        stars_CO = [star_1_CO, star_2_CO]
        fv = self.final_values

        for key in BINARYPROPERTIES:
            if POSYDON_TO_MESA['binary'][key] is None:
                setattr(self.binary, key, None)
            elif key == 'time':
                key_p = POSYDON_TO_MESA['binary'][key]
                current = getattr(self.binary, key)
                setattr(self.binary, key, current + fv[key_p])
            elif key == 'nearest_neighbour_distance':
                setattr(self.binary, key, ['None', 'None', 'None'])
            elif key in ['eccentricity', 'V_sys']:
                current = getattr(self.binary, key + '_history')[-1]
                setattr(self.binary, key, current)
            elif key in ['state', 'event', 'mass_transfer_case']:
                continue
            else:
                key_p = POSYDON_TO_MESA['binary'][key]
                setattr(self.binary, key, fv[key_p])
        for k, star in enumerate(stars):
            for key in STARPROPERTIES:
                if not stars_CO[k]:
                    if POSYDON_TO_MESA['star'][key] is None:
                        setattr(star, key, None)
                    elif key == 'mass':
                        key_p = 'star_%d_mass' % (k+1)
                        setattr(star, key, fv[key_p])
                    elif key == 'metallicity':
                        current = getattr(star, key)
                        setattr(star, key, current)
                    elif key == 'profile':
                        setattr(star, key, None)
                    elif key in ['lg_mdot', 'lg_system_mdot', 'lg_wind_mdot']:
                        key_p = POSYDON_TO_MESA['star'][key]+'_%d' % (k+1)
                        setattr(star, key, fv[key_p])
                    elif key == 'state':
                        continue
                    else:
                        key_p = 'S%d_' % (k+1)+POSYDON_TO_MESA['star'][key]
                        setattr(star, key, fv[key_p])
                else:
                    if POSYDON_TO_MESA['star'][key] is None:
                        setattr(star, key, None)
                    elif key == 'mass':
                        key_p = 'star_%d_mass' % (k+1)
                        setattr(star, key, fv[key_p])
                    elif key == 'spin':
                        current = getattr(star, 'spin')
                        setattr(star, key, current)
                    elif key == 'log_R':
                        mass = fv['star_%d_mass' % (k+1)]
                        st = infer_star_state(star_mass=mass, star_CO=True)
                        setattr(star, key, np.log10(CO_radius(mass, st)))
                    elif key == 'metallicity':
                        current = getattr(star, key)
                        setattr(star, key, current)
                    elif key == 'profile':
                        setattr(star, key, None)
                    elif key in ['lg_mdot', 'lg_system_mdot', 'lg_wind_mdot']:
                        key_p = POSYDON_TO_MESA['star'][key]+'_%d' % (k+1)
                        setattr(star, key, fv[key_p])
                    elif key == 'state':
                        continue
                    elif star.state == 'WD' and key in ['co_core_mass','he_core_mass','center_h1','center_he4','center_c12','center_n14','center_o16']:
                        continue
                    else:
                        setattr(star, key, None)
        # EXPERIMENTAL feature
        # infer stellar states
        interpolation_class = self.classes['interpolation_class']
        setattr(self.binary, f'interp_class_{self.grid_type}', interpolation_class)
        mt_history = self.classes['mt_history'] # mass transfer history (TF12 plot label)
        setattr(self.binary, f'mt_history_{self.grid_type}', mt_history)

        S1_state_inferred = cf.check_state_of_star(self.binary.star_1,
                                                   star_CO=star_1_CO)
        S2_state_inferred = cf.check_state_of_star(self.binary.star_2,
                                                   star_CO=star_2_CO)
        #S1_state_classified = self.classes['S1_state']
        #S2_state_classified = self.classes['S2_state']

        if interpolation_class != 'initial_MT':
            # DEBUG
            # if S1_state_inferred != S1_state_classified:
            #     warnings.warn('Inferred stellar state of star_1 %s is '
            #                   'different from classified state %s, note that'
            #                   'by default we use the inferred!' %
            #                   (S1_state_inferred,S1_state_classified))
            # if S2_state_inferred != S2_state_classified:
            #     warnings.warn('Inferred stellar state of star_2 %s is '
            #                   'different from classified state %s, note that'
            #                   'by default we use the inferred!' %
            #                   (S2_state_inferred,S2_state_classified))
            setattr(self.binary.star_1, 'state', S1_state_inferred)
            setattr(self.binary.star_2, 'state', S2_state_inferred)
        # else keep the current state

        binary_state, binary_event, MT_case = (
            cf.get_binary_state_and_event_and_mt_case(
                self.binary, interpolation_class, verbose=self.verbose))
        setattr(self.binary, 'state', binary_state)
        setattr(self.binary, 'event', binary_event)
        setattr(self.binary, 'mass_transfer_case', MT_case)

        if binary.state == 'initial_RLOF':
            return

        if (star_2_CO or star_1_CO):
            # Updating Bondi-Hoyle accretion
            for k, star in enumerate(stars):
                if stars_CO[k]:
                    accretor = star
                    k_bh = k
                else:
                    donor = star
            key_bh = POSYDON_TO_MESA['star']['lg_mdot']+'_%d' % (k_bh+1)
            tmp_lg_mdot = np.log10(
                10**fv[key_bh] + cf.bondi_hoyle(
                    binary, accretor, donor, idx=-1,
                    wind_disk_criteria=True, scheme='Kudritzki+2000'))
            mdot_edd = cf.eddington_limit(binary, idx=-1)[0]
            if 10**tmp_lg_mdot > mdot_edd:
                tmp_lg_mdot = np.log10(mdot_edd)
            setattr(accretor, 'lg_mdot', tmp_lg_mdot)
        # update post processed quanties
        key_post_processed = ['avg_c_in_c_core_at_He_depletion',
                              'co_core_mass_at_He_depletion',
                              'm_core_CE_1cent',
                              'm_core_CE_10cent',
                              'm_core_CE_30cent',
                              'm_core_CE_pure_He_star_10cent',
                              'r_core_CE_1cent',
                              'r_core_CE_10cent',
                              'r_core_CE_30cent',
                              'r_core_CE_pure_He_star_10cent']
        stars = [self.binary.star_1, self.binary.star_2]
        stars_CO = [star_1_CO, star_2_CO]
        for k, star in enumerate(stars):
            if not stars_CO[k]:
                for key in key_post_processed:
                    if (interpolation_class == 'unstable_MT'
                        and (key == 'avg_c_in_c_core_at_He_depletion'
                             or key == 'co_core_mass_at_He_depletion')):
                        setattr(star, key, None)
                    else:
                        setattr(star, key, fv['S%d_%s' % (k+1, key)])

        # update interpolated core collapse quantites
        if interpolation_class != 'unstable_MT':
            for MODEL_NAME in MODELS.keys():
                for i, star in enumerate(stars):
                    if (not stars_CO[i] and
                        self.classes[f'S{i+1}_{MODEL_NAME}_CO_type'] != 'None'):
                        values = {}
                        for key in ['state', 'SN_type', 'f_fb', 'mass', 'spin',
                                    'm_disk_accreted', 'm_disk_radiated']:
                            if key == "state":
                                state = self.classes[f'S{i+1}_{MODEL_NAME}_CO_type']
                                values[key] = state
                            elif key == "SN_type":
                                values[key] = self.classes[f'S{i+1}_{MODEL_NAME}_{key}']
                            else:
                                values[key] = fv[f'S{i+1}_{MODEL_NAME}_{key}']
                        setattr(star, MODEL_NAME, values)
                    else:
                        setattr(star, key, None)

    # STOPPING METHODS

    def stop_at_end(self,
                    binary,
                    property=None,
                    value=None,
                    star=None,
                    interpolate=None,
                    star_1_CO=False,
                    star_2_CO=False):
        """Update the binary event if max time has been exceeded."""
        if binary.properties.max_simulation_time - binary.time < 0.0:
            binary.event = 'MaxTime_exceeded'
            return

    def stop_at_condition(self,
                          binary,
                          property="time",
                          value=None,
                          star="star_1",
                          interpolate=False,
                          star_1_CO=False,
                          star_2_CO=False):
        if property in STARPROPERTIES:
            if star == 'star_1':
                current_property = getattr(binary.star_1, property)
                property_history = getattr(binary.star_1,
                                           property + "_history")
                np.array(property_history.append(current_property))
            elif star == 'star_2':
                current_property = getattr(binary.star_2, property)
                property_history = getattr(binary.star_2,
                                           property + "_history")
                np.array(property_history.append(current_property))
            else:
                raise ValueError(
                    'Star can only be "star_1" or "star_2", you passed {0}'.
                    format(star))
            # if property_history[-1] > value:
            #     binary.state += ' (OutsideGrid)'
            #     binary.event = 'END'
            #     return
            # if value > property_history[-1]:
            #     #t = delta_t[-1]
            #     #binary.state += ' (OutsideGrid)'
            #     #binary.event = 'END'
            #     binary.event = 'MaxTime_exceeded'
            #     return
            i = np.where(np.array(property_history) <= value)[0][-1]

        elif property in BINARYPROPERTIES:
            current_property = getattr(binary, property)
            property_history = getattr(binary, property + "_history")
            np.array(property_history.append(current_property))

            # if value > property_history[-1]:
            #     binary.event = 'MaxTime_exceeded'
            #     return
            i = np.where(np.array(property_history) <= value)[0][-1]

            # time at which to interpolate all quantities
            if property == 'time':
                t = value
        else:
            raise ValueError(
                'stop_at_condition does not support property = {0}'.format(
                    property))
        # interpolate to the desired value
        if interpolate:
            if property != 'time':
                raise ValueError('We currently support stop_at_conditon '
                                 'interpolation only for property=time!')
            # interpolate between i and i+1
            interpolated_quanties = {}
            interpolated_quanties['binary'] = {}
            interpolated_quanties['star_1'] = {}
            interpolated_quanties['star_2'] = {}
            for key in STARPROPERTIES:
                for star in ['star_1', 'star_2']:
                    star_obj = getattr(binary, star)
                    if key in ['state', 'profile', 'metallicity']:
                        interpolated_quanties[star][key] = getattr(
                            star_obj, key + "_history")[-1]
                    else:
                        current_property = getattr(star_obj, key)
                        property_history = getattr(star_obj, key + '_history')
                        np.array(property_history.append(current_property))

                        t_before = getattr(binary, 'time_history')[i]
                        t_after = getattr(binary, 'time_history')[i + 1]

                        v_before = property_history[i]
                        v_after = property_history[i + 1]
                        if v_before is None or v_after is None:
                            interpolated_quanties[star][key] = None
                            continue
                        # Debug
                        interpolated_quanties[star][
                            key] = self.interpolate_at_t(
                                t, t_before, t_after, v_before, v_after)
                        # except:
                        # # DEBUG
                        #     print('star', star, 'key', key)
                        #     print('time', t_before, t, t_after)
                        #     print('key', v_before,
                        #           interpolated_quanties[star][key], v_after)
            for key in BINARYPROPERTIES:
                if key in [
                        'state', 'event', 'mass_transfer_case'
                ]:
                # if key in [
                #         'state', 'event', 'mass_transfer_case',
                #         'nearest_neighbour_distance'
                # ]:
                    interpolated_quanties['binary'][key] = getattr(
                        binary, key + "_history")[-1]
                elif key == 'time':
                    interpolated_quanties['binary'][key] = t
                else:
                    current_property = getattr(binary, key)
                    property_history = getattr(binary, key + '_history')
                    np.array(property_history.append(current_property))

                    t_before = getattr(binary, 'time_history')[i]
                    t_after = getattr(binary, 'time_history')[i + 1]

                    v_before = property_history[i]
                    v_after = property_history[i + 1]

                    if v_before is None or v_after is None:
                        interpolated_quanties['binary'][key] = None
                        continue

                    # To deal with V_sys, nearest_neighbour, etc.
                    if hasattr(v_before, '__iter__'):
                        out = []
                        for j in range(len(v_before)):
                            out.append(self.interpolate_at_t(t,
                                                             t_before,
                                                             t_after,
                                                             v_before[j],
                                                             v_after[j]))
                        interpolated_quanties['binary'][key] = out
                    else:
                        interpolated_quanties['binary'][
                            key] = self.interpolate_at_t(t, t_before, t_after,
                                                         v_before, v_after)

        # delete the history above index i and reset all properties to index i
        for key in STARPROPERTIES:
            key_history = key + '_history'
            for star in ['star_1', 'star_2']:
                star_obj = getattr(binary, star)
                if interpolate:
                    # restore the star to interpolated state
                    setattr(star_obj, key, interpolated_quanties[star][key])
                    # detele history from i to the end
                    setattr(star_obj, key_history,
                            getattr(star_obj, key_history)[:i + 1])
                else:
                    # restore the star to state i
                    setattr(star_obj, key, getattr(star_obj, key_history)[i])
                    # detele history from i to the end
                    setattr(star_obj, key_history,
                            getattr(star_obj, key_history)[:i])

        for key in BINARYPROPERTIES:
            key_history = key + '_history'

            if interpolate:
                # restore the binary to interpolated state
                setattr(binary, key, interpolated_quanties['binary'][key])
                # detele history from i to the end
                setattr(binary, key_history,
                        getattr(binary, key_history)[:i + 1])
            else:
                # restore the binarya to state i
                setattr(binary, key, getattr(binary, key_history)[i])
                # detele history from i to the end
                setattr(binary, key_history, getattr(binary, key_history)[:i])

        # in case track_interpolation = False we will flush all the history
        if self.flush_history:
            # DEBUG
            # print('Flushing history between', self.flush_entries,
            #        len(getattr(binary, 'time_history')))
            if self.flush_entries is None:
                raise ValueError('flush_entries cannot be None!')
            for key in STARPROPERTIES:
                key_history = key + '_history'
                for star in ['star_1', 'star_2']:
                    star_obj = getattr(binary, star)
                    # detele last N entries from history
                    setattr(star_obj, key_history, getattr(
                        star_obj, key_history)[:self.flush_entries])
            for key in BINARYPROPERTIES:
                key_history = key + '_history'
                # detele last N entries from history
                setattr(binary, key_history,
                        getattr(binary, key_history)[:self.flush_entries])
            # Check star and binary state, event, MT case
            setattr(
                self.binary.star_1, 'state',
                cf.check_state_of_star(self.binary.star_1, star_CO=star_1_CO))
            setattr(
                self.binary.star_2, 'state',
                cf.check_state_of_star(self.binary.star_2, star_CO=star_2_CO))
            binary_state, binary_event, MT_case = (
                cf.get_binary_state_and_event_and_mt_case(
                    self.binary, verbose=self.verbose))
            setattr(self.binary, 'state', binary_state)
            setattr(self.binary, 'event', binary_event)
            setattr(self.binary, 'mass_transfer_case', MT_case)

    def interpolate_at_t(self, t, t_before, t_after, v_before, v_after):
        """Linear interpolation in time between two points.

        Parameters
        ----------
        t : type
            Description of parameter `t`.
        t_before : type
            Description of parameter `t_before`.
        t_after : type
            Description of parameter `t_after`.
        v_before : type
            Description of parameter `v_before`.
        v_after : type
            Description of parameter `v_after`.

        Returns
        -------
        type
            Description of returned object.

        """

        # Error handling
        if v_before == "None" or v_after == "None":
            return "None"

        slope = (v_after - v_before) / (t_after - t_before)
        v_t = (t - t_before) * slope + v_before

        return v_t


class MS_MS_step(MesaGridStep):
    """Class for performing the MESA step for a MS-MS binary."""

    def __init__(self, metallicity=1., grid_name=None, *args, **kwargs):
        """Initialize a MS_MS_step instance."""
        self.grid_type = 'HMS_HMS'
        self.interp_in_q = True
        if grid_name is None:
            metallicity = convert_metallicity_to_string(metallicity)
            grid_name = 'HMS-HMS/' + metallicity + '_Zsun.h5'
        super().__init__(metallicity=metallicity,
                         grid_name=grid_name,
                         *args, **kwargs)
        # special stuff for my step goes here
        # If nothing to do, no init necessary

        # load grid boundaries
        self.m1_min = min(self._psyTrackInterp.grid.initial_values['star_1_mass'])
        self.m1_max = max(self._psyTrackInterp.grid.initial_values['star_1_mass'])
        self.m2_min = min(self._psyTrackInterp.grid.initial_values['star_2_mass'])
        self.m2_max = max(self._psyTrackInterp.grid.initial_values['star_2_mass'])
        self.q_min = 0.05 # can be computed m2_min/m1_min
        self.q_max = 1. # note that for MESA stability we actually run q_max = 0.99
        self.p_min = min(self._psyTrackInterp.grid.initial_values['period_days'])
        self.p_max = max(self._psyTrackInterp.grid.initial_values['period_days'])

    def __call__(self, binary):
        """Apply the MS-MS step on a BinaryStar."""
        # grid set up, both stars are NOT CO
        self.star_1_CO = False
        self.star_2_CO = False
        # check binary is ready before calling the step
        self.binary = binary
        state_1 = self.binary.star_1.state
        state_2 = self.binary.star_2.state
        event = self.binary.event
        m1 = self.binary.star_1.mass
        m2 = self.binary.star_2.mass
        mass_ratio = m2/m1
        p = self.binary.orbital_period
        # check if the binary is in the grid
        if (state_1 == 'H-rich_Core_H_burning' and
            state_2 == 'H-rich_Core_H_burning' and
            event == 'ZAMS' and
            self.m1_min <= m1 <= self.m1_max and
            np.max([self.q_min, 0.5/m1]) <= mass_ratio <= self.q_max and
            self.p_min <= p <= self.p_max):
            self.flip_stars_before_step = False
            super().__call__(self.binary)
        # binary in grid but masses flipped
        elif (state_1 == 'H-rich_Core_H_burning' and
              state_2 == 'H-rich_Core_H_burning' and
              event == 'ZAMS' and
              self.m1_min <= m2 <= self.m1_max and
              np.max([self.q_min, 0.5/m1]) <= 1./mass_ratio <= self.q_max and
              self.p_min <= p <= self.p_max):
            self.flip_stars_before_step = True
            super().__call__(self.binary)
        # redirect if outside grid for period
        elif (state_1 == 'H-rich_Core_H_burning' and
              state_2 == 'H-rich_Core_H_burning' and
              event == 'ZAMS' and
              p > self.p_max):
            self.binary.event = 'redirect_from_ZAMS'
            return
        # redirect if period smaller than the minimum period
        elif (state_1 == 'H-rich_Core_H_burning' and
              state_2 == 'H-rich_Core_H_burning' and
              event == 'ZAMS' and
              p < self.p_min):
            self.binary.event = 'redirect_from_ZAMS'
            return
        # outside the mass grid for m1
        elif (state_1 == 'H-rich_Core_H_burning' and
              state_2 == 'H-rich_Core_H_burning' and
              event == 'ZAMS' and
              self.p_min <= p <= self.p_max and
              (m1 < self.m1_min or m1 > self.m1_max)):
            set_binary_to_failed(self.binary)
            raise ValueError(f'The mass of m1 ({m1}) is outside the grid,'
                             'while the period is inside the grid.')
        # outside the mass grid for m2
        # because m2_min is 0.5 Msun or from q_min, the minimum mass ratio is either
        # q_min or 0.5/m1.
        elif (state_1 == 'H-rich_Core_H_burning' and
              state_2 == 'H-rich_Core_H_burning' and
              event == 'ZAMS' and
              self.p_min <= p <= self.p_max and
              (mass_ratio < np.max([self.q_min, 0.5/m1]) or mass_ratio > self.q_max)):
            set_binary_to_failed(self.binary)
            raise ValueError(f'The mass of m2 ({m2}) is outside the grid,'
                             'while the period is inside the grid.')
        # redirect if CC1
        elif (state_1 == 'H-rich_Central_C_depletion'):
            self.binary.event = 'CC1'
            return
        # redirect if CC2
        elif (state_2 == 'H-rich_Central_C_depletion'):
            self.binary.event = 'CC2'
            return
        else:
            set_binary_to_failed(self.binary)
            raise ValueError('The star_1.state = %s, star_2.state = %s, '
                             'binary.event = %s and not H-rich_Core_H_burning '
                             '- H-rich_Core_H_burning - * - ZAMS'
                             % (state_1, state_2, event))


class CO_HMS_RLO_step(MesaGridStep):
    """Class for performing the MESA step for a CO-HMS_RLO binary."""

    def __init__(self, metallicity=1., grid_name=None, *args, **kwargs):
        """Initialize a CO_HMS_RLO_step instance."""
        self.grid_type = 'CO_HMS_RLO'
        self.interp_in_q = False
        if grid_name is None:
            metallicity = convert_metallicity_to_string(metallicity)
            grid_name = 'CO-HMS_RLO/' + metallicity + '_Zsun.h5'
        super().__init__(metallicity=metallicity,
                         grid_name=grid_name,
                         *args, **kwargs)
        # special stuff for my step goes here
        # If nothing to do, no init necessary

        # load grid boundaries
        self.m1_min = min(self._psyTrackInterp.grid.initial_values['star_1_mass'])
        self.m1_max = max(self._psyTrackInterp.grid.initial_values['star_1_mass'])
        self.m2_min = min(self._psyTrackInterp.grid.initial_values['star_2_mass'])
        self.m2_max = max(self._psyTrackInterp.grid.initial_values['star_2_mass'])
        self.p_min = min(self._psyTrackInterp.grid.initial_values['period_days'])
        self.p_max = max(self._psyTrackInterp.grid.initial_values['period_days'])

    def __call__(self, binary):
        """Evolve a binary using the MESA step."""
        # grid set up assume CO is star_2
        self.star_1_CO = False
        self.star_2_CO = True
        # check binary is ready before calling the step
        self.binary = binary
        event = self.binary.event
        state = binary.state
        state_1 = binary.star_1.state
        state_2 = binary.star_2.state
        m1 = self.binary.star_1.mass
        m2 = self.binary.star_2.mass
        p = self.binary.orbital_period
        ecc = self.binary.eccentricity

        # TODO: import states from flow_chart.py
        FOR_RLO_STATES = ["H-rich_Core_H_burning",
                          "H-rich_Shell_H_burning",
                          "H-rich_Core_He_burning",
                          "H-rich_Central_He_depleted",
                          "H-rich_Core_C_burning",
                          "stripped_He_Core_H_burning",
                          "H-rich_Central_C_depletion",  # filtered out below
                          "H-rich_non_burning"]

        # check the star states
        # TODO: import states from flow_chart.py
        if (state_2 in ["WD", "NS", "BH"]
                and (state_1 in FOR_RLO_STATES) and event == "oRLO1"):
            self.flip_stars_before_step = False
            # catch and redirect double core collapse, this happens if q=1:
            if self.binary.star_1.state == 'H-rich_Central_C_depletion':
                self.binary.event = 'CC1'
                return
        # TODO: import states from flow_chart.py
        elif (state_1 in ["WD", "NS", "BH"] and (state_2 in FOR_RLO_STATES)
                and event == "oRLO2"):
            self.flip_stars_before_step = True
            # catch and redirect double core collapse, this happens if q=1:
            if self.binary.star_2.state == 'H-rich_Central_C_depletion':
                self.binary.event = 'CC2'
                return
        else:
            raise ValueError(
                'The star_1.state = %s, star_2.state = %s, binary.state = %s, '
                'binary.event = %s and not CO - HMS - oRLO1/oRLO2!'
                % (state_1, state_2, state, event))
        # redirect if outside grids
        if ((not self.flip_stars_before_step and
            self.m1_min <= m1 <= self.m1_max and
            self.m2_min <= m2 <= self.m2_max and
            self.p_min <= p <= self.p_max and
            ecc == 0.) or (self.flip_stars_before_step and
            self.m1_min <= m2 <= self.m1_max and
            self.m2_min <= m1 <= self.m2_max and
            self.p_min <= p <= self.p_max and
            ecc == 0.)):
            super().__call__(self.binary)
        # period inside the grid, but m1 outside the grid
        elif ((not self.flip_stars_before_step and
               self.p_min <= p <= self.p_max and
               (m1 < self.m1_min or m1 > self.m1_max)
               )):
            set_binary_to_failed(self.binary)
            raise ValueError(f'The mass of m1 ({m1}) is outside the grid,'
                                'while the period is inside the grid.')

        # period inside the grid, but m2 outside the grid
        elif ((not self.flip_stars_before_step and
               self.p_min <= p <= self.p_max and
               (m2 < self.m2_min or m2 > self.m2_max)
               )):
            set_binary_to_failed(self.binary)
            raise ValueError(f'The mass of m2 ({m2}) is outside the grid,'
                                'while the period is inside the grid.')

        else:
            if len(self.binary.state_history) > 2:
                if self.binary.state_history[-2] == 'detached':
                    set_binary_to_failed(self.binary)
                    raise ValueError('CO_HMS_RLO binary outside grid and coming from detached')
            self.binary.state = "detached"
            self.binary.event = "redirect_from_CO_HMS_RLO"
            return

class CO_HeMS_RLO_step(MesaGridStep):
    """Class for performing the MESA step for a CO-HeMS_RLO binary."""

    def __init__(self, metallicity=1., grid_name=None, *args, **kwargs):
        """Initialize a CO_HeMS_RLO_step instance."""
        self.grid_type = 'CO_HeMS_RLO'
        self.interp_in_q = False
        if grid_name is None:
            metallicity = convert_metallicity_to_string(metallicity)
            grid_name = 'CO-HeMS_RLO/' + metallicity + '_Zsun.h5'
        super().__init__(metallicity=metallicity,
                         grid_name=grid_name,
                         *args, **kwargs)
        # special stuff for my step goes here
        # If nothing to do, no init necessary

        # load grid boundaries
        self.m1_min = min(self._psyTrackInterp.grid.initial_values['star_1_mass'])
        self.m1_max = max(self._psyTrackInterp.grid.initial_values['star_1_mass'])
        self.m2_min = min(self._psyTrackInterp.grid.initial_values['star_2_mass'])
        self.m2_max = max(self._psyTrackInterp.grid.initial_values['star_2_mass'])
        self.p_min = min(self._psyTrackInterp.grid.initial_values['period_days'])
        self.p_max = max(self._psyTrackInterp.grid.initial_values['period_days'])

    def __call__(self, binary):
        """Evolve a binary using the MESA step."""
        # grid set up assume CO is star_2
        self.star_1_CO = False
        self.star_2_CO = True
        # check binary is ready before calling the step
        self.binary = binary
        event = self.binary.event
        state = binary.state
        state_1 = binary.star_1.state
        state_2 = binary.star_2.state
        m1 = self.binary.star_1.mass
        m2 = self.binary.star_2.mass
        p = self.binary.orbital_period
        ecc = self.binary.eccentricity

        # TODO: import states from flow_chart.py
        CO_He_STATES = [
            'stripped_He_Core_He_burning',
            'stripped_He_Shell_He_burning',
            'stripped_He_Central_He_depleted',
            'stripped_He_Central_C_depletion',  # filtered out below
            # include systems that are on the brink of He exhaustion
            'stripped_He_non_burning',
            # include systems post CE with core_definition_H_fraction=0.1
            'H-rich_non_burning'
                        ]

        # check the star states
        # TODO: import states from flow_chart.py
        if (state_2 in ["WD", "NS", "BH"]
                and (state_1 in CO_He_STATES) and event == "oRLO1"):
            self.flip_stars_before_step = False
            # catch and redirect double core collapse, this happens if q=1:
            if self.binary.star_1.state == 'stripped_He_Central_C_depletion':
                self.binary.event = 'CC1'
                return
        # TODO: import states from flow_chart.py
        elif (state_1 in ["WD", "NS", "BH"] and (state_2 in CO_He_STATES)
                and event == "oRLO2"):
            self.flip_stars_before_step = True
            # catch and redirect double core collapse, this happens if q=1:
            if self.binary.star_2.state == 'stripped_He_Central_C_depletion':
                self.binary.event = 'CC2'
                return
        else:
            raise ValueError(
                'The star_1.state = %s, star_2.state = %s, binary.state = %s, '
                'binary.event = %s and not CO - HeMS - oRLO1/oRLO2!'
                % (state_1, state_2, state, event))
        # redirect if outside grids
        if ((not self.flip_stars_before_step and
            self.m1_min <= m1 <= self.m1_max and
            self.m2_min <= m2 <= self.m2_max and
            self.p_min <= p <= self.p_max and
            ecc == 0.) or (self.flip_stars_before_step and
            self.m1_min <= m2 <= self.m1_max and
            self.m2_min <= m1 <= self.m2_max and
            self.p_min <= p <= self.p_max and
            ecc == 0.)):
            super().__call__(self.binary)
        # period inside the grid, but m1 outside the grid
        elif ((not self.flip_stars_before_step and
               self.p_min <= p <= self.p_max and
                (m1 < self.m1_min or m1 > self.m1_max)
                )):
            set_binary_to_failed(self.binary)
            raise ValueError(f'The mass of m1 ({m1}) is outside the grid,'
                                'while the period is inside the grid.')
        # period inside the grid, but m2 outside the grid
        elif ((not self.flip_stars_before_step and
               self.p_min <= p <= self.p_max and
               (m2 < self.m2_min or m2 > self.m2_max)
               )):
            set_binary_to_failed(self.binary)
            raise ValueError(f'The mass of m2 ({m2}) is outside the grid,'
                                'while the period is inside the grid.')

        else:
            if len(self.binary.state_history) > 2:
                if self.binary.state_history[-2] == 'detached':
                    set_binary_to_failed(self.binary)
                    raise ValueError('CO_HeMS_RLO binary outside grid and coming from detached')

            self.binary.state = "detached"
            self.binary.event = "redirect_from_CO_HeMS_RLO"
            return

class CO_HeMS_step(MesaGridStep):
    """Class for performing the MESA step for a CO-HeMS binary."""

    def __init__(self, metallicity=1., grid_name=None, *args, **kwargs):
        """Initialize a CO_HeMS_step instance."""
        self.grid_type = 'CO_HeMS'
        self.interp_in_q = False
        if grid_name is None:
            metallicity = convert_metallicity_to_string(metallicity)
            grid_name = 'CO-HeMS/' + metallicity + '_Zsun.h5'
        super().__init__(metallicity=metallicity,
                         grid_name=grid_name,
                         *args, **kwargs)
        # special stuff for my step goes here
        # If nothing to do, no init necessary

        # load grid boundaries
        self.m1_min = min(self._psyTrackInterp.grid.initial_values['star_1_mass'])
        self.m1_max = max(self._psyTrackInterp.grid.initial_values['star_1_mass'])
        self.m2_min = min(self._psyTrackInterp.grid.initial_values['star_2_mass'])
        self.m2_max = max(self._psyTrackInterp.grid.initial_values['star_2_mass'])
        self.p_min = min(self._psyTrackInterp.grid.initial_values['period_days'])
        self.p_max = max(self._psyTrackInterp.grid.initial_values['period_days'])

    def __call__(self, binary):
        """Apply the CO_HeMS step to a BinaryStar object."""
        # grid set up assume CO is star_2
        self.star_1_CO = False
        self.star_2_CO = True
        # check binary is ready before calling the step
        self.binary = binary
        state_1 = self.binary.star_1.state
        state_2 = self.binary.star_2.state
        event = self.binary.event
        m1 = self.binary.star_1.mass
        m2 = self.binary.star_2.mass
        p = self.binary.orbital_period
        ecc = self.binary.eccentricity

        # TODO: import states from flow_chart.py
        CO_He_STATES = [
            'stripped_He_Core_He_burning',
            'stripped_He_Shell_He_burning',
            'stripped_He_Central_He_depleted',
            'stripped_He_Central_C_depletion',  # filtered out below
            # include systems that are on the brink of He exhaustion
            'stripped_He_non_burning',
            # include systems post CE with core_definition_H_fraction=0.1
            'H-rich_non_burning'
                        ]
        # TODO: import states from flow_chart.py
        if (state_2 in ['WD', 'NS', 'BH']
                and state_1 in CO_He_STATES and event is None):
            self.flip_stars_before_step = False
            # catch and redirect double core collapse, this happens if q=1:
            if self.binary.star_1.state == 'stripped_He_Central_C_depletion':
                self.binary.event = 'CC1'
                # REMOVED assume circularisation after first CC
                # new_separation = self.binary.separation*(
                #     1.-self.binary.eccentricity**2)
                # self.binary.separation = new_separation
                # self.binary.orbital_period = orbital_period_from_separation(
                #     new_separation, m1, m2)
                # self.binary.eccentricity = 0.
                return
        # TODO: import states from flow_chart.py
        elif (state_1 in ['WD', 'NS', 'BH']
                and state_2 in CO_He_STATES and event is None):
            self.flip_stars_before_step = True
            # catch and redirect double core collapse, this happens if q=1:
            if self.binary.star_2.state == 'stripped_He_Central_C_depletion':
                self.binary.event = 'CC2'
                return
        else:
            raise ValueError(
                'The star_1.state = %s, star_2.state = %s, binary.event = %s '
                'not supported by CO - HeMS grid!' % (state_1, state_2, event))

        # redirect if outside grids
        # remember that in MESA the CO object is star_2
        if ((not self.flip_stars_before_step and
            self.m1_min <= m1 <= self.m1_max and
            self.m2_min <= m2 <= self.m2_max and
            self.p_min <= p <= self.p_max and
            ecc == 0.) or (self.flip_stars_before_step and
            self.m1_min <= m2 <= self.m1_max and
            self.m2_min <= m1 <= self.m2_max and
            self.p_min <= p <= self.p_max and
            ecc == 0.)):
            super().__call__(binary)
        else:
            if len(self.binary.state_history) > 2:
                if self.binary.state_history[-2] == 'detached':
                    set_binary_to_failed(self.binary)
                    raise ValueError('CO_HeMS binary outside grid and coming from detached')
            self.binary.event = 'redirect_from_CO_HeMS'
            return
