"""The definition of the `Sampler` class in PSY-CRIS."""


__authors__ = [
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
    "Scott Coughlin <scottcoughlin2014@u.northwestern.edu>",
]


import sys
import time
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
from scipy.spatial.distance import pdist


class Sampler:
    """Class implementing PTMCMC and MCMC for PSY-CRIS algorith.

    Modular implementation of PTMCMC and MCMC designed to
    implement the PSY-CRIS algorithm of sampling points
    in a target distribution constructed with a Classifier
    and Regressor. After a posterior is generated, methods
    in this class are also used to downsample.
    """

    def __init__(self, classifier=None, regressor=None):
        """Initialize the sampler.

        Parameters
        ----------
        classifier : instance of <class, Classifier>
            A trained classifier object.
        regressor : instance of <class, Regressor>, optional
            A trained regressor object.

        """
        self._Classifier_ = classifier
        self._Regressor_ = regressor

        if (self._Classifier_ is not None) or (self._Regressor_ is not None):
            if self._Classifier_ is not None:
                self._TableData_ = self._Classifier_._TableData_
            else:
                self._TableData_ = self._Regressor_._TableData_
            # Find the bounds of the walker - should be a TableData attribute
            self._max_vals_ = self._TableData_._max_input_vals
            self._min_vals_ = self._TableData_._min_input_vals

        # You can save chains_history in here
        self._chain_step_hist_holder_ = OrderedDict()
        # Not fully implemented yet I'm pretty sure....
        self._MAX_APC_str_ = []

    def TD_2d_analytic(self, name, args, **kwargs):
        r"""2-dimensional analytic target distribution for testing MCMC/PTMCMC.

        The function:
        $\frac{16}{3\pi} \left( \exp\left[-\mu^2 - (9 + 4\mu^2 + 8\nu)^2\right]
        + \frac{1}{2} \exp\left[- 8 \mu^2 - 8 (\nu-2)^2\right] \right)$

        Parameters
        ----------
        name : str
            Name of algorithm to use. For this method None.
        args : array
            2D location to get the value of the function.
        **kwargs
            Kwargs for more complex target distributions.

        Returns
        -------
        array or float

        """
        mu, nu = args
        arg1 = -(mu ** 2) - (9 + 4 * mu ** 2 + 8 * nu) ** 2
        arg2 = -8 * mu ** 2 - 8 * (nu - 2) ** 2
        return (16) / (3 * np.pi) * (np.exp(arg1) + 0.5 * np.exp(arg2))

    def get_TD_classification_data(self, *args, **kwargs):
        """Get target-distribution classification data.

        Calculate terms relevant for creating target distributions
        with classification terms.

        Parameters
        ----------
        classifier_name : str
            Trained classifier name to use for predictions.
        position : array
            Position in parameter space to eval
        **kwargs
            TD_verbose : bool
                Print useful output

        Returns
        -------
        max_probs : array
            Maximum probabilities at each query point
        position : array
            Position in parameter space being queried
        cls_key : array
            Classification key predicted for each query position

        """
        classifier_name, position = args
        position = np.array(position)
        if position.ndim == 1:
            position = position[np.newaxis, :]
        TD_verbose = kwargs.get("TD_verbose", False)
        normalized_probs, where_not_nan = self._Classifier_.return_probs(
            classifier_name, position, verbose=TD_verbose
        )
        where_nan = np.array(
            [i for i in range(len(position)) if i not in where_not_nan]
        )
        max_probs = np.max(normalized_probs, axis=1)
        cls_key = np.argmax(normalized_probs, axis=1)
        position = np.array(position)
        if len(where_nan) > 0:
            max_probs = np.zeros(position.shape[0]) + 1e-16
            cls_key = [None] * position.shape[0]
            if TD_verbose:
                print('\t', max_probs, position, cls_key)
            return max_probs, position, cls_key
        else:
            max_probs = np.where(max_probs == 1, 1-1e-16, max_probs)
            if TD_verbose:
                print('\t', max_probs, position, cls_key)
            return max_probs, position, cls_key

    def TD_classification(self, classifier_name, position, **kwargs):
        r"""Target distribution using classification.

        $f(x) = 1 - max[P_{\rm class}(x)]$

        Parameters
        ----------
        classifier_name : str
            String to specify the trained classification algorithm to use.
        position : array
            Single location in parameter space for the target distribution to
            be evaluated at.
        **kwargs
            TD_BETA : float
                Exponent of target distribution - $f(x)^{\rm TD_BETA}$
                Used for smoothing or sharpening.
            TD_verbose : bool
                Extra print output every method call.

        Returns
        -------
        array
            If classification probability is Nan: f(x) = 1E-16

        """
        # TD_verbose = kwargs.get("TD_verbose", False)
        TD_BETA = kwargs.get("TD_BETA", 1.0)

        max_probs, pos, cls_keys = self.get_TD_classification_data(
            classifier_name, position, **kwargs
        )

        theoretical_max_TD_cls_term = 1 - 1 / self._TableData_.num_classes
        return ((1 - max_probs) * 1 / theoretical_max_TD_cls_term) ** (TD_BETA)

    def TD_classification_regression(self, names, args, **kwargs):
        r"""Target distribution using both classification & regression.

        Classification: $1 - max[P_{\rm class}(x)]$
        Regression: $ A_0 \log( A_1* abs( max[APC_n [loc]]) + 1 )$

        Parameters
        ----------
        names : list like
            Iterable containing the two strings specifying the classification
            and regression algorithm to use.
        args : array
            Position in parameter space to evaluate the target distribution at.
        **kwargs
            TD_A1 : float, optional
                Scaling factor inside the Log regression error term.
                (Default = 0.5)
            TD_TAU : float, optional
                Relative weight of classification to regression term.
                (Default = 0.5)
            TD_BETA : float, optional
                Exponent of the entire target distribution. Used
                for smoothing or sharpening the distribution. Default is 1.
            TD_verbose : bool, optional
                Print more diagnostic information.

        Rreturns
        --------
        array
        """
        normalized_probs, where_not_nan = self._Classifier_.return_probs(
            names[0], args, verbose=False
        )
        max_probs = np.max(normalized_probs, axis=1)
        pred_class_ids = np.argmax(normalized_probs, axis=1)
        cls_key = [self._Classifier_.class_id_mapping[i]
                   for i in pred_class_ids]

        theoretical_max_TD_cls_term = (
            1 - 1 / self._Classifier_._TableData_.num_classes)

        if max_probs == 1:
            max_probs = 1 - 1e-16
        classification_term = (1 - max_probs) * 1 / theoretical_max_TD_cls_term

        if len(where_not_nan) != len(normalized_probs):
            return 1e-16
        else:
            if isinstance(
                self._Regressor_.regr_dfs_per_class[cls_key[0]], pd.DataFrame
            ):
                max_APC, which_col_max = self._Regressor_.get_max_APC_val(
                    names[1], cls_key[0], args
                )
                self._MAX_APC_str_.append(which_col_max)

                if self._Regressor_.abs_max_APC is None:
                    raise ValueError("No max APC value found in TableData...")

                A1 = kwargs.get("TD_A1", 0.5)
                scaling_log_func = lambda A1, x: np.log10(A1 * np.abs(x) + 1)
                A0 = 1 / scaling_log_func(A1, self._Regressor_.abs_max_APC)
                regression_term = A0 * scaling_log_func(A1, max_APC)
            else:
                regression_term = 1e-16

        TD_TAU = kwargs.get("TD_TAU", 0.5)
        TD_BETA = kwargs.get("TD_BETA", 1.0)
        if kwargs.get("TD_verbose", False):
            print("TD_TAU: {0} | TD_BETA: {1}".format(TD_TAU, TD_BETA))
        return (TD_TAU * classification_term
                + (1 - TD_TAU) * regression_term) ** TD_BETA

    def save_chain_step_history(self, key, chain_step_history,
                                overwrite=False):
        """Save PTMCMC output chain_step_history inside the sampler object."""
        if not (key not in self._chain_step_hist_holder_.keys() or overwrite):
            raise Exception(
                "\nYou are about to overwrite an existing element in '{0}'\n\n"
                "\tUse the option 'overwrite=True' to reassign.".format(key))

        self._chain_step_hist_holder_[key] = chain_step_history
        print("Saved chain to '{0}'.".format(key))

    def get_saved_chain_step_history(self, key, return_all=False):
        """Return the saved chain step history."""
        if return_all:
            return self._chain_step_hist_holder_
        else:
            return self._chain_step_hist_holder_[key]

    def run_PTMCMC(self, T_max, N_tot, target_dist, classifier_name,
                   init_pos=None, N_draws_per_swap=3, c_spacing=1.2,
                   alpha=None, upper_limit_reject=1e5, verbose=False,
                   trace_plots=False, **TD_kwargs):
        """Run a Paralel Tempered MCMC with user-specified target distribution.

        Calls the method `run_MCMC`.

        Parameters
        ----------
        T_max : float
            Sets the maximum temperature MCMC in the chain.
        N_tot : int
            The total number of iterations for the PTMCMC.
        target_dist : callable
            The target distribution to sample.
            Must take arguments (method_name, location_to_eval)
            (A 2D analytic function is provided - analytic_target_dist)
        classifier_name : str, list
            A single string or list of strings specifying the interpolator to
            use for classification or classification & regression respectively.
        init_pos : array
            Initial position of walkers in each axis. Default is the median of
            the input data in TableData.
        N_draws_per_swap : int, optional
            Number of draws to perform for each MCMC before swap proposals.
        c_spacing : float, optional
            Sets the spacing of temperatures in each chain.
            T_{i+1} = T_{i}^{1/c}, range: [T_max , T=1]
        alpha : float, optional
            Sets the standard deviation of steps taken by the walkers. Default
            is 1/5 the range of training data from TableData.
        upper_limit_reject : float, optional
            Sets the upper limit of rejected points.
        verbose : bool, optional
            Useful print statements during execution.

        Returns
        -------
        chain_step_history : dict
            Hold the step history for every chain. Keys are integers that range
            from 0 (max T) to the total number of chains -1 (min T).
        T_list : array
            Array filled with the temperatures of each chain from max to min.

        Notes
        -----
        There is a zero prior on the PTMCMC outside the range of training data.

        """
        # create list of Temperatures: T_{i+1}=T_{i}^{1/c}, range: [T_max, T=1]
        T_list = [T_max]
        while T_list[-1] > 1.3:
            T_list.append(T_list[-1] ** (1 / c_spacing))
        T_list.append(1)
        T_list = np.array(T_list)

        num_chains = len(T_list)
        if verbose:
            print("Num chains: {0}\nTemperatures: {1}\n".format(num_chains,
                                                                T_list))

        # data storage
        chain_holder = OrderedDict()
        for i in range(len(T_list)):
            chain_holder[i] = []

        N_loops = int(N_tot / N_draws_per_swap)

        # Initial conditions for all links in chain
        # ADD: init_pos can be a unitless position in the range of the axes of
        # the data. This change should also be applied to alpha - user just
        # gives num(0, 1]
        if init_pos is None:
            init_pos = np.median(self._TableData_._input_.values, axis=0)
        if not isinstance(init_pos, np.ndarray):
            init_pos = np.array(init_pos)
        if init_pos.ndim > 1:
            raise ValueError(
                "init_pos has {0} dimensions, must be one dimensional.".format(
                    init_pos.ndim
                )
            )

        if alpha is None:
            alpha = [abs(max_val-min_val)/5 for min_val, max_val in
                     zip(self._min_vals_, self._max_vals_)]

        # start chains in the same position
        this_iter_step_loc = [init_pos] * num_chains

        # Accept ratio tracker
        total_acc = np.zeros(num_chains)
        total_rej = np.zeros(num_chains)
        acc_ratio_holder = np.zeros(num_chains)

        start_time = time.time()
        for counter in range(N_loops):
            # Number of draws before swap
            N_draws = N_draws_per_swap

            last_step_holder = []
            for i in range(num_chains):
                # Run MCMC as f(T) N_draw times
                step_history = [this_iter_step_loc[i]]
                steps, acc, rej = self.run_MCMC(
                    N_draws,
                    alpha,
                    step_history,
                    target_dist,
                    classifier_name,
                    T=T_list[i],
                    upper_limit_reject=upper_limit_reject,
                    **TD_kwargs
                )
                # save 'current' params for each T
                last_step_holder.append(steps[-1])
                total_acc[i] += acc
                total_rej[i] += rej
                acc_ratio_holder[i] = total_acc[i] / (total_acc[i]
                                                      + total_rej[i])
                # data storage
                chain_holder[i].append(np.array(steps))

            if verbose:
                # useful output during the PTMCMC
                num_bars = 20
                how_close = int((counter / (N_loops - 1)) * num_bars)
                progress_bar = (
                    "|"
                    + how_close * "="
                    + ">"
                    + abs(num_bars - how_close) * " "
                    + "|"
                    + "{0:.1f}%".format(counter / (N_loops - 1) * 100)
                )
                b = ("num_acc/total: Tmax={0:.4f} Tmin={1:.4f} loop #{2}, {3}".
                     format(acc_ratio_holder[0], acc_ratio_holder[-1],
                            counter, progress_bar))
                sys.stdout.write("\r" + b)

            # Calc H to see if chains SWAP
            accept = 0
            reject = 0
            for i in range(len(T_list) - 1):
                args_i = last_step_holder[i]
                args_i_1 = last_step_holder[i + 1]

                top = (target_dist(classifier_name, args_i_1, **TD_kwargs))**(
                    1 / T_list[i]
                ) * (target_dist(classifier_name, args_i, **TD_kwargs)) ** (
                    1 / T_list[i + 1]
                )
                bot = (target_dist(classifier_name, args_i, **TD_kwargs)) ** (
                    1 / T_list[i]
                ) * (target_dist(classifier_name, args_i_1, **TD_kwargs)) ** (
                    1 / T_list[i + 1]
                )

                # can get div 0 errors when using linear because of nans
                if bot == 0:
                    ratio = 0
                else:
                    ratio = top / bot

                # inter-chain transition probability
                H = min(1, ratio)

                chance = np.random.uniform(low=0, high=1)
                if chance <= H:
                    accept += 1
                    # SWAP the params between the two chains
                    last_step_holder[i] = args_i_1
                    last_step_holder[i + 1] = args_i
                else:
                    reject += 1
            # print(accept, reject); print(last_step_holder)

            # Update current params (could be swapped)
            # to read into MCMC on next iteration!!!
            this_iter_step_loc = last_step_holder

        chain_step_history = OrderedDict()
        for pos, steps in chain_holder.items():
            chain_step_history[pos] = np.concatenate(steps)

        if verbose:
            print("\nLength of chains: \n{0}".format(
                np.array([len(chain_step_history[i])
                          for i in range(num_chains)])))
            fin_time_s = time.time() - start_time
            print(
                "Finished in {0:.2f} seconds, {1:.2f} minutes.".format(
                    fin_time_s, fin_time_s / 60
                )
            )
            if trace_plots:
                self.make_trace_plot(chain_step_history, T_list,
                                     0, save_fig=False)
                self.make_trace_plot(chain_step_history, T_list,
                                     num_chains - 1, save_fig=False)

        return chain_step_history, T_list

    def run_MCMC(self, N_trials, alpha, step_history, target_dist,
                 classifier_name, T=1, upper_limit_reject=1e4, **TD_kwargs):
        """Run a Markov chain Monte Carlo given a target distribution.

        Parameters
        ----------
        N_trials : int
            Number of proposals or trial steps to take before stopping.
        alpha : float
            Related to the step size of the MCMC walker.
            Defines the standard deviation of a zero mean normal
            from which the step is randomly drawn.
        step_history : list
            Initial starting location in parameter space.
            Could contain an arbitrary number of previous steps
            but a walker will start at the last step in the list.
        targe_dist : callable
            The target distribution to sample.
            Must take arguments ( method_name, element_of_step_history )
            (A 2D analytic function is provided - TD_2d_analytic)
        classifier_name : str
            Name of interpolation technique used in the target_dist.
        T : float, optional
            Temperature of the MCMC.
        upper_limit_reject : int, optional
            Sets the maximum number of rejected steps before the MCMC
            stops walking. Avoiding a slowly converging walk with few
            accepted points.

        Returns
        -------
        step_history : array
            An array containing all accepted steps of the MCMC.
        accept : int
            Total number of accepted steps.
        reject : int
            Total number of rejected steps.

        Notes
        -----
        Assumes uniform priors and a symetric jump proposal (gaussian).
        """
        if not isinstance(step_history, np.ndarray):
            step_history = np.array(step_history)
        if step_history.ndim == 1:
            step_history = np.array([step_history])
        # We will be appending to the list
        step_history = list(step_history)

        accept = 0
        reject = 0

        while (accept + reject < N_trials
               and reject < abs(int(upper_limit_reject))):
            current_step = step_history[-1]

            # f(θ)
            val = target_dist(classifier_name, current_step, **TD_kwargs)
            # θ+Δθ
            trial_step = current_step + np.random.normal(
                0, alpha, size=len(current_step)
            )
            # f(θ+Δθ)
            trial_val = target_dist(classifier_name, trial_step, **TD_kwargs)

            # check if the trial step is in the range of data
            for i, step_in_axis in enumerate(trial_step):
                if (
                    step_in_axis <= self._max_vals_[i]
                    and step_in_axis >= self._min_vals_[i]
                ):
                    pass
                else:
                    trial_val = 0  # essential reject points outside of range

            if val == 0:  # avoid div 0 errors
                ratio = 0
            else:
                ratio = trial_val / val

            accept_prob = min(1, (ratio) ** (1 / T))

            chance = np.random.uniform(low=0, high=1)
            if chance <= accept_prob:
                accept += 1
                step_history.append(trial_step)
            else:
                reject += 1

        return np.array(step_history), accept, reject

    def normalize_step_history(self, step_history):
        """Take steps and normalize [0,1] according to min/max in each axis.

        The max and min are taken from the original data set from TableData.

        """
        normed_steps = np.copy(step_history)
        for j, steps_in_axis in enumerate(step_history.T):
            normed_steps.T[j] = (steps_in_axis - self._min_vals_[j]) / (
                self._max_vals_[j] - self._min_vals_[j]
            )
        return normed_steps

    def undo_normalize_step_history(self, normed_steps):
        """Rescale step history.

        Take normed steps from [0,1] and return their value
        in the original range of the axes based off the range in TableData.

        """
        mapped_steps = np.copy(normed_steps)
        for j, steps_in_axis in enumerate(normed_steps.T):
            mapped_steps.T[j] = (
                steps_in_axis * (self._max_vals_[j] - self._min_vals_[j])
                + self._min_vals_[j]
            )
        return mapped_steps

    def do_simple_density_logic(self, step_history, N_points, Kappa,
                                var_mult=None, add_mvns_together=False,
                                include_training_data=True, verbose=False):
        """Perform multivariate normal density logic on a given step history.

        This is a simplified version of the method 'do_density_logic'.
        It assumes that every accepted point will have the same exact MVN.

        Each proposal distribution starts with the training set from TableData
        which keeps training data from being proposed again.

        Parameters
        ----------
        step_history : ndarray
            List of points from a PTMCMC or MCMC. (posterior)
        N_points : int
            Number of points desired to be drawn from the posterior but may not
            actually be the number of points accepted. Contributes to the
            length scale of the MVN distribution of accepted points
            (along with kappa).
        Kappa : float
            Scaling factor that sets the initial size of the MVN for accepted
            points. This should be proportional to the filling factor of the
            area of iterest described by the target distribution used to
            create the posterior.
        var_mult : float, ndarray, optional
            Variance multiplier for the MVN of accepted points.
        add_mvns_together : bool, optional
            Add MVNs together when creating the accepted point distribution.
        include_training_data : bool, optional
            Include the trainind data in the target distribution before
            sampling.
        verbose : bool, optional
            Print useful diagnostic information.

        Returns
        -------
        accepted_points : ndarray
            Accepted points from the posterior to be labled by the user.
            (query points)
        rejected_points : ndarray
            Rejected points from the posterior.

        Notes
        -----
        The accepted laguage here is indicative of query points for the oracle
        to label in an active learning scheme. It is not accepted vs rejected
        normally used for MCMC.

        """
        if include_training_data:
            if self._TableData_ is None:
                raise ValueError("No TableData found in sampler. "
                                 "Set `include_training_data` to false.")
            original_training_data = self._TableData_.get_data("input")
        else:
            original_training_data = []
        how_many_training_data = len(original_training_data)
        n_dim = len(original_training_data[0])
        # approximate scaling of the variance given a set of N proposal points
        # the filling factor Kappa is not generally known a priori
        if var_mult is None:
            var_mult = 1
        sigma = Kappa * 0.5 * (N_points) ** (-1.0 / n_dim) * var_mult
        Covariance = sigma * np.identity(n_dim)
        if verbose:
            print("Covariance: \n{0}".format(Covariance))

        single_MVN = scipy.stats.multivariate_normal(np.zeros(n_dim),
                                                     Covariance)
        max_val = 1 / np.sqrt((
            2 * np.pi) ** (n_dim) * np.linalg.det(Covariance))

        # treat the training data as already accepted points
        accepted_points = list(original_training_data)
        rejected_points = []
        for step in step_history:
            if len(accepted_points) < 1:
                accepted_points.append(step)
                continue

            # distance: how far is this point from all the accepted_points
            dist_to_acc_pts = step - np.array(accepted_points)
            pdf_val_at_point = single_MVN.pdf(dist_to_acc_pts)
            if isinstance(pdf_val_at_point, float):
                pdf_val_at_point = [pdf_val_at_point]

            if add_mvns_together:
                pdf_val_at_point = [np.sum(pdf_val_at_point)]
            random_chances = np.random.uniform(
                low=0, high=max_val, size=len(pdf_val_at_point)
            )
            chance_above_distr = random_chances > pdf_val_at_point

            if chance_above_distr.all():
                accepted_points.append(step)
            else:
                rejected_points.append(step)

        # remove training data from accepted points
        only_new_accpeted_points = accepted_points[how_many_training_data:]
        return np.array(only_new_accpeted_points), np.array(rejected_points)

    def do_density_logic(self, step_history, N_points, Kappa,  shuffle=False,
                         norm_steps=False, var_mult=None,
                         add_mvns_together=False, pre_acc_points=None,
                         verbose=False):
        """Do the density based of the normal gaussian kernel on each point.

        This method automatically takes out the first 5% of steps of the MCMC
        so that the initial starting points are not chosen automatically
        (if you start in a non-ideal region). Wait for the burn in.

        Parameters
        ----------
        step_history : ndarray
        N_points : int
        Kappa : float
        shuffle : bool, optional
        norm_steps : bool, optional
        var_mult : float, optional
        add_mvns_together : bool, optional
        verbose : bool, optional

        Returns
        -------
        accepted_points : ndarray
        rejected_points : ndarray
        accepted_sigmas : ndarray

        """
        if shuffle:  # shuffle the order of the steps
            if verbose:
                print("Shuffling steps....")
            np.random.shuffle(step_history)  # returns none

        if norm_steps:  # normalize steps
            if verbose:
                print("Normalizing steps....")
            step_history = self.normalize_step_history(step_history)

        # Set the default average length scale
        num_dim = len(self._Classifier_.input_data[0])
        sigma = Kappa * 0.5 * (N_points) ** (-1.0 / num_dim)

        # We assume the covaraince is the identity - later we may pass the
        # entire array but for now we just assume you pass a
        # variance multiplier (var_mult)
        if var_mult is None:
            var_mult = np.array([1] * num_dim)
        else:
            var_mult = np.array(var_mult)
            if var_mult.ndim != num_dim:
                raise ValueError(
                    "var_mult must be the same dimensionality as input data."
                )

        if verbose:
            print("Num dims: {0}".format(num_dim))
            print("length scale sigma: {0}".format(sigma))
            print("var_mult: {0}".format(var_mult))
            print("Kappa: {0}".format(Kappa))

        # -> Forcing a few key points to always be accepted, for example
        if pre_acc_points is None:
            accepted_points = []
        elif isinstance(pre_acc_points, np.ndarray):
            accepted_points = list(pre_acc_points)
        else:
            accepted_points = []

        accepted_points = []
        accepted_sigmas = []
        max_val_holder = []
        accepted_mvn_holder = []
        rejected_points = []

        skip_steps = int(len(step_history) * 0.05)
        good_steps = step_history[
            skip_steps:
        ]  # skip first 5% of steps to get into a good region

        for i in range(len(good_steps)):
            proposed_step = good_steps[i]

            accept = False
            if len(accepted_points) == 0:
                accept = True
            else:
                # If you enter you must have accepted one point
                Sigma = accepted_sigmas[-1:]
                k = len(Sigma)
                max_val = 1 / np.sqrt((2 * np.pi)**k * np.linalg.det(Sigma))
                max_val_holder.append(max_val)
                rnd_chance = np.random.uniform(
                    low=0, high=np.max(max_val_holder), size=1
                )
                # we will choose the chance from [0, highest point in distr]

                distr_holder = []
                for point in accepted_mvn_holder:
                    eval_mvn_at_new_step = point.pdf(proposed_step)
                    distr_holder.append(eval_mvn_at_new_step)

                if add_mvns_together:
                    # instead of checking each individual point keeping all
                    # mvns seperate, we want to add them
                    # together and get upper bound
                    # IF we do this we need to change the MVN to not be
                    # normalized !!!!
                    # THE UNORMALIZED MVN IS NOT IMPLEMENTED
                    total_chance_above_distr = (
                        rnd_chance > np.sum(distr_holder))
                else:
                    total_chance_above_distr = np.sum(
                        rnd_chance > distr_holder)

                if len(accepted_points) == total_chance_above_distr:
                    accept = True
                else:
                    pass  # REJECT

            if accept:
                # https://stackoverflow.com/questions/619335
                # corner = np.random.normal(0,0.1, 1)
                # A = np.array( [ [sigma, corner], [corner, sigma] ] )
                # Sigma = np.dot(A,A.transpose())
                # Sigma = [ [sigma*var_mult[0], 0.], [0., sigma*var_mult[1]] ]
                Sigma = (sigma * np.identity(len(var_mult))
                         * np.array([var_mult]))
                mvn = scipy.stats.multivariate_normal(proposed_step, Sigma)

                accepted_mvn_holder.append(mvn)
                accepted_sigmas.append(Sigma)
                accepted_points.append(proposed_step)
            else:
                rejected_points.append(proposed_step)

        if verbose:
            print("Num accepted: {0}".format(len(accepted_points)))
            print("Num rejected: {0}".format(len(rejected_points)))

        if norm_steps:
            if verbose:
                print("Unormalizing steps....")
            accepted_points = self.undo_normalize_step_history(accepted_points)
            rejected_points = self.undo_normalize_step_history(rejected_points)

        return (
            np.array(accepted_points),
            np.array(rejected_points),
            np.array(accepted_sigmas),
        )

    def get_proposed_points(self, step_history, N_points, Kappa, shuffle=False,
                            norm_steps=False, add_mvns_together=False,
                            include_training_data=True, var_mult=None,
                            seed=None, n_repeats=1, max_iters=1e3,
                            verbose=False, **kwargs):
        """Get proposed points in parameter space given a MCMC step history.

        The desnity logic is not deterministic, so multiple iterations
        may be needed to converge on a desired number of proposed points.
        This method performs multiple calls to do_density_logic while
        changing Kappa in order to return the desired number of points. After
        n_iters instances of the correct number of N_points, the distibution
        with the largest average distance is chosen.

        Warning: This algorithm has not been tested for large N data sets and
                 may struggle to converge.

        Parameters
        ----------
        step_history : ndarray
            Posterior from which to sample new query points.
        N_points : int
            N query points to converge to.
        Kappa : float
            Multiplies the length scale of MVNs and changes such that the
            desired number of query points is found.
        shuffle : bool, optional
            Shuffle points in posterior in place before sampling.
        norm_steps : bool, optional
            Normalize steps before sampling.
        add_mvns_together : bool, optional
            Add MVNs of accepted point distribution together.
        include_training_data : bool, optional
            Include training data in the accpeted point distribution before
            sampling.
        var_mult : ndarray, optional
            Variance multiplier.
        seed : float, optional
            Random seed to use for random sampling.
        n_repeats: int, optional
            Number of times to converge to the correct number of points. Each
            iteration may be a different realization of the posterior.
        verbose : bool, optional
            Print useful information.
        **kwargs
            show_plots : bool, optional
            Show 2D plot of proposed points with step history & training data.

        Returns
        -------
        acc_pts : ndarray
            Array of proposed points to be used as initial
            conditions in new simulations.
        Kappa : float
            Scaling factor which reproduced the desired number
            of accepted points.

        Notes
        -----
        Will automatically exit if it goes through max_iters iterations
        without converging on the desired number of points.

        """
        if seed is not None:
            np.random.seed(seed=seed)
            print("Setting seed: {0}".format(seed))

        if verbose:
            print(
                "Converging to {0} points, {1} times.".format(N_points,
                                                              int(n_repeats))
            )

        enough_good_pts = False
        how_many_good_pts = int(n_repeats)
        good_n_points = []
        avg_distances = []
        good_kappas = []
        iters = 1
        start_time = time.time()
        while not (enough_good_pts) and iters < max_iters:

            acc_pts, rej_pts = self.do_simple_density_logic(
                step_history,
                N_points,
                Kappa,
                var_mult=var_mult,
                add_mvns_together=add_mvns_together,
                include_training_data=include_training_data,
                verbose=False,
            )
            if len(acc_pts) == 0:
                Kappa = Kappa * 0.5
                iters += 1
                continue

            if len(acc_pts) == N_points:
                average_dist_between_acc_points = np.mean(pdist(acc_pts))
                good_n_points.append(acc_pts)
                avg_distances.append(average_dist_between_acc_points)
                good_kappas.append(Kappa)
                if len(good_n_points) >= how_many_good_pts:
                    enough_good_pts = True

            if verbose:
                if len(acc_pts) == N_points:
                    print_str = "  *{0}*  {1:2.2f}s".format(
                        len(good_n_points), abs(start_time - time.time())
                    )
                    ending = "\n"
                else:
                    print_str = ""
                    ending = "\r"
                print(
                    "\t acc_pts: {0}, Kappa = {1:.3f}".format(len(acc_pts),
                                                              Kappa)
                    + print_str,
                    end=ending,
                )

            diff = abs(len(acc_pts) - N_points)

            change_factor = 0.01 / (max(1, np.log10(iters)))
            if len(acc_pts) > N_points:
                Kappa = Kappa * (1 + change_factor * diff)  # increase kappa
            elif len(acc_pts) < N_points:
                if (1 - change_factor * diff) < 0:
                    Kappa = Kappa * 0.1
                else:
                    Kappa = Kappa * (1 - change_factor*diff)  # decrease kappa

            iters += 1

        if iters == max_iters:
            print("Reached max iters before converging!")
        if verbose:
            print("\nFinal Kappa = {0}\nConverged in {1} iters.".format(Kappa,
                                                                        iters))

        # we want 1/r dependance to penalize closely spaced points
        where_best_distribution = np.argmax(avg_distances)
        best_acc_pts = good_n_points[where_best_distribution]
        best_Kappa = good_kappas[where_best_distribution]
        if verbose:
            print("Average Distances: \n{0}".format(np.array(avg_distances)))
            print("Kappas: \n{0}".format(np.array(good_kappas)))
            print("loc: {0}".format(where_best_distribution))

        if kwargs.get("show_plots", False):
            self.make_prop_points_plots(step_history, best_acc_pts)
        return best_acc_pts, best_Kappa

    def make_prop_points_plots(self, step_hist, prop_points, axes=(0, 1),
                               show_fig=True, save_fig=False):
        """Plot the proposed / accepted points over the step history."""
        axis1, axis2 = axes
        fig, sub = plt.subplots(nrows=1, ncols=1, figsize=(4, 4), dpi=90)
        training_data = self._Classifier_.input_data
        sub.scatter(
            training_data.T[axis1],
            training_data.T[axis2],
            color="k",
            marker="x",
            label="training data",
        )
        sub.scatter(
            step_hist.T[axis1],
            step_hist.T[axis2],
            color="C0",
            alpha=0.5,
            label="step history",
        )
        sub.scatter(
            prop_points.T[axis1],
            prop_points.T[axis2],
            color="pink",
            label="new proposed points",
        )
        sub.set_xlim(self._min_vals_[axis1], self._max_vals_[axis1])
        sub.set_ylim(self._min_vals_[axis2], self._max_vals_[axis2])
        plt.legend(loc="best", bbox_to_anchor=[1, 0, 0.22, 1])
        return fig, sub

    def make_trace_plot(self, chain_holder, T_list, Temp, save_fig=False,
                        show_fig=True):
        """Make a `step number` vs. `position of a sampler in an axis` plot.

        This function makes titles assuming you are using the data
        from the classifier.

        """
        if not show_fig and not save_fig:
            return

        which_temp = Temp  # int? index

        n_axis = len(chain_holder[which_temp].T)
        steps = chain_holder[which_temp].T
        axis_names = self._Classifier_._TableData_.get_data(
            what_data="input", return_df=True
        ).keys()

        fig, axs = plt.subplots(
            nrows=n_axis,
            ncols=2,
            figsize=(10, 2.5 * n_axis),
            dpi=100,
            gridspec_kw={"width_ratios": [1.8, 1]},
        )

        for num, ax in enumerate(axs):
            ax[0].plot(steps[num], "-", linewidth=0.5, color="C4")
            ax[0].set_title("Input: {0}".format(axis_names[num]))
            ax[0].set_ylabel("Axis {0}".format(num) + " , T={0:.2f}".
                             format(T_list[which_temp]), fontsize=12)
            ax[0].set_xlabel("N steps", fontsize=12)

            ax[1].hist(steps[num],
                       bins=50, histtype="step", density=True, color="C1")
            ax[1].set_xlabel("Axis {0}".format(num), fontsize=12)
            ax[1].set_ylabel("Posterior", fontsize=12)

        fig.subplots_adjust(hspace=0.45)
        if save_fig:
            plt.savefig("trace_plot_T{0:.0f}.pdf".format(T_list[which_temp]))
        if show_fig:
            plt.show()
        return None
