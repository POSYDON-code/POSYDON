# --------------------------------------------
# run a psy-cris sequence on synthetic data
# --------------------------------------------
import argparse
import numpy as np
import pickle
import copy
import time
import sys
import os
from posydon.active_learning.psy_cris.utils import parse_inifile
from posydon.active_learning.psy_cris.utils import do_dynamic_sampling
from posydon.active_learning.psy_cris.utils import calc_performance

parser = argparse.ArgumentParser()
parser.add_argument( "-ini", dest='inifile_path', type=str, help="Path to psy-cris inifile.")
parser.add_argument( "-dir", dest='save_data_dir', type=str, help="Path to directory to save data.")

parser.add_argument( "-id", dest='run_id', type=str, help="Integer or str used to identify runs.",
                    default = np.random.randint(low=10000, high=99999) )
parser.add_argument( "-n", dest='n_sequences', type=int, help="N sequences to run.",
                    default=1)
parser.add_argument( "-num_start", dest='n_starting_points', type=int, help="N points to start with in even grid.",
                    default=3**3)
parser.add_argument( "-num_fin", dest='n_final_points', type=int, help="Total num points to converge to.",
                    default=6**3)
parser.add_argument( "-ppi", dest='points_per_iter', type=int, help="Num points between each iteration.",
                    default=70)
parser.add_argument( "-cls", dest='performance_cls_name', type=str, help="Classification algorithm used in performance calc.",
                    default="linear")
parser.add_argument( "-regr", dest='performance_regr_name', type=str, help="Regression algorithm used in performance calc.",
                    default="rbf")
parser.add_argument( "-res", dest='performance_res', type=int, help="Resolution to compute performance.",
                    default=40)
parser.add_argument( "-length_scale_mult", dest='length_scale_multiplier', type=float, help="Small class proposal length scale mult of normal.",
                    default=0.3333)
parser.add_argument( "-percent_increase", dest='percent_increase', type=float, help="Icrease points per iter by a percent of the training set.",
                    default=-1)


args = parser.parse_args()

inifile_path = args.inifile_path
save_data_dir = args.save_data_dir
if not os.path.isfile(inifile_path):
    raise ValueError("Inifile not found. Please check that the path is correct.\nGiven: {}".format(inifile_path))
if not os.path.isdir(save_data_dir):
    raise ValueError("Directory not found. Please check that save_data_dir exists.\nGiven: {}".format(save_data_dir))

run_ID = args.run_id
n_psycris_sequences = args.n_sequences

n_starting_points = args.n_starting_points
n_final_points = args.n_final_points
new_points_per_iter = args.points_per_iter
length_scale_mult = args.length_scale_multiplier
if args.percent_increase < 0:
    percent_increase = None
else:
    percent_increase = args.percent_increase

psy_cris_kwargs_dict = parse_inifile(inifile_path, verbose=False)
dimensionality = len(psy_cris_kwargs_dict["TableData_kwargs"]["input_cols"])

print("\n\n"+">"*18 + " RUNNING PSYCRIS SEQUENCE ID: {} ".format(run_ID) + "<"*18)
print("FILE: {}".format(inifile_path))
print("SAVE DIR: {}\n".format(save_data_dir))
print("n_starting_points: {}".format(n_starting_points))
print("n_final_points: {}".format(n_final_points))
print("new_points_per_iter: {}".format(new_points_per_iter))
print("dimensionality: {}".format(dimensionality))
print("length_scale_multiplier: {}".format(length_scale_mult))
print("percent_increase: {}\n".format(percent_increase))


start_iters_time = time.time()
for i in range(n_psycris_sequences):
    print("START - do_dynamic_sampling - {0}i{1}".format(run_ID,i))
    kwargs_per_iter = copy.deepcopy(psy_cris_kwargs_dict)
    dfs_per_iter, preds_per_iter = do_dynamic_sampling(N_starting_points=n_starting_points,
                                                       N_final_points=n_final_points,
                                                       new_points_per_iter=new_points_per_iter,
                                                       verbose=True,
                                                       threshold= 1e-6,
                                                       jitter=True,
                                                       dim=dimensionality,
                                                       length_scale_mult=length_scale_mult,
                                                       percent_increase=percent_increase,
                                                       **kwargs_per_iter )
    try:
        print("\n\tSaving dfs...")
        f_name_backup = save_data_dir + "/{0}i{1}_dfs_per_iter".format(run_ID,i)
        with open(f_name_backup, "wb") as f:
            pickle.dump( dfs_per_iter, f)

        print("\tSTART - calc_performance - {0}i{1}".format(run_ID,i))
        print("\tcls: {}, regr: {}, res: {}".format(args.performance_cls_name, args.performance_regr_name, args.performance_res) )
        acc_per_iter, conf_matrix_per_iter, abs_regr_frac_diffs_per_iter = \
                                    calc_performance(dfs_per_iter,
                                                     cls_name=args.performance_cls_name,
                                                     regr_name=args.performance_regr_name,
                                                     resolution=args.performance_res, verbose=False)

        f1_path = save_data_dir + "/{0}i{1}_acc_per_iter".format(run_ID,i)
        f2_path = save_data_dir + "/{0}i{1}_conf_matrix_per_iter".format(run_ID,i)
        f3_path = save_data_dir + "/{0}i{1}_abs_regr_frac_diffs_per_iter".format(run_ID,i)
        print("\tSaving files:")
        try:
            np.save( f1_path, acc_per_iter, allow_pickle=True )
            print("\t\t{}".format(f1_path))
        except Exception as err:
            print("\t\tFAILED: {0}, err: {1}".format(f1_path, err))
        try:
            np.save( f2_path, conf_matrix_per_iter, allow_pickle=True )
            print("\t\t{}".format(f2_path))
        except Exception as err:
            print("\t\tFAILED: {0}, err: {1}".format(f2_path, err))
        try:
            np.save( f3_path, np.array(abs_regr_frac_diffs_per_iter, dtype=object), allow_pickle=True )
            print("\t\t{}".format(f3_path))
        except Exception as err:
            print("\t\tFAILED: {0}, err: {1}".format(f3_path, err))

    except Exception as exc:
        print("\tPerformance calculation Failed!\n\tErr:{}\n\n".format(exc))
print("END {0}\nTotal Time : {1:.3f} min".format(run_ID,(time.time()-start_iters_time)/60) )
