"""Provides functions to evaluate the performance of grid downsampling."""


__authors__ = [
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
]


from posydon.grids.psygrid import PSyGrid
import numpy as np
import matplotlib.pyplot as plt
import tqdm


def report_DS(path, path_DS, max_err_used=None,
              nmax=None, emax=None, emin=None):
    """Report on compression ratio and interpolation error statistics."""
    TABLES = ["binary_history", "history1", "history2",
              "final_profile1", "final_profile2"]
    STEP = 0.01
    grid = PSyGrid(path)
    grid_DS = PSyGrid(path_DS)

    x_H, x_P1, x_P2, y_H, y_P1, y_P2 = [[] for _ in range(6)]

    avg_interp_errors = {}
    max_interp_errors = {}
    med_interp_errors = {}

    def lengths(arrays1, arrays2):
        assert(len(arrays1) == len(arrays2))
        result1 = []
        result2 = []
        for arr1, arr2 in zip(arrays1, arrays2):
            if arr1 is None:
                assert(arr2 is None)
                result1.append(np.nan)
                result2.append(np.nan)
            else:
                assert(arr2 is not None)
                result1.append(len(arr1))
                result2.append(len(arr2))
        return result1, result2

    for i, (run, run_DS) in tqdm.tqdm(enumerate(zip(grid, grid_DS))):
        if nmax is not None and i >= nmax:
            break
        BH, H1, H2, P1, P2 = [run[table] for table in TABLES]
        BH_DS, H1_DS, H2_DS, P1_DS, P2_DS = [run_DS[table] for table in TABLES]

        (x_h, x_p1, x_p2), (y_h, y_p1, y_p2) = lengths(
            [BH, P1, P2], [BH_DS, P1_DS, P2_DS])
        x_H.append(x_h)
        x_P1.append(x_p1)
        x_P2.append(x_p2)
        y_H.append(y_h)
        y_P1.append(y_p1)
        y_P2.append(y_p2)

        for tables, fromwhere, independent in zip(
                [[BH, BH_DS], [H1, H1_DS], [H2, H2_DS],
                 [P1, P1_DS], [P2, P2_DS]],
                ["BH", "H1", "H2", "P1", "P2"],
                ["age", "star_age", "star_age", "mass", "mass"]
        ):
            orig, down = tables
            if orig is None or down is None:
                continue

            colnames = orig.dtype.names
            for colname in colnames:
                if colname == independent:
                    continue
                fullname = fromwhere + "." + colname

                t_orig = orig[independent]
                if len(t_orig) <= 2:
                    continue

                t_down = down[independent]
                X_orig = orig[colname]

                if colname not in down.dtype.names:
                    if colname == "model_number":
                        # ignore it... probably comparing original vs EEP grid
                        continue
                X_down = down[colname]

                if fromwhere in ["P1", "P2"]:
                    X_int = np.interp(t_orig[::-1], t_down[::-1], X_down[::-1])
                    X_int = X_int[::-1]
                else:
                    X_int = np.interp(t_orig, t_down, X_down)
                errors = X_int - X_orig

                # ignore cases where interpolation breaks anyway
                where_ok = np.ones_like(t_orig, dtype=bool)
                where_ok[:-1] = np.diff(t_orig) > 0
                where_ok[1:] &= np.diff(t_orig) > 0
                errors = errors[where_ok]

                if len(errors) == 0:
                    continue
                if not np.all(errors == 0.0):
                    errors = np.abs(errors / (np.max(X_orig) - np.min(X_orig)))

                avg_error = np.mean(errors)
                max_error = np.max(errors)
                med_error = np.median(errors)

                if fullname not in avg_interp_errors:
                    avg_interp_errors[fullname] = [avg_error]
                    med_interp_errors[fullname] = [med_error]
                    max_interp_errors[fullname] = [max_error]
                else:
                    avg_interp_errors[fullname].append(avg_error)
                    med_interp_errors[fullname].append(med_error)
                    max_interp_errors[fullname].append(max_error)

    x_H, x_P1, x_P2, y_H, y_P1, y_P2 = [
        np.array(arr) for arr in [x_H, x_P1, x_P2, y_H, y_P1, y_P2]]

    grid.close()
    grid_DS.close()

    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
    for x, y, label in zip([x_H, x_P1, x_P2],
                           [y_H, y_P1, y_P2],
                           ["History", "Profile1", "Profile2"]):
        x = x[np.isfinite(x)]
        y = y[np.isfinite(y)]
        if len(x) == 0 or len(y) == 0:
            continue
        ratios = y / x
        x_min = np.floor(min(ratios) / STEP) * STEP - STEP
        x_max = max(1.0, np.ceil(max(ratios) / STEP) * STEP) + STEP
        try:
            bins = np.arange(x_min, x_max + STEP / 2.0, STEP)
        except Exception:
            bins = "fd"
        ax1.hist(ratios, bins=bins, histtype="step",
                 label=label+" (n={})".format(len(ratios)))
        ax2.plot(x, y, ".", label=label)
    ax1.set_xlabel("Compression ratio")
    ax1.set_ylabel("Number of runs")
    ax2.set_xlabel("Original size")
    ax2.set_ylabel("Downsampled size")
    ax1.legend(loc="best")
    ax2.legend(loc="best")
    plt.tight_layout()
    plt.show()

    sorted_avg_interp_errors = {k: v for k, v in sorted(
        avg_interp_errors.items(), key=lambda item: np.mean(item[1]))}

    plt.figure(figsize=(8, 20))
    tick_labels = []
    tick_positions = []
    for i, (key, values) in enumerate(sorted_avg_interp_errors.items()):
        y = np.ones_like(values) * i
        tick_positions.append(i)
        tick_labels.append(key)
        plt.plot(max_interp_errors[key], y, "r.")
        plt.plot(values, y, "ko", mfc="none")
        plt.plot(med_interp_errors[key], y+0.2, "g.")
    plt.plot([], [], "r.", label="Maximum")
    plt.plot([], [], "ko", mfc="none", label="Average")
    plt.plot([], [], "g.", label="Median")
    plt.xscale("log")
    plt.yticks(tick_positions, tick_labels)
    plt.ylim(min(tick_positions)-1, max(tick_positions)-1)
    if max_err_used is not None:
        plt.axvline(max_err_used, label="max_err={:.4g}".format(max_err_used))
    plt.xlim(xmin=emin, xmax=emax)
    plt.grid()
    plt.xlabel("Average interpolation error")
    plt.legend(loc="upper left")
    plt.tight_layout()
    plt.show()


def compare_DS(path, path_DS, runs=None, useonly=None):
    """Compare the data in original and downsampled grids."""
    TABLES = ["binary_history", "history1", "history2",
              "final_profile1", "final_profile2"]
    grid = PSyGrid(path)
    grid_DS = PSyGrid(path_DS)

    def lengths(arrays1, arrays2):
        assert(len(arrays1) == len(arrays2))
        result1 = []
        result2 = []
        for arr1, arr2 in zip(arrays1, arrays2):
            if arr1 is None:
                assert(arr2 is None)
                result1.append(np.nan)
                result2.append(np.nan)
            else:
                assert(arr2 is not None)
                result1.append(len(arr1))
                result2.append(len(arr2))
        return result1, result2

    for i, (run, run_DS) in tqdm.tqdm(enumerate(zip(grid, grid_DS))):
        if runs is not None and i not in runs:
            continue
        BH, H1, H2, P1, P2 = [run[table] for table in TABLES]
        BH_DS, H1_DS, H2_DS, P1_DS, P2_DS = [run_DS[table] for table in TABLES]

        for tables, fromwhere, independent in zip(
                [[BH, BH_DS], [H1, H1_DS], [H2, H2_DS],
                 [P1, P1_DS], [P2, P2_DS]],
                ["BH", "H1", "H2", "P1", "P2"],
                ["age", "star_age", "star_age", "mass", "mass"]
        ):
            orig, down = tables
            if orig is None or down is None:
                continue

            colnames = orig.dtype.names
            for colname in colnames:
                if colname == independent:
                    continue
                fullname = fromwhere + "." + colname

                if useonly is not None:
                    if not any([fullname.startswith(col) for col in useonly]):
                        continue

                t_orig = orig[independent]
                t_down = down[independent]
                X_orig = orig[colname]
                X_down = down[colname]

                if fromwhere in ["P1", "P2"]:
                    X_int = np.interp(
                        t_orig[::-1], t_down[::-1], X_down[::-1])[::-1]
                else:
                    X_int = np.interp(t_orig, t_down, X_down)
                errors = X_int - X_orig
                if not np.all(errors == 0.0):
                    errors = np.abs(errors / (np.max(X_orig) - np.min(X_orig)))

                plt.figure(figsize=(6.4, 2.4))
                plt.suptitle("Run {}".format(i))
                plt.subplot(121)
                plt.plot(t_orig, X_orig, "k.--")
                plt.plot(t_down, X_down, "r.:")
                plt.ylabel(fullname)
                plt.subplot(122)
                plt.plot(t_orig, errors, "k-")
                plt.ylabel("Error")
                plt.tight_layout()
                plt.show()

    grid.close()
    grid_DS.close()
