"""Module for analyzing binary population simulation results."""


__authors__ = [
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
]


import pandas as pd


CHUNKSIZE = 100000
CHUNKSIZE_SMALL = 1000


def report_columns(path):
    """Print the columns in the `oneline` and `history` keys in an HDF5 file.

    Parameters
    ----------
    path : str
        The fullpath of the HDF5 file.
    """
    print("Population file:", path)
    for tablename in ["oneline", "history"]:
        for chunk in pd.read_hdf(path, tablename, chunksize=1):
            print("\nColumns in `{}`\n{}\n{}".format(
                tablename, "-"*20, " | ".join(chunk.columns)))
            break


def _decide_chunksize(n_binaries):
    """Decide the HDF5 chunk size when searching for n_binaries."""
    if n_binaries is None:
        return CHUNKSIZE
    return min(CHUNKSIZE, max(CHUNKSIZE_SMALL, n_binaries // 10))


def _filter_columns(columns, patterns):
    """Select a subset of columns based on a substring or array of substrings.

    Parameters
    ----------
    columns : array of str
        The original columns to be filtered.
    patterns : str or array of str
        The pattern or list of patterns, where pattern is a substring.

    Returns
    -------
    array of str
        Subset of the columns on the basis of the inclusion of substring(s).

    """
    if patterns is None:                    # return input if undefined pattern
        return list(columns)
    if isinstance(patterns, (str, bytes)):  # in the case of one pattern
        patterns = [patterns]

    selected_cols = []
    for col in list(columns):
        include_it = False
        for pattern in patterns:
            if pattern in col:
                include_it = True
                break
        if include_it:
            selected_cols.append(col)

    return selected_cols


def get_subsets(path, find_n=None, select_oneline=None, select_history=None,
                select_columns=None, verbose=True):
    """Get the `oneline` and `history` subsets for a selection of binaries.

    You select binaries from an HDF5 population by providing a function that
    operates either on `oneline` or `history` and returns a boolean mask in-
    dicating which rows correspond to the binaries of interest. E.g.,

        def find_no_end(df):
            return df["event_f"] != "END"

    can be used to find binaries that ended abruptly, based on the information
    in the `oneline` table. Therefore, we set `select_oneline=find_no_end`. If
    no selector has been set, all binaries are selected (combine this with the
    `select_columns` argument to get a column subset of the data.)

    Parameters
    ----------
    path : str
        The fullpath of the HDF5 file.
    find_n : int or None
        The number of binaries to select (at most) or None to select all.
    select_oneline : function or None
        The function to operate on chunks of the `oneline` table.
    select_history : function or None
        The function to operate on chunks of the `oneline` table.
    select_columns : None or str or array
        If None, all columns will be returned. If a string, only those columns
        with name containing this string will be returned. An array of strings
        can be used for multiple such patterns.

    Returns
    -------
    oneline, history
        The subsets of the dataframes with the data for the selected binaries.

    """
    def say(*args, **kwargs):
        """Print only if verbose is True."""
        if verbose:
            print(*args, **kwargs)

    if select_oneline is not None and select_history is not None:
        raise Exception("Cannot use both `oneline` and `history` selectors.")
    elif select_oneline is not None:
        selector = select_oneline
        tablename = "oneline"
    elif select_history is not None:
        selector = select_history
        tablename = "history"
    else:
        selector = None

    chunksize = _decide_chunksize(find_n)
    say("Chunk size:", chunksize)

    if select_columns is not None:
        say("Deciding which columns to include...")
        for chunk in pd.read_hdf(path, "oneline", chunksize=1):
            oneline_columns = _filter_columns(chunk.columns, select_columns)
            say("    Selected oneline columns:", oneline_columns)
            break
        for chunk in pd.read_hdf(path, "history", chunksize=1):
            history_columns = _filter_columns(chunk.columns, select_columns)
            say("    Selected history columns:", history_columns)
            break
    else:
        oneline_columns, history_columns = None, None

    say("Applying the selector...")
    if selector is None:
        say("    No selector defined. Selecting all.")
    else:
        all_indices, selected_indices = set(), set()
        for chunk in pd.read_hdf(path, tablename, chunksize=chunksize):
            all_indices.update(chunk.index)
            selection = selector(chunk)
            selected_indices.update(chunk[selection].index)
            say("   Selected {:7d} out of {:7d} binaries ({:.3g}%)...\t".
                format(len(selected_indices), len(all_indices),
                       100*len(selected_indices)/len(all_indices)), end="\r")
            if find_n is not None and len(selected_indices) >= find_n:
                selected_indices = set(list(selected_indices)[:find_n])
                say("\n    Found {} binaries.".format(find_n))
                break
        say()

    say("Obtaining `oneline` for the selection...")
    selected_oneline, included_indices = pd.DataFrame(), set()
    found_them = False
    for chunk in pd.read_hdf(path, "oneline",
                             chunksize=chunksize, columns=oneline_columns):
        if selector is not None:
            chunk = chunk[chunk.index.isin(selected_indices)]
            included_indices.update(chunk.index)
            found_them = len(included_indices) == len(selected_indices)
        selected_oneline = pd.concat([selected_oneline, chunk])
        if found_them:
            say("    Found the binaries.")
            break

    say("Obtaining `history` for the selection...")
    selected_history, included_indices = pd.DataFrame(), set()
    stop_at_next_chunk, found_them = False, False
    for chunk in pd.read_hdf(path, "history",
                             chunksize=chunksize, columns=history_columns):
        if selector is not None:
            chunk = chunk[chunk.index.isin(selected_indices)]
            included_indices.update(chunk.index)
            found_them = len(included_indices) == len(selected_indices)
        selected_history = pd.concat([selected_history, chunk])
        if found_them:
            if not stop_at_next_chunk:
                stop_at_next_chunk = True
            else:
                say("    Found the binaries.")
                break

    return selected_oneline, selected_history
