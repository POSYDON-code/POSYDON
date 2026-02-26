#!/usr/bin/env python3
"""
Compare evolution outcomes of test binaries saved by binaries_suite.py.

Reports three categories of differences:
  1. QUANTITATIVE: any numeric difference beyond floating-point representation
  2. QUALITATIVE: changes to categorical/string columns (states, events, step names, SN types, etc.)
  3. WARNINGS & ERRORS: changes to warnings raised or binaries that error out

By default, uses exact comparison (atol=0, rtol=0). The --loose flag enables
a small tolerance for cases where minor floating-point differences are expected.

Usage:
    python compare_runs.py baseline.h5 candidate.h5
    python compare_runs.py baseline.h5 candidate.h5 --verbose
    python compare_runs.py baseline.h5 candidate.h5 --loose

Authors: Elizabeth Teng
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

# Columns that represent qualitative (categorical) evolution properties.
# Any column matching these names will be compared as exact string matches
# and reported under "QUALITATIVE" differences.
QUALITATIVE_COLUMNS = {
    'state', 'event', 'step_names', 'S1_state', 'S2_state',
    'SN_type', 'S1_SN_type', 'S2_SN_type',
    'interp_class_HMS_HMS', 'interp_class_CO_HMS_RLO',
    'interp_class_CO_HeMS', 'interp_class_CO_HeMS_RLO',
    'mt_history_HMS_HMS', 'mt_history_CO_HMS_RLO',
    'mt_history_CO_HeMS', 'mt_history_CO_HeMS_RLO',
    'mass_transfer_case',
}


def classify_column(col, dtype):
    """Classify a column as 'qualitative' or 'quantitative'."""
    if col in QUALITATIVE_COLUMNS:
        return 'qualitative'
    if pd.api.types.is_numeric_dtype(dtype):
        return 'quantitative'
    # Catch-all: treat remaining object/string columns as qualitative
    return 'qualitative'


def compare_evolution_tables(base_df, cand_df, rtol, atol):
    """Compare two evolution DataFrames, reporting per-binary diffs.

    Returns:
        dict with keys 'quantitative', 'qualitative', 'structural'
        each mapping to a list of diff strings.
    """
    quant_diffs = []
    qual_diffs = []
    struct_diffs = []

    # Check that binary_id columns exist
    if 'binary_id' not in base_df.columns or 'binary_id' not in cand_df.columns:
        struct_diffs.append("'binary_id' column missing; cannot do per-binary comparison")
        return {'quantitative': quant_diffs, 'qualitative': qual_diffs, 'structural': struct_diffs}

    base_ids = set(base_df['binary_id'].unique())
    cand_ids = set(cand_df['binary_id'].unique())

    # Missing/extra binaries
    for bid in sorted(base_ids - cand_ids):
        struct_diffs.append(f"Binary {bid}: MISSING in candidate")
    for bid in sorted(cand_ids - base_ids):
        struct_diffs.append(f"Binary {bid}: EXTRA in candidate")

    common_ids = sorted(base_ids & cand_ids)

    for bid in common_ids:
        b = base_df[base_df['binary_id'] == bid].reset_index(drop=True)
        c = cand_df[cand_df['binary_id'] == bid].reset_index(drop=True)

        # ── Error status changes ──────────────────────────────────────
        base_failed = 'exception_type' in b.columns and b['exception_type'].notna().any()
        cand_failed = 'exception_type' in c.columns and c['exception_type'].notna().any()

        if base_failed != cand_failed:
            if cand_failed:
                exc = c['exception_type'].dropna().iloc[0] if 'exception_type' in c.columns else "unknown"
                msg = c['exception_message'].dropna().iloc[0] if 'exception_message' in c.columns else ""
                struct_diffs.append(f"Binary {bid}: NEWLY FAILING ({exc}: {msg})")
            else:
                struct_diffs.append(f"Binary {bid}: NEWLY PASSING (was failing in baseline)")
            continue

        if base_failed and cand_failed:
            b_exc = str(b['exception_type'].dropna().iloc[0]) if 'exception_type' in b.columns else ""
            c_exc = str(c['exception_type'].dropna().iloc[0]) if 'exception_type' in c.columns else ""
            b_msg = str(b['exception_message'].dropna().iloc[0]) if 'exception_message' in b.columns else ""
            c_msg = str(c['exception_message'].dropna().iloc[0]) if 'exception_message' in c.columns else ""
            if b_exc != c_exc:
                struct_diffs.append(f"Binary {bid}: error type changed ('{b_exc}' -> '{c_exc}')")
            if b_msg != c_msg:
                struct_diffs.append(f"Binary {bid}: error message changed ('{b_msg}' -> '{c_msg}')")
            continue

        # ── Step count ────────────────────────────────────────────────
        if len(b) != len(c):
            struct_diffs.append(
                f"Binary {bid}: evolution step count differs "
                f"(baseline={len(b)}, candidate={len(c)})"
            )

        # ── Column presence ───────────────────────────────────────────
        base_only_cols = set(b.columns) - set(c.columns) - {'binary_id'}
        cand_only_cols = set(c.columns) - set(b.columns) - {'binary_id'}
        if base_only_cols:
            struct_diffs.append(f"Binary {bid}: columns only in baseline: {sorted(base_only_cols)}")
        if cand_only_cols:
            struct_diffs.append(f"Binary {bid}: columns only in candidate: {sorted(cand_only_cols)}")

        # ── Per-column comparison ─────────────────────────────────────
        common_cols = sorted(set(b.columns) & set(c.columns) - {'binary_id'})
        min_rows = min(len(b), len(c))

        for col in common_cols:
            b_col = b[col].iloc[:min_rows]
            c_col = c[col].iloc[:min_rows]
            col_type = classify_column(col, b_col.dtype)

            if col_type == 'quantitative':
                b_arr = b_col.to_numpy(dtype=float)
                c_arr = c_col.to_numpy(dtype=float)

                # NaN handling: both NaN = match, one NaN = mismatch
                both_nan = np.isnan(b_arr) & np.isnan(c_arr)
                one_nan = np.isnan(b_arr) ^ np.isnan(c_arr)

                if one_nan.any():
                    nan_steps = np.where(one_nan)[0].tolist()
                    direction = []
                    for s in nan_steps[:5]:
                        bv = "NaN" if np.isnan(b_arr[s]) else f"{b_arr[s]:.6g}"
                        cv = "NaN" if np.isnan(c_arr[s]) else f"{c_arr[s]:.6g}"
                        direction.append(f"step {s}: {bv} -> {cv}")
                    quant_diffs.append(
                        f"Binary {bid}, '{col}': NaN mismatch at {len(nan_steps)} step(s): "
                        + "; ".join(direction)
                    )

                # Compare non-NaN values
                valid = ~(np.isnan(b_arr) | np.isnan(c_arr))
                if valid.any():
                    b_valid = b_arr[valid]
                    c_valid = c_arr[valid]
                    not_equal = b_valid != c_valid

                    if rtol == 0 and atol == 0:
                        # Exact comparison
                        if not_equal.any():
                            diff_indices = np.where(valid)[0][not_equal]
                            abs_diff = np.abs(b_valid[not_equal] - c_valid[not_equal])
                            worst = np.argmax(abs_diff)
                            worst_step = diff_indices[worst]
                            quant_diffs.append(
                                f"Binary {bid}, '{col}': {not_equal.sum()} value(s) differ. "
                                f"Largest abs diff = {abs_diff[worst]:.6e} "
                                f"at step {worst_step} "
                                f"(baseline={b_valid[not_equal][worst]:.15g}, "
                                f"candidate={c_valid[not_equal][worst]:.15g})"
                            )
                    else:
                        # Tolerance-based comparison
                        if not np.allclose(b_valid, c_valid, rtol=rtol, atol=atol):
                            abs_diff = np.abs(b_valid - c_valid)
                            with np.errstate(divide='ignore', invalid='ignore'):
                                denom = np.maximum(np.abs(b_valid), atol)
                                rel_diff = abs_diff / denom
                            worst = np.argmax(abs_diff)
                            worst_step = np.where(valid)[0][worst]
                            quant_diffs.append(
                                f"Binary {bid}, '{col}': numeric mismatch "
                                f"(max abs diff = {abs_diff[worst]:.6e}, "
                                f"max rel diff = {rel_diff[worst]:.6e}, "
                                f"at step {worst_step}, "
                                f"baseline={b_valid[worst]:.15g}, "
                                f"candidate={c_valid[worst]:.15g})"
                            )

            else:
                # Qualitative comparison: exact string match
                b_str = b_col.astype(str).values
                c_str = c_col.astype(str).values
                mismatches = np.where(b_str != c_str)[0]
                if len(mismatches) > 0:
                    details = []
                    for s in mismatches[:5]:
                        details.append(f"step {s}: '{b_str[s]}' -> '{c_str[s]}'")
                    qual_diffs.append(
                        f"Binary {bid}, '{col}': {len(mismatches)} step(s) differ: "
                        + "; ".join(details)
                    )

    return {'quantitative': quant_diffs, 'qualitative': qual_diffs, 'structural': struct_diffs}


def compare_warnings_tables(base_df, cand_df):
    """Compare warning tables between baseline and candidate.

    Returns list of diff strings.
    """
    diffs = []

    if base_df is None and cand_df is None:
        return diffs
    if base_df is None:
        diffs.append(f"Candidate has {len(cand_df)} warning(s), baseline has none")
        return diffs
    if cand_df is None:
        diffs.append(f"Baseline has {len(base_df)} warning(s), candidate has none")
        return diffs

    if len(base_df) != len(cand_df):
        diffs.append(f"Total warning count differs (baseline={len(base_df)}, candidate={len(cand_df)})")

    # Per-binary warning comparison
    if 'binary_id' in base_df.columns and 'binary_id' in cand_df.columns:
        base_grouped = base_df.groupby('binary_id')
        cand_grouped = cand_df.groupby('binary_id')
        all_ids = sorted(set(base_df['binary_id'].unique()) | set(cand_df['binary_id'].unique()))

        for bid in all_ids:
            b_warnings = base_grouped.get_group(bid) if bid in base_grouped.groups else pd.DataFrame()
            c_warnings = cand_grouped.get_group(bid) if bid in cand_grouped.groups else pd.DataFrame()

            b_count = len(b_warnings)
            c_count = len(c_warnings)

            if b_count == 0 and c_count > 0:
                cats = c_warnings['category'].unique().tolist() if 'category' in c_warnings.columns else ['unknown']
                diffs.append(f"Binary {bid}: {c_count} NEW warning(s) in candidate ({', '.join(str(c) for c in cats)})")
            elif b_count > 0 and c_count == 0:
                diffs.append(f"Binary {bid}: {b_count} warning(s) REMOVED in candidate")
            elif b_count != c_count:
                diffs.append(f"Binary {bid}: warning count changed ({b_count} -> {c_count})")
            elif b_count > 0:
                # Same count — check if warning categories or messages changed
                if 'category' in b_warnings.columns and 'category' in c_warnings.columns:
                    b_cats = sorted(b_warnings['category'].astype(str).tolist())
                    c_cats = sorted(c_warnings['category'].astype(str).tolist())
                    if b_cats != c_cats:
                        diffs.append(f"Binary {bid}: warning categories changed ({b_cats} -> {c_cats})")

                if 'message' in b_warnings.columns and 'message' in c_warnings.columns:
                    b_msgs = sorted(b_warnings['message'].astype(str).tolist())
                    c_msgs = sorted(c_warnings['message'].astype(str).tolist())
                    if b_msgs != c_msgs:
                        diffs.append(f"Binary {bid}: warning messages changed")

    return diffs


def read_table_safe(store, key):
    """Read a table from HDFStore, returning None if it doesn't exist."""
    try:
        if key in store:
            return store[key]
    except Exception:
        pass
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Compare baseline and candidate binary evolution HDF5 files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
By default, uses EXACT comparison (any numeric difference is reported).
Use --loose to allow small floating-point tolerances (rtol=1e-12, atol=1e-15).
        """,
    )
    parser.add_argument("baseline", help="Path to baseline HDF5 file")
    parser.add_argument("candidate", help="Path to candidate HDF5 file")
    parser.add_argument("--loose", action="store_true",
                        help="Allow small floating-point tolerance (rtol=1e-12, atol=1e-15)")
    parser.add_argument("--rtol", type=float, default=None,
                        help="Override relative tolerance (default: 0, or 1e-12 with --loose)")
    parser.add_argument("--atol", type=float, default=None,
                        help="Override absolute tolerance (default: 0, or 1e-15 with --loose)")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Print extra diagnostic info")
    args = parser.parse_args()

    # Set tolerances
    if args.loose:
        rtol = args.rtol if args.rtol is not None else 1e-12
        atol = args.atol if args.atol is not None else 1e-15
    else:
        rtol = args.rtol if args.rtol is not None else 0.0
        atol = args.atol if args.atol is not None else 0.0

    for f in [args.baseline, args.candidate]:
        if not os.path.exists(f):
            print(f"ERROR: File not found: {f}", file=sys.stderr)
            sys.exit(2)

    quant_diffs = []
    qual_diffs = []
    struct_diffs = []
    warn_diffs = []

    try:
        with pd.HDFStore(args.baseline, mode='r') as base_store, \
             pd.HDFStore(args.candidate, mode='r') as cand_store:

            base_keys = set(base_store.keys())
            cand_keys = set(cand_store.keys())

            if args.verbose:
                print(f"Baseline keys:  {sorted(base_keys)}")
                print(f"Candidate keys: {sorted(cand_keys)}")

            # ── Evolution table ───────────────────────────────────────
            base_evol = read_table_safe(base_store, '/evolution')
            cand_evol = read_table_safe(cand_store, '/evolution')

            if base_evol is None and cand_evol is None:
                struct_diffs.append("Neither file contains an 'evolution' table")
            elif base_evol is None:
                struct_diffs.append("Baseline missing 'evolution' table")
            elif cand_evol is None:
                struct_diffs.append("Candidate missing 'evolution' table")
            else:
                if args.verbose:
                    n_base = base_evol['binary_id'].nunique() if 'binary_id' in base_evol.columns else '?'
                    n_cand = cand_evol['binary_id'].nunique() if 'binary_id' in cand_evol.columns else '?'
                    print(f"Baseline:  {n_base} binaries, {len(base_evol)} total rows")
                    print(f"Candidate: {n_cand} binaries, {len(cand_evol)} total rows")

                evol_results = compare_evolution_tables(base_evol, cand_evol, rtol, atol)
                quant_diffs.extend(evol_results['quantitative'])
                qual_diffs.extend(evol_results['qualitative'])
                struct_diffs.extend(evol_results['structural'])

            # ── Warnings table ────────────────────────────────────────
            base_warn = read_table_safe(base_store, '/warnings')
            cand_warn = read_table_safe(cand_store, '/warnings')
            warn_diffs.extend(compare_warnings_tables(base_warn, cand_warn))

            # ── Extra/missing top-level keys ──────────────────────────
            for k in sorted(base_keys - cand_keys):
                if k not in ['/evolution', '/warnings', '/metadata']:
                    struct_diffs.append(f"Table '{k}' missing in candidate")
            for k in sorted(cand_keys - base_keys):
                if k not in ['/evolution', '/warnings', '/metadata']:
                    struct_diffs.append(f"Table '{k}' extra in candidate")

    except Exception as e:
        print(f"ERROR reading HDF5 files: {e}", file=sys.stderr)
        sys.exit(2)

    # ── Report ────────────────────────────────────────────────────────────
    total_diffs = len(quant_diffs) + len(qual_diffs) + len(struct_diffs) + len(warn_diffs)
    tol_label = f"rtol={rtol}, atol={atol}" if rtol > 0 or atol > 0 else "EXACT (rtol=0, atol=0)"

    print("=" * 70)
    print("POSYDON Binary Validation — Comparison Report")
    print(f"  Baseline:    {args.baseline}")
    print(f"  Candidate:   {args.candidate}")
    print(f"  Tolerances:  {tol_label}")
    print("=" * 70)

    if struct_diffs:
        print(f"\n--- STRUCTURAL ({len(struct_diffs)}) ---")
        print("  (missing/extra binaries, step count changes, newly failing/passing, errors)\n")
        for d in struct_diffs:
            print(f"  - {d}")

    if qual_diffs:
        print(f"\n--- QUALITATIVE ({len(qual_diffs)}) ---")
        print("  (state, event, step name, SN type, interpolation class changes)\n")
        for d in qual_diffs:
            print(f"  - {d}")

    if quant_diffs:
        print(f"\n--- QUANTITATIVE ({len(quant_diffs)}) ---")
        print("  (any numeric value change)\n")
        for d in quant_diffs:
            print(f"  - {d}")

    if warn_diffs:
        print(f"\n--- WARNINGS ({len(warn_diffs)}) ---")
        print("  (new, removed, or changed warnings)\n")
        for d in warn_diffs:
            print(f"  - {d}")

    print("\n" + "=" * 70)
    if total_diffs == 0:
        print("RESULT: IDENTICAL — candidate matches baseline exactly.")
        sys.exit(0)
    else:
        print(f"RESULT: {total_diffs} DIFFERENCE(S) DETECTED")
        print(f"  Structural:   {len(struct_diffs)}")
        print(f"  Qualitative:  {len(qual_diffs)}")
        print(f"  Quantitative: {len(quant_diffs)}")
        print(f"  Warnings:     {len(warn_diffs)}")
        print("=" * 70)
        sys.exit(1)


if __name__ == "__main__":
    main()
