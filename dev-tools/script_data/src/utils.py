from formatting import columns_to_show, line_length, target_rows


def print_warnings(captured_warnings):
    # Show warnings if any were captured
    if captured_warnings:
        print(f"⚠️  {len(captured_warnings)} warning(s) raised during evolution:")
        for i, warning in enumerate(captured_warnings[:3], 1):  # Show max 3 warnings
            print(f"   {i}. {warning['category']}: {warning['message']}")
        if len(captured_warnings) > 3:
            print(f"   ... and {len(captured_warnings) - 3} more warning(s)")
        elif len(captured_warnings) <= 3:
            for i in range(4-len(captured_warnings)):
                print("")
    else:
        print(f"No warning(s) raised during evolution\n\n")

def print_pop_settings(population):

    print("\nPopulation settings:")

    ignore_kwargs = ["extra_columns", "only_select_columns", "scalar_names",
                     "include_S1", "S1_kwargs", "include_S2", "S2_kwargs",
                     "population_properties", "warnings_verbose", "history_verbose",
                      "error_checking_verbose", "use_MPI", "read_samples_from_file",
                      "RANK", "size", "optimize_ram", "ram_per_cpu",
                      "dump_rate", "tqdm", "temp_directory", "breakdown_to_df"]

    for key, val in population.kwargs.items():
        if key in ignore_kwargs:
            continue
        else:
            print(f"\t {key} : {val}")

    print("\n")


def write_binary_to_screen(binary):
    """Writes a binary DataFrame prettily to the screen

    Args:
        binary: BinaryStar object with evolved data
    """
    df = binary.to_df(**{'extra_columns':{'step_names':'str'}})

    # Filter to only existing columns
    available_columns = [col for col in columns_to_show if col in df.columns]
    df_filtered = df[available_columns]

    # Reset index to use a counter instead of NaN
    df_filtered = df_filtered.reset_index(drop=True)

    print("=" * line_length)

    # Print the DataFrame
    df_string = df_filtered.to_string(index=True, float_format='%.3f')
    print(df_string)

    # Add empty lines to reach exactly 10 rows of output
    current_rows = len(df_filtered) + 1 # add one for header

    if current_rows < target_rows:
        # Calculate the width of the output to print empty lines of the same width
        lines = df_string.split('\n')
        if len(lines) > 1:
            # Use the width of the data lines (skip header)
            empty_lines_needed = target_rows - current_rows
            for i in range(empty_lines_needed):
                print("")

    print("-" * line_length)


def print_failed_binary(binary, e, max_error_lines=3):

    print("=" * line_length)
    print(f"🚨 Binary Evolution Failed!")
    print(f"Exception: {type(e).__name__}")
    print(f"Message: {e}")

    # Get the binary's current state and limit output
    try:
        df = binary.to_df(**{'extra_columns':{'step_names':'str'}})
        if len(df) > 0:
            # Select only the desired columns

            available_columns = [col for col in columns_to_show if col in df.columns]
            df_filtered = df[available_columns]

            # Reset index to use a counter instead of NaN
            df_filtered = df_filtered.reset_index(drop=True)

            # Limit to max_error_lines
            if len(df_filtered) > max_error_lines:
                df_filtered = df_filtered.tail(max_error_lines)
                print(f"\nShowing last {max_error_lines} evolution steps before failure:")
            else:
                print(f"\nEvolution steps before failure ({len(df_filtered)} steps):")

            df_string = df_filtered.to_string(index=True, float_format='%.3f')
            print(df_string)

            current_rows = len(df_filtered) + 1 + 5  # add one for header
            empty_lines_needed = target_rows - current_rows
            for i in range(empty_lines_needed):
                print("")
        else:
            print("\nNo evolution steps recorded before failure.")
    except Exception as inner_e:
        print(f"\nCould not retrieve binary state: {inner_e}")

    print("-" * line_length)
