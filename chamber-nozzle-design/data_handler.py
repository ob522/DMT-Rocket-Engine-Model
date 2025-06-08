import os
import shutil
import re
import pandas as pd

def write_arrays_to_file(file_name, *arrays):
    """
    Writes the provided arrays to a text file in the 'results' folder, with each line containing one value from each array separated by a space.
    Creates a backup of the existing file before overwriting it.

    Parameters:
        file_name (str): The name of the output file.
        *arrays: Variable-length argument list of arrays (must all have the same length).

    Returns:
        None
    """
    # Ensure the 'results' folder exists
    results_folder = "results"
    os.makedirs(results_folder, exist_ok=True)

    # Prepend the folder path to the file name
    file_path = os.path.join(results_folder, file_name)

    # Check if the file exists
    if os.path.exists(file_path):
        backup_name = file_path + ".backup"
        shutil.copy(file_path, backup_name)  # Create a backup
        print(f"Existing file backed up as '{backup_name}'")

    # Validate that all arrays have the same length
    if not all(len(arr) == len(arrays[0]) for arr in arrays):
        raise ValueError("All input arrays must have the same length.")

    # Open the file for writing
    with open(file_path, 'w') as file:
        # Iterate over the indices of the arrays
        for i in range(len(arrays[0])):
            # Extract the i-th element from each array and join them with a space
            line = ' '.join(str(arr[i]) for arr in arrays)
            # Write the line to the file
            file.write(line + '\n')


def filter_and_extract_numbers(input_path: str, output_path: str) -> None:
    """
    Reads the file at input_path, removes the first 8 lines, skips any lines starting with '#',
    extracts all integer and floating-point numbers from each remaining line, and writes
    the cleaned, number-only lines to output_path.
    """
    # Compile regex to match signed/unsigned integers and floats
    number_pattern = re.compile(r"[-+]?\d*\.\d+|\d+")

    # Read all lines from the input file
    with open(input_path, 'r') as infile:
        lines = infile.readlines()

    # Remove the first 8 header lines before any other processing
    lines = lines[8:]

    filtered_lines = []
    for line in lines:
        # Skip any comment lines beginning with '#'
        if line.lstrip().startswith('#'):
            continue
        # Extract all numbers from the line
        numbers = number_pattern.findall(line)
        # If there are numbers, join them and add a newline
        if numbers:
            filtered_lines.append(' '.join(numbers) + '\n')

    # Write the processed data to the output file
    with open(output_path, 'w') as outfile:
        outfile.writelines(filtered_lines)

def split_data_columns(file_path: str):
    """
    Reads a cleaned data file (space-separated numeric columns) and returns seven lists:
      - x_values   = column 1
      - r_values   = column 2
      - q_values   = column 6
      - Twg_values = column 7
      - Twc_values = column 9
      - Tc_values  = column 10
      - Pc_values  = column 11

    :param file_path: path to the cleaned data file
    :return: tuple of lists (x_values, r_values, q_values, Twg_values, Twc_values, Tc_values, Pc_values)
    """
    x_values = []
    r_values = []
    q_values = []
    Twg_values = []
    Twc_values = []
    Tc_values = []
    Pc_values = []

    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            # skip any malformed lines
            if len(parts) < 11:
                continue

            # map 1-indexed columns onto 0-based Python lists
            x_values.append(float(parts[0]))   # column 1
            r_values.append(float(parts[1]))   # column 2
            q_values.append(float(parts[5]))   # column 6
            Twg_values.append(float(parts[6])) # column 7
            Twc_values.append(float(parts[8])) # column 9
            Tc_values.append(float(parts[9]))  # column 10
            Pc_values.append(float(parts[10])) # column 11

    return x_values, r_values, q_values, Twg_values, Twc_values, Tc_values, Pc_values


##############################################################
# Hot Fire data
##############################################################
"""
def split_tests_by_time(
    input_path: str,
    output_dir: str,
    time_ranges: list[tuple[float, float]],
    sheet_name: str = None
):
    
    # Splits the Excel dataset into multiple test segments based on time ranges.

    # Args:
    #     input_path (str): Path to the original Excel file.
    #     output_dir (str): Directory to save the split test files.
    #     time_ranges (list of tuples): List of (start_time, end_time) for each test.
    #     sheet_name (str): Optional sheet name if Excel file has multiple sheets.
    
    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_excel(input_path, sheet_name=sheet_name)

    for i, (start, end) in enumerate(time_ranges, start=1):
        test_df = df[(df['Time(s)'] >= start) & (df['Time(s)'] <= end)]
        output_path = os.path.join(output_dir, f"test{i}_raw.xlsx")
        test_df.to_excel(output_path, index=False)
        print(f"Saved test{i} to {output_path}")
"""
def split_tests_by_time(
    input_path: str,
    output_dir: str,
    time_ranges: list[tuple[float, float]],
    sheet_name: str = None
):
    os.makedirs(output_dir, exist_ok=True)

    # Handle multiple sheets by defaulting to the first one
    if sheet_name:
        df = pd.read_excel(input_path, sheet_name=sheet_name)
    else:
        df_dict = pd.read_excel(input_path, sheet_name=None)
        first_sheet_name = next(iter(df_dict))
        df = df_dict[first_sheet_name]

    df.columns = df.columns.str.strip()

    # Find the time column robustly
    time_col = next((col for col in df.columns if "time" in col.lower()), None)
    if not time_col:
        raise ValueError("No column containing 'time' found.")

    for i, (start, end) in enumerate(time_ranges, start=1):
        test_df = df[(df[time_col] >= start) & (df[time_col] <= end)]
        output_path = os.path.join(output_dir, f"test{i}_raw.xlsx")
        test_df.to_excel(output_path, index=False)
        print(f"Saved test{i} to {output_path}")

def find_closest_index(time_array, target_time):

    closest_index = min(range(len(time_array)),
                        key=lambda i: abs(time_array[i] - target_time))
    return closest_index