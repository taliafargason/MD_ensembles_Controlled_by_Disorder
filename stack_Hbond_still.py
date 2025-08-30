import numpy as np
import pandas as pd
import os
import sys
import re
import argparse

# Argument parser for command line flags
parser = argparse.ArgumentParser()
parser.add_argument("--conformers", help="number of conformers", default=30, type=int)
parser.add_argument("--directory", help="directory", default=os.getcwd())
parser.add_argument("--folder", help="name of folder with files to stack", default="hbond_endpoints")
parser.add_argument("--o", help="name of output file", default="Step7c_endpoint_hbonds.csv")
args = parser.parse_args()

conf = args.conformers
d = args.directory
f = args.folder
outputfile = args.o

# Initialize an empty DataFrame
final_data = pd.DataFrame()

# Change to the specified directory
os.chdir(d)
# Store the original directory
original_directory = os.getcwd()
# Identify the folder with the files
folder_path = os.path.join(d, f)
os.chdir(folder_path)

# Process each conformer file
for i in range(1, conf + 1):
    filename = f"avg_Hbond_endpoint_{i}.dat"
    confname = f"MD_{i}"
    if os.path.isfile(filename):
        with open(filename, "r") as file:
            lines = file.readlines()
            # Process each line to remove whitespace and split into columns
            data = []
            for line in lines:
                # Use re.split() with a regular expression pattern to split based on whitespace
                columns = re.split(r'\s+', line.strip())
                data.append(columns)
            # Convert the list of lists to a DataFrame
            df = pd.DataFrame(data)
            # Select columns 0, 1, 3, 4, and 5
            selected_columns = df.iloc[:, [0, 1, 3, 4, 5]]
            # Add filename as the first row
            filename_row = pd.DataFrame([[confname]*selected_columns.shape[1]], columns=selected_columns.columns)
            headers_row = pd.DataFrame([selected_columns.columns], columns=selected_columns.columns)
            # Concatenate filename, headers, and data
            combined_df = pd.concat([filename_row, headers_row, selected_columns], ignore_index=True)
            # Concatenate with the main DataFrame
            final_data = pd.concat([final_data, combined_df], axis=1)

# Change back to the original directory
os.chdir(original_directory)
# Save the stacked columns to a new file
final_data.to_csv(outputfile, index=False)