import pandas as pd
import os
import sys
import re

# Initialize an empty DataFrame
final_data = pd.DataFrame()

# Determine the directory to use
if len(sys.argv) > 1:
    directory = sys.argv[1]
else:
    directory = os.getcwd()

# Store the original directory
original_directory = os.getcwd()

# Define the function to split column values and assign sub-column headers
def split_and_assign_headers(value):
    parts = re.split(r'[_@]', value)
    if len(parts) == 4:
        return parts
    else:
        return [None] * 4  # Return a list of Nones if the split doesn't produce 4 parts

# Process each folder in the directory
for folder in os.listdir(directory):
    folder_path = os.path.join(directory, folder)
    if os.path.isdir(folder_path):
        print(f"Changing directory to {folder_path}")
        try:
            x = int(folder.split("_")[1])  # Extract the conformer number
            os.system(f"cp SRSF1_{x}.prmtop {folder_path}")
        except (IndexError, ValueError):
            print(f"Skipping folder '{folder}' as it is not in the expected format.")
            continue
        
        os.chdir(folder_path)
        Hbond = f"SRSF1_{x}_Avg.hbond"
        
        with open(Hbond, 'r') as file:
            lines = file.readlines()
        # Process each line to remove whitespace and split into columns
        data = []
        for line in lines:
            # Use re.split() with a regular expression pattern to split based on whitespace
            columns = re.split(r'\s+', line.strip())
            data.append(columns)
        # Convert the list of lists to a DataFrame
        df = pd.DataFrame(data)
        
        # Split and assign headers for columns 0 and 1
        df[[0, 1]] = df[[0, 1]].applymap(split_and_assign_headers)
        
        # Expand columns 0 and 1 into separate columns
        expanded_0 = pd.DataFrame(df[0].tolist(), index=df.index)
        expanded_1 = pd.DataFrame(df[1].tolist(), index=df.index)
        
        # Concatenate expanded columns with the rest of the DataFrame
        df = pd.concat([expanded_0, expanded_1, df.iloc[:, 2:]], axis=1)
        
        # Select columns 0, 1, 3, 4, and 5
        selected_columns = df.iloc[:, [0, 1, 2, 3, 4, 5, 6, 9, 10, 11]]
        
        # Add headers for sub-columns
        headers = [
            (0, 'resn'), (0, 'resid'), (0, 'atom_n'), (0, 'atomID'),
            (1, 'resn'), (1, 'resid'), (1, 'atom_n'), (1, 'atomID'),
            (3, 'frames'), (4, 'fraction'), (5, 'avg_dist')
        ]
        selected_columns.columns = pd.MultiIndex.from_tuples(headers)
        
        # Insert a blank row after the second row for the columns starting from index 2 onwards
        blank_row = pd.DataFrame([[None] * len(selected_columns.columns)], columns=selected_columns.columns)
        selected_columns = pd.concat([selected_columns.iloc[:2], blank_row, selected_columns.iloc[2:]]).reset_index(drop=True)
        
        # Concatenate with the main DataFrame
        final_data = pd.concat([final_data, selected_columns], axis=1)

# Save the stacked columns to a new file
final_data.to_excel('Step7c_Hbonds.xlsx', index=False)
