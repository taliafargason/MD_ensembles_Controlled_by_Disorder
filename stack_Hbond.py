# -*- coding: utf-8 -*-
"""
Created on Sun May 19 22:30:58 2024

@author: tfarg
"""

import numpy as np
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

# Process each folder in the directory
for folder in os.listdir(directory):
    folder_path = os.path.join(directory, folder)
    if os.path.isdir(folder_path):
        print(f"Changing directory to {folder_path}")
        try:
            x = int(folder.split("_")[1])  # Extract the conformer number
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
        # Select columns 0, 1, 3, 4, and 5
        selected_columns = df.iloc[:, [0, 1, 3, 4, 5]]
        # Uncomment the line below if you want to directly stack the selected columns to final_data
        #final_data = pd.concat([final_data, selected_columns], axis=1)
        # Add filename as the first row
        filename_row = pd.DataFrame({col: [folder] for col in selected_columns})
        headers_row = pd.DataFrame({col: [col] for col in selected_columns})

        # Concatenate filename, headers, and data
        combined_df = pd.concat([filename_row, headers_row, selected_columns], ignore_index=True)

        # Concatenate with the main DataFrame
        final_data = pd.concat([final_data, combined_df], axis=1)
os.chdir(original_directory)
# Save the stacked columns to a new file
#final_data.to_excel('Step7c_Hbonds.csv', index=False)
final_data.to_csv('Step7c_hbond.csv', index=False)
