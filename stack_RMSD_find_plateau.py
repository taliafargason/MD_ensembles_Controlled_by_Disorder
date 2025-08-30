# -*- coding: utf-8 -*-
"""
Created on Sun May 19 22:30:58 2024

@author: tfarg
"""

import numpy as np
import pandas as pd
import os
import sys

# Calculate the slope
def calculate_slope(x, y):
    return np.diff(y) / np.diff(x)

# Identify plateaus
#def identify_plateaus(slope, threshold, min_plateau_length):
 #   plateaus = []
  #  current_plateau = None
   # for i, s in enumerate(slope):
    #    if abs(s) < threshold:
     #       if current_plateau is None:
      #          current_plateau = [i]
       #     else:
        #        current_plateau.append(i)
        #else:
         #   if current_plateau is not None and len(current_plateau) >= min_plateau_length:
          #      plateaus.append(current_plateau)
           # current_plateau = None
#    if current_plateau is not None and len(current_plateau) >= min_plateau_length:
 #       plateaus.append(current_plateau)
  #  return plateaus

# Find the final plateau
#def find_final_plateau(plateaus):
 #   if plateaus:
  #      return plateaus[-1]
   # else:
    #    return None

# Initialize an empty DataFrame
df = pd.DataFrame()

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
            os.system(f"cp SRSF1_{x}.prmtop {folder_path}")
        except (IndexError, ValueError):
            print(f"Skipping folder '{folder}' as it is not in the expected format.")
            continue
        
        os.chdir(folder_path)
        RMSDfile = f"SRSF1__1st_{x}.rmsd"
        
        try:
            # Read the RMSD file, extract the columns
            with open(RMSDfile, 'r') as file:
                headers = file.readline().strip().split()
                data = np.loadtxt(file)

            x_values = data[:, 0]
            y_values = data[:, 1]

            # Calculate slope and identify plateaus
            #slope = calculate_slope(x_values, y_values)
            #plateaus = identify_plateaus(slope, threshold=1, min_plateau_length=5)
            #final_plateau = find_final_plateau(plateaus)

            # Extract the final plateau values
            #if final_plateau is not None:
              #  final_x_values = x_values[final_plateau]
               # final_y_values = y_values[final_plateau]
            #else:
             #   final_x_values = np.array([])
              #  final_y_values = np.array([])

            # Create a DataFrame for the current file
            file_df = pd.DataFrame({
                f'{RMSDfile}_x': x_values,
                f'{RMSDfile}_y': y_values,
               # f'{RMSDfile}_plateau_x': np.concatenate([final_x_values, np.full(x_values.shape[0] - final_x_values.shape[0], np.nan)]),
            #    f'{RMSDfile}_plateau_y': np.concatenate([final_y_values, np.full(y_values.shape[0] - final_y_values.shape[0], np.nan)])
            })

            # Add filename as the first row
            filename_row = pd.DataFrame({col: [folder] for col in file_df.columns})
            headers_row = pd.DataFrame({col: [col] for col in file_df.columns})

            # Concatenate filename, headers, and data
            combined_df = pd.concat([filename_row, headers_row, file_df], ignore_index=True)

            # Concatenate with the main DataFrame
            df = pd.concat([df, combined_df], axis=1)

        except Exception as e:
            print(f"An error occurred while processing folder {folder}: {e}")
        
        finally:
            os.chdir(original_directory)  # Change back to the original directory
            print(f"Changing directory back to {original_directory}")

# Save the stacked columns to a new file
df.to_csv('Step7a_RMSD.csv', index=False, header=False, sep = ' ')
