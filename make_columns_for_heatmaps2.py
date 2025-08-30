# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 19:25:36 2024

@author: tfarg
"""

import pandas as pd
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns

# Load the CSV file
output_name='pi_pi_hypo'
file_path = 'Filtered_Step7c_pi-pi_15_frames_hypo.csv'  # Replace with your file path
data = pd.read_csv(file_path, header=None)

# Extract the header row
header = data.iloc[0]
data = data[1:]  # Remove the header row from the data

# Identify unique columns groups based on headers
unique_headers = header.unique()
column_groups = {col: header[header == col].index.tolist() for col in unique_headers}

# Initialize a dictionary to store counts of each pair
pair_counts = {pair: 0 for pair in combinations(range(1, 249), 2)}

# Function to check if a pair is in a row
def contains_pair(row, pair):
    return pair[0] in row and pair[1] in row


# Iterate through each column group and count pairs
for cols in column_groups.values():
    sub_data = data[cols].dropna(how='all')  # Drop rows where all values are NaN
    seen_pairs = set()  # Track pairs already counted in this group
    for _, row in sub_data.iterrows():
        row_numbers = set()
        for i in row:
            try:
                row_numbers.add(int(i))
            except ValueError:
                continue

        for pair in combinations(range(1, 249), 2):
            if pair[0] in row_numbers and pair[1] in row_numbers:
                if pair not in seen_pairs:  # Check if pair has already been counted
                    pair_counts[pair] += 1
                    seen_pairs.add(pair)  # Mark pair as counted

# Calculate the percentage of columns containing each pair
total_columns = len(column_groups)
pair_percentages = {pair: count / total_columns * 100 for pair, count in pair_counts.items()}

# Write the results to a text file
output_file = f'{output_name}_for_heatmap.txt'
with open(output_file, 'w') as f: 
    with open(f"parsimonius_{output_file}","w") as g:
        f.write("number1\tnumber2\tpercent\n")
        for (num1, num2), percent in pair_percentages.items():
            f.write(f"{num1}\t{num2}\t{percent:.2f}\n")
            if percent > 0:
               g.write(f"{num1}\t{num2}\t{percent:.2f}\n") 

# Load the data
data = pd.read_csv(output_file, delimiter='\t')

# Ensure the column names are correctly set
data.columns = ['number1', 'number2', 'percent']

# Pivot the data to get the format suitable for a heatmap
pivot_table = data.pivot("number2", "number1", "percent")

# Create the heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(pivot_table, 
            annot=False, 
            fmt=".2f", 
            cmap="PuOr",
            vmin = 0,
            vmax = 10,
            cbar_kws={'label': '% conformers with this interaction'},
            linewidths=0, 
            annot_kws={"size": 20})

# Add labels and title
plt.xlabel('SRSF1 Residue A')
plt.ylabel('SRSF1 Residue B')
plt.title(f'{output_name}')

# Save the plot as a PDF
plt.savefig(f'heatmap_{output_name}.pdf', format='pdf')

# Show the plot
plt.show()
