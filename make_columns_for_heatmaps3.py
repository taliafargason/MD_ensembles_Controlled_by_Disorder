# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 16:05:11 2024

@author: tfarg
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 19:25:36 2024

@author: tfarg
"""

import pandas as pd
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load the CSV file
output_name = 'arg_arg_hypo_test2'
file_path = 'Step7c_arg-arg_stills.csv'  # Replace with your file path
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
    for _, row in sub_data.iterrows():
        row_numbers = set()
        for i in row:
            try:
                row_numbers.add(int(i))
            except ValueError:
                continue

        for pair in combinations(range(1, 249), 2):
            if pair[0] in row_numbers and pair[1] in row_numbers:
                pair_counts[pair] += 1

# Calculate the percentage of columns containing each pair
total_columns = len(column_groups)
pair_percentages = {pair: count / total_columns * 100 for pair, count in pair_counts.items()}

# Write the results to a text file
output_file = f'{output_name}_for_heatmap.txt'
with open(output_file, 'w') as f:
    f.write("number1\tnumber2\tpercent\n")
    for (num1, num2), percent in pair_percentages.items():
        f.write(f"{num1}\t{num2}\t{percent:.2f}\n")

# Load the data
data = pd.read_csv(output_file, delimiter='\t')

# Ensure the column names are correctly set
data.columns = ['number1', 'number2', 'percent']

# Bin the percentage values into discrete categories (multiples of 5)
data['percent_binned'] = np.floor(data['percent'] / 5) * 5

# Create a color palette for multiples of 5% (0%, 5%, 10%, ..., 100%)
colors = sns.color_palette("YlGnBu", n_colors=21)  # 21 colors for 0%, 5%, ..., 100%
color_dict = {i*5: colors[i] for i in range(21)}

# Map the binned percent values to colors
data['color'] = data['percent_binned'].map(color_dict)

# Pivot the data to get the format suitable for a heatmap
pivot_table = data.pivot("number2", "number1", "color")

# Create the heatmap
plt.figure(figsize=(10, 8))

# Plot heatmap with custom colors
sns.heatmap(pivot_table, 
            annot=False,  # Do not display annotations
            fmt="",
            cmap=sns.color_palette(colors, as_cmap=True),
            cbar=False,  # No color bar since we are using discrete colors
            linewidths=0)

# Add labels and title
plt.xlabel('SRSF1 Residue A')
plt.ylabel('SRSF1 Residue B')
plt.title('Arg-Arg stacking between SRSF1 Residues')

# Save the plot as a PDF
plt.savefig(f'heatmap_{output_name}.pdf', format='pdf')

# Show the plot
plt.show()
