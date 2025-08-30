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
import argparse

# Argument parser for command line flags
parser = argparse.ArgumentParser(description='Process input and output file names.')
parser.add_argument("-o", "--out", help="Output file name", required=True)
parser.add_argument("-f1", "--file1", help="Input file name 1", required=True)
parser.add_argument("-f2", "--file2", help="Input file name 2", required=True)
parser.add_argument("-f3", "--file3", help="Input file name 3", required=True)

args = parser.parse_args()

# Function to process each file and return the pivot table
def process_file(file_path):
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

    # Convert the pair percentages to a DataFrame for heatmap plotting
    result_data = []
    for (num1, num2), percent in pair_percentages.items():
        result_data.append([num1, num2, percent])

    result_df = pd.DataFrame(result_data, columns=['number1', 'number2', 'percent'])

    # Pivot the data to get the format suitable for a heatmap
    pivot_table = result_df.pivot("number2", "number1", "percent")

    return pivot_table

# Reflect the matrix
def reflect_matrix(matrix):
    return np.transpose(matrix)

# Overlay heatmaps
def overlay_heatmaps(matrix1, matrix2, label1, label2, output_file):
    # Create a figure
    plt.figure(figsize=(10, 8))
    
    # Plot the first matrix
    sns.heatmap(matrix1, cmap='Blues', alpha=0.5, cbar=False, square=True)
    
    # Overlay the second matrix
    sns.heatmap(matrix2, cmap='Reds', alpha=0.5, cbar=False, square=True)

    # Add labels
    plt.xlabel('Residue')
    plt.ylabel('Residue')

    # Create custom legend
    blue_patch = plt.Line2D([0], [0], color='blue', lw=4, label=f'{label1}')
    red_patch = plt.Line2D([0], [0], color='red', lw=4, label=f'{label2}')
    plt.legend(handles=[blue_patch, red_patch], loc='upper left', bbox_to_anchor=(1.0, 1.0))

    # Annotate points with value > 1
    for i in range(matrix1.shape[0]):
        for j in range(matrix1.shape[1]):
            if matrix1.iloc[i, j] > 1:
                plt.annotate(f"({j+1},{i+1})", xy=(j+0.5, i+0.5), xytext=(j+0.5, i-0.3), 
                             color='black', fontsize='medium', ha='center')

    for i in range(matrix2.shape[0]):
        for j in range(matrix2.shape[1]):
            if matrix2.iloc[i, j] > 1:
                plt.annotate(f"({j+1},{i+1})", xy=(j+0.5, i+0.5), xytext=(j+0.5, i-0.3), 
                             color='black', fontsize='medium', ha='center')

    plt.tight_layout()
    for ext in ['pdf', 'png', 'jpg']:
        plt.savefig(f"{output_file}_overlay.{ext}", format=ext, bbox_inches='tight')
    plt.show()

# Process each file
pivot_table1 = process_file(args.file1)
pivot_table2 = process_file(args.file2)
pivot_table3 = process_file(args.file3)

# Create subplots
fig, axes = plt.subplots(1, 3, figsize=(30, 10))
sns.set(font_scale=1.2)

# Plot heatmap for the first file
sns.heatmap(pivot_table1, 
            annot=False, 
            fmt=".2f", 
            cmap="PuOr",
            vmin=0,
            vmax=10,
            cbar_kws={'label': '% conformers with this interaction'},
            linewidths=0,
            ax=axes[0])

axes[0].set_xlabel('Residue')
axes[0].set_ylabel('Residue')
axes[0].set_title(f'File 1: {args.file1}')

# Plot heatmap for the second file
sns.heatmap(pivot_table2, 
            annot=False, 
            fmt=".2f", 
            cmap="PuOr",
            vmin=0,
            vmax=10,
            cbar_kws={'label': '% conformers with this interaction'},
            linewidths=0,
            ax=axes[1])

axes[1].set_xlabel('Residue')
axes[1].set_ylabel('Residue')
axes[1].set_title(f'File 2: {args.file2}')

# Plot heatmap for the third file
sns.heatmap(pivot_table3, 
            annot=False, 
            fmt=".2f", 
            cmap="PuOr",
            vmin=0,
            vmax=10,
            cbar_kws={'label': '% conformers with this interaction'},
            linewidths=0,
            ax=axes[2])

axes[2].set_xlabel('Residue')
axes[2].set_ylabel('Residue')
axes[2].set_title(f'File 3: {args.file3}')

# Adjust layout
plt.tight_layout()

# Save the plot as a PDF, PNG, and JPG
plt.savefig(f'heatmap_{args.out}.pdf', format='pdf')
plt.savefig(f'heatmap_{args.out}.png', format='png')
plt.savefig(f'heatmap_{args.out}.jpg', format='jpg')

# Show the plot
plt.show()

# Reflect the second pivot table
reflected_pivot_table2 = reflect_matrix(pivot_table2)

# Overlay heatmaps from the first and the reflected second pivot table
overlay_heatmaps(pivot_table1, reflected_pivot_table2, args.file1, args.file2, f'heatmap_{args.out}')
