import pandas as pd
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
from matplotlib.colors import LinearSegmentedColormap

# Argument parser for command line flags
parser = argparse.ArgumentParser(description='Process input and output file names.')
parser.add_argument("-o", "--out", help="Output file name", required=True)
parser.add_argument("-f1", "--file1", help="Input file name 1", required=True)
parser.add_argument("-f2", "--file2", help="Input file name 2", required=True)
parser.add_argument("-f3", "--file3", help="Input file name 3", required=True)
parser.add_argument("-f4", "--file4", help="Input file name 4", required=True)

args = parser.parse_args()

# Function to process each file and return the pivot table
def process_file(file_path, output_csv=False, output_name=None):
    print(f"Processing file: {file_path}")
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
    print(f"Processed file: {file_path}")

    if output_csv and output_name:
        pivot_table.to_csv(output_name)
        print(f"Saved matrix to {output_name}")

    return pivot_table

# Reflect the matrix
def reflect_matrix(matrix):
    return np.transpose(matrix)

# Create a custom colormap
def create_custom_colormap():
    colors = ["white", "blue", "purple", "red", "orange"]
    cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)
    return cmap

# Overlay heatmaps
def overlay_heatmaps(matrix1, matrix2, matrix4, matrix_reflect, label1, label2, label4, output_file, cmap):
    # Create a figure
    plt.figure(figsize=(10, 8))
    
    # Plot the first matrix
    sns.heatmap(matrix1, cmap=cmap, alpha=0.5, cbar=False, square=True)
    
    # Overlay the second matrix
    sns.heatmap(matrix2, cmap=cmap, alpha=0.5, cbar=False, square=True)
    
    # Overlay the fourth matrix
    sns.heatmap(matrix4, cmap=cmap, alpha=0.5, cbar=False, square=True)
    
    # Overlay the reflected matrix without annotations
    sns.heatmap(matrix_reflect, cmap=cmap, alpha=0.5, cbar=False, square=True, annot=False)

    # Add labels
    plt.xlabel('Residue')
    plt.ylabel('Residue')

    # Create custom legend
    blue_patch = plt.Line2D([0], [0], color='blue', lw=4, label=f'{label1}')
    red_patch = plt.Line2D([0], [0], color='red', lw=4, label=f'{label2}')
    yellow_patch = plt.Line2D([0], [0], color='green', lw=4, label=f'{label4}')
    purple_patch = plt.Line2D([0], [0], color='purple', lw=4, label='Reflected')
    plt.legend(handles=[blue_patch, red_patch, yellow_patch, purple_patch], loc='upper left', bbox_to_anchor=(1.0, 1.0))

    # Annotate points with value > 1 for the first matrix
    for i in range(matrix1.shape[0]):
        for j in range(matrix1.shape[1]):
            if matrix1.iloc[i, j] > 1:
                plt.annotate(f"({j+1},{i+2})", xy=(j+0.5, i+0.5), xytext=(j+0.5, i+0.3), 
                             color='blue', fontsize='small', ha='center')

    # Annotate points with value > 1 for the second matrix
    for i in range(matrix2.shape[0]):
        for j in range(matrix2.shape[1]):
            if matrix2.iloc[i, j] > 1:
                plt.annotate(f"({j+1},{i+2})", xy=(j+0.5, i+0.5), xytext=(j+0.5, i+0.3), 
                             color='red', fontsize='small', ha='center')

    # Annotate points with value > 1 for the fourth matrix
    for i in range(matrix4.shape[0]):
        for j in range(matrix4.shape[1]):
            if matrix4.iloc[i, j] > 1:
                plt.annotate(f"({j+1},{i+2})", xy=(j+0.5, i+0.5), xytext=(j+0.5, i+0.3), 
                             color='green', fontsize='small', ha='center')

    # Add colorbar
    norm = plt.Normalize(0, 100)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, label='% conformers with this interaction')

    # Set axis limits and ticks
    plt.xlim(0, matrix1.shape[1])
    plt.ylim(0, matrix1.shape[0])
    plt.xticks(np.arange(0, matrix1.shape[1], step=10), labels=np.arange(0, matrix1.shape[1], step=10))
    plt.yticks(np.arange(0, matrix1.shape[0], step=10), labels=np.arange(0, matrix1.shape[0], step=10))

    plt.tight_layout()
    for ext in ['pdf', 'png', 'jpg']:
        plt.savefig(f"{output_file}_overlay.{ext}", format=ext, bbox_inches='tight')
        print(f"Saved {output_file}_overlay.{ext}")
    plt.show()

# Function to create group matrix and save as CSV
def create_group_matrix(file_path, output_name):
    print(f"Processing group matrix for file: {file_path}")
    data = pd.read_csv(file_path, header=None)

    # Extract the header row
    header = data.iloc[0]
    data = data[1:]  # Remove the header row from the data

    # Identify unique columns groups based on headers
    unique_headers = header.unique()
    column_groups = {col: header[header == col].index.tolist() for col in unique_headers}

    # Define the groups
    groups = {
        "1-13": range(1, 14),
        "14-90": range(14, 91),
        "91-120": range(91, 121),
        "121-195": range(121, 196),
        "196-248": range(196, 249),
        "50-51 and 85": [50, 51, 85],
        "19 and 56-58": [19,56,58],
        "87-88": range(87, 89),
        "29-39": range(29, 40),
        "65-75": range(65, 76)
    }

    # Initialize a dictionary to store counts of each group pair
    group_counts = {pair: 0 for pair in combinations(groups.keys(), 2)}

    # Iterate through each column group and count group pairs
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

            for group_pair in combinations(groups.keys(), 2):
                group1, group2 = groups[group_pair[0]], groups[group_pair[1]]
                if any(num in row_numbers for num in group1) and any(num in row_numbers for num in group2):
                    if group_pair not in seen_pairs:  # Check if pair has already been counted
                        group_counts[group_pair] += 1
                        seen_pairs.add(group_pair)  # Mark pair as counted

    # Calculate the percentage of columns containing each group pair
    total_columns = len(column_groups)
    group_percentages = {pair: count / total_columns * 100 for pair, count in group_counts.items()}

    # Convert the group percentages to a DataFrame
    result_data = []
    for (group1, group2), percent in group_percentages.items():
        result_data.append([group1, group2, percent])

    result_df = pd.DataFrame(result_data, columns=['Group1', 'Group2', 'Percent'])

    # Pivot the data to get the format suitable for a matrix
    pivot_table = result_df.pivot("Group2", "Group1", "Percent")

    # Save the group matrix to a CSV file
    pivot_table.to_csv(output_name)
    print(f"Saved group matrix to {output_name}")

# Main execution
def main():
    print("Starting script execution...")
    # Read the input files and process them
    matrix1 = process_file(args.file1, output_csv=True, output_name=f"{args.out}_{args.file1.split('.')[0]}.csv")
    matrix2 = process_file(args.file2, output_csv=True, output_name=f"{args.out}_{args.file2.split('.')[0]}.csv")
    matrix3 = process_file(args.file3, output_csv=True, output_name=f"{args.out}_{args.file3.split('.')[0]}.csv")
    matrix4 = process_file(args.file4, output_csv=True, output_name=f"{args.out}_{args.file4.split('.')[0]}.csv")

    # Reflect the third matrix
    matrix_reflect = reflect_matrix(matrix3)

    # Create a custom colormap
    cmap = create_custom_colormap()

    # Overlay the heatmaps
    overlay_heatmaps(matrix1, matrix2, matrix4, matrix_reflect, f'{args.file1}', f'{args.file2}', f'{args.file4}', args.out, cmap)

    # Create and save group matrices for each input file
    create_group_matrix(args.file1, f"{args.out}_group_matrix_{args.file1.split('.')[0]}.csv")
    create_group_matrix(args.file2, f"{args.out}_group_matrix_{args.file2.split('.')[0]}.csv")
    create_group_matrix(args.file3, f"{args.out}_group_matrix_{args.file3.split('.')[0]}.csv")
    create_group_matrix(args.file4, f"{args.out}_group_matrix_{args.file4.split('.')[0]}.csv")

    print("Script execution finished.")

if __name__ == "__main__":
    main()
