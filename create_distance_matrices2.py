import numpy as np
import pandas as pd
import os
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def main(num_conformers, num_residues, output_file, individual_output_file):
    distance_matrices = []
    sum_r_minus_6 = np.zeros((num_residues, num_residues))
    all_files_correct = True  # Flag to check if all files are correct

    # Loop to load and process distance matrices
    for i in range(1, num_conformers + 1):
        # Assuming the files are named distances_residues_1.dat, distances_residues_2.dat, ..., distances_residues_30.dat
        filename = f'distances_residues_{i}.dat'
        os.system(f'cp MD_{i}/{filename} {os.getcwd()}')
        if os.path.exists(filename):
            # Load the data
            data = pd.read_csv(filename, delim_whitespace=True, header=None)

            # Check the size of the data
            expected_size = num_residues * num_residues
            if data.size != expected_size:
                print(f"Error: File {filename} does not contain the expected number of elements ({expected_size} elements expected).")
                all_files_correct = False
                break

            # Reshape to a num_residues x num_residues matrix
            distance_matrix = data.values.reshape((num_residues, num_residues))

            # Accumulate r^-6 values
            sum_r_minus_6 += np.power(distance_matrix, -6)

            # Append to the list of distance matrices
            distance_matrices.append(distance_matrix)
        else:
            print(f"Error: File {filename} not found.")
            all_files_correct = False
            break

    if not all_files_correct:
        print("Aborting due to incorrect file sizes or missing files.")
        sys.exit(1)

    # Stack matrices into a 3D array (num_conformers, num_residues, num_residues)
    distance_matrices = np.array(distance_matrices)

    # Compute average distance matrix
    average_matrix = np.mean(distance_matrices, axis=0)

    # Compute ((sum(r^-6))^(-1/6)) matrix
    sum_r_minus_6_averaged = np.power(sum_r_minus_6 / num_conformers, -1/6)

    # Plot heatmap of the average distance matrix
    plt.figure(figsize=(12, 5))

    # Plot average distance matrix
    plt.subplot(1, 2, 1)
    sns.heatmap(average_matrix, cmap='viridis', cbar=True, square=True)
    plt.title('Average Distance Matrix')
    plt.xlabel('Residue')
    plt.ylabel('Residue')
    plt.xticks(np.arange(0, num_residues, 10), np.arange(0, num_residues, 10))
    plt.yticks(np.arange(0, num_residues, 10), np.arange(0, num_residues, 10))

    # Plot ((sum(r^-6))^(-1/6)) matrix
    plt.subplot(1, 2, 2)
    sns.heatmap(sum_r_minus_6_averaged, cmap='viridis', cbar=True, square=True)
    plt.title('((sum(r^-6))^(-1/6)) Matrix')
    plt.xlabel('Residue')
    plt.ylabel('Residue')
    plt.xticks(np.arange(0, num_residues, 10), np.arange(0, num_residues, 10))
    plt.yticks(np.arange(0, num_residues, 10), np.arange(0, num_residues, 10))

    plt.tight_layout()
    for ext in ['pdf', 'png', 'jpg']:
        plt.savefig(f"{output_file}.{ext}", format=ext)
    plt.show()

    # Create a new figure for the individual heatmaps
    num_rows = 5
    num_cols = 6
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(24, 20))
    axes = axes.flatten()

    for i, matrix in enumerate(distance_matrices):
        sns.heatmap(matrix, cmap='viridis', cbar=False, square=True, ax=axes[i])
        axes[i].set_title(f'Conformer {i + 1}')
        axes[i].set_xlabel('Residue')
        axes[i].set_ylabel('Residue')
        axes[i].set_xticks(np.arange(0, num_residues, 10))
        axes[i].set_yticks(np.arange(0, num_residues, 10))

    # Hide any unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    for ext in ['pdf', 'png', 'jpg']:
        plt.savefig(f"{individual_output_file}.{ext}", format=ext)
    plt.show()

if __name__ == "__main__":
    # Set default values
    DEFAULT_NUM_CONFORMERS = 30
    DEFAULT_NUM_RESIDUES = 248

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Process distance matrices.')
    parser.add_argument('--num_conformers', type=int, default=DEFAULT_NUM_CONFORMERS,
                        help='Number of conformers to process (default: 30)')
    parser.add_argument('--num_residues', type=int, default=DEFAULT_NUM_RESIDUES,
                        help='Number of residues per trajectory (default: 248)')
    parser.add_argument('--o', type=str, default='distance_matrices', help='Output file name for average heatmaps (default: distance_matrices)')
    parser.add_argument('--io', type=str, default='individual_heatmaps', help='Output file name for individual heatmaps (default: individual_heatmaps)')
    args = parser.parse_args()

    # Call the main function with command-line arguments
    main(args.num_conformers, args.num_residues, args.o, args.io)
