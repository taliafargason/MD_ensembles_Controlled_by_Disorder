import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Reading the data from the file
file_path = "Step7a_RMSD.csv"
df = pd.read_csv(file_path, dtype=str)  # Read all data as strings initially

# Convert all columns to numeric, forcing non-numeric values to NaN
df = df.apply(pd.to_numeric, errors='coerce')

# Fill missing values if necessary
df.fillna(method='ffill', inplace=True)

# Identifying the column clusters
column_clusters = [df.columns[i:i+4] for i in range(0, len(df.columns), 4)]

# Extract the numerical part of the title for sorting
sorted_clusters = sorted(column_clusters, key=lambda cluster: int(cluster[1].split("_")[1].split(".")[0]))

# Set up the number of rows and columns for the subplot grid
n_rows = 6
n_cols = 5

# Creating scatterplots and saving them to a PDF, PNG, and JPG
with PdfPages("RMSD_plots.pdf") as pdf:
    num_clusters = len(sorted_clusters)
    num_pages = (num_clusters + (n_rows * n_cols) - 1) // (n_rows * n_cols)  # Ceiling division

    for page in range(num_pages):
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 10))
        axes = axes.flatten()
        
        for i in range(n_rows * n_cols):
            cluster_idx = page * (n_rows * n_cols) + i
            if cluster_idx >= num_clusters:
                axes[i].axis('off')
                continue
            
            cluster = sorted_clusters[cluster_idx]
            x_col, y_col = cluster[0], cluster[1]
            conf = y_col.split("_")[1].split(".")[0]
            title = f"Conformer {conf}"
            
            axes[i].scatter(df[x_col], df[y_col])
            axes[i].set_xlabel("Simulation time (ns)")
            axes[i].set_ylabel("RMSD (\u00C5)")
            axes[i].set_title(title)
        
        plt.tight_layout()
        pdf.savefig(fig)
        fig.savefig(f"RMSD_plots_page_{page + 1}.png", format='png')
        fig.savefig(f"RMSD_plots_page_{page + 1}.jpg", format='jpg')
        plt.close(fig)

print("PDF, PNG, and JPG creation complete.")
