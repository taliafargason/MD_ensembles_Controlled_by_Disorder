import pandas as pd
import math
import argparse

# Argument parser for command line flags
parser = argparse.ArgumentParser()
parser.add_argument("--f", help="file_path", default="Step7c_hbond.csv")
args = parser.parse_args()

# Load the data from the provided file
file_path = args.f
data = pd.read_csv(file_path)

# Function to split and expand columns by both "@" and "_"
def split_columns(df):
    new_columns = {}
    
    for col in df.columns:
        # Split by "@"
        parts = df[col].str.split('@', expand=True)
        
        # Process each part further by splitting by "_"
        for i, part in enumerate(parts.columns):
            sub_parts = parts[part].str.split('_', expand=True)
            
            # Create new column names and collect them
            for j in sub_parts.columns:
                new_col_name = f'{col}_part{i+1}_sub{j+1}'
                new_columns[new_col_name] = sub_parts[j]
    
    # Concatenate all new columns into a single DataFrame
    new_df = pd.concat(new_columns, axis=1)
    
    return new_df

# Split columns by "@" and "_" and expand them into new sub-columns
split_data = split_columns(data)

# Save the processed data to a new CSV file
output_file_path = f'Processed_{file_path}'
split_data.to_csv(output_file_path, index=False)

# Print columns in split_data before filtering
print("Columns in split_data:", split_data.columns)

filtered_output_file_path = f'Filtered_{file_path}'

# Read the original file and filter its contents
with open(output_file_path, "r") as f, open(filtered_output_file_path, "w") as g:
    for line_number, line in enumerate(f):
        # Skip row 0, row 2, and row 3
        if line_number in {0, 2}:
            continue
        elif line_number ==1:
            columns = line.strip().split(",")
            filtered_columns =[]
            for j in range(len(columns)):
                offset = math.floor(j/14)*14
                conf_number = columns[offset + 1]
                if j % 14 in {1,5}:
                    filtered_columns.append(conf_number)
            g.write(",".join(filtered_columns) + "\n")
        else:
            columns = line.strip().split(",")
            filtered_columns = []
            for j in range(len(columns)):
                # Keep indexes 0, 1, 4, 5, 12, 13 in every group of 14 columns
                if j % 14 in {1, 5}:
                    filtered_columns.append(columns[j])
            g.write(",".join(filtered_columns) + "\n")

print("Filtered data has been saved to", filtered_output_file_path)


        