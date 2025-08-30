import pandas as pd
import os
import re
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pi_pi", type=str, default="pi_pi_stacks3_h.txt", help="name of pi-pi stacking input file (default = pi_pi_stacks3_h.txt)")
    parser.add_argument("--cat_pi", type=str, default="cation_pi_stacks_h.txt", help="name of cation-pi stacking input file (default = cation_pi_stacks_h.txt)")
    parser.add_argument("--arg_arg", type=str, default="arginine_stacks_h.txt", help="name of arg-arg stacking input file (default = arginine_stacks_h.txt)")
    parser.add_argument("--directory", type=str, default=os.getcwd(), help="Path to directory containing files to process (default = current directory)")
    parser.add_argument("--pi_pi_out",type=str, default="Step7c_pi-pi_stills.csv", help="name of pi-pi stacking output file (default = pi_pi_stacks3_h.txt)")
    parser.add_argument("--cat_pi_out",type=str, default="Step7c_cat-pi_stills.csv", help="name of cat-pi stacking output file (default = cat_pi_stacks3_h.txt)")
    parser.add_argument("--arg_arg_out",type=str, default="Step7c_arg_arg_stills.csv", help="name of arg-arg stacking output file (default = arg_arg_stacks3_h.txt)")
    return parser.parse_args()

def stack_stills(infile, outfile, directory):
    print(f"Executing stack_stills({infile}, {outfile})")
    final_data = pd.DataFrame()

    k = 1
    while True:
        folder = f"MD_{k}"
        folder_path = os.path.join(directory, folder)
        if not os.path.isdir(folder_path):
            break
        k += 1

        try:
            x = int(folder.split("_")[1])  # Extract the conformer number
            file_path = os.path.join(folder_path, infile)
            if os.path.exists(file_path):
                print(f"Processing file: {file_path}")
                with open(file_path, 'r') as file:
                    lines = file.readlines()

                data = []
                for line in lines[1:]:
                    columns = re.split(r'\s+', line.strip())
                    if len(columns) < 7:
                        print(f"Encountered incomplete line {columns}. Skipping.")
                        continue  # Skip incomplete lines
                    columns_stripped = [columns[0].split(".")[1], columns[1].split("-")[0], columns[2].split(".")[1], columns[3].split(":")[0], columns[4].split("=")[1], columns[6].split("=")[1]]
                    data.append(columns_stripped)
                    print(f"Stripped columns: {columns_stripped}")

                df = pd.DataFrame(data)
                if not df.empty:
                    headers = [(f'conf{x}', 'resn_a'), (f'conf{x}', 'resid_a'), (f'conf{x}', 'resn_b'), (f'conf{x}', 'resid_b'), (f'conf{x}', 'dist(Ang)'), (f'conf{x}', 'ang(Deg)')]
                    df.columns = pd.MultiIndex.from_tuples(headers)
                    final_data = pd.concat([final_data, df], axis=1)
                    print(f"{len(df)} interaction(s) found in {infile}. Added to final DataFrame.")
                else:
                    # Create a DataFrame with the correct headers and a single row of blank space
                    print(f"No interactions found in {file_path}. Creating a single-row blank DataFrame.")
                    selected_columns = pd.DataFrame([["", "", "", "", "", ""]], columns=pd.MultiIndex.from_tuples([(f'conf{x}', 'resn_a'), (f'conf{x}', 'resid_a'), (f'conf{x}', 'resn_b'), (f'conf{x}', 'resid_b'), (f'conf{x}', 'dist(Ang)'), (f'conf{x}', 'ang(Deg)')]))
                    final_data = pd.concat([final_data, selected_columns], axis=1)
            else:
                print(f"File {file_path} does not exist.")
        except (IndexError, ValueError) as e:
            print(f"Skipping folder '{folder}' due to error: {e}")
            continue

    final_data.to_csv(outfile, index=False)
    print(f"Saved final data to {outfile}")

if __name__ == "__main__":
    args = parse_arguments()
    stack_stills(args.pi_pi, args.pi_pi_out, args.directory)
    stack_stills(args.cat_pi, args.cat_pi_out, args.directory)
    stack_stills(args.arg_arg, args.arg_arg_out, args.directory)
