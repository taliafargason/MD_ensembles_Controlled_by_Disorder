import os

# Directory containing PDB files
pdb_directory = os.getcwd()

# Iterate through PDB files in the directory
for filename in os.listdir(pdb_directory):
    if filename.endswith(".pdb"):
        input_file = os.path.join(pdb_directory, filename)
        output_file = os.path.join(pdb_directory, "modified_" + filename)

        # Check if the file is empty
        if os.path.getsize(input_file) == 0:
            print(f"Skipping empty file: {input_file}")
            continue

        with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
            for line in f_in:
                # Check for "HIE" or "HIS" in the line
                contains_hie_his_thr_phe = "HIE" in line or "HIS" in line or "THR" in line or "PHE" in line or "HETATM" in line
                # Count instances of H{i}
                h_count = sum(1 for i in "ABCDEFGHIJKLMNOPQRSTUVWXYZ" if f"H{i}" in line)

                # Skip the line if:
                # - It contains "HIE" or "HIS" and there are at least 2 instances of H{i}
                # - It does not contain "HIE" or "HIS" and there is at least 1 instance of H{i}
                if (contains_hie_his_thr_phe and h_count >= 2) or (not contains_hie_his_thr_phe and h_count >= 1):
                    continue

                # Skip lines containing "TAU" or "ZP"
                if "TAU" in line or "ZP" in line or "CRYST1" in line or "CONECT" in line:
                    continue
                
                
                # Replace "OP1" with "O1P", "OP2" with "O2P", and "CYSP" with "CYS "
                modified_line = line.replace("OP1", "O1P").replace("OP2", "O2P").replace("CYSP", "CYS ").replace("HETATM", "ATOM  ")
                f_out.write(modified_line)

        print(f"Modified PDB file saved: {output_file}")
