# -*- coding: utf-8 -*-
"""
Created on Sat Apr 13 12:36:17 2024

@author: tfarg
"""

import os
import shutil

def replace_in_file(file_path, old_str, new_str):
    with open(file_path, 'r') as file:
        file_data = file.read()
    
    file_data = file_data.replace(old_str, new_str)
    
    with open(file_path, 'w') as file:
        file.write(file_data)

def process_directory(input_dir, output_dir, index):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for filename in os.listdir(input_dir):
        if os.path.isfile(os.path.join(input_dir, filename)):
            input_file_path = os.path.join(input_dir, filename)
            output_file_path = os.path.join(output_dir, filename.replace("_1", "_{}".format(index)))
            shutil.copy(input_file_path, output_file_path)
            replace_in_file(output_file_path, "_1", "_{}".format(index))

if __name__ == "__main__":
    for i in range(2, 21):
        input_directory = "R:\\CHEM_Zhang_Students\\tfarg\\Modelling\\TEMPLATES_FOR_AMBER\\ff99SB+OPC_ensemble_GPU_accelerated\\MD_1"
        output_directory = "R:\\CHEM_Zhang_Students\\tfarg\\Modelling\\TEMPLATES_FOR_AMBER\\ff99SB+OPC_ensemble_GPU_accelerated\\MD_{}".format(i)
        process_directory(input_directory, output_directory, i)

