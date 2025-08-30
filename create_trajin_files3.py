import os
import sys

# Determine the directory to use
if len(sys.argv) > 1:
    directory = sys.argv[1]
else:
    directory = os.getcwd()

# Store the original directory
original_directory = os.getcwd()
os.system(f"mkdir hbond_endpoints")
# Process each folder in the directory
for folder in os.listdir(directory):
    folder_path = os.path.join(directory, folder)
    if os.path.isdir(folder_path):
        print(f"Changing directory to {folder_path}")
        try:
            x = int(folder.split("_")[1])  # Extract the conformer number
            os.system(f"cp SRSF1_{x}.prmtop {folder_path}")
        except (IndexError, ValueError):
            print(f"Skipping folder '{folder}' as it is not in the expected format.")
            continue
        os.system(f"cp find-pi* {folder_path}")
        print(f"copying find-pi* to {folder_path}")
        
        os.chdir(folder_path)
        os.system("mkdir create_endpoint_hbond_data")
        os.system("cp ../OP1_to_O1P* create_endpoint_hbond_data")
        os.system("cp *final* create_endpoint_hbond_data")
        os.system("pwd")
        trajin_file_path = f"mk_binpos_{x}.txt"
        jobfile_path = f"trajectory_{x}.job.txt"
        RMSDfile_path = f"RMSD_analysis_{x}.txt"
        save_parm = f"{folder_path}/create_endpoint_hbond_data/save_parm_{x}.in"
        hbond_end = f"{folder_path}/create_endpoint_hbond_data/Hbond_endpoint_{x}.in"
        hbond_job = f"{folder_path}/create_endpoint_hbond_data/Hbond_endpoint_{x}.job.txt"
        try:
            with open(trajin_file_path, "w") as trajin_file:
                trajin_file.write(f"parm SRSF1-SolvI_{x}.prmtop\n")
                mdcrd = f"SRSF1-md_{x}.mdcrd"
                if os.path.isfile(mdcrd):
                    trajin_file.write(f"trajin {mdcrd}\n")
                else:
                    print(f"Difficulty adding {mdcrd} to trajin file")
                
                y = 1
                while True:
                    ext_file = f"SRSF1-ext{y}_{x}.mdcrd"
                    print(f"looking for file SRSF1-ext{y}_{x}.mdcrd")
                    if os.path.isfile(ext_file):
                        trajin_file.write(f"trajin {ext_file}\n")
                        y += 1
                        print(f"found file")
                    else:
                        print(f"did not find file")
                        break
                trajin_file.write("autoimage\n")
                trajin_file.write("rms fit :1-248\n")
                trajin_file.write("strip :WAT,Na+,Cl-\n")
               # trajin_file.write(f"trajout SRSF1_{x}_ensemble.pdb pdb nobox\n")
                trajin_file.write(f"trajout SRSF1_{x}.binpos binpos nobox\n")
            with open(RMSDfile_path, "w") as RMSD:
                RMSD.write(f"parm SRSF1_{x}.prmtop\n")
                RMSD.write(f"trajin SRSF1_{x}.binpos\n")
                #try:
                 #   os.makedirs(f"Output_{x}", exist_ok=True)
                #except FileExistsError:
                 #   pass  # Directory already exists, so no need to do anything
                #RMSD.write(f"rms first out Output_{x}/SRSF1__1st_{x}.rmsd @CA time 0.01\n")
                RMSD.write(f"rms first out SRSF1__1st_{x}.rmsd @CA time 0.01\n")
                RMSD.write("run\n")
            with open(jobfile_path, "w") as jobfile, \
                open(save_parm,"w") as parm, \
                open(hbond_end,"w") as hbond, \
                open(hbond_job,"w") as h_job:
                jobfile.write("#!/bin/bash\n")
                jobfile.write("#SBATCH --partition=express\n")
                jobfile.write(f"#SBATCH --job-name=trajectory_{x}\n")
                jobfile.write(f"#SBATCH --error=trajectory_{x}.err\n")
                jobfile.write(f"#SBATCH --output=trajectory_{x}.out\n")
                jobfile.write("#SBATCH --ntasks=1\n")
                jobfile.write("#SBATCH --time=2:00:00\n")
                jobfile.write("#SBATCH --mem-per-cpu=2G\n")
                jobfile.write("#SBATCH --mail-type=ALL\n")
                jobfile.write("#SBATCH --mail-user=tfarg@uab.edu\n")
                jobfile.write("module reset\n")
                jobfile.write("module load OpenMPI/4.1.5-GCC-12.3.0\n")
                jobfile.write("module load Amber/20.11-fosscuda-2020b-AmberTools-21.3\n")
                jobfile.write("module load CUDA/11.6.0\n\n")
                
                h_job.write("#!/bin/bash\n")
                h_job.write("#SBATCH --partition=express\n")
                h_job.write(f"#SBATCH --job-name=endpoint_hbond_{x}\n")
                h_job.write(f"#SBATCH --error=endpoint_hbond_{x}.err\n")
                h_job.write(f"#SBATCH --output=endpoint_hbond_{x}.out\n")
                h_job.write("#SBATCH --ntasks=1\n")
                h_job.write("#SBATCH --time=2:00:00\n")
                h_job.write("#SBATCH --mem-per-cpu=2G\n")
                h_job.write("#SBATCH --mail-type=ALL\n")
                h_job.write("#SBATCH --mail-user=tfarg@uab.edu\n")
                h_job.write("module reset\n")
                h_job.write("module load Amber/20.11-foss-2020b-AmberTools-21.3\n\n")
                
                parm.write("source leaprc.protein.ff14SB\n")
                parm.write("source leaprc.RNA.OL3\n")
                parm.write("source leaprc.phosaa14SB\n")
                parm.write("source leaprc.water.tip3p\n\n")
                
                if y > 1:
                    jobfile.write(f"ambpdb -c SRSF1-ext{y-1}_{x}.rst -p SRSF1-SolvI_{x}.prmtop > SRSF1-ext{y-1}_{x}_final.pdb\n")
                    final = f"SRSF1-ext{y-1}_{x}_final"

                else:
                    jobfile.write(f"ambpdb -c SRSF1-md_{x}.rst -p SRSF1-SolvI_{x}.prmtop > SRSF1-md_{x}_final.pdb\n")
                    final = f"SRSF1-md_{x}_final"
                    
                h_job.write(f"cp {folder_path}/{final}.pdb {original_directory}/create_endpoint_hbond_data\n")
                parm.write(f'mol = loadpdb "modified_{final}.pdb"\n')
                parm.write(f'saveamberparm mol {final}.prmtop {final}.inpcrd\n')
                parm.write("quit")
                
                hbond.write(f"parm modified_modified_{final}.pdb\n")
                hbond.write(f"trajin modified_modified_{final}.pdb\n")
                hbond.write(f"hbond out Hbond_endpoint_{x}.dat avgout avg_Hbond_endpoint_{x}.dat series printatomnum\n")
                hbond.write("run")
                
                h_job.write("module load Anaconda3\nsource activate Talias_modules\n\n")
                h_job.write("python3 OP1_to_O1P.py\n")
                h_job.write(f"tleap -f {save_parm}\n")
                h_job.write(f"ambpdb -c {final}.inpcrd -p {final}.prmtop > modified_modified_{final}.pdb\n")
                h_job.write(f"cpptraj -i {hbond_end}\n")
                h_job.write(f"cp Hbond_endpoint_{x}.dat {original_directory}/hbond_endpoints\n")
                h_job.write(f"cp avg_Hbond_endpoint_{x}.dat {original_directory}/hbond_endpoints\n")
                
                jobfile.write(f"cpptraj -i {trajin_file_path}\n")
                jobfile.write(f"cpptraj -i {RMSDfile_path}\n")
                jobfile.write("module purge\n")
                jobfile.write("module reset\n")
                jobfile.write("module load Anaconda3\n")
                jobfile.write("source activate Talias_modules\n")
                if y > 1:
                   jobfile.write(f"python3 find-pi4.py SRSF1-ext{y-1}_{x}_final.pdb -r arginine_stacks_i.txt -p pi_pi_stacks3_i.txt -cp cation_pi_stacks_i.txt >find-pi3_i.out\n")
                else:
                    jobfile.write(f"python3 find-pi4.py SRSF1-md_{x}_final.pdb -r arginine_stacks_i.txt -p pi_pi_stacks3_i.txt -cp cation_pi_stacks_i.txt >find-pi3_i.out\n") 
                jobfile.write(f"\ncd {folder_path}/create_endpoint_hbond_data\nsbatch {hbond_job}")
            os.system(f"dos2unix {jobfile_path}")
            os.system(f"sbatch {jobfile_path}")
        except Exception as e:
            print(f"An error occurred while processing folder {folder}: {e}")
        finally:
            os.chdir(original_directory)  # Change back to the original directory
            print(f"Changing directory to {original_directory}")
