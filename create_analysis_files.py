import os
import sys

#### Note: This script requires a space-delimited file called "equilibrium_regions.txt"
#### in which column1=conformer#, column2=eq_start, column3=eq_end

# Determine the directory to use
if len(sys.argv) > 1:
    directory = sys.argv[1]
else:
    directory = os.getcwd()

# Store the original directory
original_directory = os.getcwd()

# Open a list of equilibrium regions in which column1=conformer#, column2=eq_start, column3=eq_end
with open("unified_equilibrium_regions.txt", "r") as equilibrium_regions:
    eq_reg = [line.strip().split() for line in equilibrium_regions]

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
        os.chdir(folder_path)
        
        trajin_file_path = f"mk_binpos_{x}.txt"
        jobfile_path = f"eq_analysis_{x}.job.txt"
        RMSDfile_path = f"RMSD_analysis_{x}.txt"
        basic = f"BasicAnalysis_{x}.in"
        eq = f"EquilibratedAnalysis_{x}.in"
        eps = f"find_epsilon_{x}.in"
        cluster = f"Cluster_{x}.in"
        
        try:
            # Find the start and end equilibrium regions from the eq_reg list
            eq_start, eq_end = None, None
            for line in eq_reg:
                conformer = line[0]
                if float(conformer) == float(x):
                    eq_start = line[1]
                    eq_end = line[2]
                    break

            if eq_start is None or eq_end is None:
                print(f"No equilibrium region found for conformer {x}")
                continue

            with open(basic, "w") as Basic:
                Basic.write("##Trajectory Building##\n")
                Basic.write(f"parm SRSF1_{x}.prmtop\n")
                Basic.write(f"trajin SRSF1_{x}.binpos\n\n")
                Basic.write(f"trajout SRSF1_{x}-eq.binpos onlyframes {eq_start}-{eq_end}\n")
                Basic.write(f"trajout SRSF1_{x}-eq.pdb pdb onlyframes {eq_start}-{eq_end}\n")
                Basic.write(f"rms first out SRSF1_1st_{x}b.rmsd @CA time 0.01\n")
                Basic.write("average crdset SRSF1_Avg\n")
                Basic.write("createcrd SRSF1_Traj\nrun\n\n")
                Basic.write(f"crdout SRSF1_Traj SRSF1-movie.binpos crdframes {eq_start},{eq_end}\nrun\n\n")
                Basic.write("##Basic Analysis##\n#dynamic correlation of each residue#\n")
                Basic.write("crdaction SRSF1_Traj atomicfluct \\")
                Basic.write(f"\n  SRSF1_{x}.rmsf byres \\")
                Basic.write(f"\n  crdframes {eq_start},{eq_end} @CA\nrun\n\n")
                Basic.write("#average distance between any two residues (remove 'byres' to get atoms instead)#\n")
                Basic.write("crdaction SRSF1_Traj matrix \\")
                Basic.write(f"\n  SRSF1_{x}.Dist byres \\")
                Basic.write(f"\n  crdframes {eq_start},{eq_end} @CA\nrun\n\n")
                Basic.write("#distance between residues 220 and 72 over time#\n")
                Basic.write("crdaction SRSF1_Traj distance \\")
                Basic.write(f"\n  :124 :149 out SRSF1-wing-tip.agr \\")
                Basic.write(f"\nrun\n\n")
                Basic.write("#solvent accessible surface area\n")
                Basic.write(f"crdaction SRSF1_Traj molsurf out SRSF1_{x}.SASA crdframes {eq_start},{eq_end}\nrun\n\n")
                Basic.write("##Principal Component Analysis##\n")
                Basic.write("#generate covariance matrix#\n")
                Basic.write("crdaction SRSF1_Traj matrix covar \\")
                Basic.write(f"\n  name SRSF1_Covar \\")
                Basic.write(f"\n  crdframes {eq_start},{eq_end} !:WAT&!@H=\n")
                Basic.write("#diagonalize to generate first 10 PCs#\n")
                Basic.write(f"runanalysis diagmatrix SRSF1_Covar out SRSF1_{x}_EigVecs.dat \\")
                Basic.write(f"\n  vecs 50 name SRSF1_Eigvecs \\")
                Basic.write(f"\n nmwiz nmwizvecs 50 nmwizfile SRSF1_{x}.nmd nmwizmask !:WAT&!@H=\n")
                Basic.write("#project the pcs onto our trajectory and generate histograms of the first three PCs#\n")
                Basic.write("crdaction SRSF1_Traj projection SRSF1 modes SRSF1_EigVecs \\")
                Basic.write(f"\n  beg 1 end 50 crdframes {eq_start},{eq_end} !:WAT&!@H=\n")
                for i in range(1, 11):
                    Basic.write(f"hist SRSF1:{i} bins 100 out SRSF1-Hist.agr norm name SRSF1-{i}\n")
                Basic.write("run\nclear all\n\n")
                Basic.write("##PC Visualization##\n")
                Basic.write("#new input files\n")
                Basic.write(f"readdata SRSF1_{x}_EigVecs.dat\n")
                Basic.write(f"parm SRSF1_{x}.prmtop\n")
                Basic.write("parmstrip !(!:WAT&!@H=)\n")
                Basic.write(f"parmwrite out SRSF1_{x}_bkbone.prmtop\n")
                Basic.write("#build new PC trajectories for each.  These can be viewed as NetCDF files in VMD\n")
                Basic.write(f"runanalysis modes name EigVecs trajout SRSF1_{x}_mode1.nc pcmin -100 pcmax 100 tmode 1 trajout mask !:WAT&!@H= trajoutfmt netcdf")

            with open(eq, "w") as Eq:
                Eq.write(f"parm SRSF1_{x}.prmtop\n")
                Eq.write(f"trajin SRSF1_{x}-eq.binpos\n")
                Eq.write(f"#Find hydrogen bonds#\n")
                Eq.write(f"hbond HB @N,C,O,H out SRSF1.hbond avgout SRSF1_{x}_Avg.hbond series printatomnum\n")
                Eq.write(f"#summary of fraction of time spent in each structure and matrix of structure over time for each residue\n")
                Eq.write(f"secstruct out SRSF1_{x}.StrMat sumout SRSF1_{x}.SecStr\n")
                Eq.write("#calculate pairwise distances\n")
                Eq.write(f"matrix dist @CA @CA out distances_residues_{x}.dat\nrun\n")
            with open(eps, "w") as Eps:
                Eps.write(f"parm SRSF1_{x}.prmtop\n")
                Eps.write(f"trajin SRSF1_{x}-eq.binpos\n")
                Eps.write("# Build kdist plot to help select epsilon and minpoint values for dbscan clustering\n")
                Eps.write(f"cluster C0 dbscan kdist 5 rms sieve 10 kfile SRSF1_{x}_Cluster\n")
                Eps.write("run")

            with open(cluster, "w") as clst:
                clst.write(f"parm SRSF1_{x}.prmtop\n")
                clst.write(f"trajin SRSF1_{x}-eq.binpos\n")
                clst.write("cluster dbscan minpoints 5 epsilon 1.696 sievetoframe \\\n")
                clst.write("\trms \\\n")
                clst.write("\tsieve 10 random \\\n")
                clst.write(f"\tout SRSF1_{x}.clust \\\n")
                clst.write("\tsil Sil \\\n")
                clst.write(f"\tsummary SRSF1_{x}-summary.clust \\\n")
                clst.write(f"\tinfo SRSF1_{x}-info.clust \\\n")
                clst.write(f"\tcpopvtime SRSF1_{x}-time.clust normframe \\\n")
                clst.write("\trepout clust repfmt pdb \\\n")
                clst.write(f"\tsinglerepout SRSF1_{x}-2DrepClust.nc singlerepfmt netcdf \\\n")
                clst.write("\tavgout Avg avgfmt restart\n")
                clst.write("run\n")
            with open(jobfile_path, "w") as jobfile:
                jobfile.write("#!/bin/bash\n")
                jobfile.write("#SBATCH --partition=short\n")
                jobfile.write(f"#SBATCH --job-name=trajectory_{x}\n")
                jobfile.write(f"#SBATCH --error=trajectory_{x}.err\n")
                jobfile.write(f"#SBATCH --output=trajectory_{x}.out\n")
                jobfile.write("#SBATCH --ntasks=1\n")
                jobfile.write("#SBATCH --time=12:00:00\n")
                jobfile.write("#SBATCH --mem-per-cpu=4G\n")
                jobfile.write("#SBATCH --mail-type=ALL\n")
                jobfile.write("#SBATCH --mail-user=tfarg@uab.edu\n")
                jobfile.write("module reset\n")
                jobfile.write("module load OpenMPI/4.1.5-GCC-12.3.0\n")
                #jobfile.write("module load Amber/20.11-fosscuda-2020b-AmberTools-21.3\n")
                jobfile.write("module load Amber/20.11-foss-2020b-AmberTools-21.3\n")
                jobfile.write("module load CUDA/11.6.0\n")
                jobfile.write("module load Grace\n\n")
                
                jobfile.write("rm *ensemble*pdb\n")
                jobfile.write(f"dos2unix {basic}\n")
                jobfile.write(f"cpptraj -i {basic}\n")
                jobfile.write(f"dos2unix {eq}\n")
                jobfile.write(f"cpptraj -i {eq}\n")
                jobfile.write(f"dos2unix {eps}\n")
                jobfile.write(f"cpptraj -i {eps}\n")
                jobfile.write("xmgrace Kdist.5.dat\n")
                jobfile.write(f"dos2unix {cluster}\n")
                jobfile.write(f"cpptraj -i {cluster}\n")
                jobfile.write("module purge\n")
                jobfile.write("module reset\n")
                jobfile.write("module load Anaconda3\n")
                jobfile.write("source activate Talias_modules\n")
                jobfile.write('> find-pi3_ensemble.out; for file in *-eq.pdb; do python3 find-pi4.py "${file}" -r arginine_stacks_ensemble.txt -p pi_pi_stacks3_ensemble.txt -cp cation_pi_stacks_ensemble.txt >> find-pi3_ensemble.out; done')
                
            os.system(f"dos2unix {jobfile_path}\n")
            os.system(f"sbatch {jobfile_path}\n")
            os.system("mkdir ../pdb_files")
            os.system("cp *pdb ../pdb_files")
        except Exception as e:
            print(f"An error occurred while processing folder {folder}: {e}")
        finally:
            os.chdir(original_directory)  # Change back to the original directory
            print(f"Changing directory to {original_directory}")
#change directory to the one containing the pdb files
os.chdir(f"{original_directory}/pdb_files")
#create a list of equilibrated analysis files
ensemblepdbs = []
for file in os.listdir(f"{original_directory}/pdb_files"):
    if "-eq.pdb" in file:
        ensemblepdbs.append(file)
#convert into a space delimited list
spacedelimited_list = " ".join(ensemblepdbs)
#create a text file containing the line to copy/paste into vmd
with open("vmd_input.txt","w") as file:
    file.write(f"vmd {spacedelimited_list}")