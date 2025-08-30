import numpy as np
import itertools
import math
import argparse

# Set default thresholds and output files
#vers = 0.1
dthresh = 4.5
angthresh = 20
arg_dist_thresh = 4.5

# Argument parser for command line flags
parser = argparse.ArgumentParser(description="Analyze stacking interactions in PDB files.")
parser.add_argument("testfile", help="Input PDB file")
parser.add_argument("-r", "--arginine", help="Output file for arginine stacking interactions", default="arginine_stacks.txt")
parser.add_argument("-p", "--pi_pi", help="Output file for pi-pi stacking interactions", default="pi_pi_stacks.txt")
parser.add_argument("-cp", "--cation_pi", help="Output file for cation-pi stacking interactions", default="cation_pi_stacks.txt")
parser.add_argument("--headless", help="Run in headless mode", action="store_true")

args = parser.parse_args()

testfile = args.testfile
output_arg_stacking = args.arginine
output_pi_stacking = args.pi_pi
output_cation_pi_stacking = args.cation_pi
headless = args.headless

def read_pdb_find_aromatics_and_args(file):
    '''find all the aromatics and arginines, get the coordinates of relevant atoms'''
    aromatics = {}      # dic {name:{CG:array(coords),CD1:array(coords), and etc...}} 
    arginines = {}      # dic {name:{CZ:array(coords),NH1:array(coords),NH2:array(coords)}}
    data = open(file,'r').readlines()
    for line in data:
        if line[0:4] == 'ATOM':
            if line[17:20] in ['PHE','TYR'] and line[13:16] in ['CD1','CD2','CE1','CE2','CZ']:
                aaname ='{0}.{1}{2}'.format(line[21],line[17:20],line[22:26])
                if aaname not in aromatics:
                    aromatics[aaname] = {}
                aromatics[aaname][line[13:16].strip(' ')] = np.array([float(line[31:38]),float(line[39:46]),float(line[47:54])])
            elif line[17:20] == 'TRP' and line[13:16] in ['CD2','CE3','CZ3','CH2','CZ2','CE2']:
                aaname ='{0}.{1}{2}'.format(line[21],line[17:20],line[22:26])
                if aaname not in aromatics:
                    aromatics[aaname] = {}
                aromatics[aaname][line[13:16].strip(' ')] = np.array([float(line[31:38]),float(line[39:46]),float(line[47:54])])
            elif line[17:20] in ['HIS','HIE'] and line[13:16] in [ 'ND1', 'CD2', 'CE1', 'NE2']:
                aaname ='{0}.{1}{2}'.format(line[21],line[17:20],line[22:26])
                if aaname not in aromatics:
                    aromatics[aaname] = {}
                aromatics[aaname][line[13:16].strip(' ')] = np.array([float(line[31:38]),float(line[39:46]),float(line[47:54])])
            elif line[17:20] == 'ARG' and line[13:16] in ['CZ', 'NH1', 'NH2', 'NE']:
                aaname ='{0}.{1}{2}'.format(line[21],line[17:20],line[22:26])
                if aaname not in arginines:
                    arginines[aaname] = {}
                arginines[aaname][line[13:16].strip(' ')] = np.array([float(line[31:38]),float(line[39:46]),float(line[47:54])])
    
    return aromatics, arginines

def find_ring_center(ring):
    '''find the center of an aromatic ring input dic entry in format aaname:{CG:array(coords),CD1:array(coords), and etc...}'''

    # Define atom sets for different aromatic residues
    atom_sets = {
        'PHE': [ 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'TYR': [ 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'TRP': ['CD2', 'CE3', 'CZ3', 'CH2', 'CZ2', 'CE2'],
        'HIS': [ 'ND1', 'CD2', 'CE1', 'NE2'],
        'HIE': [ 'ND1', 'CD2', 'CE1', 'NE2']
    }


    # Determine which atom set to use based on the atoms present in the ring
    residue_type = None
    for res, atoms in atom_sets.items():
        if all(atom in ring for atom in atoms):
            residue_type = res
            required_atoms = atoms
            break

    if residue_type is None:
        print(f"residue type is none {[atom for atom in ring]} ")
        return None, None  # Skip this residue if it does not contain necessary atoms

    # Calculate the center and face normal of the ring
    center = np.mean([ring[atom] for atom in required_atoms], axis=0)
    vec1 = ring[required_atoms[0]] - center
    vec2 = ring[required_atoms[1]] - center
    face_norm = np.cross(vec1, vec2)
    return center, face_norm


def find_arg_center(arg):
    '''find the center of the guanidinium group in an arginine sidechain'''
    required_atoms = ['CZ', 'NH1', 'NH2','NE']
    for atom in required_atoms:
        if atom not in arg:
            return None, None  # Skip this arginine if it does not contain necessary atoms
    center = (arg['CZ'] + arg['NH1'] + arg['NH2'] + arg['NE']) / 4
    vec1 = arg['NH1'] - center
    vec2 = arg['NH2'] - center
    face_norm = np.cross(vec1, vec2)
    return center, face_norm

def unit_vector(vector):
    """ Returns the unit vector of the vector. """
    if np.all(vector == 0):
        vector = np.array([0.001,0.001,0.001])
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Return the angle in degrees between two vectors"""
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = math.degrees(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))
    if angle == 0:
        angle = 0.000001
    if angle < 90:
        return angle
    else:
        return 180.0 - angle

#def point_plane_distance(point, plane_point, plane_normal):
 #   """ Returns the distance from a point to a plane defined by a point and a normal vector """
  #  return np.dot(plane_normal, point - plane_point) / np.linalg.norm(plane_normal)

def is_pointing_to_plane(point, plane_point, plane_normal):
    """ Checks if the vector from the point to the plane is in the same direction as the plane normal """
    vector_to_plane = plane_point - point
    return np.dot(vector_to_plane, plane_normal) > 0

try:
    aromatics, arginines = read_pdb_find_aromatics_and_args(testfile)

    # Print the input data structures
    print("Aromatics:", aromatics)
    print("Arginines:", arginines)

    # get the centers and face normal vectors of all aromatics and arginines
    aromatic_centers = {} # center {coords:aaname}
    aromatic_face_norms = {}
    for i in aromatics:
        center, face_norm = find_ring_center(aromatics[i])
        if center is not None and face_norm is not None:
            aromatic_centers[i] = center
            aromatic_face_norms[i] = face_norm
        else:
            print(f"Warning: Skipping aromatic residue {i} due to missing atoms")

    arginine_centers = {} # center {coords:aaname}
    arginine_face_norms = {}
    for i in arginines:
        center, face_norm = find_arg_center(arginines[i])
        if center is not None and face_norm is not None:
            arginine_centers[i] = center
            arginine_face_norms[i] = face_norm
        else:
            print(f"Warning: Skipping arginine residue {i} due to missing atoms")

    # Print the calculated centers and face normals
    print("Aromatic Centers:", aromatic_centers)
    print("Aromatic Face Norms:", aromatic_face_norms)
    print("Arginine Centers:", arginine_centers)
    print("Arginine Face Norms:", arginine_face_norms)

    # find the different possible combinations of aromatics and calculate their distances
    aromatic_combinations = list(itertools.product(aromatic_centers, repeat=2))

    aromatic_distances = {}
    for i in aromatic_combinations:
        if i[0] != i[1]:
            combo = '{0}-{1}'.format(i[0],i[1])
            aromatic_distances[combo] = np.linalg.norm(aromatic_centers[i[0]] - aromatic_centers[i[1]])

    # find the different possible combinations of arginines and calculate their distances
    arginine_combinations = list(itertools.product(arginine_centers, repeat=2))

    arginine_distances = {}
    for i in arginine_combinations:
        if i[0] != i[1]:
            combo = '{0}-{1}'.format(i[0],i[1])
            arginine_distances[combo] = np.linalg.norm(arginine_centers[i[0]] - arginine_centers[i[1]])

    # find the different possible combinations of arginines and aromatics and calculate their distances
    arg_aro_combinations = list(itertools.product(arginine_centers, aromatic_centers))

    arg_aro_distances = {}
    for i in arg_aro_combinations:
        combo = '{0}-{1}'.format(i[0],i[1])
        arg_aro_distances[combo] = np.linalg.norm(arginine_centers[i[0]] - aromatic_centers[i[1]])

    # Print the distance calculations
    print("Aromatic Distances:", aromatic_distances)
    print("Arginine Distances:", arginine_distances)
    print("Arginine-Aromatic Distances:", arg_aro_distances)

    ## Analysis for pi-pi stacking
    # find the distances less than threshold and find chains of multiple within threshold pairs
    threshed = []
    threshed_pairs = []
    #threshed_aas = set()
    for i in aromatic_distances:
        if aromatic_distances[i] < dthresh:
            
            aro1,aro2 = i.split('-')
            angle = angle_between(aromatic_face_norms[aro1], aromatic_face_norms[aro2]) 
            if angle < angthresh:
                # Ensure the norm of the arginine plane intersects with the plane of the aromatic residue
                if is_pointing_to_plane(aromatic_centers[aro1], aromatic_centers[aro2], aromatic_face_norms[aro2]):
                    threshed.append(i)
                    forward = (aro1,aro2,aromatic_distances[i],angle)
                    reverse = (aro2,aro1,aromatic_distances[i],angle)
                    if forward not in threshed_pairs and reverse not in threshed_pairs:
                        threshed_pairs.append((aro1,aro2,aromatic_distances[i],angle))
                    #threshed_pairs2.append(i.split('-'))
                        #threshed_aas.add(i.split('-')[0])
                        #threshed_aas.add(i.split('-')[1])
    with open(output_pi_stacking,'w') as f:
        f.write("pi-pi stacking interactions:\n")
        f.write("\n".join([f"{pair[0]}-{pair[1]}: Distance={pair[2]:.2f} Å, Angle={pair[3]:.2f} degrees" for pair in threshed_pairs]))

   # find the sets of three chains with each residue in a pair of the others
    #chains = []
    #chains.append(list(threshed_aas))
    #for i in threshed_aas:
     #   for j in threshed_aas:
      #      if i != j and set([i,j]) not in chains and set([j,i]) not in chains:
       #         if [i,j] in threshed_pairs2 or [j,i] in threshed_pairs2:
        #            chains.append([i,j])

    # Print the identified pi-pi interactions
   # print("Pi-Pi Stacking Interactions:", chains)

   # with open(output_pi_stacking, 'w') as f:
    #    f.write("pi-pi stacking interactions:\n")
     #   f.write("\n".join([str(chain) for chain in chains]))
  
    ## Analysis for cation-pi stacking
    # find the arg-aro distances and angles lower than threshold
    # Analysis for cation-pi stacking
    cation_pi_pairs = []
    for i in arg_aro_distances:
        if arg_aro_distances[i] < dthresh:
            arg, aro = i.split('-')
            angle = angle_between(arginine_face_norms[arg], aromatic_face_norms[aro])
            if angle < angthresh:
                # Try using NE, CZ, NH1, and NH2 as the point
                points_to_try = [arginines[arg].get(atom) for atom in ['NE', 'CZ', 'NH1', 'NH2'] if atom in arginines[arg]]
                for point in points_to_try:
                    if is_pointing_to_plane(point, aromatic_centers[aro], aromatic_face_norms[aro]):
                        cation_pi_pairs.append((arg, aro, arg_aro_distances[i], angle))
                        break  # If any point gives a positive dot product, append the pair and stop checking other points
    
    # Print the identified cation-pi interactions
    print("Cation-Pi Stacking Interactions:", cation_pi_pairs)
    
    with open(output_cation_pi_stacking, 'w') as f:
        f.write("cation-pi stacking interactions:\n")
        f.write("\n".join([f"{pair[0]}-{pair[1]}: Distance={pair[2]:.2f} Å, Angle={pair[3]:.2f} degrees" for pair in cation_pi_pairs]))

    ## Analysis for arginine-arginine stacking
    # find the arg-arg distances and angles less than threshold
    arg_arg_pairs = []
    for i in arginine_distances:
        if arginine_distances[i] < arg_dist_thresh:
            arg1, arg2 = i.split('-')
            angle = angle_between(arginine_face_norms[arg1], arginine_face_norms[arg2])
            if angle < angthresh:
                forward = (arg1, arg2, arginine_distances[i],angle)
                reverse = (arg2, arg1, arginine_distances[i],angle)
                if forward not in arg_arg_pairs and reverse not in arg_arg_pairs:
                    arg_arg_pairs.append((arg1, arg2, arginine_distances[i], angle))

    # Print the identified arginine-arginine interactions
    print("Arginine-Arginine Stacking Interactions:", arg_arg_pairs)

    with open(output_arg_stacking, 'w') as f:
        f.write("arginine-arginine stacking interactions:\n")
        f.write("\n".join([f"{pair[0]}-{pair[1]}: Distance={pair[2]:.2f} Å, Angle={pair[3]:.2f} degrees" for pair in arg_arg_pairs]))

except KeyError as e:
    if 'aromatics' in str(e):
        residue, atom = e.args[0].split("'")[1].split('.')
        print(f"Error: Missing atom '{atom}' in aromatic residue '{residue}'")
    elif 'arginines' in str(e):
        residue, atom = e.args[0].split("'")[1].split('.')
        print(f"Error: Missing atom '{atom}' in arginine residue '{residue}'")
    else:
        print(f"Error: Missing atom in residue input file: {e}")

except Exception as e:
    print(f"An error occurred: {e}")
    #sys.exit(1)