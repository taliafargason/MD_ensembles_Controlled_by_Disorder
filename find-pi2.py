import sys
import numpy as np
import itertools
import math
import argparse

# Set default thresholds and output files
vers = 0.1
dthresh = 6
angthresh = 20
arg_dist_thresh = 5

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
            if line[17:20] in ['PHE','TYR'] and line[13:16] in ['CG ','CD1','CD2','CE1','CE2','CZ ']:
                aaname ='{0}.{1}{2}'.format(line[21],line[17:20],line[22:26])
                if aaname not in aromatics:
                    aromatics[aaname] = {}
                aromatics[aaname][line[13:16].strip(' ')] = np.array([float(line[31:38]),float(line[39:46]),float(line[47:54])])
            elif line[17:20] == 'TRP' and line[13:16] in ['CD2','CE3','CZ3','CH2','CZ2','CE2']:
                aaname ='{0}.{1}{2}'.format(line[21],line[17:20],line[22:26])
                if aaname not in aromatics:
                    aromatics[aaname] = {}
                aromatics[aaname][line[13:16].strip(' ')] = np.array([float(line[31:38]),float(line[39:46]),float(line[47:54])])
            elif line[17:20] == 'ARG' and line[13:16] in ['CZ ', 'NH1', 'NH2']:
                aaname ='{0}.{1}{2}'.format(line[21],line[17:20],line[22:26])
                if aaname not in arginines:
                    arginines[aaname] = {}
                arginines[aaname][line[13:16].strip(' ')] = np.array([float(line[31:38]),float(line[39:46]),float(line[47:54])])
    
    # Filter out incomplete aromatic rings
    baddies = [i for i in aromatics if (len(aromatics[i]) != 6 and 'TRP' not in i) or (len(aromatics[i]) != 6 and 'TRP' in i)]
    for i in baddies:
        if not headless:
            print('dropped {0}: missing carbon(s)'.format(i))
        del aromatics[i]

    # Filter out incomplete arginine side chains
    baddies = [i for i in arginines if len(arginines[i]) != 3]
    for i in baddies:
        if not headless:
            print('dropped {0}: missing atom(s)'.format(i))
        del arginines[i]

    return aromatics, arginines

def find_ring_center(ring):
    '''find the center of an aromatic ring input dic entry in format aaname:{CG:array(coords),CD1:array(coords), and etc...}'''
    if 'CG' in ring:  # PHE and TYR
        center = (ring['CG'] + ring['CD1'] + ring['CD2'] + ring['CE1'] + ring['CE2'] + ring['CZ']) / 6
        vec1 = ring['CG'] - center
        vec2 = ring['CD1'] - center
    elif 'CD2' in ring:  # TRP
        center = (ring['CD2'] + ring['CE3'] + ring['CZ3'] + ring['CH2'] + ring['CZ2'] + ring['CE2']) / 6
        vec1 = ring['CD2'] - center
        vec2 = ring['CE3'] - center
    else:
        raise ValueError("Unknown ring structure")
    face_norm = np.cross(vec1, vec2)
    return center, face_norm

def find_arg_center(arg):
    '''find the center of the guanidinium group in an arginine sidechain'''
    required_atoms = ['CZ', 'NH1', 'NH2']
    for atom in required_atoms:
        if atom not in arg:
            raise KeyError(f'Missing {atom} in arginine: {arg}')
    center = (arg['CZ'] + arg['NH1'] + arg['NH2']) / 3
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

def point_plane_distance(point, plane_point, plane_normal):
    """ Returns the distance from a point to a plane defined by a point and a normal vector """
    return np.dot(plane_normal, point - plane_point) / np.linalg.norm(plane_normal)

def is_pointing_to_plane(point, plane_point, plane_normal):
    """ Checks if the vector from the point to the plane is in the same direction as the plane normal """
    vector_to_plane = plane_point - point
    return np.dot(vector_to_plane, plane_normal) > 0

try:
    aromatics, arginines = read_pdb_find_aromatics_and_args(testfile)

    # get the centers and face normal vectors of all aromatics and arginines
    aromatic_centers = {} # center {coords:aaname}
    aromatic_face_norms = {}
    for i in aromatics:
        center, face_norm = find_ring_center(aromatics[i])
        aromatic_centers[i] = center
        aromatic_face_norms[i] = face_norm

    arginine_centers = {} # center {coords:aaname}
    arginine_face_norms = {}
    for i in arginines:
        center, face_norm = find_arg_center(arginines[i])
        arginine_centers[i] = center
        arginine_face_norms[i] = face_norm

    # find the different possible combinations of aromatics
    aromatic_combinations = list(itertools.combinations(aromatic_centers, 2))

    # calculate the distances between all of the aromatic centers
    aromatic_distances = {}
    for i in aromatic_combinations:
        combo = '{0}-{1}'.format(i[0],i[1])
        aromatic_distances[combo] = np.linalg.norm(aromatic_centers[i[0]] - aromatic_centers[i[1]])

    # find the different possible combinations of arginines
    arginine_combinations = list(itertools.combinations(arginine_centers, 2))

    # calculate the distances between all of the arginine centers
    arginine_distances = {}
    for i in arginine_combinations:
        combo = '{0}-{1}'.format(i[0],i[1])
        arginine_distances[combo] = np.linalg.norm(arginine_centers[i[0]] - arginine_centers[i[1]])

    # find the different possible combinations of arginines and aromatics
    arg_aro_combinations = list(itertools.product(arginine_centers, aromatic_centers))

    # calculate the distances between all of the arginine and aromatic centers
    arg_aro_distances = {}
    for i in arg_aro_combinations:
        combo = '{0}-{1}'.format(i[0],i[1])
        arg_aro_distances[combo] = np.linalg.norm(arginine_centers[i[0]] - aromatic_centers[i[1]])

    ## Analysis for pi-pi stacking
    # find the aromatic pairs within distance and angle threshold
    threshed = []
    threshed_pairs = []
    threshed_aas = set()
    for i in aromatic_distances:
        if aromatic_distances[i] < dthresh:
            threshed.append(i)
            threshed_pairs.append(i.split('-'))
            threshed_aas.add(i.split('-')[0])
            threshed_aas.add(i.split('-')[1])

    # find the sets of three chains with each residue in a pair of the others
    chains = []
    chains.append(list(threshed_aas))
    for i in threshed_aas:
        for j in threshed_aas:
            if i != j and set([i,j]) not in chains and set([j,i]) not in chains:
                if [i,j] in threshed_pairs or [j,i] in threshed_pairs:
                    chains.append([i,j])

    with open(output_pi_stacking, 'w') as f:
        f.write("pi-pi stacking interactions:\n")
        f.write("\n".join([str(chain) for chain in chains]))

    ## Analysis for cation-pi stacking
    # find the arg-aro distances and angles less than threshold
    cation_pi_pairs = []
    for i in arg_aro_distances:
        if arg_aro_distances[i] < dthresh:
            arg, aro = i.split('-')
            angle = angle_between(arginine_face_norms[arg], aromatic_face_norms[aro])
            if angle < angthresh:
                # Ensure the norm of the arginine plane intersects with the plane of the aromatic residue
                if is_pointing_to_plane(arginine_centers[arg], aromatic_centers[aro], aromatic_face_norms[aro]):
                    cation_pi_pairs.append((arg, aro, arg_aro_distances[i], angle))

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
                arg_arg_pairs.append((arg1, arg2, arginine_distances[i], angle))

    with open(output_arg_stacking, 'w') as f:
        f.write("arginine-arginine stacking interactions:\n")
        f.write("\n".join([f"{pair[0]}-{pair[1]}: Distance={pair[2]:.2f} Å, Angle={pair[3]:.2f} degrees" for pair in arg_arg_pairs]))

except KeyError as e:
    print(f"Error: Missing atom in input file: {e}")
    sys.exit(1)
except Exception as e:
    print(f"An error occurred: {e}")
    sys.exit(1)
