
import numpy as np
from scipy.spatial import distance_matrix

def coords_fetch(pdbfile,chain):
    #pdbfile,chain = pdbfilechain.split('_')
    f = open(pdbfile,'r')
    rows = f.readlines()
    f.close()
    coords_CA = []
    coords_N = []
    coords_C = []
    coords_O = []
    coords_CB = []
    coords_CG1 = []
    coords_CG2 = []
    for row in rows:
        row = row.strip()
        if row[:4] == 'ATOM' and row[21] in chain:
            if row[13:16] =='CA ':
                coords_CA.append([ float(row[27:38]),float(row[38:46]),float(row[46:54])])
            if row[13:16] =='N  ':
                coords_N.append([ float(row[27:38]),float(row[38:46]),float(row[46:54])])
            if row[13:16] =='O ':
                coords_O.append([ float(row[27:38]),float(row[38:46]),float(row[46:54])])
            if row[13:16] =='C  ':
                coords_C.append([ float(row[27:38]),float(row[38:46]),float(row[46:54])])
            if row[13:16] =='CB ':
                coords_CB.append([ float(row[27:38]),float(row[38:46]),float(row[46:54])])
            if row[13:16] =='CG1':
                coords_CG1.append([ float(row[27:38]),float(row[38:46]),float(row[46:54])])
            if row[13:16] =='CG2':
                coords_CG2.append([ float(row[27:38]),float(row[38:46]),float(row[46:54])])
    atom_coords = dict()
    atom_coords['CA'] = coords_CA
    atom_coords['C'] = coords_C
    atom_coords['N'] = coords_N
    atom_coords['O'] = coords_N
    atom_coords['CB'] = coords_CB
    atom_coords['CG1'] = coords_CG1
    atom_coords['CG2'] = coords_CG2

    dismap = []
    for atom1 in ['CA','N','C','O','CB','CG1','CG2']:
        for atom2 in ['CA','N','C','O','CB','CG1','CG2']:
            dismap.append( distance_matrix( atom_coords[atom1],atom_coords[atom2] )  )
    return np.stack(dismap,axis = -1)
