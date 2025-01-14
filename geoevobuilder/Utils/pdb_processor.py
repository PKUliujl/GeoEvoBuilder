# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 16:17:36 2020

@author: Administrator
"""
import pandas as pd
import Bio
import Bio.PDB
from Bio.PDB.DSSP import DSSP
import numpy as np
from scipy.spatial import distance_matrix
from math import sin,cos

import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore')


standard_aa_names = [
    "ALA",
    "CYS",
    "ASP",
    "GLU",
    "PHE",
    "GLY",
    "HIS",
    "ILE",
    "LYS",
    "LEU",
    "MET",
    "ASN",
    "PRO",
    "GLN",
    "ARG",
    "SER",
    "THR",
    "VAL",
    "TRP",
    "TYR",
    ]

## 20 standard anino acids

AA = {'A':0,
      'C':1,
      'D':2,
      'E':3,
      'F':4,
      'G':5,
      'H':6,
      'I':7,
      'K':8,
      'L':9,
      'M':10,
      'N':11,
      'P':12,
      'Q':13,
      'R':14,
      'S':15,
      'T':16,
      'V':17,
      'W':18,
      'Y':19}

AA_reverse = dict( [val,key ] for key,val in AA.items() )

Rotamer_number = {0:0,
                1:1,
                2:2,
                3:3,
                4:2,
                5:0,
                6:2,
                7:2,
                8:4,
                9:2,
                10:3,
                11:2,
                12:1,
                13:3,
                14:4,
                15:1,
                16:1,
                17:1,
                18:2,
                19:2
                }


def eight2three(ss):
    if ss=='H':
        return [1,0,0,0,0,0,0,0]
    elif ss=='G':
        return [0,1,0,0,0,0,0,0]
    elif ss=='I':
        return [0,0,1,0,0,0,0,0]
    elif ss=='E':
        return [0,0,0,1,0,0,0,0]
    elif ss=='B':
        return [0,0,0,0,1,0,0,0]
    elif ss=='T':
        return [0,0,0,0,0,1,0,0]
    elif ss=='S':
        return [0,0,0,0,0,0,1,0]
    else:
        return [0,0,0,0,0,0,0,1]

def fetchseq( files, chainID):
    structure = Bio.PDB.PDBParser(PERMISSIVE=1).get_structure(files, files)
    assert len(structure)==1,'Only 1 model is supported for this version'
    for model in structure:
        for chain in model:
            if chain.id == chainID:
                poly=Bio.PDB.Polypeptide.Polypeptide(chain)
                sequence_list=poly.get_sequence()
                print('Original sequence: ', sequence_list)
                return sequence_list, [AA[i] for i in sequence_list]





def loadpdbs(files, chainID, interval_b, interval_e):
    structure = Bio.PDB.PDBParser(PERMISSIVE=1).get_structure('tmp', files)

    #delete irregular residue atoms
    for model in structure:
        for chain in model:
         if chain.id == chainID:
            delete_REScluster=[]
            for residue in chain:
                #print(residue.id)
                if residue.get_resname() not in standard_aa_names:
                    delete_REScluster.append(residue.id)
            if delete_REScluster!=[]:
                for delete_res in delete_REScluster:
                    chain.detach_child(delete_res)

    #neighbor search
    atoms  = Bio.PDB.Selection.unfold_entities(structure, 'A')
    ns = Bio.PDB.NeighborSearch(atoms)

    for model in structure:
            #print(residue.id)
            if residue.get_resname() not in standard_aa_names:
                delete_REScluster.append(residue.id)
            if delete_REScluster!=[]:
                #print(delete_REScluster)
                for delete_res in delete_REScluster:
                    chain.detach_child(delete_res)

    #neighbor search
    atoms  = Bio.PDB.Selection.unfold_entities(structure, 'A')
    ns = Bio.PDB.NeighborSearch(atoms)

    for model in structure:
     dssp = DSSP(model, files, dssp='dssp' )

     for chain in model:
      if chain.id == chainID:  
        #print('\n-------\n%s %s:\n'%(structure.id,chain.id))
        poly=Bio.PDB.Polypeptide.Polypeptide(chain)
        phi_psi_list = poly.get_phi_psi_list()
        #print("phi_psi_list: ",phi_psi_list)  #dtype: list
        sequence_list=poly.get_sequence()
        #print("sequence list: ",sequence_list)  #dtype
        residue_list=[residue.id for residue in chain]
        #print("residue list: ",residue_list) #dtype
        matrix = [[i,chain.id,j,l,dssp[(chain.id,chain[int(l)].id)][2],n,p] for i,(j,(k,l,m),(n,p)) in enumerate(zip(sequence_list,residue_list,phi_psi_list))  ]
        #print(matrix)
        matrix = pd.DataFrame(matrix)
        node_feature = []
        node_category = []
        edge_attributes = []
        edge_index =[]
        position_encoding=[]
        absolutepostion = int(matrix.iloc[0,3])  ##first node
        for i in range( int(interval_b), int(interval_e)): # len(matrix)
            each_node_feature = []
            ####first four columns are sin *cos (phi psi)
            if pd.isna(matrix.iloc[i,-2]):
                each_node_feature.extend( [ 0,0,sin(matrix.iloc[i,-1]),cos(matrix.iloc[i,-1]),   0,0,0,0,0,0,sin(matrix.iloc[i,-1]*2),sin(matrix.iloc[i,-1]*3),sin(matrix.iloc[i,-1]*4),cos(matrix.iloc[i,-1]*2),cos(matrix.iloc[i,-1]*3),cos(matrix.iloc[i,-1]*4) ]     )
            elif pd.isna(matrix.iloc[i,-1]):
                each_node_feature.extend( [ sin(matrix.iloc[i,-2]),cos(matrix.iloc[i,-2]),0,0,   sin(matrix.iloc[i,-2]*2),sin(matrix.iloc[i,-2]*3),sin(matrix.iloc[i,-2]*4),cos(matrix.iloc[i,-2]*2),cos(matrix.iloc[i,-2]*3),cos(matrix.iloc[i,-2]*4),0,0,0,0,0,0 ]      )
            else:
                each_node_feature.extend( [ sin(matrix.iloc[i,-2]),cos(matrix.iloc[i,-2]),sin(matrix.iloc[i,-1]),cos(matrix.iloc[i,-1]),       sin(matrix.iloc[i,-2]*2),sin(matrix.iloc[i,-2]*3),sin(matrix.iloc[i,-2]*4),cos(matrix.iloc[i,-2]*2),cos(matrix.iloc[i,-2]*3),cos(matrix.iloc[i,-2]*4),       sin(matrix.iloc[i,-1]*2),sin(matrix.iloc[i,-1]*3),sin(matrix.iloc[i,-1]*4),cos(matrix.iloc[i,-1]*2),cos(matrix.iloc[i,-1]*3),cos(matrix.iloc[i,-1]*4)  ]     )
            ###The next three columns are 3 states of DSSP
            each_node_feature.extend( eight2three( matrix.iloc[i,4] )  )
            node_category.append(AA[matrix.iloc[i,2]])
            ##local coordinatate 
            CA_coord = chain[int(matrix .iloc[i,3])]['CA'].coord
            C_coord = chain[int(matrix .iloc[i,3])]['C'].coord
            N_coord = chain[int(matrix .iloc[i,3])]['N'].coord
            O_coord = chain[int(matrix .iloc[i,3])]['O'].coord
            CA_C = C_coord - CA_coord    ## norm   np.linalg.norm( CA_C )
            CA_N = N_coord - CA_coord 
            CA_O = O_coord - CA_coord
            orthognal = np.cross(CA_C, CA_N)
            #each_node_feature.extend(  CA_C.tolist() )
            #each_node_feature.extend( CA_N.tolist() )
            #each_node_feature.extend( orthognal.tolist())
            #node_feature.append(each_node_feature )
            NA = ns.search(chain[int(matrix .iloc[i,3])]['CA'].coord, 12)  ##Neighbor Atoms
            #print(i,matrix.iloc[i,3],[atom.get_parent().id[1] for atom in NA if atom.id=='CA'])
            All_CA = [atom for atom in NA if atom.id=='CA']   # atom.get_parent().id[1]
            if len(All_CA) ==1:
                each_node_feature.extend([len(All_CA)-1, 10000])
            else:
                each_node_feature.extend([len(All_CA)-1, 1/(len(All_CA)-1)])
            node_feature.append(each_node_feature )
            each_residue_position=[]
            for ik in range(0,13):
                each_residue_position.extend([sin( (matrix.iloc[i,3]-absolutepostion)/1000**(ik/13) ), cos( (matrix.iloc[i,3]-absolutepostion)/1000**(ik/13) )])
            position_encoding.append(each_residue_position)
            for CA in All_CA: 
                if (CA.coord != chain[int(matrix .iloc[i,3])]['CA'].coord).any() and CA.get_parent().get_parent().id == chain.id :
                    each_edge_attributes = []
                    each_edge_index = []
                    #print(chain[int(matrix .iloc[i,3])].id[1], CA.get_parent().id[1] )
                    CA_CA_dis = np.linalg.norm( CA.coord - CA_coord )  #distance between two CA
                    CA_CA_orientation = [np.dot( (CA.coord - CA_coord), CA_C )/np.linalg.norm( CA_C ), np.dot( (CA.coord - CA_coord), CA_N )/np.linalg.norm( CA_N ), np.dot( (CA.coord - CA_coord), orthognal )/np.linalg.norm( orthognal )]
                    ## CA -> C
                    CA_C_dis = np.linalg.norm( CA.get_parent()['C'].coord - CA_coord )
                    CA_C_orientation = [np.dot( CA.get_parent()['C'].coord - CA_coord , CA_C )/np.linalg.norm( CA_C ), np.dot( CA.get_parent()['C'].coord - CA_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['C'].coord - CA_coord , orthognal )/np.linalg.norm( orthognal )]
                    ## CA -> N
                    CA_N_dis = np.linalg.norm( CA.get_parent()['N'].coord - CA_coord )
                    CA_N_orientation = [np.dot( CA.get_parent()['N'].coord - CA_coord , CA_C )/np.linalg.norm( CA_C ), np.dot( CA.get_parent()['N'].coord - CA_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['N'].coord - CA_coord , orthognal )/np.linalg.norm( orthognal )]
                    ## CA_O 
                    CA_O_dis = np.linalg.norm( CA.get_parent()['O'].coord - CA_coord )
                    CA_O_orientation = [np.dot( CA.get_parent()['O'].coord - CA_coord , CA_C )/np.linalg.norm( CA_C ), np.dot( CA.get_parent()['O'].coord - CA_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['O'].coord - CA_coord , orthognal )/np.linalg.norm( orthognal )]
                    ##CA_CB
                    if CA.get_parent().resname !='GLY':
                        CA_CB_dis = np.linalg.norm( CA.get_parent()['CB'].coord - CA_coord )
                        CA_CB_orientation = [np.dot( CA.get_parent()['CB'].coord - CA_coord , CA_C )/np.linalg.norm( CA_C ), np.dot( CA.get_parent()['CB'].coord - CA_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['CB'].coord - CA_coord , orthognal )/np.linalg.norm( orthognal )]
                    else:
                        CA_CB_dis = CA_CA_dis
                        CA_CB_orientation = CA_CA_orientation
                    CA_CG1_dis = np.linalg.norm( CA.get_parent()['CG1'].coord - CA_coord )
                    CA_CG1_orientation = [np.dot( CA.get_parent()['CG1'].coord - CA_coord , CA_C )/np.linalg.norm( CA_C ), np.dot( CA.get_parent()['CG1'].coord - CA_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['CG1'].coord - CA_coord , orthognal )/np.linalg.norm( orthognal )]
                    CA_CG2_dis = np.linalg.norm( CA.get_parent()['CG2'].coord - CA_coord )
                    CA_CG2_orientation = [np.dot( CA.get_parent()['CG2'].coord - CA_coord , CA_C )/np.linalg.norm( CA_C ), np.dot( CA.get_parent()['CG2'].coord - CA_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['CG2'].coord - CA_coord , orthognal )/np.linalg.norm( orthognal )]
                    
                    ## O -> N, where O is from the central residue and N is from the neighbor residue
                    O_N_dis = np.linalg.norm( CA.get_parent()['N'].coord - O_coord )
                    O_N_orientation = [np.dot( CA.get_parent()['N'].coord - O_coord , CA_C )/np.linalg.norm( CA_C ), np.dot( CA.get_parent()['N'].coord - O_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['N'].coord - O_coord , orthognal )/np.linalg.norm( orthognal )]
                    ## O ->C
                    O_C_dis = np.linalg.norm( CA.get_parent()['C'].coord - O_coord )
                    O_C_orientation = [np.dot( CA.get_parent()['C'].coord - O_coord , CA_C )/np.linalg.norm( CA_C ), np.dot( CA.get_parent()['C'].coord - O_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['C'].coord - O_coord , orthognal )/np.linalg.norm( orthognal )]
                    ## O ->O
                    O_O_dis = np.linalg.norm( CA.get_parent()['O'].coord - O_coord )
                    O_O_orientation = [np.dot( CA.get_parent()['O'].coord - O_coord , CA_C )/np.linalg.norm( CA_C ), np.dot( CA.get_parent()['O'].coord - O_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['O'].coord - O_coord , orthognal )/np.linalg.norm( orthognal )]
                    ## O ->CA
                    O_CA_dis = np.linalg.norm( CA.coord - O_coord )
                    O_CA_orientation = [np.dot( CA.coord - O_coord , CA_C )/np.linalg.norm( CA_C ), np.dot( CA.coord - O_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.coord - O_coord , orthognal )/np.linalg.norm( orthognal )]
                    ## o ->CB
                    if CA.get_parent().resname !='GLY':
                        O_CB_dis = np.linalg.norm( CA.get_parent()['CB'].coord - O_coord )
                        O_CB_orientation = [np.dot( CA.get_parent()['CB'].coord - O_coord , CA_C )/np.linalg.norm( CA_C ), np.dot( CA.get_parent()['CB'].coord - O_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['CB'].coord - O_coord , orthognal )/np.linalg.norm( orthognal )]
                    else:
                        O_CB_dis = O_CA_dis
                        O_CB_orientation = O_CA_orientation
                    O_CG1_dis = np.linalg.norm( CA.get_parent()['CG1'].coord - O_coord )
                    O_CG1_orientation = [np.dot( CA.get_parent()['CG1'].coord - O_coord , CA_C )/np.linalg.norm( CA_C ), np.dot( CA.get_parent()['CG1'].coord - O_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['CG1'].coord - O_coord , orthognal )/np.linalg.norm( orthognal )]
                    O_CG2_dis = np.linalg.norm( CA.get_parent()['CG2'].coord - O_coord )
                    O_CG2_orientation = [np.dot( CA.get_parent()['CG2'].coord - O_coord , CA_C )/np.linalg.norm( CA_C ), np.dot( CA.get_parent()['CG2'].coord - O_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['CG2'].coord - O_coord , orthognal )/np.linalg.norm( orthognal )]
                    
                    ## N -> O, where N is from the central residue and O is from the neighbor residue
                    N_O_dis = np.linalg.norm( CA.get_parent()['O'].coord - N_coord )
                    N_O_orientation = [np.dot( CA.get_parent()['O'].coord - N_coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['O'].coord - N_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['O'].coord - N_coord , orthognal )/np.linalg.norm( orthognal )]
                    ## N -> C
                    N_C_dis = np.linalg.norm( CA.get_parent()['C'].coord - N_coord )
                    N_C_orientation = [np.dot( CA.get_parent()['C'].coord - N_coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['C'].coord - N_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['C'].coord - N_coord , orthognal )/np.linalg.norm( orthognal )]
                    ## N -> N
                    N_N_dis = np.linalg.norm( CA.get_parent()['N'].coord - N_coord )
                    N_N_orientation = [np.dot( CA.get_parent()['N'].coord - N_coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['N'].coord - N_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['N'].coord - N_coord , orthognal )/np.linalg.norm( orthognal )]
                    ## N -> CA
                    N_CA_dis = np.linalg.norm( CA.coord - N_coord )
                    N_CA_orientation = [np.dot( CA.coord - N_coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.coord - N_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.coord - N_coord , orthognal )/np.linalg.norm( orthognal )]
                    ## N ->CB
                    if CA.get_parent().resname !='GLY':
                        N_CB_dis = np.linalg.norm( CA.get_parent()['CB'].coord - N_coord )
                        N_CB_orientation= [np.dot( CA.get_parent()['CB'].coord - N_coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['CB'].coord - N_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['CB'].coord - N_coord , orthognal )/np.linalg.norm( orthognal )]
                    else:
                        N_CB_dis = N_CA_dis
                        N_CB_orientation = N_CA_orientation
                    N_CG1_dis = np.linalg.norm( CA.get_parent()['CG1'].coord - N_coord )
                    N_CG1_orientation= [np.dot( CA.get_parent()['CG1'].coord - N_coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['CG1'].coord - N_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['CG1'].coord - N_coord , orthognal )/np.linalg.norm( orthognal )]
                    N_CG2_dis = np.linalg.norm( CA.get_parent()['CG2'].coord - N_coord )
                    N_CG2_orientation= [np.dot( CA.get_parent()['CG2'].coord - N_coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['CG2'].coord - N_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['CG2'].coord - N_coord , orthognal )/np.linalg.norm( orthognal )]    
                    
                    ## C -> C, 
                    C_C_dis = np.linalg.norm( CA.get_parent()['C'].coord - C_coord )
                    C_C_orientation = [np.dot( CA.get_parent()['C'].coord - C_coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['C'].coord - C_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['C'].coord - C_coord , orthognal )/np.linalg.norm( orthognal )]
                    ## C -> N
                    C_N_dis = np.linalg.norm( CA.get_parent()['N'].coord - C_coord )
                    C_N_orientation = [np.dot( CA.get_parent()['N'].coord - C_coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['N'].coord - C_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['N'].coord - C_coord , orthognal )/np.linalg.norm( orthognal )]
                    ## C -> O
                    C_O_dis = np.linalg.norm( CA.get_parent()['O'].coord - C_coord )
                    C_O_orientation = [np.dot( CA.get_parent()['O'].coord - C_coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['O'].coord - C_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['O'].coord - C_coord , orthognal )/np.linalg.norm( orthognal )]
                    ## C -> CA
                    C_CA_dis = np.linalg.norm( CA.coord - C_coord )
                    C_CA_orientation = [np.dot( CA.coord - C_coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.coord - C_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.coord - C_coord , orthognal )/np.linalg.norm( orthognal )]
                    ## C ->CB
                    if CA.get_parent().resname !='GLY':
                        C_CB_dis = np.linalg.norm( CA.get_parent()['CB'].coord - C_coord )
                        C_CB_orientation = [np.dot( CA.get_parent()['CB'].coord - C_coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['CB'].coord - C_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['CB'].coord - C_coord , orthognal )/np.linalg.norm( orthognal )]
                    else:
                        C_CB_dis = C_CA_dis
                        C_CB_orientation = C_CA_orientation
                    C_CG1_dis = np.linalg.norm( CA.get_parent()['CG1'].coord - C_coord )
                    C_CG1_orientation = [np.dot( CA.get_parent()['CG1'].coord - C_coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['CG1'].coord - C_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['CG1'].coord - C_coord , orthognal )/np.linalg.norm( orthognal )]
                    C_CG2_dis = np.linalg.norm( CA.get_parent()['CG2'].coord - C_coord )
                    C_CG2_orientation = [np.dot( CA.get_parent()['CG2'].coord - C_coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['CG2'].coord - C_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['CG2'].coord - C_coord , orthognal )/np.linalg.norm( orthognal )]
                    
                    ## CB ->CB
                    if chain[int(matrix .iloc[i,3])].resname !='GLY':
                        CB_C_dis = np.linalg.norm( CA.get_parent()['C'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord )
                        CB_C_orientation = [np.dot( CA.get_parent()['C'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['C'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['C'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord , orthognal )/np.linalg.norm( orthognal )]
                        CB_N_dis = np.linalg.norm( CA.get_parent()['N'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord )
                        CB_N_orientation =  [np.dot( CA.get_parent()['N'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['N'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['N'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord , orthognal )/np.linalg.norm( orthognal )]
                        CB_O_dis = np.linalg.norm( CA.get_parent()['O'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord )
                        CB_O_orientation = [np.dot( CA.get_parent()['O'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['O'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['O'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord , orthognal )/np.linalg.norm( orthognal )]
                        CB_CA_dis = np.linalg.norm( CA.get_parent()['CA'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord )
                        CB_CA_orientation = [np.dot( CA.get_parent()['CA'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['CA'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['CA'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord , orthognal )/np.linalg.norm( orthognal )]
                        if CA.get_parent().resname !='GLY':
                            CB_CB_dis = np.linalg.norm( CA.get_parent()['CB'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord )
                            CB_CB_orientation = [np.dot( CA.get_parent()['CB'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['CB'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['CB'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord , orthognal )/np.linalg.norm( orthognal )]
                        else:
                            CB_CB_dis = CB_CA_dis
                            CB_CB_orientation = CB_CA_orientation
                    else:
                        CB_C_dis = CA_C_dis #np.linalg.norm( CA.get_parent()['C'].coord - CA_coord )
                        CB_C_orientation = CA_C_orientation #[np.dot( (CA.get_parent()['C'].coord - CA_coord), CA_C )/np.linalg.norm( CA_C ), np.dot( (CA.get_parent()['C'].coord - CA_coord), CA_N )/np.linalg.norm( CA_N ), np.dot( (CA.get_parent()['C'].coord - CA_coord), orthognal )/np.linalg.norm( orthognal )]
                        CB_N_dis = CA_N_dis #np.linalg.norm( CA.get_parent()['N'].coord - CA_coord )
                        CB_N_orientation = CA_N_orientation #[np.dot( (CA.get_parent()['N'].coord - CA_coord), CA_C )/np.linalg.norm( CA_C ), np.dot( (CA.get_parent()['N'].coord - CA_coord), CA_N )/np.linalg.norm( CA_N ), np.dot( (CA.get_parent()['N'].coord - CA_coord), orthognal )/np.linalg.norm( orthognal )]    
                        CB_O_dis = CA_O_dis #np.linalg.norm( CA.get_parent()['O'].coord - CA_coord )
                        CB_O_orientation = CA_O_orientation #[np.dot( (CA.get_parent()['O'].coord - CA_coord), CA_C )/np.linalg.norm( CA_C ), np.dot( (CA.get_parent()['O'].coord - CA_coord), CA_N )/np.linalg.norm( CA_N ), np.dot( (CA.get_parent()['O'].coord - CA_coord), orthognal )/np.linalg.norm( orthognal )]    
                        CB_CA_dis = CA_CA_dis
                        CB_CA_orientation = CA_CA_orientation
                        if CA.get_parent().resname !='GLY':
                            CB_CB_dis = np.linalg.norm( CA.get_parent()['CB'].coord - CA_coord )
                            CB_CB_orientation = [np.dot( CA.get_parent()['CB'].coord - CA_coord , CA_C )/np.linalg.norm( CA_C ), np.dot( CA.get_parent()['CB'].coord - CA_coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['CB'].coord - CA_coord , orthognal )/np.linalg.norm( orthognal )]
                        else:
                            CB_CB_dis = CA_CA_dis
                            CB_CB_orientation = CA_CA_orientation
                    CB_CG1_dis = np.linalg.norm( CA.get_parent()['CG1'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord )
                    CB_CG1_orientation = [np.dot( CA.get_parent()['CG1'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['CG1'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['CG1'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord , orthognal )/np.linalg.norm( orthognal )]
                    CB_CG2_dis = np.linalg.norm( CA.get_parent()['CG2'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord )
                    CB_CG2_orientation = [np.dot( CA.get_parent()['CG2'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['CG2'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['CG2'].coord - chain[int(matrix .iloc[i,3])]['CB'].coord , orthognal )/np.linalg.norm( orthognal )]        
                    
                    ##CG1 --> C
                    CG1_C_dis = np.linalg.norm( CA.get_parent()['C'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord )
                    CG1_C_orientation = [np.dot( CA.get_parent()['C'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['C'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['C'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord , orthognal )/np.linalg.norm( orthognal )]
                    CG1_N_dis = np.linalg.norm( CA.get_parent()['N'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord )
                    CG1_N_orientation = [np.dot( CA.get_parent()['N'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['N'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['N'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord , orthognal )/np.linalg.norm( orthognal )]
                    CG1_O_dis = np.linalg.norm( CA.get_parent()['O'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord )
                    CG1_O_orientation = [np.dot( CA.get_parent()['O'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['O'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['O'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord , orthognal )/np.linalg.norm( orthognal )]
                    CG1_CA_dis = np.linalg.norm( CA.get_parent()['CA'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord )
                    CG1_CA_orientation = [np.dot( CA.get_parent()['CA'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['CA'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['CA'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord , orthognal )/np.linalg.norm( orthognal )]
                    CG1_CB_dis = np.linalg.norm( CA.get_parent()['CB'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord )
                    CG1_CB_orientation = [np.dot( CA.get_parent()['CB'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['CB'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['CB'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord , orthognal )/np.linalg.norm( orthognal )]
                    CG1_CG1_dis = np.linalg.norm( CA.get_parent()['CG1'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord )
                    CG1_CG1_orientation = [np.dot( CA.get_parent()['CG1'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['CG1'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['CG1'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord , orthognal )/np.linalg.norm( orthognal )]
                    CG1_CG2_dis = np.linalg.norm( CA.get_parent()['CG2'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord )
                    CG1_CG2_orientation = [np.dot( CA.get_parent()['CG2'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['CG2'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['CG2'].coord - chain[int(matrix .iloc[i,3])]['CG1'].coord , orthognal )/np.linalg.norm( orthognal )]

                    ##CG2-->
                    CG2_C_dis = np.linalg.norm( CA.get_parent()['C'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord )
                    CG2_C_orientation = [np.dot( CA.get_parent()['C'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['C'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['C'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord , orthognal )/np.linalg.norm( orthognal )]
                    CG2_N_dis = np.linalg.norm( CA.get_parent()['N'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord )
                    CG2_N_orientation = [np.dot( CA.get_parent()['N'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['N'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['N'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord , orthognal )/np.linalg.norm( orthognal )]
                    CG2_O_dis = np.linalg.norm( CA.get_parent()['O'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord )
                    CG2_O_orientation = [np.dot( CA.get_parent()['O'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['O'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['O'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord , orthognal )/np.linalg.norm( orthognal )]
                    CG2_CA_dis = np.linalg.norm( CA.get_parent()['CA'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord )
                    CG2_CA_orientation = [np.dot( CA.get_parent()['CA'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['CA'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['CA'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord , orthognal )/np.linalg.norm( orthognal )]
                    CG2_CB_dis = np.linalg.norm( CA.get_parent()['CB'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord )
                    CG2_CB_orientation = [np.dot( CA.get_parent()['CB'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['CB'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['CB'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord , orthognal )/np.linalg.norm( orthognal )]
                    CG2_CG1_dis = np.linalg.norm( CA.get_parent()['CG1'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord )
                    CG2_CG1_orientation = [np.dot( CA.get_parent()['CG1'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['CG1'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['CG1'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord , orthognal )/np.linalg.norm( orthognal )]
                    CG2_CG2_dis = np.linalg.norm( CA.get_parent()['CG2'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord )
                    CG2_CG2_orientation = [np.dot( CA.get_parent()['CG2'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord , CA_C )/np.linalg.norm( CA_C ),  np.dot( CA.get_parent()['CG2'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord , CA_N )/np.linalg.norm( CA_N ), np.dot( CA.get_parent()['CG2'].coord - chain[int(matrix .iloc[i,3])]['CG2'].coord , orthognal )/np.linalg.norm( orthognal )]

                    each_edge_attributes.extend([CA_CA_dis,1/CA_CA_dis**2,1/CA_CA_dis**4,1/CA_CA_dis**6])
                    each_edge_attributes.extend(CA_CA_orientation)
                    each_edge_attributes.extend([CA_C_dis,1/CA_C_dis**2,1/CA_C_dis**4,1/CA_C_dis**6])
                    each_edge_attributes.extend(CA_C_orientation)
                    each_edge_attributes.extend([CA_O_dis,1/CA_O_dis**2,1/CA_O_dis**4,1/CA_O_dis**6])
                    each_edge_attributes.extend(CA_O_orientation)
                    each_edge_attributes.extend([CA_N_dis,1/CA_N_dis**2,1/CA_N_dis**4,1/CA_N_dis**6])
                    each_edge_attributes.extend(CA_N_orientation)
                    each_edge_attributes.extend([CA_CB_dis,1/CA_CB_dis**2,1/CA_CB_dis**4,1/CA_CB_dis**6])
                    each_edge_attributes.extend(CA_CB_orientation)
                    each_edge_attributes.extend([CA_CG1_dis,1/CA_CG1_dis**2,1/CA_CG1_dis**4,1/CA_CG1_dis**6])
                    each_edge_attributes.extend(CA_CG1_orientation)
                    each_edge_attributes.extend([CA_CG2_dis,1/CA_CG2_dis**2,1/CA_CG2_dis**4,1/CA_CG2_dis**6])
                    each_edge_attributes.extend(CA_CG2_orientation)
                    
                    each_edge_attributes.extend([O_N_dis,1/O_N_dis**2,1/O_N_dis**4,1/O_N_dis**6])
                    each_edge_attributes.extend(O_N_orientation)
                    each_edge_attributes.extend([O_C_dis,1/O_C_dis**2,1/O_C_dis**4,1/O_C_dis**6])
                    each_edge_attributes.extend(O_C_orientation)
                    each_edge_attributes.extend([O_CA_dis,1/O_CA_dis**2,1/O_CA_dis**4,1/O_CA_dis**6])
                    each_edge_attributes.extend(O_CA_orientation)
                    each_edge_attributes.extend([O_O_dis,1/O_O_dis**2,1/O_O_dis**4,1/O_O_dis**6])
                    each_edge_attributes.extend(O_O_orientation)
                    each_edge_attributes.extend([O_CB_dis,1/O_CB_dis**2,1/O_CB_dis**4,1/O_CB_dis**6])
                    each_edge_attributes.extend(O_CB_orientation)
                    each_edge_attributes.extend([O_CG1_dis,1/O_CG1_dis**2,1/O_CG1_dis**4,1/O_CG1_dis**6])
                    each_edge_attributes.extend(O_CG1_orientation)
                    each_edge_attributes.extend([O_CG2_dis,1/O_CG2_dis**2,1/O_CG2_dis**4,1/O_CG2_dis**6])
                    each_edge_attributes.extend(O_CG2_orientation)

                    each_edge_attributes.extend([N_O_dis,1/N_O_dis**2,1/N_O_dis**4,1/N_O_dis**6])
                    each_edge_attributes.extend(N_O_orientation)
                    each_edge_attributes.extend([N_C_dis,1/N_C_dis**2,1/N_C_dis**4,1/N_C_dis**6])
                    each_edge_attributes.extend(N_C_orientation)
                    each_edge_attributes.extend([N_CA_dis,1/N_CA_dis**2,1/N_CA_dis**4,1/N_CA_dis**6])
                    each_edge_attributes.extend(N_CA_orientation)
                    each_edge_attributes.extend([N_N_dis,1/N_N_dis**2,1/N_N_dis**4,1/N_N_dis**6])
                    each_edge_attributes.extend(N_N_orientation)
                    each_edge_attributes.extend([N_CB_dis,1/N_CB_dis**2,1/N_CB_dis**4,1/N_CB_dis**6])
                    each_edge_attributes.extend(N_CB_orientation)
                    each_edge_attributes.extend([N_CG1_dis,1/N_CG1_dis**2,1/N_CG1_dis**4,1/N_CG1_dis**6])
                    each_edge_attributes.extend(N_CG1_orientation)
                    each_edge_attributes.extend([N_CG2_dis,1/N_CG2_dis**2,1/N_CG2_dis**4,1/N_CG2_dis**6])
                    each_edge_attributes.extend(N_CG2_orientation)

                    each_edge_attributes.extend([C_O_dis,1/C_O_dis**2,1/C_O_dis**4,1/C_O_dis**6])
                    each_edge_attributes.extend(C_O_orientation)
                    each_edge_attributes.extend([C_N_dis,1/C_N_dis**2,1/C_N_dis**4,1/C_N_dis**6])
                    each_edge_attributes.extend(C_N_orientation)
                    each_edge_attributes.extend([C_C_dis,1/C_C_dis**2,1/C_C_dis**4,1/C_C_dis**6])
                    each_edge_attributes.extend(C_C_orientation)
                    each_edge_attributes.extend([C_CA_dis,1/C_CA_dis**2,1/C_CA_dis**4,1/C_CA_dis**6])
                    each_edge_attributes.extend(C_CA_orientation)
                    each_edge_attributes.extend([C_CB_dis,1/C_CB_dis**2,1/C_CB_dis**4,1/C_CB_dis**6])
                    each_edge_attributes.extend(C_CB_orientation)
                    each_edge_attributes.extend([C_CG1_dis,1/C_CG1_dis**2,1/C_CG1_dis**4,1/C_CG1_dis**6])
                    each_edge_attributes.extend(C_CG1_orientation)
                    each_edge_attributes.extend([C_CG2_dis,1/C_CG2_dis**2,1/C_CG2_dis**4,1/C_CG2_dis**6])
                    each_edge_attributes.extend(C_CG2_orientation)

                    each_edge_attributes.extend([CB_O_dis,1/CB_O_dis**2,1/CB_O_dis**4,1/CB_O_dis**6])
                    each_edge_attributes.extend(CB_O_orientation)
                    each_edge_attributes.extend([CB_N_dis,1/CB_N_dis**2,1/CB_N_dis**4,1/CB_N_dis**6])
                    each_edge_attributes.extend(CB_N_orientation)
                    each_edge_attributes.extend([CB_C_dis,1/CB_C_dis**2,1/CB_C_dis**4,1/CB_C_dis**6])
                    each_edge_attributes.extend(CB_C_orientation)
                    each_edge_attributes.extend([CB_CA_dis,1/CB_CA_dis**2,1/CB_CA_dis**4,1/CB_CA_dis**6])
                    each_edge_attributes.extend(CB_CA_orientation)
                    each_edge_attributes.extend([CB_CB_dis,1/CB_CB_dis**2,1/CB_CB_dis**4,1/CB_CB_dis**6])
                    each_edge_attributes.extend(CB_CB_orientation)
                    each_edge_attributes.extend([CB_CG1_dis,1/CB_CG1_dis**2,1/CB_CG1_dis**4,1/CB_CG1_dis**6])
                    each_edge_attributes.extend(CB_CG1_orientation)
                    each_edge_attributes.extend([CB_CG2_dis,1/CB_CG2_dis**2,1/CB_CG2_dis**4,1/CB_CG2_dis**6])
                    each_edge_attributes.extend(CB_CG2_orientation)

                    each_edge_attributes.extend([CG1_O_dis,1/CG1_O_dis**2,1/CG1_O_dis**4,1/CG1_O_dis**6])
                    each_edge_attributes.extend(CG1_O_orientation)
                    each_edge_attributes.extend([CG1_N_dis,1/CG1_N_dis**2,1/CG1_N_dis**4,1/CG1_N_dis**6])
                    each_edge_attributes.extend(CG1_N_orientation)
                    each_edge_attributes.extend([CG1_C_dis,1/CG1_C_dis**2,1/CG1_C_dis**4,1/CG1_C_dis**6])
                    each_edge_attributes.extend(CG1_C_orientation)
                    each_edge_attributes.extend([CG1_CA_dis,1/CG1_CA_dis**2,1/CG1_CA_dis**4,1/CG1_CA_dis**6])
                    each_edge_attributes.extend(CG1_CA_orientation)
                    each_edge_attributes.extend([CG1_CB_dis,1/CG1_CB_dis**2,1/CG1_CB_dis**4,1/CG1_CB_dis**6])
                    each_edge_attributes.extend(CG1_CB_orientation)
                    each_edge_attributes.extend([CG1_CG1_dis,1/CG1_CG1_dis**2,1/CG1_CG1_dis**4,1/CG1_CG1_dis**6])
                    each_edge_attributes.extend(CG1_CG1_orientation)
                    each_edge_attributes.extend([CG1_CG2_dis,1/CG1_CG2_dis**2,1/CG1_CG2_dis**4,1/CG1_CG2_dis**6])
                    each_edge_attributes.extend(CG1_CG2_orientation)

                    each_edge_attributes.extend([CG2_O_dis,1/CG2_O_dis**2,1/CG2_O_dis**4,1/CG2_O_dis**6])
                    each_edge_attributes.extend(CG2_O_orientation)
                    each_edge_attributes.extend([CG2_N_dis,1/CG2_N_dis**2,1/CG2_N_dis**4,1/CG2_N_dis**6])
                    each_edge_attributes.extend(CG2_N_orientation)
                    each_edge_attributes.extend([CG2_C_dis,1/CG2_C_dis**2,1/CG2_C_dis**4,1/CG2_C_dis**6])
                    each_edge_attributes.extend(CG2_C_orientation)
                    each_edge_attributes.extend([CG2_CA_dis,1/CG2_CA_dis**2,1/CG2_CA_dis**4,1/CG2_CA_dis**6])
                    each_edge_attributes.extend(CG2_CA_orientation)
                    each_edge_attributes.extend([CG2_CB_dis,1/CG2_CB_dis**2,1/CG2_CB_dis**4,1/CG2_CB_dis**6])
                    each_edge_attributes.extend(CG2_CB_orientation)
                    each_edge_attributes.extend([CG2_CG1_dis,1/CG2_CG1_dis**2,1/CG2_CG1_dis**4,1/CG2_CG1_dis**6])
                    each_edge_attributes.extend(CG2_CG1_orientation)
                    each_edge_attributes.extend([CG2_CG2_dis,1/CG2_CG2_dis**2,1/CG2_CG2_dis**4,1/CG2_CG2_dis**6])
                    each_edge_attributes.extend(CG2_CG2_orientation)


                    #edge_attributes.append(each_edge_attributes)
                    ##beginnode to endnode 
                    endnode = matrix[ matrix.iloc[:,3] == CA.get_parent().id[1] ].iloc[0,0]
                    
                    #each_edge_attributes.extend([sin(endnode- i ),cos(endnode-i)])
                    each_edge_attributes.extend([sin(matrix.iloc[endnode,3] - matrix.iloc[i,3] ),cos(matrix.iloc[endnode,3] - matrix.iloc[i,3])])
                    edge_attributes.append(each_edge_attributes)

                    each_edge_index.extend([i,endnode ])
                    edge_index.append(each_edge_index)
        
        assert len(node_feature)==len(position_encoding),'Wrong dimensional number when processing this pdbfile, please carefully check the pdbfile'

        return {'Node_features': np.array(node_feature), 'edge_attrs':np.array(edge_attributes), 'edge_index': np.array(edge_index), 'position_encoding': np.array(position_encoding)}


def triangularexpression(pdbfiles, chainID,edge_index, interval_b, interval_e):
    structure = Bio.PDB.PDBParser(PERMISSIVE=1).get_structure('tmp', pdbfiles)
    atoms  = Bio.PDB.Selection.unfold_entities(structure, 'A')
    ns = Bio.PDB.NeighborSearch(atoms)
    for model in structure:
        dssp = DSSP(model, pdbfiles, dssp='dssp')
        for chain in model:
            if chain.id == chainID:
                poly=Bio.PDB.Polypeptide.Polypeptide(chain)
                phi_psi_list = poly.get_phi_psi_list()
                sequence_list=poly.get_sequence()
                residue_list=[residue.id for residue in chain]
                matrix = [[i,chain.id,j,l,dssp[(chain.id,chain[int(l)].id)][2],n,p] for i,(j,(k,l,m),(n,p)) in enumerate(zip(sequence_list,residue_list,phi_psi_list))  ]
                matrix = pd.DataFrame(matrix)
                triangle_index_total = []
                triangle2line_total = []
                for i in range(int(interval_b), int(interval_e)): # len(matrix)
                    ## surrounding residues
                    NA = ns.search(chain[int(matrix .iloc[i,3])]['CA'].coord, 12)  ##Neighbor Atoms
                    #print(i,matrix.iloc[i,3],[atom.get_parent().id[1] for atom in NA if atom.id=='CA'])
                    All_CA = [atom for atom in NA if atom.id=='CA' and (atom.coord != chain[int(matrix .iloc[i,3])]['CA'].coord).any() and atom.get_parent().get_parent().id == chain.id]   # atom.get_parent().id[1]
                    CA_coords = []
                    endnode = []
                    for CA in All_CA:
                        endnode.append(matrix[ matrix.iloc[:,3] == CA.get_parent().id[1] ].iloc[0,0])
                        CA_coords.append(CA.coord.tolist())
                    dis_map = distance_matrix(CA_coords,CA_coords)
                    dis_map = np.triu(dis_map)  ## 上三角
                    index_less12 = np.argwhere((dis_map<12) & (dis_map>0))
                    triangle_index = []
                    for idx_i,idx_j in index_less12:
                        triangle_index.append([endnode[idx_i],i,endnode[idx_j]])
                    triangle2line = []
                    for tr2l_i, tr2l_j,tr2l_k in triangle_index:
                        #print(structure.id,chain.id,'ij',np.argwhere( (edge_index[:,0]==tr2l_i) & (edge_index[:,1]==tr2l_j)),end ='\t')
                        ij_index = np.argwhere( (edge_index[:,0]==tr2l_i) & (edge_index[:,1]==tr2l_j)).item()
                        #print(structure.id,chain.id,'kj',np.argwhere( (edge_index[:,0]==tr2l_k) & (edge_index[:,1]==tr2l_j)),end ='\t')
                        kj_index = np.argwhere( (edge_index[:,0]==tr2l_k) & (edge_index[:,1]==tr2l_j)).item()
                        #print(structure.id,chain.id,'ik',np.argwhere( (edge_index[:,0]==tr2l_i) & (edge_index[:,1]==tr2l_k)),end ='\t')
                        ik_index = np.argwhere( (edge_index[:,0]==tr2l_i) & (edge_index[:,1]==tr2l_k)).item()
                        #print(structure.id,chain.id,'ki',np.argwhere( (edge_index[:,0]==tr2l_k) & (edge_index[:,1]==tr2l_i)),end ='\t')
                        ki_index = np.argwhere( (edge_index[:,0]==tr2l_k) & (edge_index[:,1]==tr2l_i)).item()
                        triangle2line.append([ij_index,kj_index,kj_index,ki_index])
                    triangle_index_total.extend(triangle_index)
                    triangle2line_total.extend(triangle2line)
                
                return {'triangleindex': np.array(triangle_index_total), 'triangle2line':np.array(triangle2line_total)}

from multiprocessing import *
num_processes = 6
def batch_fea( files,chainID,seqlength ):
    interval = [int(i) for i in np.linspace(0,seqlength, num_processes) ]
    pools = Pool(num_processes)
    results = []
    #loadpdbs(files, chainID, interval_b, interval_e)
    for i in range(len(interval)-1):
        subresult = pools.apply_async( loadpdbs,args=(files,chainID,interval[i],interval[i+1],))
        results.append(subresult)
    pools.close()
    pools.join()
    Node_features = np.empty([0,26])
    edge_attrs = np.empty([0,345])
    edge_index = np.empty([0,2])
    
    for result in results:
        result = result.get()
        Node_features = np.concatenate((Node_features,result['Node_features']),axis=0)
        edge_attrs = np.concatenate( (edge_attrs,result['edge_attrs']),axis=0)
        edge_index = np.concatenate((edge_index,result['edge_index']),axis=0)
    assert Node_features.shape[0]== seqlength,'Sequence length does not match the shape of protein nodes'
    return {'Node_features': Node_features, 'edge_attrs': edge_attrs, 'edge_index': edge_index }

def batch_tri( files,chainID,seqlength, edgeindex):
    interval = [int(i) for i in np.linspace(0,seqlength, num_processes) ]
    pools = Pool(num_processes)
    results = []
    #triangularexpression(pdbfiles, chainID,edge_index, interval_b, interval_e)
    for i in range(len(interval)-1):
        subresult = pools.apply_async( triangularexpression,args=(files,chainID, edgeindex, interval[i],interval[i+1],))
        results.append(subresult)
    pools.close()
    pools.join()

    triangleindex = np.empty((0,3))
    triangle2line = np.empty((0,4))
    for result in results:
        result = result.get()
        print(result['triangleindex'].shape,result['triangle2line'].shape)
        triangleindex = np.concatenate((triangleindex,result['triangleindex']),axis=0)
        triangle2line = np.concatenate((triangle2line,result['triangle2line']),axis=0)
    return {'triangleindex': triangleindex, 'triangle2line':triangle2line}




