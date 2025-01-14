
import Bio
import Bio.PDB
from Bio.PDB.DSSP import DSSP

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


def loadpdbs(files, chainID):
    structure = Bio.PDB.PDBParser(PERMISSIVE=1).get_structure('tmp', files)

    #delete irregular residue atoms
    for model in structure:
        for chain in model:
        if chain.id == argv[3]:
            delete_REScluster=[]
            for residue in chain:
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
        dssp = DSSP(model, files,dssp='dssp')

        for chain in model:
            if chain.id == chainID:
                poly=Bio.PDB.Polypeptide.Polypeptide(chain)
                phi_psi_list = poly.get_phi_psi_list()
                print("phi_psi_list: ",phi_psi_list)  #dtype: list
                sequence_list=poly.get_sequence()
                print("sequence list: ",sequence_list)  #dtype
                residue_list=[residue.id for residue in chain]
                print("residue list: ",residue_list) #dtype
                matrix = [[i,chain.id,j,l,dssp[(chain.id,chain[int(l)].id)][2],n,p] for i,(j,(k,l,m),(n,p)) in enumerate(zip(sequence_list,residue_list,phi_psi_list))  ]

                matrix = pd.DataFrame(matrix)
                node_feature = []
                edge_index =[]
                edge_attributes = []
                for i in range(len(matrix)):
                    each_node_feature = []
                    if pd.isna(matrix.iloc[i,-2]):
                        each_node_feature.extend( [ 0,0,sin(matrix.iloc[i,-1]),cos(matrix.iloc[i,-1]),   0,0,0,0,0,0,sin(matrix.iloc[i,-1]*2),sin(matrix.iloc[i,-1]*3),sin(matrix.iloc[i,-1]*4),cos(matrix.iloc[i,-1]*2),cos(matrix.iloc[i,-1]*3),cos(matrix.iloc[i,-1]*4) ]     )
                    elif pd.isna(matrix.iloc[i,-1]):
                        each_node_feature.extend( [ sin(matrix.iloc[i,-2]),cos(matrix.iloc[i,-2]),0,0,   sin(matrix.iloc[i,-2]*2),sin(matrix.iloc[i,-2]*3),sin(matrix.iloc[i,-2]*4),cos(matrix.iloc[i,-2]*2),cos(matrix.iloc[i,-2]*3),cos(matrix.iloc[i,-2]*4),0,0,0,0,0,0 ]      )
                    else:
                        each_node_feature.extend( [ sin(matrix.iloc[i,-2]),cos(matrix.iloc[i,-2]),sin(matrix.iloc[i,-1]),cos(matrix.iloc[i,-1]),       sin(matrix.iloc[i,-2]*2),sin(matrix.iloc[i,-2]*3),sin(matrix.iloc[i,-2]*4),cos(matrix.iloc[i,-2]*2),cos(matrix.iloc[i,-2]*3),cos(matrix.iloc[i,-2]*4),       sin(matrix.iloc[i,-1]*2),sin(matrix.iloc[i,-1]*3),sin(matrix.iloc[i,-1]*4),cos(matrix.iloc[i,-1]*2),cos(matrix.iloc[i,-1]*3),cos(matrix.iloc[i,-1]*4)  ]     )
                    each_node_feature.extend( eight2three( matrix.iloc[i,4] )  )

                    NA = ns.search(chain[int(matrix .iloc[i,3])]['CA'].coord, 12)  ##Neighbor Atoms
                    All_CA = [atom for atom in NA if atom.id=='CA']
                    if len(All_CA) ==1:
                        each_node_feature.extend([len(All_CA)-1, 10000])
                    else:
                        each_node_feature.extend([len(All_CA)-1, 1/(len(All_CA)-1)])
                    
                    node_feature.append(each_node_feature )

                    for CA in All_CA:
                        if (CA.coord != chain[int(matrix .iloc[i,3])]['CA'].coord).any() and CA.get_parent().get_parent().id == chain.id :
                            each_edge_attributes = []
                            each_edge_index =[]
                            ##beginnode to endnode
                            endnode = matrix[ matrix.iloc[:,3] == CA.get_parent().id[1] ].iloc[0,0]
                            #each_edge_attributes.extend([sin(endnode- i ),cos(endnode-i)])
                            each_edge_attributes.extend([sin(matrix.iloc[endnode,3] - matrix.iloc[i,3] ),cos(matrix.iloc[endnode,3] - matrix.iloc[i,3])])
                            edge_attributes.append(each_edge_attributes)

                            each_edge_index.extend([i,endnode ])
                            edge_index.append(each_edge_index)
                
                return {'Node_features': np.array(node_feature), 'edge_attrs':np.array(edge_attributes), 'edge_index': np.array(edge_index) }

