# -*- coding: utf-8 -*-

###################################################################################
#### This code is a part from OPUS-Rota4, https://doi.org/10.1093/bib/bbab529.  ###
#### If you find use, please cite it.                                           ###
###################################################################################
"""
Created on Wed Apr 27 13:27:07 2022
 
@author: liujl
"""

num_dihedrals_dict = {"G":0, "A":0, "S":1, "C":1, "V":1, "I":2, "L":2, "T":1, "R":4, "K":4,
                      "D":2, "N":2, "E":3, "Q":3, "M":3, "H":2, "P":1, "F":2, "Y":2, "W":2}
num_side_chain_atoms_dict = {"G":0, "A":0, "S":1, "C":1, "V":2, "I":3, "L":3, "T":2, "R":6, "K":4,
                             "D":3, "N":3, "E":4, "Q":4, "M":3, "H":5, "P":2, "F":6, "Y":7, "W":9}
class Residue:
    def __init__(self, resid, resname, chain1):
        self.resid = resid       
        self.resname = resname
        self.chain1 = chain1
        self.atoms = {}
        
        self.main_chain_atoms_matrixid = {}
        self.side_chain_atoms_matrixid = {}
        
        self.num_dihedrals = num_dihedrals_dict[resname]
        self.num_side_chain_atoms = num_side_chain_atoms_dict[resname]
        if resname == 'G':
            self.num_atoms = self.num_side_chain_atoms + 4
        else:
            self.num_atoms = self.num_side_chain_atoms + 5

        self.resname_tri = triResname(resname)
        
def getResidueData(atomsData):
    
    residuesData = []
    last_resid = None
    for atom in atomsData:
        
        if atom.resname == "I" and atom.name1 == "CD":
            atom.name1 = "CD1"
            
        if(atom == atomsData[0]):           
            residue = Residue(atom.resid,atom.resname,atom.chain1) 
            residue.atoms[atom.name1] = atom
        elif(atom == atomsData[-1]):
            if(last_resid == atom.resid):                
                residue.atoms[atom.name1] = atom  
                residuesData.append(residue)
            else:                
                residue = Residue(atom.resid,atom.resname,atom.chain1) 
                residue.atoms[atom.name1] = atom
                residuesData.append(residue)
        else:
            if(last_resid == atom.resid):                
                residue.atoms[atom.name1] = atom          
            else:
                residuesData.append(residue)
                residue = Residue(atom.resid,atom.resname,atom.chain1) 
                residue.atoms[atom.name1] = atom                  
        
        last_resid = atom.resid
    
    return residuesData

def singleResname(AA):
    if(len(AA) == 1):
        return AA
    else:
        if(AA in ['GLY','AGLY']):
            return "G"
        elif(AA in ['ALA','AALA']):
            return "A"
        elif(AA in ['SER','ASER']):
            return "S"
        elif(AA in ['CYS','ACYS']):
            return "C"
        elif(AA in ['VAL','AVAL']):
            return "V"
        elif(AA in ['ILE','AILE']):
            return "I"
        elif(AA in ['LEU','ALEU']):
            return "L"
        elif(AA in ['THR','ATHR']):
            return "T"
        elif(AA in ['ARG','AARG']):
            return "R"
        elif(AA in ['LYS','ALYS']):
            return "K"
        elif(AA in ['ASP','AASP']):
            return "D"
        elif(AA in ['GLU','AGLU']):
            return "E"
        elif(AA in ['ASN','AASN']):
            return "N"
        elif(AA in ['GLN','AGLN']):
            return "Q"
        elif(AA in ['MET','AMET']):
            return "M"
        elif(AA in ['HIS','AHIS','HSD']):
            return "H"
        elif(AA in ['PRO','APRO']):
            return "P"
        elif(AA in ['PHE','APHE']):
            return "F"
        elif(AA in ['TYR','ATYR']):
            return "Y"
        elif(AA in ['TRP','ATRP']):
            return "W"
        else:
            return None

def triResname(AA):
    if(len(AA) == 3):
        return AA
    else:
        if(AA == "G"):
            return "GLY"
        elif(AA == "A"):
            return "ALA"
        elif(AA == "S"):
            return "SER"
        elif(AA == "C"):
            return "CYS"
        elif(AA == "V"):
            return "VAL"
        elif(AA == "I"):
            return "ILE"
        elif(AA == "L"):
            return "LEU"
        elif(AA == "T"):
            return "THR"
        elif(AA == "R"):
            return "ARG"
        elif(AA == "K"):
            return "LYS"
        elif(AA == "D"):
            return "ASP"
        elif(AA == "E"):
            return "GLU"
        elif(AA == "N"):
            return "ASN"
        elif(AA == "Q"):
            return "GLN"
        elif(AA == "M"):
            return "MET"
        elif(AA == "H"):
            return "HIS"
        elif(AA == "P"):
            return "PRO"
        elif(AA == "F"):
            return "PHE"
        elif(AA == "Y"):
            return "TYR"
        elif(AA == "W"):
            return "TRP"
        else:
            return None
