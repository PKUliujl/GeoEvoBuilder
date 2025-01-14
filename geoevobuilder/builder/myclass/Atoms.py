# -*- coding: utf-8 -*-
###################################################################################
#### This code is a part from OPUS-Rota4, https://doi.org/10.1093/bib/bbab529.  ###
#### If you find use, please cite it.                                           ###
###################################################################################
"""
Created on Wed Apr 27 13:29:44 2022

@author: liujl
"""

from geoevobuilder.builder.myclass import Residues

class Atom:
    def __init__(self, atomid, name1, resname, chain1, resid, position):
        self.atomid = atomid
        self.name1 = name1
        self.chain1 = chain1
        self.resname = Residues.singleResname(resname)
        self.resid = resid
        self.position = position
        
        if self.name1 in ['N','CA','C','O','CB']:
            self.ismainchain = True
        else:
            self.ismainchain = False 
