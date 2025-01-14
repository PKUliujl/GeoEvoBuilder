'''This module is part of the PeptideBuilder library,
written by Matthew Z. Tien, Dariya K. Sydykova,
Austin G. Meyer, and Claus O. Wilke.

The Geometry module contains the default geometries of
all 20 amino acids. The main function to be used is the
geometry() function, which returns the default geometry
for the requested amino acid.

This file is provided to you under the GNU General Public
License, version 2.0 or later.'''

class Geo():
    '''Geometry base class'''
    def __repr__(self):
        repr = ""
        for var in dir(self):
            if var in self.__dict__: # exclude member functions, only print member variables
                repr += "%s = %s\n" % ( var, self.__dict__[var] )
        return repr


class GlyGeo(Geo):
    '''Geometry of Glycine'''
    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 110.8914

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.5117
        self.N_CA_C_O_diangle = 180.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.residue_name= 'G'
    
class AlaGeo(Geo):
    '''Geometry of Alanin'''
    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 111.068

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.5
        self.N_CA_C_O_diangle = -60.5

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277
    

        self.CA_CB_length = 1.53
        self.C_CA_CB_angle = 110.5
        self.N_C_CA_CB_diangle = -122.5

        self.residue_name= 'A'

class SerGeo(Geo):
    '''Geometry of Serine'''
    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 111.2812

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.5
        self.N_CA_C_O_diangle = -60.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277
    

        self.CA_CB_length = 1.53
        self.C_CA_CB_angle = 110.5
        self.N_C_CA_CB_diangle = -122.5

        self.CB_OG_length = 1.417
        self.CA_CB_OG_angle = 110.8
        self.N_CA_CB_OG_diangle=-63.3

        self.residue_name= 'S'

    def inputRotamers(self, rotamers):
        self.N_CA_CB_OG_diangle=rotamers[0]
        
class CysGeo(Geo):                                        
    '''Geometry of Cystine'''
    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 110.8856

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.5
        self.N_CA_C_O_diangle = -60.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277
    

        self.CA_CB_length = 1.53
        self.C_CA_CB_angle = 110.5
        self.N_C_CA_CB_diangle = -122.5

        self.CB_SG_length = 1.807
        self.CA_CB_SG_angle = 114
        self.N_CA_CB_SG_diangle=-62.2

        self.residue_name= 'C'

    def inputRotamers(self, rotamers):
        self.N_CA_CB_SG_diangle=rotamers[0]
            
class ValGeo(Geo):                            #x    cg1cg2
    '''Geometry of Valine'''
    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 109.7698

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.5686
        self.N_CA_C_O_diangle = -60.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.546
        self.C_CA_CB_angle = 111.5
        self.N_C_CA_CB_diangle = -122.5

        self.CB_CG1_length = 1.521
        self.CA_CB_CG1_angle = 110.5
        self.N_CA_CB_CG1_diangle=177.2

        self.CB_CG2_length = 1.521
        self.CG1_CB_CG2_angle = 110.5 #CT-CT-CT
        self.CG1_CA_CB_CG2_diangle = 122.6#est

        self.residue_name= 'V'

    def inputRotamers(self, rotamers):
        self.N_CA_CB_CG1_diangle=rotamers[0]

#modify x2
class IleGeo(Geo):                                          #x cg1cg2
    '''Geometry of Isoleucine'''
    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 109.7202

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.5403
        self.N_CA_C_O_diangle = -60.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.546
        self.C_CA_CB_angle = 111.5
        self.N_C_CA_CB_diangle = -122.5

        self.CB_CG1_length = 1.53
        self.CA_CB_CG1_angle = 110.3
        self.N_CA_CB_CG1_diangle=59.7

        self.CB_CG2_length = 1.521
        self.CG1_CB_CG2_angle = 110.5 #CT-CT-CT
        self.CG1_CA_CB_CG2_diangle = -122.6#est

        self.CG1_CD1_length = 1.516
        self.CB_CG1_CD1_angle = 114
        self.CA_CB_CG1_CD1_diangle=169.8

        self.residue_name= 'I'

    def inputRotamers(self, rotamers):
        self.N_CA_CB_CG1_diangle=rotamers[0]
        self.CA_CB_CG1_CD1_diangle=rotamers[1]
                 
class LeuGeo(Geo):                              #x cd1cd2
    '''Geometry of Leucine'''
    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 110.8652

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.4647
        self.N_CA_C_O_diangle = 120.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.53
        self.C_CA_CB_angle = 110.5
        self.N_C_CA_CB_diangle = -122.5

        self.CB_CG_length = 1.53
        self.CA_CB_CG_angle = 116.3
        self.N_CA_CB_CG_diangle=-60.1

        self.CG_CD1_length = 1.521
        self.CB_CG_CD1_angle = 110.5
        self.CA_CB_CG_CD1_diangle=174.9

        self.CG_CD2_length = 1.521
        self.CD1_CG_CD2_angle = 110.5 #CT-CT-CT
        self.CD1_CB_CG_CD2_diangle = 122.6#est

        self.residue_name= 'L'

    def inputRotamers(self, rotamers):
        self.N_CA_CB_CG_diangle=rotamers[0]
        self.CA_CB_CG_CD1_diangle=rotamers[1]
            
class ThrGeo(Geo):                                 #x og1 cg2
    '''Geometry of Threonine'''
    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 110.7014

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.5359
        self.N_CA_C_O_diangle = 120.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.542
        self.C_CA_CB_angle = 111.5
        self.N_C_CA_CB_diangle = -122

        self.CB_OG1_length = 1.433
        self.CA_CB_OG1_angle = 109.5
        self.N_CA_CB_OG1_diangle=60.0

        self.CB_CG2_length = 1.521
        self.OG1_CB_CG2_angle = 110.5 #CT-CT-OH
        self.OG1_CA_CB_CG2_diangle = -120#est

        self.residue_name= 'T'

    def inputRotamers(self, rotamers):
        self.N_CA_CB_OG1_diangle=rotamers[0]

class ArgGeo(Geo):
    '''Geometry of Arginine'''
    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 110.98

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.54
        self.N_CA_C_O_diangle = 120.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.53
        self.C_CA_CB_angle = 110.5
        self.N_C_CA_CB_diangle = -122.5

        self.CB_CG_length = 1.52
        self.CA_CB_CG_angle = 114.1
        self.N_CA_CB_CG_diangle=-65.2

        self.CG_CD_length = 1.52
        self.CB_CG_CD_angle = 111.5
        self.CA_CB_CG_CD_diangle=-179.2

        self.CD_NE_length = 1.461
        self.CG_CD_NE_angle = 112
        self.CB_CG_CD_NE_diangle=-179.3

        self.NE_CZ_length = 1.33
        self.CD_NE_CZ_angle = 124.5
        self.CG_CD_NE_CZ_diangle=-178.7

        self.CZ_NH1_length = 1.326
        self.NE_CZ_NH1_angle = 120
        self.CD_NE_CZ_NH1_diangle = 0.0

        self.NH1_NH2_length = 1.326
        self.NE_CZ_NH2_angle = 120
        self.NH1_NE_CZ_NH2_diangle=180.0

        self.residue_name= 'R'

    def inputRotamers(self, rotamers):
        self.N_CA_CB_CG_diangle = rotamers[0]
        self.CA_CB_CG_CD_diangle = rotamers[1]
        self.CB_CG_CD_NE_diangle = rotamers[2]
        self.CG_CD_NE_CZ_diangle = rotamers[3]

class LysGeo(Geo):
    '''Geometry of Lysine'''
    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 111.08

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.54
        self.N_CA_C_O_diangle = 120.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.53
        self.C_CA_CB_angle = 110.5
        self.N_C_CA_CB_diangle = -122.5

        self.CB_CG_length = 1.52
        self.CA_CB_CG_angle = 114.1
        self.N_CA_CB_CG_diangle=-64.5

        self.CG_CD_length = 1.52
        self.CB_CG_CD_angle = 111.5
        self.CA_CB_CG_CD_diangle=-178.1

        self.CD_CE_length = 1.52
        self.CG_CD_CE_angle = 111.5
        self.CB_CG_CD_CE_diangle=-179.6

        self.CE_NZ_length = 1.489
        self.CD_CE_NZ_angle = 112
        self.CG_CD_CE_NZ_diangle=179.6

        self.residue_name= 'K'

    def inputRotamers(self, rotamers):
        self.N_CA_CB_CG_diangle=rotamers[0]
        self.CA_CB_CG_CD_diangle=rotamers[1]
        self.CB_CG_CD_CE_diangle=rotamers[2]
        self.CG_CD_CE_NZ_diangle=rotamers[3]
    
class AspGeo(Geo):
    '''Geometry of Aspartic Acid'''
    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 111.03

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.51
        self.N_CA_C_O_diangle = 120.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277
    
        self.CA_CB_length = 1.53
        self.C_CA_CB_angle = 110.5
        self.N_C_CA_CB_diangle = -122.5

        self.CB_CG_length = 1.516
        self.CA_CB_CG_angle = 112.7
        self.N_CA_CB_CG_diangle=-66.4

        self.CG_OD1_length = 1.25
        self.CB_CG_OD1_angle = 118.5
        self.CA_CB_CG_OD1_diangle=-46.7

        self.CG_OD2_length = 1.25
        self.CB_CG_OD2_angle = 118.5
        self.OD1_CB_CG_OD2_diangle = 180.0

        self.residue_name= 'D'

    def inputRotamers(self, rotamers):
        self.N_CA_CB_CG_diangle=rotamers[0]
        self.CA_CB_CG_OD1_diangle=rotamers[1]

class AsnGeo(Geo):
    '''Geometry of Asparagine'''
    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 111.5

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.4826
        self.N_CA_C_O_diangle = -60.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277
        
        self.CA_CB_length = 1.53
        self.C_CA_CB_angle = 110.5
        self.N_C_CA_CB_diangle = -122.5

        self.CB_CG_length = 1.516
        self.CA_CB_CG_angle = 112.7
        self.N_CA_CB_CG_diangle=-65.5

        self.CG_OD1_length = 1.231
        self.CB_CG_OD1_angle = 120.8
        self.CA_CB_CG_OD1_diangle=-58.3

        self.CG_ND2_length = 1.328
        self.CB_CG_ND2_angle = 116.5
        self.OD1_CB_CG_ND2_diangle = 180.0

        self.residue_name= 'N'

    def inputRotamers(self, rotamers):
        self.N_CA_CB_CG_diangle=rotamers[0]
        self.CA_CB_CG_OD1_diangle=rotamers[1]

class GluGeo(Geo):
    '''Geometry of Glutamic Acid'''
    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 111.1703

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.511
        self.N_CA_C_O_diangle = 120.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.53
        self.C_CA_CB_angle = 110.5
        self.N_C_CA_CB_diangle = -122.5

        self.CB_CG_length = 1.52
        self.CA_CB_CG_angle = 114.1
        self.N_CA_CB_CG_diangle=-63.8

        self.CG_CD_length = 1.516
        self.CB_CG_CD_angle = 112.7
        self.CA_CB_CG_CD_diangle=-179.8

        self.CD_OE1_length = 1.25
        self.CG_CD_OE1_angle = 118.5
        self.CB_CG_CD_OE1_diangle=-6.2

        self.CD_OE2_length = 1.25
        self.CG_CD_OE2_angle = 118.5    
        self.OE1_CG_CD_OE2_diangle = 180.0

        self.residue_name= 'E'

    def inputRotamers(self, rotamers):
        self.N_CA_CB_CG_diangle=rotamers[0]
        self.CA_CB_CG_CD_diangle=rotamers[1]
        self.CB_CG_CD_OE1_diangle=rotamers[2]
        
class GlnGeo(Geo):                                
    '''Geometry of Glutamine'''
    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 111.0849

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.5029
        self.N_CA_C_O_diangle = 120.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.53
        self.C_CA_CB_angle = 110.5
        self.N_C_CA_CB_diangle = -122.5

        self.CB_CG_length = 1.52
        self.CA_CB_CG_angle = 114.1
        self.N_CA_CB_CG_diangle=-60.2

        self.CG_CD_length = 1.516
        self.CB_CG_CD_angle = 112.7
        self.CA_CB_CG_CD_diangle=-69.6

        self.CD_OE1_length = 1.231
        self.CG_CD_OE1_angle = 120.8
        self.CB_CG_CD_OE1_diangle=-50.5

        self.CD_NE2_length = 1.328
        self.CG_CD_NE2_angle = 116.5
        self.OE1_CG_CD_NE2_diangle = 180.0

        self.residue_name= 'Q'

    def inputRotamers(self, rotamers):
        self.N_CA_CB_CG_diangle=rotamers[0]
        self.CA_CB_CG_CD_diangle=rotamers[1]
        self.CB_CG_CD_OE1_diangle=rotamers[2]
            
class MetGeo(Geo):                                    
    '''Geometry of Methionine'''
    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 110.9416

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.4816
        self.N_CA_C_O_diangle = 120.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.53
        self.C_CA_CB_angle = 110.5
        self.N_C_CA_CB_diangle = -122.5

        self.CB_CG_length = 1.52
        self.CA_CB_CG_angle = 114.1
        self.N_CA_CB_CG_diangle=-64.4

        self.CG_SD_length = 1.807
        self.CB_CG_SD_angle = 112.7
        self.CA_CB_CG_SD_diangle=-179.6

        self.SD_CE_length = 1.789
        self.CG_SD_CE_angle = 100.8
        self.CB_CG_SD_CE_diangle=70.1

        self.residue_name= 'M'
    def inputRotamers(self, rotamer):
        self.N_CA_CB_CG_diangle=rotamer[0]
        self.CA_CB_CG_SD_diangle=rotamer[1]
        self.CB_CG_SD_CE_diangle=rotamer[2]
    

class HisGeo(Geo):                               
    '''Geometry of Histidine'''
    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 111.0859

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.4732
        self.N_CA_C_O_diangle = 120.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.53
        self.C_CA_CB_angle = 110.5
        self.N_C_CA_CB_diangle = -122.5

        self.CB_CG_length = 1.5
        self.CA_CB_CG_angle = 113.8
        self.N_CA_CB_CG_diangle=-63.2
        
        self.CG_ND1_length = 1.378
        self.CB_CG_ND1_angle = 122.7
        self.CA_CB_CG_ND1_diangle=-75.7          

        self.CG_CD2_length = 1.354
        self.CB_CG_CD2_angle = 131
        self.ND1_CB_CG_CD2_diangle = 180.0

        self.ND1_CE1_length = 1.32
        self.CG_ND1_CE1_angle = 109.2
        self.CB_CG_ND1_CE1_diangle = 180.0

        self.CD2_NE2_length = 1.374
        self.CG_CD2_NE2_angle = 107.2
        self.CB_CG_CD2_NE2_diangle = 180.0

        self.residue_name= 'H'

    def inputRotamers(self, rotamers):
        self.N_CA_CB_CG_diangle=rotamers[0]
        self.CA_CB_CG_ND1_diangle=rotamers[1]

class ProGeo(Geo):
    '''Geometry of Proline'''
    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 112.7499

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.2945
        self.N_CA_C_O_diangle = -45.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.53
        self.C_CA_CB_angle = 103.2
        self.N_C_CA_CB_diangle = -120

        self.CB_CG_length = 1.495
        self.CA_CB_CG_angle = 104.5
        self.N_CA_CB_CG_diangle=29.6

        self.CG_CD_length = 1.507
        self.CB_CG_CD_angle = 105.5
        self.CA_CB_CG_CD_diangle=-34.8

        self.residue_name= 'P'
    
    def inputRotamers(self, rotamer):
        self.N_CA_CB_CG_diangle=rotamer[0]
        self.CA_CB_CG_CD_diangle=rotamer[1]
            
class PheGeo(Geo):
    '''Geometry of Phenylalanine'''
    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 110.7528

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.5316
        self.N_CA_C_O_diangle = 120.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.53
        self.C_CA_CB_angle = 110.5
        self.N_C_CA_CB_diangle = -122.5

        self.CB_CG_length = 1.50
        self.CA_CB_CG_angle = 113.8
        self.N_CA_CB_CG_diangle=-64.7

        self.CG_CD1_length = 1.391
        self.CB_CG_CD1_angle = 120.7
        self.CA_CB_CG_CD1_diangle=93.3

        self.CG_CD2_length = 1.391
        self.CB_CG_CD2_angle = 120.7
        self.CD1_CB_CG_CD2_diangle = 180.0

        self.CD1_CE1_length = 1.393
        self.CG_CD1_CE1_angle = 120.7
        self.CB_CG_CD1_CE1_diangle = 180.0

        self.CD2_CE2_length = 1.393
        self.CG_CD2_CE2_angle = 120.7
        self.CB_CG_CD2_CE2_diangle = 180.0

        self.CE1_CZ_length = 1.39
        self.CD1_CE1_CZ_angle = 120.0
        self.CG_CD1_CE1_CZ_diangle = 0.0

        self.residue_name= 'F'

    def inputRotamers(self, rotamers):
        self.N_CA_CB_CG_diangle=rotamers[0]
        self.CA_CB_CG_CD1_diangle=rotamers[1]

class TyrGeo(Geo):                                             
    '''Geometry of Tyrosine'''
    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 110.9288

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.5434
        self.N_CA_C_O_diangle = 120.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.53
        self.C_CA_CB_angle = 110.5
        self.N_C_CA_CB_diangle = -122.5

        self.CB_CG_length = 1.511
        self.CA_CB_CG_angle = 113.8
        self.N_CA_CB_CG_diangle=-64.3

        self.CG_CD1_length = 1.394
        self.CB_CG_CD1_angle = 120.8
        self.CA_CB_CG_CD1_diangle=93.1

        self.CG_CD2_length = 1.394
        self.CB_CG_CD2_angle = 120.8
        self.CD1_CB_CG_CD2_diangle = 180

        self.CD1_CE1_length = 1.392
        self.CG_CD1_CE1_angle = 121.1
        self.CB_CG_CD1_CE1_diangle = 180.0

        self.CD2_CE2_length = 1.392
        self.CG_CD2_CE2_angle = 121.1
        self.CB_CG_CD2_CE2_diangle = 180.0

        self.CE1_CZ_length = 1.385
        self.CD1_CE1_CZ_angle = 119.5
        self.CG_CD1_CE1_CZ_diangle = 0.0

        self.CZ_OH_length = 1.376
        self.CE1_CZ_OH_angle = 119.7
        self.CD1_CE1_CZ_OH_diangle = 180.0

        self.residue_name= 'Y'

    def inputRotamers(self, rotamers):
        self.N_CA_CB_CG_diangle=rotamers[0]
        self.CA_CB_CG_CD1_diangle=rotamers[1]
            
class TrpGeo(Geo):                              
    '''Geometry of Tryptophan'''
    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 110.8914

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.5117
        self.N_CA_C_O_diangle = 120.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.53
        self.C_CA_CB_angle = 110.5
        self.N_C_CA_CB_diangle = -122.5

        self.CB_CG_length = 1.50
        self.CA_CB_CG_angle = 113.8
        self.N_CA_CB_CG_diangle=-66.4

        self.CG_CD1_length = 1.365
        self.CB_CG_CD1_angle = 126.9
        self.CA_CB_CG_CD1_diangle=96.3

        self.CG_CD2_length = 1.433
        self.CB_CG_CD2_angle = 126.7
        self.CD1_CB_CG_CD2_diangle = 180

        self.CD1_NE1_length = 1.375
        self.CG_CD1_NE1_angle = 110.2
        self.CB_CG_CD1_NE1_diangle = 180.0

        self.CD2_CE2_length = 1.413
        self.CG_CD2_CE2_angle = 107.2
        self.CB_CG_CD2_CE2_diangle = 180.0

        self.CD2_CE3_length = 1.4
        self.CG_CD2_CE3_angle = 133.9
        self.CB_CG_CD2_CE3_diangle = 0.0

        self.CE2_CZ2_length = 1.399
        self.CD2_CE2_CZ2_angle = 122.4
        self.CG_CD2_CE2_CZ2_diangle = 180.0

        self.CE3_CZ3_length = 1.392
        self.CD2_CE3_CZ3_angle = 118.7
        self.CG_CD2_CE3_CZ3_diangle=180.0

        self.CZ2_CH2_length = 1.372
        self.CE2_CZ2_CH2_angle = 117.5
        self.CD2_CE2_CZ2_CH2_diangle = 0.0

        self.residue_name= 'W'

    def inputRotamers(self, rotamers):
        self.N_CA_CB_CG_diangle=rotamers[0]
        self.CA_CB_CG_CD1_diangle=rotamers[1]

def geometry(AA):
    '''Generates the geometry of the requested amino acid.
    The amino acid needs to be specified by its single-letter
    code. If an invalid code is specified, the function
    returns the geometry of Glycine.'''
    if(AA=='G'):
        return GlyGeo()
    elif(AA=='A'):
        return AlaGeo()
    elif(AA=='S'):
        return SerGeo()
    elif(AA=='C'):
        return CysGeo()
    elif(AA=='V'):
        return ValGeo()
    elif(AA=='I'):
        return IleGeo()
    elif(AA=='L'):
        return LeuGeo()
    elif(AA=='T'):
        return ThrGeo()
    elif(AA=='R'):
        return ArgGeo()
    elif(AA=='K'):
        return LysGeo()
    elif(AA=='D'):
        return AspGeo()
    elif(AA=='E'):
        return GluGeo()
    elif(AA=='N'):
        return AsnGeo()
    elif(AA=='Q'):
        return GlnGeo()
    elif(AA=='M'):
        return MetGeo()
    elif(AA=='H'):
        return HisGeo()
    elif(AA=='P'):
        return ProGeo()
    elif(AA=='F'):
        return PheGeo()
    elif(AA=='Y'):
        return TyrGeo()
    elif(AA=='W'):
        return TrpGeo()
    else:
        print ("Geometry.geometry() wrong")

