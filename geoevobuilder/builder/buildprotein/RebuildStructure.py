##################################################################################
# This code is a part from OPUS-Rota4, https://doi.org/10.1093/bib/bbab529. ###
# If you find use, please cite it.                                          ###
##################################################################################


from geoevobuilder.builder.buildprotein import Geometry, PeptideBuilder


def extractmc(atomsData_real, seq_tobe_designed=None):
    if seq_tobe_designed:
        """
        f = open(seq_tobe_designed,'r')
        rows = f.readlines()
        seq = rows[1].strip()
        f.close()
        """
        seq = seq_tobe_designed
    # main_chain_only
    atomsData_main_chain_only = []
    if seq_tobe_designed:
        atoms_n = 0
        for atom in atomsData_real:
            if atom.name1 in ["N", "CA", "C", "O"]:
                atom.resname = seq[int(atoms_n / 4)]
                atomsData_main_chain_only.append(atom)
                atoms_n += 1
        # print( atoms_n, len(seq))
        assert atoms_n == 4 * len(seq), (
            "Total atoms number of backbone does not match the length of designed sequence."
        )
    else:
        for atom in atomsData_real:
            if atom.name1 in ["N", "CA", "C", "O"]:
                atomsData_main_chain_only.append(atom)

    return atomsData_main_chain_only


def getGeosData(residuesData):

    geosData = []
    for residue in residuesData:
        geo = Geometry.geometry(residue.resname)
        geosData.append(geo)

    return geosData


def rebuild_cb(residuesData, geosData):

    for residue, geo in zip(residuesData, geosData):
        if residue.resname == "G":
            geo = Geometry.geometry("A")
            residue.atoms["CB"] = PeptideBuilder.get_cb(geo, residue)
        else:
            residue.atoms["CB"] = PeptideBuilder.get_cb(geo, residue)

    return residuesData


def rebuild(rotamers, residuesData, geosData, init_atoms_matrix):

    count = 0
    atoms_matrix = []
    atoms_matrix_name = []
    for idx, (residue, geo) in enumerate(zip(residuesData, geosData)):
        num_rotamers = residue.num_dihedrals
        rotamer = rotamers[count : count + num_rotamers]

        atoms, atoms_names = PeptideBuilder.get_coordinate(rotamer, residue, geo, init_atoms_matrix)
        atoms_matrix.extend(atoms)
        atoms_matrix_name.append(atoms_names)
        count += num_rotamers

    assert count == len(rotamers)

    return atoms_matrix, atoms_matrix_name


def rebuild_CBonly(residuesData, geosData, init_atoms_matrix):

    count = 0
    atoms_matrix = []
    atoms_matrix_name = []
    for idx, (residue, geo) in enumerate(zip(residuesData, geosData)):
        # num_rotamers = residue.num_dihedrals
        # rotamer = rotamers[count: count+num_rotamers]

        atoms, atoms_names = PeptideBuilder.get_coordinate_CB(residue, geo, init_atoms_matrix)
        atoms_matrix.extend(atoms)
        atoms_matrix_name.append(atoms_names)
    # count += num_rotamers

    # assert count == len(rotamers)

    return atoms_matrix, atoms_matrix_name


def rebuild_Val(residuesData, geosData, init_atoms_matrix):

    count = 0
    atoms_matrix = []
    atoms_matrix_name = []
    for idx, (residue, geo) in enumerate(zip(residuesData, geosData)):
        num_rotamers = residue.num_dihedrals
        # rotamer = rotamers[count: count+num_rotamers]
        rotamer = [175.0]
        atoms, atoms_names = PeptideBuilder.get_coordinate(rotamer, residue, geo, init_atoms_matrix)
        atoms_matrix.extend(atoms)
        atoms_matrix_name.append(atoms_names)
    # count += num_rotamers

    # assert count == len(rotamers)

    return atoms_matrix, atoms_matrix_name


def rebuild_Ile(residuesData, geosData, init_atoms_matrix):

    count = 0
    atoms_matrix = []
    atoms_matrix_name = []
    for idx, (residue, geo) in enumerate(zip(residuesData, geosData)):
        num_rotamers = residue.num_dihedrals
        # rotamer = rotamers[count: count+num_rotamers]
        rotamer = [(-65.0), 170.0]
        atoms, atoms_names = PeptideBuilder.get_coordinate(rotamer, residue, geo, init_atoms_matrix)
        atoms_matrix.extend(atoms)
        atoms_matrix_name.append(atoms_names)
    # count += num_rotamers

    # assert count == len(rotamers)

    return atoms_matrix, atoms_matrix_name


def make_atoms_matrix(residuesData, atoms_matrix):

    counter = 0
    for idx, residue in enumerate(residuesData):
        atoms_matrix[counter] = residue.atoms["N"].position
        residue.main_chain_atoms_matrixid["N"] = counter
        counter += 1
        atoms_matrix[counter] = residue.atoms["CA"].position
        residue.main_chain_atoms_matrixid["CA"] = counter
        counter += 1
        atoms_matrix[counter] = residue.atoms["C"].position
        residue.main_chain_atoms_matrixid["C"] = counter
        counter += 1
        atoms_matrix[counter] = residue.atoms["O"].position
        residue.main_chain_atoms_matrixid["O"] = counter
        counter += 1

        if residue.resname == "G":
            atoms_matrix[counter] = residue.atoms["CB"].position
            counter += 1
            continue

        atoms_matrix[counter] = residue.atoms["CB"].position
        residue.main_chain_atoms_matrixid["CB"] = counter
        counter += 1

        if residue.resname == "S":
            residue.side_chain_atoms_matrixid["OG"] = counter
            counter += 1
        elif residue.resname == "C":
            residue.side_chain_atoms_matrixid["SG"] = counter
            counter += 1
        elif residue.resname == "V":
            residue.side_chain_atoms_matrixid["CG1"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CG2"] = counter
            counter += 1
        elif residue.resname == "I":
            residue.side_chain_atoms_matrixid["CG1"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CG2"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CD1"] = counter
            counter += 1
        elif residue.resname == "L":
            residue.side_chain_atoms_matrixid["CG"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CD1"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CD2"] = counter
            counter += 1
        elif residue.resname == "T":
            residue.side_chain_atoms_matrixid["OG1"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CG2"] = counter
            counter += 1
        elif residue.resname == "R":
            residue.side_chain_atoms_matrixid["CG"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CD"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["NE"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CZ"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["NH1"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["NH2"] = counter
            counter += 1
        elif residue.resname == "K":
            residue.side_chain_atoms_matrixid["CG"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CD"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CE"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["NZ"] = counter
            counter += 1
        elif residue.resname == "D":
            residue.side_chain_atoms_matrixid["CG"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["OD1"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["OD2"] = counter
            counter += 1
        elif residue.resname == "N":
            residue.side_chain_atoms_matrixid["CG"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["OD1"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["ND2"] = counter
            counter += 1
        elif residue.resname == "E":
            residue.side_chain_atoms_matrixid["CG"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CD"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["OE1"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["OE2"] = counter
            counter += 1
        elif residue.resname == "Q":
            residue.side_chain_atoms_matrixid["CG"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CD"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["OE1"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["NE2"] = counter
            counter += 1
        elif residue.resname == "M":
            residue.side_chain_atoms_matrixid["CG"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["SD"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CE"] = counter
            counter += 1
        elif residue.resname == "H":
            residue.side_chain_atoms_matrixid["CG"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["ND1"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CD2"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CE1"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["NE2"] = counter
            counter += 1
        elif residue.resname == "P":
            residue.side_chain_atoms_matrixid["CG"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CD"] = counter
            counter += 1
        elif residue.resname == "F":
            residue.side_chain_atoms_matrixid["CG"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CD1"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CD2"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CE1"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CE2"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CZ"] = counter
            counter += 1
        elif residue.resname == "Y":
            residue.side_chain_atoms_matrixid["CG"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CD1"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CD2"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CE1"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CE2"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CZ"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["OH"] = counter
            counter += 1
        elif residue.resname == "W":
            residue.side_chain_atoms_matrixid["CG"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CD1"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CD2"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["NE1"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CE2"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CE3"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CZ2"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CZ3"] = counter
            counter += 1
            residue.side_chain_atoms_matrixid["CH2"] = counter
            counter += 1

    assert counter == atoms_matrix.shape[0]

    return atoms_matrix
