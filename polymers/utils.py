from pymatgen.io.babel import BabelMolAdaptor

import openbabel


__author__ = "Kiran Mathew"


def get_topology(molecule):
    """
    Return the molecular topology obtained via openabel

    Args:
        molecule (Molecule): pymatgen Molecule object

    Returns:
        atoms (list): list of atoms and their force field mapping
        bonds (list): [[i,j, bond_type], ...] where bond_type is
            a sorted tuple of the force field names of atoms i and j.
        angles (list): [[i,j, k, angle_type], ...] where angle_type is
            a sorted tuple of the force field names of atoms i, j and k.
        dihedrals (list): [[i,j, k, l, dihedral_type], ...] where
            dihedral_type is a sorted tuple of the force field names of
            atoms i, j, k and l.
    """
    bma = BabelMolAdaptor(molecule)
    obmol = bma.openbabel_mol
    #print obmol.NumAtoms(), obmol.NumBonds()
    atoms = [x.GetIdx()-1 for x in openbabel.OBMolAtomIter(obmol)]
    bonds = [[x.GetBeginAtomIdx()-1, x.GetEndAtomIdx()-1] for x in openbabel.OBMolBondIter(obmol)]
    angles = [list(x) for x in openbabel.OBMolAngleIter(obmol)]
    dihedrals = [list(x) for x in openbabel.OBMolTorsionIter(obmol)]
    #print len(atoms), len(bonds), len(angles), len(dihedrals)
    atoms = [tuple([str(molecule[x].specie),
                    molecule[x].ff_map]) for x in atoms ]
    bonds = [x+[tuple((molecule[x[0]].ff_map,
                       molecule[x[1]].ff_map))] for x in bonds]
    angles = [x+[tuple(((molecule[x[0]].ff_map,
                         molecule[x[1]].ff_map,
                         molecule[x[2]].ff_map)))] for x in angles]
    dihedrals =  [x+[tuple((molecule[x[0]].ff_map,
                            molecule[x[1]].ff_map,
                            molecule[x[2]].ff_map,
                            molecule[x[3]].ff_map))] for x in dihedrals]
    return atoms, bonds, angles, dihedrals
