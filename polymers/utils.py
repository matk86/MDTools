import openbabel

from pymatgen.io.babel import BabelMolAdaptor


__author__ = "Kiran Mathew"


def get_topology(molecule):
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
    bonds = [x+[tuple(sorted((molecule[x[0]].ff_map,
                              molecule[x[1]].ff_map)))] for x in bonds ]
    angles = [x+[tuple(((molecule[x[0]].ff_map,
                         molecule[x[1]].ff_map,
                         molecule[x[2]].ff_map)))] for x in angles ]
    dihedrals =  [x+[tuple(sorted((molecule[x[0]].ff_map,
                                   molecule[x[1]].ff_map,
                                   molecule[x[2]].ff_map,
                                   molecule[x[3]].ff_map)))] for x in dihedrals ]
    return atoms, bonds, angles, dihedrals