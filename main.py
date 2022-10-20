"""
@author Zack Meyer
"""

import time

import openbabel
import pybel
import glob
import re
import LigandExtract
from tqdm import tqdm

def find_nitrogen_multi_bond(atom):
    num_bonds = 0
    for _ in openbabel.OBAtomBondIter(atom.OBAtom):
        num_bonds += 1

    if num_bonds < 3:
        return True
    return False


class MonoFinder:
    def __init__(self, mol):
        self.multi_bond = False
        self.metal_hit = 0
        # comment line tells us the number of bonds on the metal
        self.num_metal_bonds = int(re.findall("MND = (\d+)", str(mol))[0])
        self.mol = mol

        self.find_nearest_atoms()
        self.set_bond_id()

    def find_ligands(self, atom):
        a = atom.OBAtom
        # print("Looking at bonds for atom with atomic number: ", a.GetAtomicNum())
        # print("Atom Index num: ", a.GetIndex()+1)
        for bond in openbabel.OBAtomBondIter(a):
            bond_str = str(a.GetIndex()+1) + "-" + str(self.mol.atoms[bond.GetNbrAtomIdx(a) - 1].idx)
            # seeing if bond has been visited
            if bond.GetId() == 0:
                # print(f"Bond {bond_str} not visited, investigating...")
                # marking visited bond
                bond.SetId(1)

                # if the atom on the other side of the bond is H we can ignore it
                if self.mol.atoms[bond.GetNbrAtomIdx(a) - 1].atomicnum == 1:
                    # print("Skipping Hydrogen ", bond_str)
                    continue
                # if we're looking at the bond to the metal we ignore it
                elif self.mol.atoms[bond.GetNbrAtomIdx(a) - 1].OBAtom.IsMetal():
                    # print("Metal bond identified on bond ", bond_str)
                    self.metal_hit += 1
                    continue
                # otherwise we recursively call find ligands on the atom
                else:
                    # print(f'Traversing across bond {bond_str}')
                    self.find_ligands(self.mol.atoms[bond.GetNbrAtomIdx(a) - 1])
            else:
                # print(f"Bond {bond_str} already visited, continuing...")
                pass

    def start(self, metal):
        bond_iter = 0
        total_metal_hits = 0
        m = metal.OBAtom
        # Loop on each bond attached to the our metal atom
        for bond in openbabel.OBAtomBondIter(m):
            # checking to see if the bond has been visited before
            print("Bond iter: ", bond_iter, " with bond ID ", bond.GetId())
            bond_str = str(metal.idx) + "-" + str(self.mol.atoms[bond.GetNbrAtomIdx(m) - 1].idx)
            print(f'Looking at bond {bond_str}')
            if bond.GetId() == 0:
                # marking the visited bond

                bond.SetId(1)
                # bond.GetNmbAtomIdx takes in a specified atom and returns its neighbor
                # in this case we pass in our metal base and we get back the 
                # atom on the other side which serves as a start point to the ligand
                ligand_start = self.mol.atoms[bond.GetNbrAtomIdx(m) - 1]
                

                # Nitrogen bond special case
                if ligand_start.atomicnum == 7:
                    if find_nitrogen_multi_bond(ligand_start):
                        self.metal_hit += 1  # If the start of the ligand is a nitrogen with a double or triple bond
                        # the metal, then add one to the metal count to ignore this ligand

                if not do_ligand_filter:
                    # identify the rest of the ligand, using our identified start point
                    self.find_ligands(ligand_start)
                elif do_ligand_filter and ligand_start.atomicnum == ligand_filter:
                    self.find_ligands(ligand_start)
                else:
                    # skipping a filtered ligand
                    bond_iter += 1
                    continue

                # TODO: add logic for when our metal hit count is more than 0 (multidentate)
                print("Metal hit count after searching ligands: ", self.metal_hit)
                total_metal_hits += self.metal_hit
                if self.metal_hit == 0:
                    copy_molecule = pybel.Molecule(openbabel.OBMol(self.mol.OBMol))
                    sub = LigandExtract.LigandExtractor(copy_molecule, bond_iter, mol_num)
                    sub.extract_ligand()
                    # return  # This causes the first iteration to be the only iteration.  If we want to do all the monodentates, then remove this.
                else:
                    print("Attempting to extract polydentate!")
                    copy_molecule = pybel.Molecule(openbabel.OBMol(self.mol.OBMol))
                    sub = LigandExtract.LigandExtractor(copy_molecule, bond_iter, mol_num)
                    sub.extract_ligand()
                    self.metal_hit = 0

            bond_iter += 1
        print("Total metal hits: ", total_metal_hits)

    def find_nearest_atoms(self):
        metal_bonded_list = []
        for atom in self.mol:
            if atom.OBAtom.IsMetal():
                for otherAtom in self.mol:
                    if not otherAtom.OBAtom.IsMetal():
                        metal_bonded_list.append((otherAtom, atom.OBAtom.GetDistance(otherAtom.OBAtom)))

        new_metal_bonded_list = sorted(metal_bonded_list, key=lambda x: x[1])

        for atom in self.mol:
            if atom.OBAtom.IsMetal():
                for j in range(self.num_metal_bonds):
                    if self.mol.OBMol.GetBond(atom.OBAtom, new_metal_bonded_list[j][0].OBAtom) is None:
                        if new_metal_bonded_list[j][0].atomicnum == 1:
                            i = 0
                            for _ in openbabel.OBAtomBondIter(new_metal_bonded_list[j][0].OBAtom):
                                i += 1
                            if i > 0:
                                break
                            else:
                                bond = openbabel.OBBond()
                                bond.SetBegin(atom.OBAtom)
                                bond.SetEnd(new_metal_bonded_list[j][0].OBAtom)
                                self.mol.OBMol.AddBond(bond)
                        else:
                            bond = openbabel.OBBond()
                            bond.SetBegin(atom.OBAtom)
                            bond.SetEnd(new_metal_bonded_list[j][0].OBAtom)
                            self.mol.OBMol.AddBond(bond)

    def set_bond_id(self):
        for atom in self.mol:
            a = atom.OBAtom
            for bond in openbabel.OBAtomBondIter(a):
                bond.SetId(0)


filter_input = input("Do you want to only replace ligands for a certain metal? ").lower()
do_metal_filter = filter_input == 'y' or filter_input == 'yes'
if do_metal_filter:
    metal_filter = int(input("Please enter the atomic number of the metal you want: "))

filter_input = input("Do you want to only replace ligands with a certain starting atom? ").lower()
do_ligand_filter = filter_input == 'y' or filter_input == 'yes'
if do_ligand_filter:
    ligand_filter = int(input("Please enter the atomic number of the ligand's starting atom: "))

start_time = time.time()
mol_num = 0
for file in glob.glob("*.xyz"):
    print("Entering file")
    for molecule in pybel.readfile("xyz", file):
        print("looking at molecule...")
        mol_num += 1
        finder = MonoFinder(molecule)
        print("Identified metal bonds = ", finder.num_metal_bonds)
        # look for metal atom to start search
        for metal_atom in molecule:
            if metal_atom.OBAtom.IsMetal() and not do_metal_filter:
                finder.start(metal_atom)
            elif metal_atom.OBAtom.IsMetal() and metal_atom.atomicnum == metal_filter:
                finder.start(metal_atom)
            else:
                break
end_time = time.time()
elapsed_time = end_time - start_time
minutes = elapsed_time // 60
seconds = elapsed_time % 60
print(f"Ran through {mol_num} molecules in {minutes:.0f} minutes and {seconds:.2f} seconds.")
