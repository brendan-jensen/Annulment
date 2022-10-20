"""
@author Zack Meyer
"""

import os
import pybel
import openbabel

import HeaderReplacer
import LigandChargeFinder
import ChangeCharge


class LigandExtractor:
    def __init__(self, mol, bond_iter, num_atom):
        self.mol = mol
        self.bond_iter = bond_iter
        self.num_atom = num_atom
        self.metal_ind = 0
        self.new_carbon_ind = 0
        self.start_index = 0
        for atom in self.mol:
            atom.OBAtom.SetId(0)
            if atom.OBAtom.IsMetal():
                self.type = atom.type
        self.set_bond_id()
        self.new_charge = self.modify_charge()
        self.set_bond_id()

    def modify_charge(self):
        for atom in self.mol:
            if atom.OBAtom.IsMetal():
                counter = 0
                for bond in openbabel.OBAtomBondIter(atom.OBAtom):
                    if self.bond_iter == counter:
                        ligand_atom = self.mol.atoms[bond.GetNbrAtomIdx(atom.OBAtom) - 1]
                        charge_finder = LigandChargeFinder.LigandChargeFinder(self.mol)
                        new_charge = charge_finder.change_charge(ligand_atom)

                        # This chunk of code sorts the files into folders based on what their charges are - kept here for future
                        # testing needs
                        # if new_charge == -1:
                        #     if not os.path.exists('Charges'):
                        #         os.makedirs('Charges')
                        #     if not os.path.exists('Charges/MinusOne'):
                        #         os.makedirs('Charges/MinusOne')
                        #     self.mol.write("xyz", f"Charges/MinusOne/{self.type}{self.num_atom}-{self.bond_iter}-Before.xyz", True)
                        # elif new_charge == 1:
                        #     if not os.path.exists('Charges'):
                        #         os.makedirs('Charges')
                        #     if not os.path.exists('Charges/PlusOne'):
                        #         os.makedirs('Charges/PlusOne')
                        #     self.mol.write("xyz", f"Charges/PlusOne/{self.type}{self.num_atom}-{self.bond_iter}-Before.xyz",
                        #                    True)
                        # elif new_charge == 0:
                        #     if not os.path.exists('Charges'):
                        #         os.makedirs('Charges')
                        #     if not os.path.exists('Charges/Zero'):
                        #         os.makedirs('Charges/Zero')
                        #     self.mol.write("xyz", f"Charges/Zero/{self.type}{self.num_atom}-{self.bond_iter}-Before.xyz",
                        #                    True)
                        # else:
                        #     if not os.path.exists('Charges'):
                        #         os.makedirs('Charges')
                        #     if not os.path.exists('Charges/Other'):
                        #         os.makedirs('Charges/Other')
                        #     self.mol.write("xyz", f"Charges/Other/{self.type}{self.num_atom}-{self.bond_iter}-Before.xyz",
                        #                    True)
                        return new_charge

                    counter += 1

    def set_bond_id(self):
        for atom in self.mol:
            a = atom.OBAtom
            for bond in openbabel.OBAtomBondIter(a):
                bond.SetId(0)

    def find_ligands(self, atom):
        a = atom.OBAtom
        # marking visited atom
        a.SetId(2)
        # iterate on each bond of the atom
        for bond in openbabel.OBAtomBondIter(a):
            # if bond hasn't been visited mark, and call find_ligands on other side of it
            # Here LigandExtract assumes we are looking at a monodentate, but we need to deal with bonds that will connect again
            # during the process
            bond_str = str(a.GetIndex()+1) + "-" + str(self.mol.atoms[bond.GetNbrAtomIdx(a) - 1].idx)
            # print("LigExtract: ", bond_str)
            if bond.GetId() == 0:
                # print("Unvisited")
                bond.SetId(1)

                if self.mol.atoms[bond.GetNbrAtomIdx(a) - 1].OBAtom.IsMetal():
                    # print("Metal found, skipping...")
                    continue
                else:
                    self.find_ligands(self.mol.atoms[bond.GetNbrAtomIdx(a) - 1])
            else:
                # print("already visited")
                pass
            
                

    def extract_ligand(self):
        print("Enter extract_ligand...")
        # for each atom in the molecule
        for atom in self.mol:
            # looks for metal
            if atom.OBAtom.IsMetal():
                counter = 0
                # iterates over each bond of the metal
                for bond in openbabel.OBAtomBondIter(atom.OBAtom):
                    # marks visited bond
                    bond.SetId(1)
                    # if the specified bond iter equals the counter (basically have we arrived at the ligand we want to extract)
                    if self.bond_iter == counter:
                        print("LigExtract calling find_ligands...")
                        self.find_ligands(self.mol.atoms[bond.GetNbrAtomIdx(atom.OBAtom) - 1])
                        # Sets the first atom in the ligand's ID to one so it doesn't get deleted
                        self.mol.atoms[bond.GetNbrAtomIdx(atom.OBAtom) - 1].OBAtom.SetId(3)

                        # break


                        # Change the atom to a carbon
                        # self.mol.atoms[bond.GetNbrAtomIdx(atom.OBAtom) - 1].OBAtom.SetAtomicNum(6)
                        #
                        # # Set the length to 2.1 so it's close enough to be counted as bonded to the metal
                        # bond.SetLength(2.1)

                    counter += 1

        for atom in self.mol:
            print(atom.OBAtom.GetId())
            if atom.OBAtom.GetId() < 2:
                self.mol.OBMol.DeleteAtom(atom.OBAtom)
                print("Atom removed")
            # if atom.OBAtom.GetId() == 0:
            #     self.mol.OBMol.DeleteAtom(atom.OBAtom)
            #     print("Atom removed")

        # This is after the ligand gets deleted because when they get deleted the indexes of the atoms change

        atom_contact = 0
        for atom in self.mol:
            if atom.OBAtom.GetId() == 3:
                self.start_index = atom.OBAtom.GetIndex() + 1
                atom_contact = atom.atomicnum
                print(atom_contact)

        # This chunk of code sorts the files into folders based on what their charges are - kept here for future
        # testing needs
        # if self.new_charge == -1:
        #     if not os.path.exists('Charges'):
        #         os.makedirs('Charges')
        #     if not os.path.exists('Charges/MinusOne'):
        #         os.makedirs('Charges/MinusOne')
        #     self.mol.write("xyz", f"Charges/MinusOne/{self.type}{self.num_atom}-{self.bond_iter}-After.xyz", True)
        # elif self.new_charge == 1:
        #     if not os.path.exists('Charges'):
        #         os.makedirs('Charges')
        #     if not os.path.exists('Charges/PlusOne'):
        #         os.makedirs('Charges/PlusOne')
        #     self.mol.write("xyz", f"Charges/PlusOne/{self.type}{self.num_atom}-{self.bond_iter}-After.xyz",
        #                    True)
        # elif self.new_charge == 0:
        #     if not os.path.exists('Charges'):
        #         os.makedirs('Charges')
        #     if not os.path.exists('Charges/Zero'):
        #         os.makedirs('Charges/Zero')
        #     self.mol.write("xyz", f"Charges/Zero/{self.type}{self.num_atom}-{self.bond_iter}-After.xyz",
        #                    True)
        # else:
        #     if not os.path.exists('Charges'):
        #         os.makedirs('Charges')
        #     if not os.path.exists('Charges/Other'):
        #         os.makedirs('Charges/Other')
        #     self.mol.write("xyz", f"Charges/Other/{self.type}{self.num_atom}-{self.bond_iter}-After.xyz",
        #                    True)


        print("Writing ligand to file type xyz...")

        if not os.path.exists('ExtractedLigand'):
            os.makedirs('ExtractedLigand')
        if not os.path.exists(f'ExtractedLigand/Elem {atom_contact}'):
            os.makedirs(f'ExtractedLigand/Elem {atom_contact}')
        self.mol.write("xyz", f"ExtractedLigand/Elem {atom_contact}/{self.type}{self.num_atom}-{self.bond_iter}.xyz", True)
        charge_changer = ChangeCharge.ChargeChanger(f"ExtractedLigand/Elem {atom_contact}/{self.type}{self.num_atom}-{self.bond_iter}.xyz")
        charge_changer.change(self.new_charge)

        self.mol.OBMol.SetTotalCharge(self.new_charge)
        # TODO Find out how to figure out which atom has the charge, then give it its formal charge using OBAtom.SetFormalCharge()
        # openbabel.OBChargeModel.ComputeCharges(self.mol.OBMol)
        # new_mol = openbabel.OBMol(self.mol.OBMol)
        # openbabel.OBChargeModel.ComputeCharges(new_mol)

        # if not os.path.exists('UpdatedCharge'):
        #     os.makedirs('UpdatedCharge')
        # if not os.path.exists('UpdatedCharge/smile'):
        #     os.makedirs('UpdatedCharge/smile')
        # if not os.path.exists('UpdatedCharge/mol'):
        #     os.makedirs('UpdatedCharge/mol')

        # file_mol = f"UpdatedCharge/mol/{self.type}{self.num_atom}-{self.bond_iter}.mol"
        # file_smile = f"UpdatedCharge/smile/{self.type}{self.num_atom}-{self.bond_iter}.smi"
        # print("Writing ligand to file type mol...")
        # self.mol.write("mol", file_mol, True)
        # print("Writing ligand to file type smile...")
        # self.mol.write("smi", file_smile, True)
        # print("succesfully wrote to file...")

        # header_replacer = HeaderReplacer.HeaderReplacer(file_smile)
        # header_replacer.replace_header(self.start_index)
       
        # header_replacer = HeaderReplacer.HeaderReplacer(file_mol)
        # header_replacer.replace_header(self.start_index)
