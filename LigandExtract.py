"""
@author Zack Meyer
"""

import os
from openbabel import pybel
from openbabel import openbabel

import HeaderReplacer
import LigandChargeFinder
import ChangeCharge


dentate_map = {1:"Monodentate", 2:"Bidentate", 3:"Tridentate"}
atomic_nums_to_elem = {1:"H", 2:"He", 3:"Li", 4:"Be", 5:"B", 6:"C", 7:"N", 8:"O", 9:"F", 11:"Na",
                       12:"Mg", 13:"Al", 14:"Si", 15:"P", 16:"S", 17:"Cl", 19:"K", 20:"Ca", 21:"Sc",  
                       22:"Ti", 23:"V", 24:"Cr", 25:"Mn", 26:"Fe", 27:"Co", 28:"Ni", 29:"Cu", 
                       30:"Zi", 31:"Ga", 32:"Ge", 33:"As", 34:"Se", 35:"Br", 37:"Rb", 38:"Sr", 
                       39:"Y", 40:"Zr", 41:"Nb", 42:"Mo", 43:"Tc", 44:"Ru", 45:"Rh", 46:"Pd", 
                       47:"Ag", 48:"Cd", 49:"In", 50:"Sn", 51:"Sb", 52:"Te", 53:"I", 55:"Cs", 56:"Ba",
                       57:"La", 72:"Hf", 73:"Ta", 74:"W", 75:"Re", 76:"Os", 77:"Ir", 78:"Pt", 79:"Au",
                       80:"Hg", 104:"Rf", 105:"Db", 106:"Sg", 107:"Bh", 108:"Hs"}

class LigandExtractor:
    def __init__(self, mol, bond_iter, num_atom, denticity, center_metal, ligand_num):
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
        self.denticity = denticity
        self.atom_connections = []
        self.connecting_indices = []
        self.center_metal = center_metal
        self.ligand_num = ligand_num

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

    def convert_connections(self):
        self.atom_connections.sort()
        output = ""
        for i in range(len(self.atom_connections)):
            output += atomic_nums_to_elem[self.atom_connections[i]] 
        return output

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
            if bond.GetId() == 0:
                bond.SetId(1)

                if self.mol.atoms[bond.GetNbrAtomIdx(a) - 1].OBAtom.IsMetal():
                    # if atom is connected to metal AND has not been marked as a metal connection yet
                    if a.GetId() < 3:
                        a.SetId(3)
                        self.atom_connections.append(a.GetAtomicNum())
                    
                    continue
                else:
                    self.find_ligands(self.mol.atoms[bond.GetNbrAtomIdx(a) - 1])
            else:
                pass
            
    def getCSDid(self, file_path):
        with open(file_path, 'r') as file:
            content = file.read()
            index = content.find('CSD_code = ') + 11 
            csd_id = content[index:index+6]
            return csd_id

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
                        self.atom_connections.append(self.mol.atoms[bond.GetNbrAtomIdx(atom.OBAtom) - 1].atomicnum)

                        # break
                    
                        # Change the atom to a carbon
                        # self.mol.atoms[bond.GetNbrAtomIdx(atom.OBAtom) - 1].OBAtom.SetAtomicNum(6)
                        #
                        # # Set the length to 2.1 so it's close enough to be counted as bonded to the metal
                        # bond.SetLength(2.1)

                    counter += 1

        for atom in self.mol:
            if atom.OBAtom.GetId() < 2:
                self.mol.OBMol.DeleteAtom(atom.OBAtom)

        # This is after the ligand gets deleted because when they get deleted the indexes of the atoms change

        atom_contact = 0
        for atom in self.mol:
            if atom.OBAtom.GetId() == 3:
                # self.start_index = atom.OBAtom.GetIndex() + 1
                self.connecting_indices.append(atom.OBAtom.GetIndex() + 1)

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

        connecting_atoms = self.convert_connections()

        if not os.path.exists('ExtractedLigand'):
            os.makedirs('ExtractedLigand')

        if not self.denticity in dentate_map.keys():
            self.dentate_name = "TetrasAndGreater"
        else:
            self.dentate_name = dentate_map[self.denticity]

        if not os.path.exists(f'ExtractedLigand/{self.dentate_name}'):
            os.makedirs(f'ExtractedLigand/{self.dentate_name}')
        
        if not os.path.exists(f'ExtractedLigand/{self.dentate_name}/{connecting_atoms}'):
            os.makedirs(f'ExtractedLigand/{self.dentate_name}/{connecting_atoms}')


        # self.mol.write("xyz", f"ExtractedLigand/{self.dentate_name}/{connecting_atoms}/{self.type}{self.num_atom}-{self.bond_iter}.xyz", True)
        # charge_changer = ChangeCharge.ChargeChanger(f"ExtractedLigand/{self.dentate_name}/{connecting_atoms}/{self.type}{self.num_atom}-{self.bond_iter}.xyz")
        # charge_changer.change(self.new_charge)

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

        print("Writing ligand to file type mol...")
        mol_file_path = f"ExtractedLigand/{self.dentate_name}/{connecting_atoms}/{self.type}{self.num_atom}-{self.bond_iter}.mol2"

        self.mol.write("mol2", mol_file_path, True)
        print("Wrote successfully")

        # rename file using new format: 
        # get CSDid
        csd_id = self.getCSDid(mol_file_path)
        # get center metal
        center_metal_symbol = atomic_nums_to_elem[self.center_metal]
        # get ligand num as str() from self
        ligand_num = str(self.ligand_num)
        # get coordination num from len() of connecting_atoms
        coordination_num = str(len(self.atom_connections))
        # coordination atoms we already have
        new_file_name = f'ExtractedLigand/{self.dentate_name}/{connecting_atoms}/'
        new_file_name += csd_id + '_' + center_metal_symbol + '_' + \
                        ligand_num + '_' + coordination_num + '_'+ connecting_atoms + '.mol2'
        
        print("renaming file to ", new_file_name)
        os.rename(mol_file_path, new_file_name)

        # file_smile = f"UpdatedCharge/smile/{self.type}{self.num_atom}-{self.bond_iter}.smi"
        # self.mol.write("smi", file_smile, True)

        # header_replacer = HeaderReplacer.HeaderReplacer(file_smile)
        # header_replacer.replace_header(self.start_index)
       
        # header_replacer = HeaderReplacer.HeaderReplacer(mol_file_path)
        # header_replacer.replace_header(self.connecting_indices)
