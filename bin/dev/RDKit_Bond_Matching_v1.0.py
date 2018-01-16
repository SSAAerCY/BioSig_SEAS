
"""
@author: Janusz Petkowski

RDkit Bond Matching v1.0. 

1.Allows for checking if only a Chemical Bond is present in a larger molecule - Natural Product (NP).
2.Uses the most general substructure search - ignores valence, hybridization etc. 
3.Compares only identity of two atoms and a bond between them
4.Returns NPs where Chemical Motif was found - results_bonds.txt file
5.Returns how many times Chemical Motif was found in a NPs databse - freq_bonds.txt file

"""

#from rdkit.Chem import AllChem as Chem
#from rdkit import rdBase
#from rdkit import RDConfig

"""
##### Importing required modules #####

"""

from rdkit import Chem
from rdkit.Chem import rdFMCS

from rdkit import Chem
from rdkit.Chem import MCS


"""
##### opening the Chemical Bonds txt file and NP databse txt file #####

"""


def data_reading(filename):
    # need guard against empty
    f = open(filename,"r")
    data = f.read().split("\n")
    return data
bonds_smiles_filename = "C:\Users\Janusz\Dropbox\Personal_projects\MIT Projects\Halocarbons biosignatures\List_of_all_bonds_CONSPhal.txt"
NP_smiles_filename = "C:\Users\Janusz\Dropbox\Personal_projects\MIT Projects\Halocarbons biosignatures\mASM_4.3_SMILES.txt"
    
bonds_smiles = data_reading(bonds_smiles_filename)
db_smiles = data_reading(NP_smiles_filename)
 

"""

##### The Main function in RDKit Bond Matching v1.0 #####

1. Converts bond and NP smiles to RDKit Mol type
2. Checks if the converted Mol type is valid if not ignores it and procceds further
3. Checks if bond (m2) is present as a general substructure in NP (m1)
4. Calculates frequency of each bond in the NP database - how many NPs in the database contain a given bond
5. Prints the final output in results_bonds.txt and freq_bonds.txt files

"""


def rdkit_bond_matching():

    with open("freq_bonds_mASM_4.3_all_bonds.txt", "w") as file_freq:
        with open("results_bonds_mASM_4.3_all_bonds.txt", "w") as out_:
            for mot_smiles in bonds_smiles:
                print mot_smiles
                m2 = Chem.MolFromSmiles(mot_smiles)
                frequency = 0
                for smiles in db_smiles:
                    m1 = Chem.MolFromSmiles(smiles)
                    if m1 == None:
                        continue
                    is_present = m1.HasSubstructMatch(m2)
                    if is_present == True:
                        print >> out_, "1", smiles, mot_smiles
                        frequency += 1
                    else:
                        continue
                        #print >> out_, "0", smiles, mot_smiles
                #Final result: The frequency of each bond in the NP database 
                        
                print >> file_freq, mot_smiles, frequency

           
        return "Done"


if __name__ == "__main__":
    print rdkit_bond_matching()    



























    