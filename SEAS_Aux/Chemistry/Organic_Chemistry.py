
"""

This is the chemistry module


it requires rdkit to run.
For mac computers, please install rdkit via conda install.





"""




import os
import sys

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))


import rdkit
from rdkit import Chem



def molecule_bond(smiles,bond):
    
    m1 = Chem.MolFromSmiles(smiles)
    m2 = Chem.MolFromSmiles(bond)
    
    if m1==m2:
        return False
    return m1.HasSubstructMatch(m2)


def molecule_matching():
    
    m = Chem.MolFromSmiles('OC(C#C)=O')
    patt = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")

    print m.HasSubstructMatch(patt)
    print m.GetSubstructMatch(patt)
    
    
    n = Chem.MolFromSmiles('C=C')
    patt1 = Chem.MolFromSmarts("[$([CX2](=C)=C)]")
    patt2 = Chem.MolFromSmarts("[$([CX3]=[CX3])]")
    
    print n.HasSubstructMatch(patt1),n.HasSubstructMatch(patt2)


def feature_profile(feature):
    """
    Each feature have the following Profile:
        Amplitude: Strong (90%), Medium (60%), Weak (30%)
        Shape: Broad (round top), Narrow (sharp peak)
        Span: 
    
    """
    
    shape = feature["shape"]
    amplitude = feature["amplitude"]
    
    
    
    profile = []
    
    
    
    
    return profile
    
    

def predict_spectroscopy():
    """
    Figure out a way to predict the spectra feature based on existing bonds in the molecule
    """
    
    has_OH = True
    
    feature_band = []
    
    if has_OH:
        band = [3000,3500]
        feature_band.append(band)
        
        
        
    return feature_band
    
    
def check_CH_bonds():
    
    
    m = Chem.MolFromSmiles("CC=C")
    patt = Chem.MolFromSmarts("H")
    
    print m.HasSubstructMatch(patt)
    

if __name__ == "__main__":
    
    patt = "[$([CX2](=C)=C)]"
    
    molecule_matching()





