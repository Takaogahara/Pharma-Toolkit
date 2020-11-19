from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import numpy as np
import pandas as pd


class Utils():

    def getDescriptors(smiles:str):
        RO5_violations = 0

        mol = Chem.MolFromSmiles(smiles)

        # ! Descriptors
        desc_LogP = Descriptors.MolLogP(mol)
        desc_MW = Descriptors.MolWt(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        
        desc_NOCount = Lipinski.NOCount(mol)
        desc_NHOHCount = Lipinski.NHOHCount(mol)
        desc_RotatableBonds = Lipinski.NumRotatableBonds(mol)
        desc_HeavyAtomCount = Descriptors.HeavyAtomCount(mol)
        desc_HeavyAtomMW = Descriptors.HeavyAtomMolWt(mol)

        desc_RingCount = Lipinski.RingCount(mol)
        desc_NumAromaticRings = Lipinski.NumAromaticRings(mol)
        desc_AromaticAtoms = Utils.AromaticAtoms(mol)

        # ! Lipinski rule of 5
        if desc_LogP > 5:
            RO5_violations = RO5_violations + 1
        if desc_MW > 500:
            RO5_violations = RO5_violations + 1
        if desc_NumHAcceptors > 10:
            RO5_violations = RO5_violations + 1
        if desc_NumHDonors > 5:
            RO5_violations = RO5_violations + 1

        # ! Export
        values = np.array([desc_LogP,
                            desc_MW,
                            desc_NumHAcceptors,
                            desc_NumHDonors,
                            desc_NOCount,
                            desc_NHOHCount,
                            desc_RotatableBonds,
                            desc_HeavyAtomCount,
                            desc_HeavyAtomMW,
                            desc_RingCount,
                            desc_NumAromaticRings,
                            desc_AromaticAtoms,
                            RO5_violations])

        names = ['LogP',
                    'MW',
                    'Num H Acceptors',
                    'Num H Donors',
                    'NO Count',
                    'NH OH Count',
                    'Rotatable Bonds',
                    'Heavy Atom Count',
                    'Heavy Atom MW',
                    'Ring Count',
                    'Num Aromatic Rings',
                    'Aromatic Atoms',
                    'Num Rule of 5 Violations']

        return values, names

 # ---------------------------------------------------------------------------------

    def AromaticAtoms(mol):
        aromatic_atoms = [mol.GetAtomWithIdx(i).GetIsAromatic() for i in range(mol.GetNumAtoms())]
        AA_count = []
        sum_AA_count = 0

        for elem in aromatic_atoms:
            if elem==True:
                AA_count.append(1)
                sum_AA_count = sum(AA_count)

        return sum_AA_count

 # ---------------------------------------------------------------------------------