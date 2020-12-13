from rdkit import Chem
from rdkit.Chem import Descriptors
import numpy as np
import pandas as pd

from backend.bioisosteres import Bioisosteres
from backend.images.imageGen import GenImage
from backend.utils import Utils

class showDescriptors():

    def getSingle(smiles:str):
        mol = Chem.MolFromSmiles(smiles)

        # * Get IMG
        image = GenImage.getSingle(mol)

        # * Get descriptors
        data, dataName = Utils.getDescriptors(mol)
        descriptors = pd.DataFrame(data=data, columns=['Value'], index=dataName)

        # * Get prediction
        _, pred_LogS = Utils.getPredLogS(smiles, data)

        # # * Get Bioisosteres
        # isostereReactions = []
        # isosteres = Bioisosteres.getBioisosteresList():
        # for iso in isosteres:
        #     isostere = Bioisosteres.buildIsostereReaction(mol, iso)
        #     isostereReactions.append(isostere)


        # isostereReactions =[ for x in isosteres]
        # # image = Utils.getBioisosteres(mol)


        return descriptors, image, pred_LogS

 # ---------------------------------------------------------------------------------

    def getMultiple(smiles):
        mol_name = []
        desc_values = []

        smiles = smiles.split("\n")

        for elem in smiles:
            mol_name.append(elem)

            data, dataName = Utils.getDescriptors(elem)
            desc_values.append(data)

        desc_values = np.reshape(np.array(desc_values), (-1,len(dataName)))

        descriptors = pd.DataFrame(data=desc_values, columns=dataName, index=mol_name)

        return descriptors