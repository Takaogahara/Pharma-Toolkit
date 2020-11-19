from rdkit import Chem
from rdkit.Chem import Descriptors
import numpy as np
import pandas as pd

from backend.images.imageGen import GenImage
from backend.utils import Utils

class showDescriptors():

    def getSingle(smiles:str):

        image = GenImage.getSingle(smiles)
        data, dataName = Utils.getDescriptors(smiles)

        descriptors = pd.DataFrame(data=data, columns=['Value'], index=dataName)

        return descriptors, image

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