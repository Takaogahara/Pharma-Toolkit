from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt

class GenImage():

    def getSingle(mol):

        fig = Draw.MolToMPL(mol)
        plt.axis('on')
        fig.savefig('./backend/images/singleMol.jpeg', bbox_inches='tight')
        # plt.show()

        img = plt.imread('./backend/images/singleMol.jpeg')
        return img

    def getGrid(mol):

        fig = Draw.MolsToGridImage(mol)
        plt.axis('on')
        fig.savefig('./backend/images/gidMol.jpeg', bbox_inches='tight')
        # plt.show()

        img = plt.imread('./backend/images/gidMol.jpeg')
        return img   