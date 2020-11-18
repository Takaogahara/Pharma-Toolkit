from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt

class GenImage():

    def getSingle(smiles):
        mol = Chem.MolFromSmiles(smiles)

        fig = Draw.MolToMPL(mol)
        plt.axis('on')
        fig.savefig('./core/images/singleMol.jpeg', bbox_inches='tight')
        # plt.show()

        img = plt.imread('./core/images/singleMol.jpeg')
        return img