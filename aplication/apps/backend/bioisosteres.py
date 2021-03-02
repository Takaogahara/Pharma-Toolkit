from rdkit import Chem

from backend.images.imageGen import GenImage

class Bioisosteres():

    def getBioisosteresList():

        isosteres = []

        isostere_list = ('[*:1]S(=O)(=O)NC[*:2]',
                                '[*:1]C(C(F)(F)(F))NC[*:2]',
                                '[*:1]C1CC1NC[*:2]',
                                '[*:1]C(F)=CC[*:2]',
                                '[*:1]C1COC1NC[*:2]',
                                '[*:1]C1COC1OC[*:2]',
                                '[*:1]C1=NN=C([*:2])N1',
                                '[*:1]C1=NN=C([*:2])O1',
                                '[*:1]N1N=NC([*:2])=C1',
                                '[*:1]N1N=NC([*:2])=N1',
                                '[*:1]C1=NOC([*:2])=N1')
        
        for smile in isostere_list:
            isosteres.append(Chem.MolFromSmiles(smile))

        return isosteres

    def buildIsostereReaction(start, replacement):
        qps = Chem.AdjustQueryParameters()
        qps.adjustDegree = False
        qps.adjustHeavyDegree = False
        qps.adjustRingCount = False
        qps.aromatizeIfPossible = False
        qps.makeAtomsGeneric = False
        qps.makeBondsGeneric = False
        qps.makeDummiesQueries = True

        start = Chem.AdjustQueryProperties(start,qps)
        replacement = Chem.AdjustQueryProperties(replacement,qps)

        product = AllChem.ChemicalReaction()
        product.AddReactantTemplate(start)
        product.AddProductTemplate(replacement)
        product.Initialize()

        return product