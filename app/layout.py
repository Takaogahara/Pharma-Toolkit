import streamlit as st

from backend.descriptors import showDescriptors
from backend.utils import Utils

class Run:

    def SingleSmiles(smiles):
        st.write("""### Single smiles""")

        column1, column2 = st.beta_columns(2)

        if smiles:
            descriptors, img = showDescriptors.getSingle(smiles)
            try:
                pass
            except:
                column1.write('ERROR: Check input smiles')
                column2.write('ERROR: Check input smiles')

            
            column1.dataframe(data=descriptors, height=1500)
            column2.image(img, caption=smiles, width=400)

 #----------------------------------------------------------------------------

    def MultipleSmiles(smiles):
        st.write("""
        ### Multiple smiles
        """)

        if smiles:
            try:
                descriptors = showDescriptors.getMultiple(smiles)
            except:
                st.write('ERROR: Check input smiles')

            with st.beta_container():
                st.dataframe(data=descriptors, width=2000, height=2000)