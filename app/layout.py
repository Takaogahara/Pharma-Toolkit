import streamlit as st

from backend.descriptors import showDescriptors
from backend.utils import Utils

class Run:

    def SingleSmiles(smiles):
        st.write("""### **Descriptors**""")
        col_desc = st.beta_columns(2)

        st.write("""### **Solubility**""")
        col_sol = st.beta_columns(2)

        st.write("""### **Bioisosteres**""")
        col_bioiso = st.beta_columns(2)

        if smiles:

            try:
                descriptors, img, pred_LogS = showDescriptors.getSingle(smiles)

                with st.beta_container():
                    col_desc[0].dataframe(data=descriptors, height=1500)
                    col_desc[1].image(img, caption=smiles, width=400)
                
                with st.beta_container():
                    col_sol[0].write(pred_LogS)
                    col_sol[1].write('Implemented from Delaney JS. 2004 J. Chem. Inf. Model')
                    
                with st.beta_container():
                    col_bioiso[0].write('Bioisosteres')
                    # isostereReactions =[buildIsostereReaction(amide,x) for x in isosteres]
                    # image = Utils.getBioisosteres(smiles)

            except:
                col_desc[0].write('ERROR: Check input smiles')
                col_desc[1].write('ERROR: Check input smiles')

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