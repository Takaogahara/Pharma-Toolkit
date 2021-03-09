import streamlit as st

def app():
    st.markdown("""## Home""")
    with st.beta_container():
        st.markdown("""### About """)

    with st.beta_expander('Dataset Filter'):
        st.markdown("""
        ### Dataset Filter

        Placeholder
        """)


    with st.beta_expander('Molecular Fingerprint Calculator'):
        st.markdown("""
        ### Molecular Fingerprint Calculator

        This module allows the calculation of **molecular fingerprints** used in computational drug discovery projects such as for the construction of quantitative structure-activity/property relationship (QSAR/QSPR) models.  

        There are 12 **molecular fingerprints** avaliable:
        """)

        col1_fp, col2_fp, col3_fp = st.beta_columns(3)

        col1_fp.markdown("""
        - `AtomPairs2D`
        - `AtomPairs2DCount`
        - `CDK`
        - `CDKextended`
        """)

        col2_fp.markdown("""
        - `CDKgraphonly`
        - `EState`
        - `KlekotaRoth`
        - `KlekotaRothCount`
        """)

        col3_fp.markdown("""
        - `MACCS`
        - `PubChem`
        - `Substructure`
        - `SubstructureCount`
        """)

        st.markdown("""
        **Credits**  
        - App built inspired in [Chanin Nantasenamat](https://medium.com/@chanin.nantasenamat) (aka [Data Professor](http://youtube.com/dataprofessor)) application
        - Descriptor calculated using [PaDEL-Descriptor](http://www.yapcwsoft.com/dd/padeldescriptor/) software.  
        Yap CW. [PaDEL‐descriptor: An open source software to calculate molecular descriptors and fingerprints](https://doi.org/10.1002/jcc.21707). ***J Comput Chem*** 32 (2011) 1466-1474.
        
        ---
        """)