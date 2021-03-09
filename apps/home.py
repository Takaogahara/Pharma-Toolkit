import streamlit as st

def app():
    st.markdown("""## Home""")
    with st.beta_container():
        st.markdown("""### About""")

        st.markdown("""

        Toolbox created with the objective of facilitating the work with the tools and operations commonly used in Drug Discovery and making the process more productive.
        Developed in Python language allowing a high integration with scripts and API.

        The Toolkit is divided into modules with different functionalities. Currently there are 2 modules: `Dataset Filter` and` Molecular Fingerprint Calculator`.
        New modules with extra features will be added over time.

        Read the description below for each module to understand how it works.
        """)

        st.markdown("""####""")
        with st.beta_expander('About the author'):
            st.markdown("""
            Placeholder # TODO - Provide text
            """)

        st.markdown("""####""")
        with st.beta_expander('Credits'):
            st.markdown("""
            **Credits**  
            - App built inspired in [Chanin Nantasenamat](https://medium.com/@chanin.nantasenamat) (aka [Data Professor](http://youtube.com/dataprofessor)) applications
            - Molecule handling done using [RDKit: Open-source cheminformatics.](https://www.rdkit.org/docs/index.html).
            - Descriptor calculated using [PaDEL-Descriptor](http://www.yapcwsoft.com/dd/padeldescriptor/) software.  
            Yap CW. [PaDEL‐descriptor: An open source software to calculate molecular descriptors and fingerprints](https://doi.org/10.1002/jcc.21707). ***J Comput Chem*** 32 (2011) 1466-1474.
            """)

    st.markdown("""##""")
    with st.beta_expander('Dataset Filter'):
        st.markdown("""
        ### Dataset Filter

        This module allows several operations related to filtering and processing the dataset.

        These operations are:
        """)

        col1_df, col2_df = st.beta_columns(2)

        col1_df.markdown("""
        - Removal of NaN values
        - Filtering of molecules by element
        - Filtering by target organism
        - Convert standard unit and determine activity
        """)

        col2_df.markdown("""
        - Removal of duplicate entries
        - Removal of dataset columns
        - Shuffle dataset
        """)

        st.markdown("""
        Remarks:

        - Some operations require a specific order of column selection.
        - Currently filtering of molecules by element filters molecules containing **ONLY**: `C`, `O`, `N`, `S`, `P`, `F`, `I`, `Br`, `Cl`.
        - Converting standard unit and determining activity works only with the **INPUT** values: `ug.mL-1`,` uM`, `nM`.
        - Converting standard unit and determining activity works only with the **OUTPUT** values: ` uM`, `nM`.

        """)
        

    st.markdown("""##""")
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