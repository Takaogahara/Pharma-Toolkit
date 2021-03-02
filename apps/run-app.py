from layout import Run
import streamlit as st

# * Header
st.write("""
# Pharma Toolkit
test version
""")

# * Sidebar
st.sidebar.header('User Input Parameters')
option = st.sidebar.selectbox('Type',
                    ('Single smile', 'Multiple smile'))

# * Page
if option == 'Single smile':
    st.sidebar.subheader('Single smiles')
    single_smiles = str(st.sidebar.text_input('Enter smiles', 'CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O'))
    Run.SingleSmiles(single_smiles)

elif option == 'Multiple smile':
    st.sidebar.subheader('Multiple smiles')
    single_smiles = str(st.sidebar.text_area('Enter smiles', 'C1CCCC1\nCC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O\nO=C(C)Oc1ccccc1C(=O)O'))
    Run.MultipleSmiles(single_smiles)