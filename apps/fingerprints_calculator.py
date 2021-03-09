import streamlit as st
import pandas as pd
import subprocess
import os
import base64

def app():
    def _fileDownload(df, dname, text):
        '''Provide dataframe for download'''

        csv = df.to_csv(index=False)
        b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
        href = f'<a href="data:file/csv;base64,{b64}" download="{dname}.csv">{text}</a>'
        return href

    def _runPaDEL():
        '''Performs the descriptor calculation'''

        bashCommand = f'java -Xms2G -Xmx2G -Djava.awt.headless=true -jar ./apps/backend/PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./apps/backend/PaDEL-Descriptor/{str(selected_fp)} -dir ./apps/backend/PaDEL-Descriptor/molecule.smi -file ./apps/backend/PaDEL-Descriptor/descriptors_output.csv'
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        os.remove('./apps/backend/PaDEL-Descriptor/molecule.smi')
        # ./apps/backend/PaDEL-Descriptor/

        # Read in calculated descriptors 
        st.subheader('Calculated molecular descriptors')
        desc = pd.read_csv('./apps/backend/PaDEL-Descriptor/descriptors_output.csv')
        nmol = desc.shape[0]
        ndesc = desc.shape[1]
        
        # Write the data dimension and display the dataframe
        st.success(f'Selected fingerprint: {user_fp}')
        st.success(f'Number of molecules: {str(nmol)}')
        st.success(f'Number of descriptors: {str(ndesc-1)}')
        st.write(desc)
        st.markdown(_fileDownload(desc, f'descriptor_{user_fp}', 'Download CSV File'), unsafe_allow_html=True)

    # ----------------------------------------------------------------------------------------------------------
    # Page title
    st.markdown("""## Molecular Fingerprint Calculator""")
    execute = False

    # ----------------------------------------------------------------------------------------------------------
    # Sidebar
    
    # CSV file uploader
    with st.sidebar.header('1. Upload your CSV data'):
        uploaded_file = st.sidebar.file_uploader("Upload your input CSV file", type=["csv"])
        exemple_file = pd.read_csv('./apps/backend/PaDEL-Descriptor/example_CSV.csv')
        st.sidebar.markdown(_fileDownload(exemple_file, 'example_csv', 'Download example CSV file'), unsafe_allow_html=True)
    
    if uploaded_file is not None:

        

        # Select column names and delimiter
        with st.sidebar.header('2. Enter column names for **Molecule ID** and **SMILES** and privide the CSV **delimiter**'):

            delimiter_dict = {',':',', ';':';'}
            user_delimiter = st.sidebar.selectbox('Choose CSV file delimiter', list(delimiter_dict.keys()))
            selected_delimiter = delimiter_dict[user_delimiter]

            raw_df = pd.read_csv(uploaded_file, delimiter=selected_delimiter)

            index_columns = list(raw_df.columns)
            moleculeColumn = st.sidebar.selectbox('Molecule ID' , index_columns)
            smilesColumn = st.sidebar.selectbox('SMILES' , index_columns)

        # Select fingerprint type and DF size
        with st.sidebar.header('3. Set parameters'):
            fp_dict = {'AtomPairs2D':'AtomPairs2DFingerprinter.xml',
                        'AtomPairs2DCount':'AtomPairs2DFingerprintCount.xml',
                        'CDK':'Fingerprinter.xml',
                        'CDKextended':'ExtendedFingerprinter.xml',
                        'CDKgraphonly':'GraphOnlyFingerprinter.xml',
                        'EState':'EStateFingerprinter.xml',
                        'KlekotaRoth':'KlekotaRothFingerprinter.xml',
                        'KlekotaRothCount':'KlekotaRothFingerprintCount.xml',
                        'MACCS':'MACCSFingerprinter.xml',
                        'PubChem':'PubchemFingerprinter.xml',
                        'Substructure':'SubstructureFingerprinter.xml',
                        'SubstructureCount':'SubstructureFingerprintCount.xml'}
            user_fp = st.sidebar.selectbox('Choose fingerprint to calculate', list(fp_dict.keys()))
            selected_fp = fp_dict[user_fp]

            # Select number of molecules to compute
            number2calc = st.sidebar.slider('How many molecules to compute?', min_value=10, max_value=raw_df.shape[0], value=10, step=10)

            # Run fingerprint
            with st.sidebar.header('4. Run fingerprints'):
                if st.sidebar.button('Calculate'):
                    execute = not execute

    # ----------------------------------------------------------------------------------------------------------

    if uploaded_file is not None:
        # Read CSV data and display dataframe
        df = raw_df.iloc[:number2calc,1:]
        st.subheader('Initial data from CSV file')
        st.write(df)

        # Select columns to make SMI file and display dataframe
        try:
            fingerprint_df = pd.concat([df[smilesColumn], df[moleculeColumn]], axis=1)
            st.subheader('Formatted as PADEL input file')
            st.write(fingerprint_df)
        except:
            st.error('Please provide a valid column name.')

        if execute == True:
            execute == False
            # Write SMI data
            fingerprint_df.to_csv('./apps/backend/PaDEL-Descriptor/molecule.smi', sep = '\t', header = False, index = False)
            # Run PaDEL-Descriptor
            with st.spinner('Calculating descriptors...'):
                _runPaDEL()

    else:
        st.info('Awaiting for CSV file to be uploaded.')