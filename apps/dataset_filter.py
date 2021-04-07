import base64
import os

import pandas as pd
import streamlit as st

# from apps.backend.dataset_filter.filters import Filters
from backend.dataset_filter.filters import Filters #! FIXME: DELETE


def _fileDownload(df, dname, text):
    '''Provide dataframe for download'''

    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="{dname}.csv">{text}</a>'
    return href

def app():
    #@ ----------------------------------------------------------------------------------------------------------
    #@ Page title
    st.markdown("""# **Dataset Filter**""")

    selected_df = False

    #@ ----------------------------------------------------------------------------------------------------------
    #@ Sidebar

    #@ --------------------------------------------------------------
    #@ CSV file uploader

    with st.sidebar.header('Upload your CSV data'):
        uploaded_file = st.sidebar.file_uploader("Upload your input CSV file", type=["csv"])
        exemple_file = pd.read_csv('./backend/dataset_filter/example_csv.csv', delimiter=';')
        st.sidebar.markdown(_fileDownload(exemple_file, 'example_csv', 'Download example CSV file'), unsafe_allow_html=True)

    if uploaded_file is None:
        st.info('Awaiting for CSV file to be uploaded.')

    #@ --------------------------------------------------------------
    #@ Sidebar options

    if uploaded_file is not None:

        #@ Delimiter selection
        delimiter_dict = {',':',', ';':';'}
        user_delimiter = st.sidebar.selectbox('Choose CSV file delimiter', list(delimiter_dict.keys()))
        selected_delimiter = delimiter_dict[user_delimiter]

        try:
            #@ Select number of molecules to compute
            raw_df = pd.read_csv(uploaded_file, delimiter=selected_delimiter)
            number2calc = st.sidebar.slider('How many entries to compute?', min_value=10, max_value=raw_df.shape[0], value=raw_df.shape[0], step=10)

            #@ Columns selection
            index = list(raw_df.columns)
            selected_columns = st.sidebar.multiselect('Select columns', index, None)

            #@ Display original dataset
            df = raw_df.iloc[:number2calc]
            df = df[selected_columns]
            st.markdown("""## **Selected data from CSV file**""")
            st.write(df)
            st.markdown(_fileDownload(df, 'selection_csv', 'Download CSV File'), unsafe_allow_html=True)

            #@ Save selections
            if st.sidebar.button('Done!'):
                df.to_csv('./backend/dataset_filter/selected_columns.csv', index=False)

            st.markdown(""" --- """)

        except:
            st.error('Error in provided configurations.')


    #@ ----------------------------------------------------------------------------------------------------------
    try:
        selected_df = pd.read_csv('./backend/dataset_filter/selected_columns.csv')
    except:
        st.info('Awaiting column selection.')

    if (uploaded_file is not None) and (selected_df  is not None):
        st.markdown("""## **Filtering Data**""")

        #@ -------------------------------
        #@ Remove NaN
        st.markdown("""##""")
        with st.beta_expander('Remove NaN'):
            st.markdown("""### Remove NaN""")
            if st.button('Remove NaN'):
                selected_df = pd.read_csv('./backend/dataset_filter/selected_columns.csv')

                selected_df = Filters.removeNaN(selected_df)
                selected_df.to_csv('./backend/dataset_filter/selected_columns.csv', index=False)

        #@ -------------------------------
        #@ Filter elements
        st.markdown("""##""")
        with st.beta_expander('Filter Elements'):
            st.markdown("""### Filter Elements""")
            index_selection = list(selected_df.columns)
            smilesColumn = st.selectbox('SMILES Column', index_selection)

            if st.button('Filter Elements'):
                try:
                    selected_df = pd.read_csv('./backend/dataset_filter/selected_columns.csv')

                    selected_df = Filters.removeElements(selected_df, smilesColumn)
                    selected_df.to_csv('./backend/dataset_filter/selected_columns.csv', index=False)
                except:
                    st.error('Selected column is incompatiple')

        #@ -------------------------------
        #@ Filter organism
        st.markdown("""##""")
        with st.beta_expander('Filter Organism'):
            st.markdown("""### Filter Organism""")
            index_organism = list(selected_df.columns)
            organismColumn = st.selectbox('Organism Column', index_selection)

            unique_organism = selected_df[str(organismColumn)].unique()
            selected_organism = st.multiselect('Select target organisms', unique_organism, None)

            if st.button('Filter Organism'):
                selected_df = pd.read_csv('./backend/dataset_filter/selected_columns.csv')

                selected_df = Filters.removeStrain(selected_df, selected_organism, organismColumn)
                selected_df.to_csv('./backend/dataset_filter/selected_columns.csv', index=False)     

        #@ -------------------------------
        #@ Set threshold
        st.markdown("""##""")
        with st.beta_expander('Select threshold'):
            st.markdown("""### Select threshold""")
            index_threshold = list(selected_df.columns)
            selected_threshold = st.multiselect('Select: Value, Units, Molecular Weight (IN THIS ORDER)', index_threshold, None)

            threshold_unit_dict = {'uM':'uM', 'nM':'nM'}
            user_threshold_unit = st.selectbox('Choose unit for threshold', list(threshold_unit_dict.keys()))
            selected_threshold_unit = threshold_unit_dict[user_threshold_unit]

            try:
                threshold_value = float(st.text_input('Threshold value', '10'))
            except:
                st.error('Please provide a valid value. (Decimal separator is ".")')

            if st.button('Convert values'):
                selected_df = pd.read_csv('./backend/dataset_filter/selected_columns.csv')

                selected_df = Filters.convertThreshold(selected_df, threshold_value, selected_threshold_unit, selected_threshold)
                selected_df.to_csv('./backend/dataset_filter/selected_columns.csv', index=False)
        
        #@ -------------------------------
        #@ Check duplicates
        st.markdown("""##""")
        with st.beta_expander('Check duplicates'):
            st.markdown("""### Check duplicates""")
            index_duplicate = list(selected_df.columns)
            selected_duplicates = st.multiselect('Select: ID, Activity, Value (IN THIS ORDER)', index_threshold, None)

            if st.button('Check duplicates'):
                selected_df = pd.read_csv('./backend/dataset_filter/selected_columns.csv')

                selected_df = Filters.checkDuplicates(selected_df, selected_duplicates)
                selected_df.to_csv('./backend/dataset_filter/selected_columns.csv', index=False)


        #@ -------------------------------
        #@ Drop column
        st.markdown("""##""")
        with st.beta_expander('Drop columns'):
            st.markdown("""### Drop columns""")
            index_remove = list(selected_df.columns)
            selected_remove = st.multiselect('Select columns to removal', index_remove, None)

            if st.button('Drop columns'):
                selected_df = pd.read_csv('./backend/dataset_filter/selected_columns.csv')

                selected_df = Filters.dropColumns(selected_df, selected_remove)
                selected_df.to_csv('./backend/dataset_filter/selected_columns.csv', index=False)

        #@ -------------------------------
        #@ Shuffle dataset
        st.markdown("""##""")
        with st.beta_expander('Shuffle dataset'):
            st.markdown("""### Shuffle dataset""")
            if st.button('Shuffle dataset'):
                selected_df = pd.read_csv('./backend/dataset_filter/selected_columns.csv')

                selected_df = Filters.shuffleRows(selected_df)
                selected_df.to_csv('./backend/dataset_filter/selected_columns.csv', index=False)


        #@ -------------------------------
        #@ Display

        with st.beta_container():

            #@ Check number of entries, NaN values
            st.markdown(""" --- """)
            st.markdown("""### Data""")
            selected_df.to_csv('./backend/dataset_filter/selected_columns.csv', index=False)
            df_num = selected_df.shape[0]
            nan_num = sum(selected_df.isna().sum())
            st.info(f'Number of molecules: {df_num}')
            st.info(f'NaN values: {nan_num}')

            #@ Display selected dataset
            selected_df = pd.read_csv('./backend/dataset_filter/selected_columns.csv')
            st.write(selected_df)
            st.markdown(_fileDownload(selected_df, 'filtered_csv', 'Download CSV File'), unsafe_allow_html=True)

app() #! FIXME: DELETE

#! FIXME: ./backend/ -> ./apps/backend/