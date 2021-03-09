import streamlit as st
from apps import home, fingerprints_calculator, dataset_filter

# ----------------------------------------------------------------

class MultiApp:
    '''Framework for combining multiple streamlit applications.'''
    def __init__(self):
        self.apps = []

    def add_app(self, title, func):
        self.apps.append({
                'title': title,
                'function': func})

    def run(self):
        # app = st.sidebar.radio(
        app = st.sidebar.selectbox(
            'Navigation',
            self.apps,
            format_func=lambda app: app['title'])

        app['function']()

# ----------------------------------------------------------------

app = MultiApp()

st.markdown('''# Pharma Toolkit''')

app.add_app('Home', home.app)
app.add_app('Dataset Filter', dataset_filter.app)
app.add_app('Fingerprints Calculator', fingerprints_calculator.app)
# app.add_app('Calculate descriptors', home.app)
# app.add_app('Regression / Classification', home.app)
# app.add_app('SMILES to image', home.app)

app.run()