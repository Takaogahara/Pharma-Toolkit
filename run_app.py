import streamlit as st
from apps import home, calc_fingerprints

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
# app.add_app('Dataset Validation', home.app)
app.add_app('Calculate fingerprints', calc_fingerprints.app)
# app.add_app('Calculate descriptors', home.app)
# app.add_app('Regression / Classification', home.app)
# app.add_app('SMILES to image', home.app)

# The main app
app.run()