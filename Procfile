wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.2-Linux-x86_64.sh
chmod +x Miniconda3-py37_4.8.2-Linux-x86_64.sh
bash ./Miniconda3-py37_4.8.2-Linux-x86_64.sh -b -f -p /usr/local
conda install -c rdkit rdkit -y

web: sh setup.sh && streamlit run frontend/run-app.py
