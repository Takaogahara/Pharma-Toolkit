FROM heroku/miniconda

RUN conda create -c rdkit -n rdkit-env rdkit
RUN conda activate rdkit-env

ADD ./requirements.txt /tmp/requirements.txt
RUN pip install -qr /tmp/requirements.txt
