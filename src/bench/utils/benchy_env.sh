#!/bin/bash


# Create a base environment for Python dependencies and spacer-containment
micromamba clean --all # remove cache etc
micromamba create -n base_env python=3.10 -c conda-forge 
micromamba activate base_env
micromamba install polars -y
micromamba install hyperfine -y
micromamba install pyfastx -y
micromamba install needletail -y
micromamba install matplotlib -y
micromamba install seaborn -y
micromamba install altair -y


wget https://github.com/apcamargo/spacer-containment/releases/download/v1.0.0/spacer-containment-1.0.0-x86_64.tar.gz
tar -xvzf spacer-containment-1.0.0-x86_64.tar.gz
mv spacer-containment-1.0.0-x86_64 $CONDA_PREFIX/bin/spacer-containment
chmod +x $CONDA_PREFIX/bin/spacer-containment

# build the python package with the rust module
# assume we are in the git repo folder
#build-rust 
cd src/rust_simulator
maturin develop --release
cd ../../
#build-python 
# hatch version micro 
hatch build --clean
pip install -e .
