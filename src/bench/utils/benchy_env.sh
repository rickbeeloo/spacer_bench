#!/bin/bash


# Create a base environment for Python dependencies and spacer-containment
mamba clean --all # remove cache etc
mamba create -n base_env python=3.10 -c conda-forge 
mamba activate base_env
mamba install polars -y
mamba install hyperfine -y
mamba install pyfastx -y
mamba install needletail -y
mamba install matplotlib -y
mamba install seaborn -y
mamba install altair -y


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
