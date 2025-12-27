#!/bin/bash

set -euo pipefail

conda create -n geoevobuilder python=3.11 -y
conda activate geoevobuilder

# conda install -c conda-forge dssp -y
conda install -c salilab dssp==3.0.0 -y
conda install -c anaconda "libboost=1.73.*" -y

pip install uv

uv pip install -e .
