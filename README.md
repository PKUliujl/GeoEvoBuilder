# GeoEvoBuilder

GeoEvoBuilder is an efficient **zero-shot deep-learning** method that enhances **both protein thermal stability and activity**. This approach combines structure-based de novo sequence design with a protein language model (ESM2) to enable diverse functional protein design.

![Alt text](https://github.com/PKUliujl/GeoEvoBuilder/blob/main/image/flow.png)

## INSTALLATION

We recommend using conda for environment management.

```bash
# Create and activate environment
conda create -n geoevobuilder python=3.11 -y
conda activate geoevobuilder

# Install system dependencies
conda install -c salilab dssp==3.0.0 -y
conda install -c anaconda "libboost=1.73.*" -y

# Install Python dependencies using uv
pip install uv
uv pip install -e .[dev]
```

## MODEL & WEIGHT

Before running their tasks, users must download:

  1. GeoEvoBuilder package via `git clone https://github.com/PKUliujl/GeoEvoBuilder.git`
  2. model weights using the [link](https://disk.pku.edu.cn/link/AABFDFF8FB729743A8A27FEB5855B31EE0) with access code `xx7W`.

## USAGE

```bash
To better use GeoEvoBuilder toolkit for functional protein sequence design, please add some of these parameters

optional arguments:
  -h, --help            show this help message and exit
  --purpose {0,1}       0 is used for functional sequence design. Only activated at 1 when using additional MSA information (--MSA is given). default: 0
  --inputPATH INPUTPATH, -iP INPUTPATH
                        The directory path should contain the pdb file
  --chainID CHAINID     The given chain for design or repacking
  --inputfile INPUTFILE, -i INPUTFILE
                        A pdb file under inputPATH, eg. T1024_A.pdb.
  --ST ST               Sampling temperature (scaling factor) belongs to (0,1). Larger values increase the number of interation steps for convegence. default: 0.1
  --n N                 The number of designed sequences. default: 10
  --SM SM               Starting mode with different sequence initializations (how many residues are randomly masked). default = 100 (%) masked, i.e. random
                        initialization)
  --Fixed FIXED         Those residues are expected to remain unchanged with names following the format either {residuenumber_chainid_resname,residuenumber_chainid_resname},
                        e.g. 1_A_W,10_A_C, or path to a file with one residue per line. Please note that the residuenumber is renumbered begining at 0.
  --outputPATH OUTPUTPATH, -oP OUTPUTPATH
                        the directory path of the outputfile. default: inputPATH
  --verbose {0,1}       Display the intermediate sequences
  --MSA MSA             The distribution of sequence profile from MSA
```

For general sequence design:

```bash
python run_GeoEvoBuilder.py -iP examples/ -i 3mpc_A.pdb --chainID A --SM 100
```

![Alt text](https://github.com/PKUliujl/GeoEvoBuilder/blob/main/image/DHFR_results.png)

or see the design for dihydrofolate reductase (DHFR):

```bash
python run_GeoEvoBuilder.py -iP examples/ -i 3D80_A.pdb --chainID A --SM 50 --Fixed examples/3D80_A_fixed_residues
```

For sequence design with additional [processed MSA information](https://github.com/PKUliujl/GeoEvoBuilder/tree/main/MSA_Processing):

```bash
python run_GeoEvoBuilder.py -iP examples/ -i 3mpc_A.pdb --chainID A --MSA examples/3mpc_A.pt
```

## CITATION

If you find GeoEvoBuilder useful in your research, please cite it:

```bibtex
@article{liu2025geo,
  author = {Jiale Liu  and Hantian You  and Zheng Guo  and Qin Xu  and Changsheng Zhang  and Luhua Lai },
  title = {GeoEvoBuilder: A deep learning framework for efficient functional and thermostable protein design},
  journal = {Proceedings of the National Academy of Sciences},
  volume = {122},
  number = {41},
  pages = {e2504117122},
  year = {2025},
  doi = {10.1073/pnas.2504117122},
  URL = {https://www.pnas.org/doi/abs/10.1073/pnas.2504117122}
}
```
