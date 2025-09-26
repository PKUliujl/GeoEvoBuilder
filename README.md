# GeoEvoBuilder
GeoEvoBuilder, an efficient zero-shot learning method, is capable to improve both protein thermal stability and activity. This approach integrates structure-based de novo sequence design with protein language model for diverse functional protein design.

<!-- ![Alt text](https://github.com/PKUliujl/GeoEvoBuilder/blob/main/image/flow.jpg) -->
![Alt text](https://github.com/PKUliujl/GeoEvoBuilder/blob/main/image/flow.png)

### INSTALLATION
<!-- ====================== -->
1. GeoEvoBuilder relies on [ESM2](https://github.com/facebookresearch/esm), please install it using `pip install fair-esm`.
2. GeoEvoBuilder relies on the environment used in [GeoSeqBuilder](https://github.com/PKUliujl/GeoSeqBuilder/).
Alternatively, users can install the environment directly using `conda env create -f py38env.yml`.

### MODEL & WEIGHT
Before running their tasks, users must download: 
  1. GeoEvoBuilder package via `git clone https://github.com/PKUliujl/GeoEvoBuilder.git`
  2. model weights using the [link](https://disk.pku.edu.cn/link/AABFDFF8FB729743A8A27FEB5855B31EE0) with access code `xx7W`.

### USAGE
<!-- ====================== 
We are currently in the process of organizing the code and anticipate releasing the software soon. -->
```
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
```
python run_GeoEvoBuilder.py -iP examples/ -i 3mpc_A.pdb --chainID A --SM 100
```
or see the design for dihydrofolate reductase (DHFR):
```
python run_GeoEvoBuilder.py -iP examples/ -i 3D80_A.pdb --chainID A --SM 50 --Fixed examples/3D80_A_fixed_residues
```

For sequence design with additional [processed MSA information](https://github.com/PKUliujl/GeoEvoBuilder/tree/main/MSA_Processing):
```
python run_GeoEvoBuilder.py -iP examples/ -i 3mpc_A.pdb --chainID A --MSA examples/3mpc_A.pt
```

### CITATION
