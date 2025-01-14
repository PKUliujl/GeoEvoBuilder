# GeoEvoBuilder
GeoEvoBuilder, an efficient zero-shot learning method, is capable to improve both protein thermal stability and activity. This approach integrates structure-based de novo sequence design with protein language model for diverse functional protein design.

<!-- ![Alt text](https://github.com/PKUliujl/GeoEvoBuilder/blob/main/image/flow.jpg) -->

### INSTALLATION
<!-- ====================== -->
1. GeoEvoBuilder relies on [ESM2](https://github.com/facebookresearch/esm), please install it using `pip install fair-esm`.
2. GeoEvoBuilder relies on the environment used in [GeoSeqBuilder](https://github.com/PKUliujl/GeoSeqBuilder/).

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
  --Fixed FIXED         Those residues are expected to remain unchanged with names following the format {residuenumber_chainid_resname,residuenumber_chainid_resname},
                        e.g. 1_A_W,10_A_C Please note that the residuenumber is renumbered begining at 0.
  --outputPATH OUTPUTPATH, -oP OUTPUTPATH
                        the directory path of the outputfile. default: inputPATH
  --verbose {0,1}       Display the intermediate sequences
  --MSA MSA             The distribution of sequence profile from MSA
```
