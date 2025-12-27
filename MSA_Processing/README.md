# MSA Processing

## USAGE

```bash
Transform the MSA information to the sequence profile distribution

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input file (aln format)
  -s SCALE, --scale SCALE
                        Scaling factor. Utilizing smaller values leads to a sharper distribution
  -o OUTPUT, --output OUTPUT
                        Output file
```

Here is an example:

```bash
python run_sequenceprofile.py -i 3mpc_A.aln -o 3mpc_A.pt
```
