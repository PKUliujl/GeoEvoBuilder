import argparse
import os

from pathlib import Path


def parse_fixed_res(value):
    if value is None:
        return []

    try:
        if Path(value).is_file():
            with open(value, "r") as f:
                return [line.strip() for line in f if line.strip()]
    except (OSError, TypeError):
        pass

    return [item.strip() for item in value.split(",") if item.strip()]


def run_inputparameters():
    description = "To better use GeoEvoBuilder toolkit for functional protein sequence design, please add some of these parameters"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--purpose",
        type=int,
        choices=[0, 1],
        help="0 is used for functional sequence design. Only activated at 1 when using additional MSA information (--MSA is given). default: 0",
        default=0,
    )
    parser.add_argument(
        "--inputPATH",
        "-iP",
        type=str,
        help="The directory path should contain the pdb file",
        required=True,
    )
    parser.add_argument(
        "--chainID", type=str, help="The given chain for design or repacking", required=True
    )
    parser.add_argument(
        "--inputfile",
        "-i",
        type=str,
        help="A pdb file under inputPATH, eg. T1024_A.pdb.",
        required=True,
    )
    parser.add_argument(
        "--ST",
        type=float,
        help="Sampling temperature (scaling factor) belongs to (0,1). Larger values increase the number of interation steps for convegence. default: 0.1",
        default=0.1,
    )
    parser.add_argument(
        "--n", type=int, help="The number of designed sequences. default: 10", default=10
    )
    parser.add_argument(
        "--SM",
        type=int,
        help="Starting mode with different sequence initializations (how many residues are randomly masked). default = 100 (%%) masked, i.e. random initialization)",
        default=100,
    )
    parser.add_argument(
        "--Fixed",
        type=str,
        help="Those residues are expected to remain unchanged with names following the format either {residuenumber_chainid_resname,residuenumber_chainid_resname}, e.g. 1_A_W,10_A_C, or path to a file with one residue per line. Please note that the residuenumber is renumbered begining at 0.",
        default=None,
    )
    parser.add_argument(
        "--outputPATH",
        "-oP",
        type=str,
        help="the directory path of the outputfile. default: inputPATH",
    )
    parser.add_argument(
        "--verbose", type=int, choices=[0, 1], default=0, help="Display the intermediate sequences"
    )
    # parser.add_argument('--noCYS',type=int,choices={0,1},help='Do NOT generate Cys as far as possible', default = 0)
    parser.add_argument(
        "--MSA", type=str, help="The distribution of sequence profile from MSA", default=None
    )

    args = parser.parse_args()

    args.Fixed = parse_fixed_res(args.Fixed)

    if args.MSA is None:
        if args.purpose == 1:
            parser.error("With --purpose 1, please provide --MSA information")
    if args.MSA is not None:
        args.purpose = 1
    if args.ST < 0 or args.ST > 1:
        parser.error("--ST, Sampling temperature(scaling factor) belongs to (0,1)")
    if args.SM < 0 or args.SM > 100:
        parser.error("--SM, random masking rate for sequence initialization belings to (0,100)")

    if args.outputPATH is None:
        args.outputPATH = args.inputPATH

    return args


if __name__ == "__main__":
    run_inputparameters()
