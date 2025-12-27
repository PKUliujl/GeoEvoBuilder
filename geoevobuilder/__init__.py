"""A deep learning framework for efficient functional and thermostable protein design"""

# from geoevobuilder.src import Sc
# from geoevobuilder.sampling import samPling
from geoevobuilder.builder.test import Val_builder, builder
from geoevobuilder.src import Se
from geoevobuilder.Utils.pdb_processor import (
    AA,
    AA_reverse,
    Rotamer_number,
    batch_fea,
    batch_tri,
    fetchseq,
)

sequence_design = Se.ESM2finetuning_Net()
sequence_design_withMSA = Se.ESM2finetuning_Net("MSA")

__version__ = "0.1.0"
