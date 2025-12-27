"""A deep learning framework for efficient functional and thermostable protein design"""

#from geoevobuilder.src import Sc
from geoevobuilder.src import Se
#from geoevobuilder.sampling import samPling
from geoevobuilder.builder.test import builder,Val_builder
from geoevobuilder.Utils.pdb_processor import fetchseq,batch_fea,batch_tri,Rotamer_number,AA,AA_reverse

sequence_design = Se.ESM2finetuning_Net()
sequence_design_withMSA = Se.ESM2finetuning_Net('MSA')

__version__ = "0.1.0"
