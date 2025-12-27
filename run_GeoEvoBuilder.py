
import torch
from geoevobuilder.common.run_argparse import *
from geoevobuilder.sampling import samPling

def GeoEvoBuilder(args,models):
    samPling(args,models)
    return

if __name__ == '__main__':
    args = run_inputparameters()
    #device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')  
    models = torch.load("geoevobuilder/params/Se.pt", weights_only=False)#.to(device)
    GeoEvoBuilder(args,models)
    print(' '*13+'##############################################\n\
             ##      Thanks for using GeoEvoBuilder      ##\n\
             ##    More details see mdl.ipc.pku.edu.cn   ##\n\
             ##############################################\n')
