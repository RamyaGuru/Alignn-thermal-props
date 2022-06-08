#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 11:37:20 2022

@author: rlg3

Run model on full dft_3d dataset of structures
"""

import torch
from jarvis.core.atoms import Atoms
from alignn.pretrained import get_multiple_predictions
from jarvis.db.figshare import data
from alignn.models.alignn import ALIGNN,ALIGNNConfig
device = "cpu"
alignn_layers=6 #4
max_neighbors=20 #12
gcn_layers=6
if torch.cuda.is_available():
    device = torch.device("cuda") # might have to change to "cpu"
    
run = 'run21'

#For phononDos make output_features 200, also add gcn_layers etc. parameters if different from default ALIGNN parameters
model = ALIGNN(ALIGNNConfig(name="alignn", alignn_layers=alignn_layers, gcn_layers=gcn_layers, output_features=66))
model.load_state_dict(torch.load('../../{}/checkpoint_600.pt'.format(run), map_location=device)["model"])
model.to(device)
model.eval()
atoms_array=[]

cod=data('cod') #dft_3d, cod
for ii,i in enumerate(cod):
  atoms_array.append(Atoms.from_dict(i['atoms']))

get_multiple_predictions(atoms_array=atoms_array,model=model,max_neighbors=max_neighbors, output_features = 66)
#import pandas as pd
#df=pd.read_json('{}_pred_data.json'.format(run))
# df.sort_values('pred',ascending=False)
# print(df)



