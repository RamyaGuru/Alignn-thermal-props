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

#For phononDos make output_features 200, also add gcn_layers etc. parameters if different from default ALIGNN parameters
model = ALIGNN(ALIGNNConfig(name="alignn", alignn_layers=alignn_layers, gcn_layers=gcn_layers, output_features=201))
model.load_state_dict(torch.load('../../run11/temp/checkpoint_600.pt', map_location=device)["model"])
model.to(device)
model.eval()
atoms_array=[]

dft_3d=data('dft_3d') #dft_3d, cod
for ii,i in enumerate(dft_3d):
  atoms_array.append(Atoms.from_dict(i['atoms']))
  if ii==5:
    break
get_multiple_predictions(atoms_array=atoms_array,model=model,max_neighbors=max_neighbors, output_features = 201)
import pandas as pd
df=pd.read_json('pred_data.json')
df.sort_values('pred',ascending=False)
print(df)



