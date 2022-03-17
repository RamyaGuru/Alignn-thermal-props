#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 17:14:01 2022

@author: rlg3

Run Thermal Property Calculation on the pred_data.json file

Potential issue:
    1. Training set was scaled by maximum intensity for each spectrum 
    2. For the new 
"""

import pdos_integration as pint
import json
import matplotlib.pyplot as plt
import numpy as np

with open('output_files/pred_data.json') as json_file:
    pred = json.load(json_file)
    
    
#Plot each of the spectra

# for p in pred:
#     freq = np.linspace(0, 1000, 201)
#     plt.figure()
#     plt.plot(freq, p['pred'])
    
    
#


