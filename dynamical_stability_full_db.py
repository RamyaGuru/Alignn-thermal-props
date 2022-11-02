#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 12:10:14 2022

@author: rlg3

Dynamical Stability Labelling of DFT-3D JARVIS PDOS database
"""

from jarvis.db.figshare import data as jdata
import pdos_integration as pint
import numpy as np
import json
import matplotlib.pyplot as plt
import matplotlib as mpl
from jarvis.core.atoms import Atoms


# pdos = jdata("edos_pdos")

# dft_3d = jdata("dft_3d")
# print(pdos[0])
# print(dft_3d[0])

run = "run21"
dos_file = "../../{}/PhononData-Freq-300_1000_20.json".format(run)

with open(dos_file) as json_file:
    dos_dict = json.load(json_file)

jid_list = []
stable = []
unstable = []

freq = np.linspace(-300, 1000, len(dos_dict[0]['pdos_elast']))
        
def label_dynamic_stability(freq, dos_dict):
    stable_list = []
    unstable_list = []
    zero_indx = np.where(freq > 0)[0][0] - 1
    neg_freq = np.linspace(-300, 0, zero_indx)
    for d in dos_dict:
        intdos_neg_target = np.trapz(neg_freq, d['pdos_elast'][:zero_indx])
        intdos_target = np.trapz(freq, d['pdos_elast'])
        if intdos_neg_target / intdos_target > 0.1:
            unstable_list.append(d)
        else:
            stable_list.append(d)
    return stable_list, unstable_list  
        

stable_list, unstable_list = label_dynamic_stability(freq, dos_dict)

print(len(stable_list))

print(len(unstable_list))

# Histograms comparing stable versus unstable compounds

for uc in unstable_list:
    atoms = Atoms.from_dict(uc["atoms"])
    chem_formula = atoms.composition.reduced_formula
    print(chem_formula)






