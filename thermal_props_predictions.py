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
import pandas as pd
from jarvis.core.atoms import Atoms
from jarvis.core.specie import Specie
from jarvis.db.figshare import data as jdata

run = 'run21'

with open('../../{}/pred_data.json'.format(run)) as json_file:
    dos_data = json.load(json_file)
    
dft_3d = jdata("dft_3d")
    

full_freq = np.linspace(-300, 1000, len(dos_data[0]['pred']))
zero_indx = np.where(full_freq > 0)[0][0] - 1

stable_dos = []
for d in dos_data:
    if any(np.array(d['pred'][:zero_indx]) > 0.1):
        continue
    else:
        d['pred'] = d['pred'][zero_indx:]
        stable_dos.append(d)
        

#Calculate the integral of each prediction

Cv = []
Cv_kg = []
Svib = []
Svib_kg = []
iso_tau = []
id_list = []
jid_list = []
mol_mass = []
form_unit = []
dos = []

freq = np.linspace(0, 1000, len(stable_dos[0]['pred']))
for p in stable_dos:
    atoms = Atoms.from_dict(p['atoms'])
    # I should calculate molar mass myself
    #mm = atoms.composition.weight * 1e-3
    mm = 0
    for k,v in atoms.composition.reduce()[0].items():
        el_sub = Specie(k)
        mm_sub = el_sub.atomic_mass * v
        mm = mm + mm_sub
    fu = atoms.composition.reduced_formula
    form_unit.append(fu)
    mol_mass.append(mm)
    id_list.append(int(p['id']))
    jid_list.append(dft_3d[int(p['id'])]['jid'])
    int_dos = np.trapz(p['pred'], freq)
    fu_num = pint.get_natoms_form_unit(p)
    scale = (int_dos / fu_num) / 3.0
    p['pred'] = np.array(p['pred']) / scale
    dos.append(p['pred'])
    Cv_mol = pint.heat_capacity(freq, p['pred'], T = 300)
    Cv.append(Cv_mol)
    Cv_kg.append(Cv_mol / mm * 1e3)
    Svib_mol = pint.vibrational_entropy(freq, p['pred'], T = 300)
    Svib.append(Svib_mol)
    Svib_kg.append(Svib_mol / mm * 1e3)
    try:
        iso_tau.append(pint.isotopic_tau(p, freq, p['pred']))
    except:
        iso_tau.append(0)


np.save('molar_mass_dft3d', mol_mass)
np.save('form_unit', form_unit)

output = { 'id' : id_list,
        'jid' : jid_list,
          'molar_mass' : mol_mass,
          'form_unit' : form_unit,
          'isotope_scatt' : iso_tau,
          'Svib_kg' : Svib_kg,
          'Cv_kg' : Cv_kg}

df = pd.DataFrame(output)
df.to_csv('output_files/{}_thermal_props_dft3d.csv'.format(run))

# for p in pred:
#     freq = np.linspace(0, 1000, 201)
#     plt.figure()
#     plt.plot(freq, p['pred'])
    
    



