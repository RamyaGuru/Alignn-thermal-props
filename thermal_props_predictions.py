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


with open('output_files/pred_data_dft3d.json') as json_file:
    pred = json.load(json_file)
    

freq = np.linspace(0, 1000, 201)

#Calculate the integral of each prediction

Cv_mol = []
Cv_kg = []
Svib = []
Svib_kg = []
iso_tau = []
id_list = []
mol_mass = []
form_unit = []

for p in pred:
    atoms = Atoms.from_dict(p['atoms'])
    mm = atoms.composition.weight * 1e-3
    fu = atoms.composition.reduced_formula
    form_unit.append(fu)
    #mol_mass.append(mm)
    # id_list.append(p['id'])
    # int_dos = np.trapz(p['pred'], freq)
    # form_unit = pint.get_natoms_form_unit(p)
    # scale = (int_dos / form_unit) / 3.0
    # p['pred'] = np.array(p['pred']) / scale
    # Cv_mol = pint.heat_capacity(freq, p['pred'], T = 300)
    # Cv_mol.append(Cv_mol)
    # Cv_kg.append(Cv_mol / mol_mass)
    # Svib_mol = pint.vibrational_entropy(freq, p['pred'], T = 300)
    # Svib_mol.append(Svib_mol)
    # Svib_kg.append(Svib_mol / mol_mass)
    # try:
    #     iso_tau.append(pint.isotopic_tau(p, freq, p['pred']))
    # except:
    #     iso_tau.append(0)


#np.save('molar_mass_dft3d', mol_mass)
np.save('form_unit', form_unit)

# output = {'id' : id_list,
#           'isotope_scatt' : iso_tau,
#           'S_vib' : Svib,
#           'Cv_mol' : Cv_mol}


# df = pd.DataFrame(output)

# df.to_csv('thermal_props_dft3d.csv')

# for p in pred:
#     freq = np.linspace(0, 1000, 201)
#     plt.figure()
#     plt.plot(freq, p['pred'])
    
    
#


