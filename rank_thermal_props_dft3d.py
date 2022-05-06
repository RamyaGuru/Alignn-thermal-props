#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 13:48:23 2022

@author: rlg3

Multiply by molar mass and rank
"""

import pandas as pd
import numpy as np 

run = 'run21'

dft3d_df = pd.read_csv('output_files/{}_thermal_props_dft3d.csv'.format(run))

mol_mass = np.load('molar_mass_dft3d.npy')

form_unit = np.load('form_unit.npy')


Svib_kg = dft3d_df['Svib_kg'] 

Cv_kg = dft3d_df['Cv_kg'] 

Cv_max_indx = np.argsort(Cv_kg)[-50:]
Cv_min_indx = np.argsort(Cv_kg)[:50]

print('Maximum heat capacities and compositions')
Cv_max = Cv_kg[Cv_max_indx]
print(np.array(Cv_max))
print(form_unit[Cv_max_indx])

print('Minimum heat capacities and compositions')
Cv_min = Cv_kg[Cv_min_indx]
print(np.array(Cv_min))
print(form_unit[Cv_min_indx])

print("Maximum vibratonal entropies and compositions")
Svib_max_indx = np.argsort(Svib_kg)[-50:]
Svib_min_indx = np.argsort(Svib_kg)[:50]

Svib_max = Svib_kg[Svib_max_indx]
print(np.array(Svib_max))
print(form_unit[Svib_max_indx])

#Fix this to be greater than 0?
print("Minimum vibratonal entropies and compositions")
Svib_min = Svib_kg[Svib_min_indx]
print(np.array(Svib_min))
print(form_unit[Svib_min_indx])



dft3d_df_mod = dft3d_df[dft3d_df['isotope_scatt'] > 0]
form_unit_mod = form_unit[dft3d_df['isotope_scatt'] > 0]

print("Minimum isotope scattering rate")
iso_scatt_min_indx = np.argsort(dft3d_df_mod['isotope_scatt'])[:50]
print(np.array(dft3d_df_mod['isotope_scatt'])[iso_scatt_min_indx])
print(form_unit_mod[iso_scatt_min_indx])

iso_scatt_max_indx = np.argsort(dft3d_df_mod['isotope_scatt'])[-50:]
print(np.array(dft3d_df_mod['isotope_scatt'])[iso_scatt_max_indx])
print(form_unit_mod[iso_scatt_max_indx])
#iso_scatt_min = np.argpartition(dft3d_df_mod['isotope_scatt'], 100)[:100].values