#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 13:48:23 2022

@author: rlg3

Multiply by molar mass and rank
"""

import pandas as pd
import numpy as np 


dft3d_df = pd.read_csv('thermal_props_dft3d.csv')

mol_mass = np.load('molar_mass_dft3d.npy')

form_unit = np.load('form_unit.npy')


Svib_kg = dft3d_df['S_vib'] / mol_mass

Cv_kg = dft3d_df['Cv'] / mol_mass


Cv_max = np.argpartition(Cv_kg, -100)[-100:].values

Svib_max = np.argpartition(Svib_kg, -100)[-100:].values

dft3d_df_mod = dft3d_df[dft3d_df['isotope_scatt'] > 0]
iso_scatt_min = np.argpartition(dft3d_df_mod['isotope_scatt'], 100)[:100].values