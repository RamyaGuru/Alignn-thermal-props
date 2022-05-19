#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 13:48:23 2022

@author: rlg3

Multiply by molar mass and rank
"""

import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import json

pd.set_option('display.float_format', '{:.2E}'.format)

run = 'run21'

dft3d_df = pd.read_csv('output_files/{}_thermal_props_dft3d.csv'.format(run))

mol_mass = np.load('output_files/molar_mass_dft3d.npy')

form_unit = np.load('output_files/form_unit.npy')

with open('../../{}/pred_data.json'.format(run)) as json_file:
    dos_dict = json.load(json_file)

Svib_sorted = dft3d_df.sort_values(by = ['Svib_kg'])
Cv_sorted = dft3d_df.sort_values(by = ['Cv_kg'])
iso_scatt = dft3d_df.sort_values(by = ['isotope_scatt'])

iso_scatt = iso_scatt[iso_scatt['isotope_scatt'] > 0]

'''
Vibraitonal Entropy
'''

print('Minimum Vibrational Entropies')

print(Svib_sorted[:10])


print('Maximum vibrational Entropies')

print(Svib_sorted[-10:])

Svib_latex_dict = {'Low JID' : np.array(Svib_sorted[:10]['id']),\
                    'Low Composition' : np.array(Svib_sorted[:10]['form_unit']),\
                        'Low Svib' : np.array(Svib_sorted[:10]['Svib_kg']),\
                'High JID' : np.array(Svib_sorted[-10:]['id']),\
                    'High Composition' : np.array(Svib_sorted[-10:]['form_unit']),\
                        'High Svib' : np.array(Svib_sorted[-10:]['Svib_kg'])}

Svib_df = pd.DataFrame(Svib_latex_dict)

Svib_df = Svib_df.round(2)

Svib_str = Svib_df.to_latex(index = False)

print('Minimum Heat Capacities')

print(Cv_sorted[:10])

print('Maximum Heat Capacities')

print(Cv_sorted[-10:])

Cv_latex_dict = {'Low JID' : np.array(Cv_sorted[:10]['id']),\
                    'Low Composition' : np.array(Cv_sorted[:10]['form_unit']),\
                        'Low Cv' : np.array(Cv_sorted[:10]['Cv_kg']),\
                'High JID' : np.array(Cv_sorted[-10:]['id']),\
                    'High Composition' : np.array(Cv_sorted[-10:]['form_unit']),\
                        'High Cv' : np.array(Cv_sorted[-10:]['Cv_kg'])}

Cv_df = pd.DataFrame(Cv_latex_dict)

Cv_str = Cv_df.to_latex(index = False)


#Plot the highest CV spectrum

freq = np.linspace(-300, 1000, len(dos_dict[0]['pred']))

highest_Cv_id = np.array(Cv_sorted[-10:]['id'])[-1]

plt.figure()
plt.plot(freq, dos_dict[highest_Cv_id]['pred'])
plt.xlabel(r'Frequency cm$^{-1}$')
plt.ylabel('DOS')
plt.savefig('alignn_highest_Cv_spec.pdf', bbox_inches = 'tight')

print('Minimum Isotope Scattering Rates')

print(iso_scatt[:10])

print('Maximum Isotope Scattering Rates')

print(iso_scatt[-10:])

iso_latex_dict = {'Low JID' : np.array(iso_scatt[:10]['id']),\
                    'Low Composition' : np.array(iso_scatt[:10]['form_unit']),\
                        'Low tau' : np.array(iso_scatt[:10]['isotope_scatt']) / 1e9,\
                'High JID' : np.array(iso_scatt[-10:]['id']),\
                    'High Composition' : np.array(iso_scatt[-10:]['form_unit']),\
                        'High tau' : np.array(iso_scatt[-10:]['isotope_scatt']) / 1e9}
    
iso_df = pd.DataFrame(iso_latex_dict)



iso_str = iso_df.to_latex(index = False)


# Svib_kg = dft3d_df['Svib_kg'] 

# Cv_kg = dft3d_df['Cv_kg'] 

# Cv_max_indx = np.argsort(Cv_kg)[-50:]
# Cv_min_indx = np.argsort(Cv_kg)[:50]

# print('Maximum heat capacities and compositions')
# Cv_max = Cv_kg[Cv_max_indx]
# print(np.array(Cv_max))
# print(form_unit[Cv_max_indx])

# print('Minimum heat capacities and compositions')
# Cv_min = Cv_kg[Cv_min_indx]
# print(np.array(Cv_min))
# print(form_unit[Cv_min_indx])

# print("Maximum vibratonal entropies and compositions")
# Svib_max_indx = np.argsort(Svib_kg)[-50:]
# Svib_min_indx = np.argsort(Svib_kg)[:50]

# Svib_max = Svib_kg[Svib_max_indx]
# print(np.array(Svib_max))
# print(form_unit[Svib_max_indx])

# #Fix this to be greater than 0?
# print("Minimum vibratonal entropies and compositions")
# Svib_min = Svib_kg[Svib_min_indx]
# print(np.array(Svib_min))
# print(form_unit[Svib_min_indx])



# dft3d_df_mod = dft3d_df[dft3d_df['isotope_scatt'] > 0]
# form_unit_mod = form_unit[dft3d_df['isotope_scatt'] > 0]

# print("Minimum isotope scattering rate")
# iso_scatt_min_indx = np.argsort(dft3d_df_mod['isotope_scatt'])[:50]
# print(np.array(dft3d_df_mod['isotope_scatt'])[iso_scatt_min_indx])
# print(form_unit_mod[iso_scatt_min_indx])

# iso_scatt_max_indx = np.argsort(dft3d_df_mod['isotope_scatt'])[-50:]
# print(np.array(dft3d_df_mod['isotope_scatt'])[iso_scatt_max_indx])
# print(form_unit_mod[iso_scatt_max_indx])
#iso_scatt_min = np.argpartition(dft3d_df_mod['isotope_scatt'], 100)[:100].values