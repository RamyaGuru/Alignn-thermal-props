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
from pyvalem.formula import Formula
import json

pd.reset_option('display.float_format')#, '{:.2E}'.format)

run = 'run21'

json_file = 'output_files/{}_thermal_props_dft3d_mod.json'.format(run)

with open(json_file) as infile:
    dft3d_dict = json.load(infile)
    
dft3d_df = pd.DataFrame.from_dict(dft3d_dict)

#dft3d_df = pd.read_csv('output_files/{}_thermal_props_dft3d_mod.csv'.format(run))

mol_mass = np.load('output_files/molar_mass_dft3d.npy')

form_unit = np.load('output_files/form_unit.npy')


'''
Exclude n_elem < 2 and > 16. Exclude radioactive elements?
'''

radio_list = ['Tc', 'Po', 'Ac', 'Rn', 'Fr', 'Ra', 'Rf', 'Db', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']


with open('../../{}/pred_data.json'.format(run)) as json_file:
    dos_dict = json.load(json_file)

Svib_sorted = dft3d_df.sort_values(by = ['Svib_norm'])
Cv_sorted = dft3d_df.sort_values(by = ['Cv_norm'])
iso_scatt = dft3d_df.sort_values(by = ['isotope_scatt'])

iso_scatt = iso_scatt[iso_scatt['isotope_scatt'] > 0]

'''
Vibraitonal Entropy
'''

print('Minimum Vibrational Entropies')

print(Svib_sorted[:10])


print('Maximum vibrational Entropies')

print(Svib_sorted[-10:])

Svib_latex_dict = {'Low JID' : np.array(Svib_sorted[:10]['jid']),\
                    'Low Composition' : np.array(Svib_sorted[:10]['form_unit']),\
                        'Low Svib' : np.array(Svib_sorted[:10]['Svib_kg']),\
                'High JID' : np.array(Svib_sorted[-10:]['jid']),\
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



'''
Plot the lowest vibrational entropy spectra
'''

fig, ax = plt.subplots(1, 3, figsize = (9, 2))
fig.tight_layout(h_pad = 2)

freq = np.linspace(0, 1000, 51)
n = 1
for indx in range(3):
    plt.subplot(1, 3, indx+1)
    plt.plot(freq, np.array(Svib_sorted['scaled_dos'][indx:(indx + 1)])[0],\
             color = 'xkcd:black')
    cf = Formula(np.array(Svib_sorted['form_unit'][indx:(indx + 1)])[0])
    jid = np.array(Svib_sorted['jid'][indx:(indx + 1)])[0]
    title_str = '$' + cf.latex + '$' + '; JID:' + jid
    plt.title(title_str)
    plt.ylim([-0.001, 0.01])
    plt.yticks([0, 0.005, 0.01])
fig.text(-0.03, 0.5, 'Scaled DOS (a.u.)', va='center', rotation='vertical', fontsize = 12)
fig.text(0.5, -0.02, r'Frequency (cm$^{-1}$)', ha='center', fontsize = 12)
plt.savefig('lowest_svib_spectra.pdf', bbox_inches = 'tight')
    
    
    
'''
Plot the highest vibrational entropy spectra
'''

fig, ax = plt.subplots(1, 3, figsize = (9, 2))
fig.tight_layout(h_pad = 2)

freq = np.linspace(0, 1000, 51)
n = 1
for indx in range(3):
    plt.subplot(1, 3, indx+1)
    if indx == 0:
        plt.plot(freq, np.array(Svib_sorted['scaled_dos'][-1:])[0],\
                 color = 'xkcd:black')
        cf = Formula(np.array(Svib_sorted['form_unit'][-1:])[0])
        jid = np.array(Svib_sorted['jid'][-1:])[0]
    else:
        plt.plot(freq, np.array(Svib_sorted['scaled_dos'][-(indx+1):-(indx)])[0],\
                 color = 'xkcd:black')        
        cf = Formula(np.array(Svib_sorted['form_unit'][-(indx+1):-(indx)])[0])
        jid = np.array(Svib_sorted['jid'][-(indx+1):-(indx)])[0]
    title_str = '$' + cf.latex + '$' + '; JID:' + jid
    plt.title(title_str)
fig.text(-0.03, 0.5, 'Scaled DOS (a.u.)', va='center', rotation='vertical', fontsize = 12)
fig.text(0.5, -0.02, r'Frequency (cm$^{-1}$)', ha='center', fontsize = 12)
plt.savefig('highest_svib_spectra.pdf', bbox_inches = 'tight') 


'''
Plot the lowest phonon-isotope scattering
'''

fig, ax = plt.subplots(1, 3, figsize = (9, 2))
fig.tight_layout(h_pad = 2)

freq = np.linspace(0, 1000, 51)
n = 1
for indx in range(3):
    plt.subplot(1, 3, indx+1)
    plt.plot(freq, np.array(iso_scatt['scaled_dos'][indx:(indx + 1)])[0],\
             color = 'xkcd:black')
    cf = Formula(np.array(iso_scatt['form_unit'][indx:(indx + 1)])[0])
    jid = np.array(iso_scatt['jid'][indx:(indx + 1)])[0]
    title_str = '$' + cf.latex + '$' + '; JID:' + jid
    plt.title(title_str)
fig.text(-0.03, 0.5, 'Scaled DOS (a.u.)', va='center', rotation='vertical', fontsize = 12)
fig.text(0.5, -0.02, r'Frequency (cm$^{-1}$)', ha='center', fontsize = 12)
plt.savefig('lowest_tau_spectra.pdf', bbox_inches = 'tight')
    
    
    
'''
Plot the highest phonon-isotope scattering
'''

fig, ax = plt.subplots(1, 3, figsize = (9, 2))
fig.tight_layout(h_pad = 2)

freq = np.linspace(0, 1000, 51)
n = 1
for indx in range(3):
    plt.subplot(1, 3, indx+1)
    if indx == 0:
        plt.plot(freq, np.array(Cv_sorted['scaled_dos'][-1:])[0],\
                 color = 'xkcd:black')
        cf = Formula(np.array(Cv_sorted['form_unit'][-1:])[0])
        jid = np.array(iso_scatt['jid'][-1:])[0]
    else:
        plt.plot(freq, np.array(Cv_sorted['scaled_dos'][-(indx+1):-(indx)])[0],\
                 color = 'xkcd:black')        
        cf = Formula(np.array(Cv_sorted['form_unit'][-(indx+1):-(indx)])[0])
        jid = np.array(Cv_sorted['jid'][-(indx+1):-(indx)])[0]
    title_str = '$' + cf.latex + '$' + '; JID:' + jid
    plt.title(title_str)
fig.text(-0.03, 0.5, 'Scaled DOS (a.u.)', va='center', rotation='vertical', fontsize = 12)
fig.text(0.5, -0.02, r'Frequency (cm$^{-1}$)', ha='center', fontsize = 12)
plt.savefig('highest_tau_spectra.pdf', bbox_inches = 'tight') 


'''
Plot the lowest heat capacity
'''

fig, ax = plt.subplots(1, 3, figsize = (9, 2))
fig.tight_layout(h_pad = 2)

freq = np.linspace(0, 1000, 51)
n = 1
for indx in range(3):
    plt.subplot(1, 3, indx+1)
    plt.plot(freq, np.array(Cv_sorted['scaled_dos'][indx:(indx + 1)])[0],\
             color = 'xkcd:black')
    cf = Formula(np.array(Cv_sorted['form_unit'][indx:(indx + 1)])[0])
    jid = np.array(Cv_sorted['jid'][indx:(indx + 1)])[0]
    title_str = '$' + cf.latex + '$' + '; JID:' + jid
    plt.title(title_str)
fig.text(-0.03, 0.5, 'Scaled DOS (a.u.)', va='center', rotation='vertical', fontsize = 12)
fig.text(0.5, -0.02, r'Frequency (cm$^{-1}$)', ha='center', fontsize = 12)
plt.savefig('lowest_cv_spectra.pdf', bbox_inches = 'tight')
    
    
    
'''
Plot the highest heat capacity
'''

fig, ax = plt.subplots(1, 3, figsize = (9, 2))
fig.tight_layout(h_pad = 2)

freq = np.linspace(0, 1000, 51)
n = 1
for indx in range(3):
    plt.subplot(1, 3, indx+1)
    if indx == 0:
        plt.plot(freq, np.array(Cv_sorted['scaled_dos'][-1:])[0],\
                 color = 'xkcd:black')
        cf = Formula(np.array(Cv_sorted['form_unit'][-1:])[0])
        jid = np.array(iso_scatt['jid'][-1:])[0]
    else:
        plt.plot(freq, np.array(Cv_sorted['scaled_dos'][-(indx+1):-(indx)])[0],\
                 color = 'xkcd:black')        
        cf = Formula(np.array(Cv_sorted['form_unit'][-(indx+1):-(indx)])[0])
        jid = np.array(Cv_sorted['jid'][-(indx+1):-(indx)])[0]
    title_str = '$' + cf.latex + '$' + '; JID:' + jid
    plt.title(title_str)
fig.text(-0.03, 0.5, 'Scaled DOS (a.u.)', va='center', rotation='vertical', fontsize = 12)
fig.text(0.5, -0.02, r'Frequency (cm$^{-1}$)', ha='center', fontsize = 12)
plt.savefig('highest_cv_spectra.pdf', bbox_inches = 'tight') 



    
# for n in range(3,6):
#     plt.subplot(2, 3, n+1)
#     indx = n - 3
#     if indx == 0:
#         plt.plot(freq, np.array(Svib_sorted['scaled_dos'][-1:])[0],\
#                  color = 'xkcd:black')
#         cf = Formula(np.array(Svib_sorted['form_unit'][-1:])[0])
#         jid = np.array(Svib_sorted['jid'][-1:])[0]
#     else:
#         plt.plot(freq, np.array(Svib_sorted['scaled_dos'][-(indx+1):-(indx)])[0],\
#                  color = 'xkcd:black')        
#         cf = Formula(np.array(Svib_sorted['form_unit'][-(indx+1):-(indx)])[0])
#         jid = np.array(Svib_sorted['jid'][-(indx+1):-(indx)])[0]
#     title_str = '$' + cf.latex + '$' + '; JID:' + jid
#     plt.title(title_str)        
# plt.savefig('highest_svib_spectra.pdf', bbox_inches = 'tight')
    
    
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