#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 19 07:53:26 2022

@author: rlg3


DFT DOS Thermal Proeprties
"""

import numpy as np
import pdos_integration as pint
from jarvis.core.atoms import Atoms
from jarvis.db.figshare import data as jdata
from jarvis.core.spectrum import Spectrum
import  matplotlib.pyplot as plt
from jarvis.core.specie import Specie
import os
import json

#Change to directory with the phonon predictions
os.chdir('../../vasp_phonon_runs')



def transform_dos(match, dos_array):
    zero_indx = np.where(dos_array[:,0] > 0)[0][0] - 1
    min_indx = np.where(dos_array[:,0] > -300)[0][0] - 1
    max_indx = np.where(dos_array[:,0] > 1000)[0][0] - 1
    #Convert from THz to inv cm just so it works with methods
    freq = dos_array[zero_indx:max_indx, 0]
    dos = dos_array[zero_indx:max_indx, 1]
    #Scale to 3N
    int_dos = pint.integrate_dos(freq, dos)
    form_unit = pint.get_natoms_form_unit(match)
    scale = (int_dos / form_unit) / 3.0
    dos = dos / scale 
    return np.array([freq, dos])


def get_molar_mass(atoms):
    mm = 0
    for k,v in atoms.composition.reduce()[0].items():
        el_sub = Specie(k)
        mm_sub = el_sub.atomic_mass * v
        mm = mm + mm_sub
    return mm

jid_list = ['JVASP-116263',\
       'JVASP-116267', 'JVASP-116273',\
           'JVASP-117472', 'JVASP-112289', 'JVASP-62656']
    
id_list = ["49215", "49228", "49250", "53311", "52067", "5660"]

dft_3d = jdata("dft_3d")

fig, ax = plt.subplots(6, 2, figsize = (6, 15))
fig.tight_layout(h_pad = 4, w_pad = 2)

i = 1
for jid in jid_list:
    print(jid)
    os.chdir('{}_PBEBO_FD_ELAST/MAIN-ELASTIC-bulk@{}'.format(jid,jid))  
    
    match = next(d for d in dft_3d if d["jid"] == jid)
    
    atoms = Atoms.from_dict(match['atoms'])
    
    dos_file = 'total_dos.dat'
    
    icm_to_thz = 2.99792458e-2
    
    
    dos_array_full = np.genfromtxt(dos_file, skip_header = 1)
    
    
    dos_array = transform_dos(match, dos_array_full)
    
    # zero_indx = np.where(dos_array[0,:] > 0)[0][0] - 1
    
    # min_indx = np.where(dos_array[:,0] > -300)[0][0] - 1
    
    # max_indx = np.where(dos_array[:,0] > 1000)[0][0] - 1
    
    Cv = pint.heat_capacity(dos_array[0], dos_array[1], T = 300)
    Svib = pint.vibrational_entropy(dos_array[0], dos_array[1], T = 300)
    tau_iso = pint.isotopic_tau(match, dos_array[0], dos_array[1])
    
    mm = get_molar_mass(atoms)
    
    Cv_kg = Cv / mm * 1e3
    
    Svib_kg = Svib / mm * 1e3
    
    #Write the thermal properties to file
    with open("therm_props.txt", 'w') as f:
        f.write("Specific Heat Capacity: {}\n".format(Cv_kg))
        f.write("Specific vibraitonal entropy: {}\n".format(Svib_kg))
        f.write("Phonon-isotope scattering: {}\n".format(tau_iso))
        
    
    #Abridged DOS Plot
    plt.subplot(6,2,i)
    plt.plot(dos_array[0], dos_array[1])
    plt.xlabel(r'Frequency (cm$^{-1}$)')
    plt.ylabel('DOS')
    plt.title(jid)
        
    #plt.savefig('test_dos_abridged.pdf', bbox_inches = 'tight')
    
    #Full DOS plot
    plt.subplot(6,2,i+1)
    plt.plot(dos_array_full[:, 0], dos_array_full[:, 1])
    plt.xlabel(r'Frequency (cm$^{-1}$)')
    plt.ylabel('DOS')
    plt.title(jid)
        
    #plt.savefig('test_dos_full.pdf', bbox_inches = 'tight')
    i = i+2
    
    os.chdir("../../")
    
plt.savefig('vasp_dos_plots_multpanel.pdf', bbox_inches = 'tight')




'''
Plot the ALIGNN spectra
'''

run = 'run21'

#Load the stable DOS from the JSON file

dos_file = '../{}/pred_data.json'.format(run)

with open(dos_file) as json_file:
    dos_dict = json.load(json_file)
    
    
fig, ax = plt.subplots(2, 3, figsize = [10, 5])
fig.tight_layout(h_pad = 4, w_pad = 2)

freq = np.linspace(-300,1000,66)

i = 1
for jid in id_list:
    match = [d for d in dos_dict if d["id"] == jid]
    plt.subplot(2,3, i)
    plt.plot(freq, match[0]['pred'])
    plt.xlabel(r'Frequency (cm$^{-1}$)')
    plt.ylabel('DOS')
    plt.title(jid_list[i-1])
    i = i+1

plt.savefig('ALIGNN_spectra_pred.pdf', bbox_inches = 'tight')
    
    
    
    
    
    