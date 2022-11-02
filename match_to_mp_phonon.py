#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 11:25:59 2022

@author: rlg3

Parse MP Phonon Data


DOS Frequency: inverse cm
"""

import os
import json
from jarvis.db.figshare import data
import numpy as np
import pdos_integration as pint
from jarvis.core.atoms import Atoms
from jarvis.core.specie import Specie
import pandas as panda
from jarvis.core.spectrum import Spectrum
from sklearn.metrics import mean_absolute_error, r2_score

phonon_db = "../../phonon/"


'''
Materials Project Data
'''

total_ph_dict = {}

for file in os.listdir(phonon_db):
    fpath = phonon_db + file
    with open(fpath) as jfile:
        ph_dict = json.load(jfile)   
    total_ph_dict[ph_dict['metadata']['material_id']] = {'phonon' :\
                                                         ph_dict['phonon'],\
                                                             'metadata': ph_dict['metadata'],
                                                             'flags': ph_dict['flags']}


mp_mpid_list = list(total_ph_dict.keys())


'''
Jarvis data
'''

run = 'run21'
pred_file = '../../{}/pred_data.json'.format(run)

#Load data from the prediction set
with open(pred_file) as json_file:
    pred_data = json.load(json_file)

dft_3d = data("dft_3d")

pred_mpids = []

for pd in pred_data:
    jarvis_entry = dft_3d[int(pd["id"])]
    mpid_pd = jarvis_entry['reference']
    pd["mpid"] = mpid_pd
    pred_mpids.append(mpid_pd)

new_pred_file = '../../{}/pred_data_augmented.json'.format(run)

with open(new_pred_file, 'w') as outfile:
    json.dump(pred_data, outfile)

overlap = set(pred_mpids) & set(mp_mpid_list)

'''
Compute Thermal Properties of the MP Phonon DOS
'''
mpid_list = []
mp_unstable_list = []
dos_scale_list = []
freq_list = []
Cv_mol_list = []
Cv_kg_list = []
Svib_mol_list = []
Svib_kg_list = []
mol_mass = []
form_unit_list = []
mol_mass_list = []
iso_tau_list = []

for mpid in overlap:    
    mp_unstable = total_ph_dict[mpid]['flags']['has_neg_fr']
    if mp_unstable == False:
        mpid_list.append(mpid)
        full_freq = np.array(total_ph_dict[mpid]['phonon']['dos_frequencies'])
        find_zero = np.where(full_freq > 0)[0][0] - 1
        try:
            find_1000 = np.where(full_freq > 1000)[0][0] - 1
        except:
            find_1000 = None
        freq = full_freq[find_zero:find_1000]
        freq_list.append(list(freq))
        dos_full = np.array(total_ph_dict[mpid]['phonon']['ph_dos'])
        dos = dos_full[find_zero:find_1000]
        s = Spectrum(x= freq, y=dos)
        hist, bins = pint.histo_making(1000, 0, 20, s)
        dos = hist
        freq = bins
        int_dos = np.trapz(dos, freq)
        match = next(i for i in pred_data if i["mpid"] == mpid)
        atoms = Atoms.from_dict(match['atoms'])
        mm = 0
        for k,v in atoms.composition.reduce()[0].items():
            el_sub = Specie(k)
            mm_sub = el_sub.atomic_mass * v
            mm = mm + mm_sub
        fu = atoms.composition.reduced_formula
        form_unit_list.append(fu)
        mol_mass_list.append(mm)
        fu_num = pint.get_natoms_form_unit(match)
        scale = (int_dos / fu_num) / 3.0
        dos_scale = np.array(dos) / scale
        dos_scale_list.append(list(dos_scale))
        Cv_mol = pint.heat_capacity(freq, dos_scale, T = 300)
        Cv_mol_list.append(Cv_mol)
        Cv_kg_list.append(Cv_mol / mm * 1e3)
        Svib_mol = pint.vibrational_entropy(freq, dos_scale, T = 300)
        Svib_mol_list.append(Svib_mol)
        Svib_kg_list.append(Svib_mol / mm * 1e3)
        try:
            iso_tau_list.append(pint.isotopic_tau(match, freq, dos_scale))
        except:
            iso_tau_list.append(0)

output = { 'mpid' : mpid_list,
          'molar_mass' : mol_mass_list,
          'scaled_dos' : dos_scale_list,
          'freq' : freq_list,
          'form_unit' : form_unit_list,
          'isotope_scatt' : iso_tau_list,
          'Svib_mol' : Svib_mol_list,
          'Cv_mol' : Cv_mol_list,
          'Svib_kg' : Svib_kg_list,
          'Cv_kg' : Cv_kg_list}

# df = panda.DataFrame(output)
# df.to_csv('output_files/{}_thermal_props_mp_bins.csv'.format(run))


# Going to write to a JSON file instead
json_out = 'output_files/{}_thermal_props_mp_bins_mod.json'.format(run)
with open(json_out, 'w') as outfile:
    json.dump(output, outfile)


'''
Compare ALIGNN Predictions and MP Data
'''

