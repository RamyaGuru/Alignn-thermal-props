#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 15:13:09 2022

@author: rlg3


DOS-derived Thermal Properties

integrated DOS fixed to 3N
"""
import pdos_integration as pint
import numpy as np
import json
import pandas as pd
from jarvis.db.figshare import data as jdata
from jarvis.analysis.elastic.tensor import ElasticTensor
from jarvis.core.atoms import Atoms
from jarvis.core.specie import Specie

from sklearn.metrics import mean_absolute_error

import matplotlib.pyplot as plt

run = 'run21'
input_file = '../../{}/true_stable_predictions.json'.format(run)



with open(input_file) as json_file:
    dos_dict = json.load(json_file)
    
full_freq = np.linspace(-300, 1000, len(dos_dict[0]['target']))
zero_indx = np.where(full_freq > 0)[0][0] - 1
print(zero_indx)

edos_pdos = jdata("edos_pdos")
dft_3d = jdata("dft_3d")
# x = []
# for i in dos_dict:
#     if i['predictions'] != 'na':
#         x.append(i)

jid_list = []

iso_t = []
S_vib_t = []
Cv_t = []

iso_p = []
S_vib_p = []
Cv_p = []

Cv_debT_t = []
Cv_debT_p = []

debyeT_list = []
for i in dos_dict:
    # Get the JARVIS ID
    jid = i["id"]
    start = jid.find('JVASP')
    end = jid.find('.vasp')
    jid = jid[start:end]
    jid_list.append(jid)
    match = next(d for d in dft_3d if d["jid"] == jid)
    # Get number of atoms in formula unit
    #atoms = Atoms.from_dict(p['atoms'])
    form_unit = pint.get_natoms_form_unit(match)
    target = np.array(i['target'])[zero_indx:]
    pred = np.array(i['predictions'])[zero_indx:]
    freq = np.linspace(0, 1000, len(target))
    #Store target thermal properties
    intDOS_t = pint.integrate_dos(freq, target)
    scale = (intDOS_t / form_unit) / 3.0
    target = target / scale
    try:
        iso_t.append(pint.isotopic_tau(match, freq, target))
    except:
        iso_t.append(0)
    S_vib_t.append(pint.vibrational_entropy(freq, target))
    Cv_t.append(pint.heat_capacity(freq, target))
    #Store predicted thermal properties
    intDOS_p = pint.integrate_dos(freq, pred)
    scale = (intDOS_p / form_unit) / 3.0
    pred = pred / scale
    try:
        iso_p.append(pint.isotopic_tau(match, freq, pred))
    except:
        iso_p.append(0)
    S_vib_p.append(pint.vibrational_entropy(freq, pred))
    Cv_p.append(pint.heat_capacity(freq, pred))
    '''
    Equal fraction of the Debye Temperature
    '''
    atoms = Atoms.from_dict(match['atoms'])
    et = ElasticTensor(match['elastic_tensor'])
    debyeT = et.debye_temperature_toberer(atoms)
    debyeT_list.append(debyeT)
    try:
        Cv_debT_t.append(pint.heat_capacity(freq, target, debyeT / 2))
        Cv_debT_p.append(pint.heat_capacity(freq, pred, debyeT / 2))
    except:
        Cv_debT_t.append(0)
        Cv_debT_p.append(0)        

therm_output = {'JID' : jid_list,
          'iso_scattering (Hz)' : {'target' : iso_t, 'prediction' : iso_p},
          'S_vib (J/mol/K)' : {'target' : S_vib_t, 'prediction' : S_vib_p},
          'Cp (J/mol/K)' : {'target' : Cv_t, 'prediction' : Cv_p},
          'Cp (J/mol/K) DebT' : {'target': Cv_debT_t , 'prediction': Cv_debT_p}}


with open('output_files/{}_thermal_props_scale_3N.json'.format(run), 'w') as out_file:
    json.dump(therm_output, out_file)
    

#df.to_csv(run + label + '_thermal_props_orig.csv')




