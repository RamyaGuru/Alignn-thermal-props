#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 23:56:51 2022

@author: rlg3

Calculate the C_V at half of the Debye temperature
"""

import pdos_integration as pint
import numpy as np
import json
import pandas as pd
from jarvis.db.figshare import data as jdata
from jarvis.analysis.elastic.tensor import ElasticTensor
from jarvis.core.atoms import Atoms
from jarvis.core.specie import Specie

from sklearn.metrics import mean_absolute_error, r2_score
from math import isnan
from math import pi

run = 'run11'
label = 'target'

hbar = 1.0545718E-34
kB = 1.38064852E-23


def mass_difference(p):
    elements = p['atoms']['elements']
    mass = []
    for sp in elements:
        mass.append(Specie(sp).atomic_mass)
    mdiff = max(mass) - min(mass)
    return mdiff


with open('../../' + run + '/temp/multi_out_predictions.json') as json_file:
    dos_dict = json.load(json_file)

edos_pdos = jdata("edos_pdos")
dft_3d = jdata("dft_3d")

norm_dos_dict = pint.transform_normalized_dos(edos_pdos, dos_dict, dos_label = 'target')

norm_dos_dict = pint.transform_normalized_dos(edos_pdos, norm_dos_dict, dos_label = 'predictions')

jid_list = []
debyeT_list = []
mdiff_list = []
Cp_target = []
Cp_pred = []

gamma_iso_target = []
gamma_iso_pred = []

'''
Largest mass difference (heaviest atom - lightest atom)
'''


'''
Speed of sound -- might not be so different from Debye T?
'''

for i in dos_dict:
    jid = i["id"]
    start = jid.find('JVASP')
    end = jid.find('.vasp')
    jid = jid[start:end]
    jid_list.append(jid)
    target = np.array(i['target'])
    prediction = np.array(i['predictions'])
    freq = np.linspace(0, 1000, len(target))
    match = next(i for i in dft_3d if i["jid"] == jid)
    atoms = Atoms.from_dict(match['atoms'])
    et = ElasticTensor(match['elastic_tensor'])
    debyeT = et.debye_temperature_toberer(atoms)
    debyeT_list.append(debyeT)
    mdiff = mass_difference(match)
    mdiff_list.append(mdiff)
    if not isnan(debyeT):
        Cp_target.append(pint.heat_capacity(freq, target, debyeT / 2))
        Cp_pred.append(pint.heat_capacity(freq, prediction, debyeT / 2))
    else:
        Cp_target.append(pint.heat_capacity(freq, target,300))
        Cp_pred.append(pint.heat_capacity(freq, prediction, 300))
    try:
        gamma_iso_target.append(pint.isotopic_tau(match, freq, target))
        gamma_iso_pred.append(pint.isotopic_tau(match, freq, prediction))
    except:
        gamma_iso_target.append(0)
        gamma_iso_pred.append(0)


output = {'JID' : jid_list, 'Cp Target' : Cp_target, 'Cp Prediction' : Cp_pred,\
          'DebyeT' : debyeT_list, 'Mass Difference' : mdiff_list, 'Isotope Gamma Target' : gamma_iso_target,\
              'Isotope Gamma Prediction' : gamma_iso_pred}

output_df = pd.DataFrame(output)

output_df.to_csv(run + 'Cp_at_half_debyeT.csv')





