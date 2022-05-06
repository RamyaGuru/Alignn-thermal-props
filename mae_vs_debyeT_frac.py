#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 13:42:47 2022

@author: rlg3


Run MAE for mutiple fractions of the Debye Temperature
"""

import json
import numpy as np
import pdos_integration as pint
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error

from jarvis.analysis.elastic.tensor import ElasticTensor
from jarvis.core.atoms import Atoms
from jarvis.core.specie import Specie

from jarvis.db.figshare import data as jdata

#Debye Temperature Fractions
frac = np.arange(0.1, 1.1, 0.1)


input_file = '/multi_out_predictions.json'

run = 'run21'

with open('../../' + run + input_file) as json_file:
    dos_dict = json.load(json_file)

dft_3d = jdata("dft_3d")

mae_list = []
mae_svib_list = []
debyeT_list = []
for f in frac:
    print(f)
    temp_t_cv = []
    temp_p_cv = []
    temp_t_svib = []
    temp_p_svib = []
    for i in dos_dict:
        jid = i["id"]
        start = jid.find('JVASP')
        end = jid.find('.vasp')
        jid = jid[start:end]
        match = next(d for d in dft_3d if d["jid"] == jid)
        # Get number of atoms in formula unit
        #atoms = Atoms.from_dict(p['atoms'])
        form_unit = pint.get_natoms_form_unit(match)
        target = np.array(i['target'])
        pred = np.array(i['predictions'])
        freq = np.linspace(0, 1000, len(target))
        intDOS_t = pint.integrate_dos(freq, target)
        scale = (intDOS_t / form_unit) / 3.0
        target = target / scale
        intDOS_p = pint.integrate_dos(freq, pred)
        scale = (intDOS_p / form_unit) / 3.0
        pred = pred / scale
        atoms = Atoms.from_dict(match['atoms'])
        et = ElasticTensor(match['elastic_tensor'])
        debyeT = et.debye_temperature_toberer(atoms)
        debyeT_list.append(debyeT)
        if not np.isnan(debyeT): 
            temp_t_cv.append(pint.heat_capacity(freq, target, debyeT * f))
            temp_p_cv.append(pint.heat_capacity(freq, pred, debyeT * f))
            temp_t_svib.append(pint.vibrational_entropy(freq, target, debyeT * f))
            temp_p_svib.append(pint.vibrational_entropy(freq,pred, debyeT * f))
    mae_list.append(mean_absolute_error(np.array(temp_t_cv)[~np.isnan(temp_t_cv)], np.array(temp_p_cv)[~np.isnan(temp_t_cv)]))
    mae_svib_list.append(mean_absolute_error(np.array(temp_t_svib)[~np.isnan(temp_t_svib)], np.array(temp_p_svib)[~np.isnan(temp_t_svib)]))

np.save('output_files/{}_debyeT_list.npy'.format(run), debyeT_list)
np.save('output_files/{}_mae_cv_debyeT_list.npy'.format(run), mae_list)
np.save('output_files/{}_mae_svib_debyeT_list.npy'.format(run), mae_svib_list)
