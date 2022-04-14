#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 13:22:12 2022

@author: rlg3


Thermal Properties from analytic DOS
"""

import json
import numpy as np
from jarvis.db.figshare import data as jdata
import pdos_integration as pint

#Open file just to get the JVASP IDs

input_file = '/temp/multi_out_predictions.json'

run = 'run11'

with open('../../' + run + input_file) as json_file:
    dos_dict = json.load(json_file)

dft_3d = jdata("dft_3d")

freq = np.linspace(0, 1000, 200)

jid_list = []
Cv_debye_list = []
Svib_debye_list = []

Cv_bvk_list = []
Svib_bvk_list = []

for d in dos_dict:
    jid = d["id"]
    start = jid.find('JVASP')
    end = jid.find('.vasp')
    jid = jid[start:end]
    jid_list.append(jid)
    match = next(d for d in dft_3d if d["jid"] == jid)
    #Debye Approximation
    debye_dos = pint.debye_DOS(match, freq)
    Cv_debye_list.append(pint.heat_capacity(freq, debye_dos))
    Svib_debye_list.append(pint.vibrational_entropy(freq, debye_dos))
    #Born von Karman Approximation
    try:
        bvk_dos = pint.BvK_DOS(match, freq)
        Cv_bvk_list.append(pint.heat_capacity(freq, bvk_dos))
        Svib_bvk_list.append(pint.vibrational_entropy(freq, bvk_dos))
    except:
        Cv_bvk_list.append(0)
        Svib_bvk_list.append(0)
    
therm_output = {'JID' : jid_list,\
                'Cv' : Cv_debye_list,\
                    'Svib' : Svib_debye_list,\
                        'Cv_bvk' : Cv_bvk_list,\
                            'Svib_bvk' : Svib_bvk_list}
    
with open('output_files/thermal_props_debye.json', 'w') as out_file:
    json.dump(therm_output, out_file)