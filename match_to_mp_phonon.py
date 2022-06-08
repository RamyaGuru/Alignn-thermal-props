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
                                                             'metadata': ph_dict['metadata']}


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
    pred_mpids.append(mpid_pd)



overlap = set(pred_mpids) & set(mp_mpid_list)




