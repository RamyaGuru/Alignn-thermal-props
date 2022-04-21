#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 17:22:30 2022

@author: rlg3


ALIGNN Stability Prediction

Current Criteria for dynamic instability:
    Presence of Scaled DOS intensity greater than 0.1
"""

import json
import numpy as np
from sklearn.metrics import confusion_matrix

run = 'run21'

in_file = '../../{}/predictions_augmented.json'.format(run)

with open(in_file) as json_file:
    dos_dict = json.load(json_file)

freq = np.linspace(-300, 1000, len(dos_dict[0]['target']))


def label_dynamic_stability(freq, dos_dict):
    target_list = []
    pred_list = []
    zero_indx = np.where(freq > 0)[0][0] - 1
    for d in dos_dict:
        if any(np.array(d['target'][:zero_indx]) > 0.1):
            d['target_stable'] = 0
        else:
            d['target_stable'] = 1
        if any(np.array(d['predictions'][:zero_indx]) > 0.1):
            d['pred_stable'] = 0
        else:
            d['pred_stable'] = 1
        target_list.append(d['target_stable'])
        pred_list.append(d['pred_stable'])
    return target_list, pred_list
        
    

target_labels, pred_labels = label_dynamic_stability(freq, dos_dict)

C = confusion_matrix(target_labels, pred_labels)
