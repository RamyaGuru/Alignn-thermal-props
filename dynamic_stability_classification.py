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
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import matplotlib.pyplot as plt
import matplotlib as mpl

run = 'run21'

in_file = '../../{}/predictions_augmented.json'.format(run)

with open(in_file) as json_file:
    dos_dict = json.load(json_file)

freq = np.linspace(-300, 1000, len(dos_dict[0]['target']))


def label_dynamic_stability(freq, dos_dict):
    target_list = []
    pred_list = []
    true_stable_list = []
    zero_indx = np.where(freq > 0)[0][0] - 1
    neg_freq = np.linspace(-300, 0, zero_indx)
    for d in dos_dict:
        intdos_neg_target = np.trapz(neg_freq, d['target'][:zero_indx])
        intdos_target = np.trapz(freq, d['target'])
        intdos_neg_pred = np.trapz(neg_freq, d['predictions'][:zero_indx])
        intdos_pred = np.trapz(freq, d['predictions'])
        d['target_stable'] = 1
        d['pred_stable'] = 1
        if intdos_neg_target / intdos_target > 0.1:
            d['target_stable'] = 0
        if intdos_neg_pred / intdos_pred > 0.1:
            d['pred_stable'] = 0
        # if any(np.array(d['target'][:zero_indx]) > 0.1):
        #     d['target_stable'] = 0
        # else:
        #     d['target_stable'] = 1
        # if any(np.array(d['predictions'][:zero_indx]) > 0.1):
        #     d['pred_stable'] = 0
        # else:
        #     d['pred_stable'] = 1
        target_list.append(d['target_stable'])
        pred_list.append(d['pred_stable'])
        if d['target_stable'] == 1 and d['pred_stable'] == 1:
            true_stable_list.append(d)
    return target_list, pred_list, true_stable_list
        
    

target_labels, pred_labels, true_stable_list = label_dynamic_stability(freq, dos_dict)

'''
Write True Stable List to json
'''

out_file = '../../{}/true_stable_predictions.json'.format(run)

with open(out_file, 'w') as json_out:
    json.dump(true_stable_list, json_out)

C = confusion_matrix(target_labels, pred_labels)


disp = ConfusionMatrixDisplay(confusion_matrix = C, display_labels = ['Unstable', 'Stable'])
disp.plot(cmap = 'Greys')

plt.savefig('figures/{}_stability_confusion_matrix.pdf'.format(run), bbox_inches = 'tight')


'''
Plot Sample Spectra from each category
'''

for pairs in [[0,0], [1,1], [0, 1], [1, 0]]:
    n = 1
    plt.figure()
    fig, ax = plt.subplots(1, 2, figsize = (4, 2))
    fig.tight_layout(h_pad = 2)
    for d in dos_dict:
        if d['target_stable'] == pairs[0] and d['pred_stable'] == pairs[1]:
            plt.subplot(1, 2, n)
            plt.plot(freq, d['target'],color = 'xkcd:black', linewidth = 2)
            plt.plot(freq, d['predictions'], color = 'xkcd:bluish purple', alpha = 0.6, linewidth = 2)
            n = n+1
        if n == 3:
            break
    #fig.text(-0.02, 0.5, 'Scaled DOS (a.u.)', va='center', rotation='vertical', fontsize = 10)
    fig.text(0.5, -0.02, r'Frequency (cm$^{-1}$)', ha='center', fontsize = 10)
    plt.savefig('figures/{}_example{}{}.pdf'.format(run, pairs[0], pairs[1]), bbox_inches = 'tight')



