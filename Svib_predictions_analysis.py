#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 09:34:28 2022

@author: rlg3


Plot Svib spectra with large residuals
"""
import json
import numpy as np
import matplotlib.pyplot as plt
from pyvalem.formula import Formula
from jarvis.core.atoms import Atoms

run = 'run21'

with open('output_files/{}_thermal_props_scale_3N_mod.json'.format(run)) as therm_file:
    therm_dict = json.load(therm_file)
    
    
Svib = therm_dict['Cp (J/mol/K)']


resid = np.abs(np.array(Svib['target']) - np.array(Svib['prediction']))


resid_max = np.argsort(resid)[-6:]

dos_file = '../../{}/true_stable_predictions.json'.format(run)



with open(dos_file) as json_file:
    dos_dict = json.load(json_file)
    
freq = np.linspace(-300, 1000, len(dos_dict[0]['target']))    

fig, ax = plt.subplots(2, 3, figsize = (9, 4))
fig.tight_layout(h_pad = 2)

n = 1
for indx in resid_max:
    plt.subplot(2, 3, n)
    plt.plot(freq, dos_dict[indx]['target'], color = 'xkcd:black')
    plt.plot(freq, dos_dict[indx]['predictions'], color = 'xkcd:purpley blue', alpha = 0.6)
    f = Formula(dos_dict[indx]['composition'])
    jid = dos_dict[indx]["id"]
    start = jid.find('JVASP')
    end = jid.find('.vasp')
    jid = jid[start:end]
    title_str = '$' + f.latex + '$' + '; JID:' + jid
    plt.title(title_str)    
    n = n+1
    
fig.text(-0.02, 0.5, 'Scaled DOS (a.u.)', va='center', rotation='vertical', fontsize = 14)
fig.text(0.5, -0.02, r'Frequency (cm$^{-1}$)', ha='center', fontsize = 14)

plt.savefig('figures/{}_cv_high_residual.pdf'.format(run), bbox_inches = 'tight')




