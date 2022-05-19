#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 10:24:10 2022

@author: rlg3
"""

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
import numpy as np
import json
from pyvalem.formula import Formula
from jarvis.core.atoms import Atoms


run = '21'

direct_Cv_test = pd.read_csv('../../run{}/run_{}_Cv_output/prediction_results_test_set.csv'.format(run, run))

direct_Svib_test = pd.read_csv('../../run{}/run_{}_Svib_output/prediction_results_test_set.csv'.format(run, run))


# plt.figure()
# plt.scatter(direct_Cv_test['target'], direct_Cv_test['prediction'])

Cv_r2_score = r2_score(direct_Cv_test['target'], direct_Cv_test['prediction'])


label = 'Cp (J/mol/K)'

plt.figure()

histo_cv = np.histogram2d(direct_Cv_test['target'], direct_Cv_test['prediction'], bins = 100)

cx = np.digitize(direct_Cv_test['target'], histo_cv[1])
cy = np.digitize(direct_Cv_test['prediction'], histo_cv[2])

cx = np.where(cx==101, 100, cx) - 1
cy = np.where(cy==101, 100, cy) - 1
pairs = [(x,y) for x,y in zip(cx, cy)]

c_cv = []

for p in pairs:
    c_cv.append(histo_cv[0][p])

plt.scatter(direct_Cv_test['target'], direct_Cv_test['prediction'], s = 5, c= c_cv, alpha=0.5, cmap = 'Spectral_r', vmax = 50)

x = np.linspace(0,300,10)
y = np.linspace(0,300,10)

plt.plot(x, y, linewidth=1, color='xkcd:black')

plt.plot(x, y - (np.quantile(direct_Cv_test['target'], 0.5) - np.quantile(direct_Cv_test['target'], 0.25)), linewidth=1, linestyle='--', color='xkcd:black')

plt.plot(x, y + (np.quantile(direct_Cv_test['target'], 0.75) - np.quantile(direct_Cv_test['target'], 0.5)), linewidth=1, linestyle='--', color='xkcd:black')

plt.ylim([0, 300])
plt.xlim([0,300])


plt.ylabel(r'Predicted RT C$_{\mathrm{V}}$ (J/mol/K)', fontsize = 12)
plt.xlabel(r'Target RT C$_{\mathrm{V}}$ (J/mol/K)', fontsize = 12)

plt.colorbar(pad=0.01, label = 'Sample Count')

plt.savefig('direct_alignn_cv.pdf', bbox_inches = 'tight')



plt.figure()
plt.hist(direct_Cv_test['target'], 50, color = 'xkcd:bluey green')
plt.ylabel('Counts', fontsize = 12)
plt.xlabel(r'RT C$_{\mathrm{V}}$ (J/mol/K)', fontsize = 12)

plt.savefig('direct_Cv_target.pdf', bbox_inches = 'tight')

plt.figure()
plt.hist(direct_Cv_test['prediction'], 50, color = 'xkcd:bluey green')
plt.ylabel('Counts', fontsize = 12)
plt.xlabel(r'RT C$_{\mathrm{V}}$ (J/mol/K)', fontsize = 12)

plt.savefig('direct_Cv_prediction.pdf', bbox_inches = 'tight')


plt.figure()
#plt.scatter(direct_Svib_test['target'], direct_Svib_test['prediction'])
Svib_r2_score = r2_score(direct_Svib_test['target'], direct_Svib_test['prediction'])

label = 'S_vib (J/mol/K)'
  

histo_svib = np.histogram2d(direct_Svib_test['target'], direct_Svib_test['prediction'], bins = 100)

cx = np.digitize(direct_Svib_test['target'], histo_svib[1])
cy = np.digitize(direct_Svib_test['prediction'], histo_svib[2])

cx = np.where(cx==101, 100, cx) - 1
cy = np.where(cy==101, 100, cy) - 1
pairs = [(x,y) for x,y in zip(cx, cy)]

c_svib = []

for p in pairs:
    c_svib.append(histo_svib[0][p])
    

plt.scatter(direct_Svib_test['target'], direct_Svib_test['prediction'], s = 3, c= c_svib, alpha=0.5, cmap = 'Spectral_r', vmax = 80)

x = np.linspace(0,4000,10)
y = np.linspace(0,4000,10)


plt.plot(x, y, linewidth=1, color='xkcd:black')

plt.plot(x, y - (np.quantile(direct_Svib_test['target'], 0.5) - np.quantile(direct_Svib_test['target'], 0.25)), linewidth=1, linestyle='--', color='xkcd:black')

plt.plot(x, y + (np.quantile(direct_Svib_test['target'], 0.75) - np.quantile(direct_Svib_test['target'], 0.5)), linewidth=1, linestyle='--', color='xkcd:black')


plt.xlim([0, 4000])
plt.ylim([0, 4000])

plt.ylabel(r'Predicted RT S$_{\mathrm{vib}}$ (J/mol/K)', fontsize = 9)
plt.xlabel(r'Target RT S$_{\mathrm{vib}}$ (J/mol/K)', fontsize = 9)
plt.colorbar(pad=0.01)

plt.savefig('direct_alignn_svib.pdf', bbox_inches = 'tight')


plt.figure()
plt.hist(direct_Svib_test['target'], 50, color = 'xkcd:bluey green')
plt.ylabel('Counts', fontsize = 12)
plt.xlabel(r'RT S$_{\mathrm{vib}}$ (J/mol/K)', fontsize = 12)

plt.savefig('direct_Svib_target.pdf', bbox_inches = 'tight')

plt.figure()
plt.hist(direct_Svib_test['prediction'], 50, color = 'xkcd:bluey green')
plt.ylabel('Counts', fontsize = 12)
plt.xlabel(r'RT S$_{\mathrm{vib}}$ (J/mol/K)', fontsize = 12)

plt.savefig('direct_Svib_prediction.pdf', bbox_inches = 'tight')



'''
Plot High-Residual Spectra for Heat Capacity
'''

Cv_resids = np.abs(np.array(direct_Cv_test['target']) - np.array(direct_Cv_test['prediction']))

Cv_resid_max = np.argsort(Cv_resids)[-6:]


dos_file = '../../run{}/PhononData-Freq-300_1000_20.json'.format(run)

with open(dos_file) as in_file:
    dos_dict = json.load(in_file)



freq = np.linspace(-300, 1000, len(dos_dict[0]['pdos_elast']))    

fig, ax = plt.subplots(2, 3, figsize = (9, 4))
fig.tight_layout(h_pad = 2)

n = 1
for indx in Cv_resid_max:
    plt.subplot(2, 3, n)
    jid = direct_Cv_test.iloc[indx]['id']
    start = jid.find('JVASP')
    end = jid.find('.vasp')
    jid = jid[start:end]
    match = next(d for d in dos_dict if d["jid"] == jid)
    plt.plot(freq, match['pdos_elast'], color = 'xkcd:black')
#    plt.plot(freq, match['predictions'], color = 'xkcd:purpley blue', alpha = 0.6)
    atoms = Atoms.from_dict(match['atoms'])
    composition = atoms.composition.reduced_formula
    f = Formula(composition)
    title_str = '$' + f.latex + '$' + '; JID:' + jid
    plt.title(title_str)    
    n = n+1
    
fig.text(-0.02, 0.5, 'Scaled DOS (a.u.)', va='center', rotation='vertical', fontsize = 14)
fig.text(0.5, -0.02, r'Frequency (cm$^{-1}$)', ha='center', fontsize = 14)

plt.savefig('figures/{}_direct_Cv_high_residual.pdf'.format(run), bbox_inches = 'tight')

