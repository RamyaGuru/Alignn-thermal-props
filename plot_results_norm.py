#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 16:51:51 2022

@author: rlg3

1) Error Statistics

2) Plot Results after the normalization to integrated DOS = 3N


***Need to remove zero values before computing error metrics***
"""

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pdos_integration as pint
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from math import sqrt
import json



with open('output_files/thermal_props_scale_3N.json') as therm_file:
    therm_dict = json.load(therm_file)

prop_list = ['S_vib (J/mol/K)', 'Cp (J/mol/K)']

def mean_absolute_deviation(data, axis=None):
    return np.mean(np.absolute(data - np.mean(data, axis)), axis)


'''
Target distribution MAD values
'''
mad_list = []
mae_list = []
r2_list = []

for prop in prop_list:
    mad_list.append(mean_absolute_deviation(therm_dict[prop]['target']))
    mae_list.append(mean_absolute_error(therm_dict[prop]['target'], therm_dict[prop]['prediction']))
    r2_list.append(r2_score(therm_dict[prop]['target'], therm_dict[prop]['prediction']))
    
'''
Remove zeros from 0.5 Debye T heat capacities and get error
'''

prop_list_2 = ['iso_scattering (Hz)', 'Cp (J/mol/K) DebT']

for prop2 in prop_list_2:
    target = np.array(therm_dict[prop2]['target'])#[therm_dict[prop2]['target'] != 0]
    target = target[target != 0]
    target = target[~np.isnan(target)]
    prediction = np.array(therm_dict[prop2]['prediction'])#[therm_dict[prop2]['prediction'] != 0]
    prediction = prediction[prediction != 0]
    prediction = prediction[~np.isnan(prediction)]
    mad_list.append(mean_absolute_deviation(target))
    mae_list.append(mean_absolute_error(target, prediction))
    r2_list.append(r2_score(target, prediction))    
    

'''
Plot of room temperature thermal properties
'''
mpl.rcdefaults()

plt.figure()


fig, ax = plt.subplots(2,2, figsize = (6,5))

plt.tight_layout(h_pad = 2, w_pad = 2)

'''
VIBRATIONAL ENTROPY
'''
label = 'S_vib (J/mol/K)'

plt.subplot(2, 2, 1)    

histo_svib = np.histogram2d(therm_dict[label]['target'], therm_dict[label]['prediction'], bins = 100)

cx = np.digitize(therm_dict[label]['target'], histo_svib[1])
cy = np.digitize(therm_dict[label]['prediction'], histo_svib[2])

cx = np.where(cx==101, 100, cx) - 1
cy = np.where(cy==101, 100, cy) - 1
pairs = [(x,y) for x,y in zip(cx, cy)]

c_svib = []

for p in pairs:
    c_svib.append(histo_svib[0][p])
    

plt.scatter(therm_dict[label]['target'], therm_dict[label]['prediction'], s = 3, c= c_svib, alpha=0.5, cmap = 'Spectral_r', vmax = 100)

x = np.linspace(0,4000,10)
y = np.linspace(0,4000,10)


plt.plot(x, y, linewidth=2, color='xkcd:black')

plt.plot(x, y - (np.quantile(therm_dict[label]['target'], 0.5) - np.quantile(therm_dict[label]['target'], 0.25)), linestyle='--', color='xkcd:black')

plt.plot(x, y + (np.quantile(therm_dict[label]['target'], 0.75) - np.quantile(therm_dict[label]['target'], 0.5)), linestyle='--', color='xkcd:black')


plt.xlim([0, 4000])
plt.ylim([0, 4000])

plt.ylabel(r'Predicted S$_{\mathrm{vib}}$ (J/mol/K)', fontsize = 9)
plt.xlabel(r'Target S$_{\mathrm{vib}}$ (J/mol/K)', fontsize = 9)
plt.colorbar(pad=0.01)


'''
HEAT CAPACITY -- Room Temperature
'''

plt.subplot(2,2,2)

label = 'Cp (J/mol/K)'

histo_cv = np.histogram2d(therm_dict[label]['target'], therm_dict[label]['prediction'], bins = 100)

cx = np.digitize(therm_dict[label]['target'], histo_cv[1])
cy = np.digitize(therm_dict[label]['prediction'], histo_cv[2])

cx = np.where(cx==101, 100, cx) - 1
cy = np.where(cy==101, 100, cy) - 1
pairs = [(x,y) for x,y in zip(cx, cy)]

c_cv = []

for p in pairs:
    c_cv.append(histo_cv[0][p])

plt.scatter(therm_dict[label]['target'], therm_dict[label]['prediction'], s = 3, c= c_cv, alpha=0.5, cmap = 'Spectral_r', vmax = 100)

x = np.linspace(0,300,10)
y = np.linspace(0,300,10)

plt.plot(x, y, linewidth=2, color='xkcd:black')

plt.plot(x, y - (np.quantile(therm_dict[label]['target'], 0.5) - np.quantile(therm_dict[label]['target'], 0.25)), linestyle='--', color='xkcd:black')

plt.plot(x, y + (np.quantile(therm_dict[label]['target'], 0.75) - np.quantile(therm_dict[label]['target'], 0.5)), linestyle='--', color='xkcd:black')

plt.ylim([0, 300])
plt.xlim([0,300])


plt.ylabel(r'Predicted RT C$_{\mathrm{V}}$ (J/mol/K)', fontsize = 9)
plt.xlabel(r'Target RT C$_{\mathrm{V}}$ (J/mol/K)', fontsize = 9)

plt.colorbar(pad=0.01)


'''
HEAT CAPACITY -- Half Debye Temperature
'''

plt.subplot(2,2,4)

label = 'Cp (J/mol/K) DebT'

target = np.array(therm_dict[label]['target'])#[therm_dict[prop2]['target'] != 0]
target = target[target != 0]
target = target[~np.isnan(target)]    

prediction = np.array(therm_dict[label]['prediction'])#[therm_dict[prop2]['prediction'] != 0]
prediction = prediction[prediction != 0]
prediction = prediction[~np.isnan(prediction)]

histo_cv_dT = np.histogram2d(target, prediction, bins = 100)

cx = np.digitize(target, histo_cv_dT[1])
cy = np.digitize(prediction, histo_cv_dT[2])

cx = np.where(cx==101, 100, cx) - 1
cy = np.where(cy==101, 100, cy) - 1
pairs = [(x,y) for x,y in zip(cx, cy)]

c_cv_dT = []

for p in pairs:
    c_cv_dT.append(histo_cv_dT[0][p])


plt.scatter(target, prediction, s = 3, c= c_cv_dT, alpha=0.5, cmap = 'Spectral_r', vmax = 100)

x = np.linspace(0,300,10)
y = np.linspace(0,300,10)

plt.plot(x, y, linewidth=2, color='xkcd:black')

plt.plot(x, y - (np.quantile(target, 0.5) - np.quantile(target, 0.25)), linestyle='--', color='xkcd:black')

plt.plot(x, y + (np.quantile(target, 0.75) - np.quantile(target, 0.5)), linestyle='--', color='xkcd:black')

plt.ylim([0, 300])
plt.xlim([0,300])


plt.ylabel(r'Predicted 0.5$\theta_{\mathrm{D}}$ C$_{\mathrm{V}}$ (J/mol/K)', fontsize = 9)
plt.xlabel(r'Target 0.5$\theta_{\mathrm{D}}$ C$_{\mathrm{V}}$ (J/mol/K)', fontsize = 9)
#cax = plt.axes([0.85, 0.1, 0.075, 0.8])
plt.colorbar(label = 'Sample Count', pad=0.01)


#plt.savefig('figures/new_multipanel_mockup.pdf', bbox_inches = 'tight')
                                                                          

'''
Isotope Scattering
'''               

plt.subplot(2,2,3)

label = 'iso_scattering (Hz)'

target = np.array(therm_dict[label]['target'])#[therm_dict[prop2]['target'] != 0]
target = target[target != 0]
target = target[~np.isnan(target)] / 1e9

prediction = np.array(therm_dict[label]['prediction'])#[therm_dict[prop2]['prediction'] != 0]
prediction = prediction[prediction != 0]
prediction = prediction[~np.isnan(prediction)] / 1e9

histo_iso= np.histogram2d(target, prediction, bins = 100)

cx = np.digitize(target, histo_iso[1])
cy = np.digitize(prediction, histo_iso[2])

cx = np.where(cx==101, 100, cx) - 1
cy = np.where(cy==101, 100, cy) - 1
pairs = [(x,y) for x,y in zip(cx, cy)]

c_iso = []

for p in pairs:
    c_iso.append(histo_iso[0][p])


plt.scatter(target, prediction, s = 3, c= c_iso, alpha=0.5, cmap = 'Spectral_r', vmax = 100)

x = np.linspace(0,100,10)
y = np.linspace(0,100,10)

plt.plot(x, y, linewidth=2, color='xkcd:black')

plt.plot(x, y - (np.quantile(target, 0.5) - np.quantile(target, 0.25)), linestyle='--', color='xkcd:black')

plt.plot(x, y + (np.quantile(target, 0.75) - np.quantile(target, 0.5)), linestyle='--', color='xkcd:black')

plt.ylim([0, 100])
plt.xlim([0,100])


plt.ylabel(r'Predicted $\tau^{-1}_{\mathrm{i}}$ (GHz)', fontsize = 9)
plt.xlabel(r'Target $\tau^{-1}_{\mathrm{i}}$ (GHz)', fontsize = 9)

plt.colorbar(pad=0.01)
plt.savefig('figures/therm_props_multipanel_mockup.pdf', bbox_inches = 'tight')

'''
Histograms
'''
plt.figure()

fig, ax = plt.subplots(1, 4, constrained_layout = True, figsize = (8, 2))
fig.tight_layout()

plt.subplot(1,4,1)
#fig.text(0.21, 0.25, '(e)', va='center', fontsize = 18)

plt.hist(therm_dict['S_vib (J/mol/K)']['target'], 50, color = 'xkcd:bluey green')
plt.ylabel('Counts', fontsize = 12)
plt.xlabel(r'RT S$_{\mathrm{vib}}$ (J/mol/K)', fontsize = 12)

plt.subplot(1,4,2)
plt.hist(np.array(therm_dict['iso_scattering (Hz)']['target']) / 1e9, 50, color = 'xkcd:bluey green')
#plt.ylabel('Counts', fontsize = 12)
plt.xlabel(r'$\tau^{-1}_{\mathrm{i}}$ (GHz)', fontsize = 12)

plt.subplot(1,4,3)
plt.hist(therm_dict['Cp (J/mol/K)']['target'], 50, color = 'xkcd:bluey green')
plt.xlabel(r'RT C$_{\mathrm{V}}$ (J/mol/K)', fontsize = 12)
#fig.text(0.45, 0.25, '(f)', va='center', fontsize = 18)

plt.subplot(1,4,4)
plt.hist(therm_dict['Cp (J/mol/K) DebT']['target'], 50, color = 'xkcd:bluey green')
plt.xlabel(r'0.5$\theta_{\mathrm{D}}$ C$_{\mathrm{V}}$ (J/mol/K)', fontsize = 12)


plt.savefig('figures/therm_props_histograms.pdf', bbox_inches = 'tight')


'''
Fractions of Debye Temperature
'''

mae_cv_list = np.load('output_files/mae_debyeT_list.npy')
mae_svib_list = np.load('output_files/mae_svib_debyeT_list.npy')
debT_list = np.load('output_files/debyeT_list.npy')

frac = np.arange(0.1, 1.1, 0.1)

plt.figure()

fig, ax = plt.subplots(1, 2, gridspec_kw={'width_ratios' : [3,1]}, figsize = (8, 2))
fig.tight_layout(w_pad = 7.5)

plt.subplot(1,2,1)

ax2 = ax[0].twinx()

ax[0].scatter(frac, mae_cv_list, color = 'xkcd:black', s = 10, clip_on = False)
ax[0].plot(frac, mae_cv_list, color = 'xkcd:black')
plt.xlim([0,1])
ax[0].set_ylim([0, 4])

plt.xlabel(r'Fraction of $\theta_{\mathrm{D}}$')
ax[0].set_ylabel(r'MAE in C$_{\mathrm{V}}$ (J/mol/K)')

ax2.scatter(frac, mae_svib_list, color = 'xkcd:blood red', s = 10, clip_on = False)
ax2.plot(frac, mae_svib_list, color = 'xkcd:blood red')
ax2.set_ylim([0, 150])


ax[0].set_xlabel(r'Fraction of $\theta_{\mathrm{D}}$')
ax2.set_ylabel(r'MAE in S$_{\mathrm{vib}}$ (J/mol/K)', color = 'xkcd:blood red')
ax2.tick_params(axis='y', labelcolor='xkcd:blood red')

plt.subplot(1,2,2)

plt.hist(debT_list, 50, color = 'xkcd:bluey green')
plt.xlabel(r'$\theta_{\mathrm{D}}$ (K)')
plt.ylabel('Counts')

plt.savefig('Cv_MAE_T.pdf', bbox_inches = 'tight')




