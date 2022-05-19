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

run = 'run21'

#Open Target and Predicted Direct ALIGNN Predicted Thermal Properties


#Open the Target and Predicted DOS thermal properties
with open('output_files/{}_thermal_props_scale_3N.json'.format(run)) as therm_file:
    therm_dict = json.load(therm_file)

prop_list = ['S_vib (J/mol/K)', 'Cp (J/mol/K)']

def mean_absolute_deviation(data, axis=None):
    return np.mean(np.absolute(data - np.mean(data, axis)), axis)


#Open the Debye model DOS properties
with open('output_files/thermal_props_debye.json') as therm_debye_file:
    therm_debye_dict = json.load(therm_debye_file)


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
Compare target distribution to the debye approximation
'''    

bvK_Cv = np.array(therm_debye_dict['Cv_bvk'])
debye_Cv = np.array(therm_debye_dict['Cv'])
target_Cv = np.array(therm_dict['Cp (J/mol/K)']['target'])

bvK_Cv = bvK_Cv[~np.isnan(debye_Cv)]
target_Cv = target_Cv[~np.isnan(debye_Cv)]
debye_Cv = debye_Cv[~np.isnan(debye_Cv)]

mae_debye_Cv = mean_absolute_error(target_Cv, debye_Cv)
mae_bvK_Cv = mean_absolute_error(target_Cv, bvK_Cv)

r2_debye_Cv = r2_score(target_Cv, debye_Cv)
r2_bvK_Cv = r2_score(target_Cv, bvK_Cv)

mad_debye_Cv = mean_absolute_deviation(target_Cv)
mad_bvK_Cv = mean_absolute_deviation(target_Cv)

bvK_Svib = np.array(therm_debye_dict['Svib_bvk'])
debye_Svib = np.array(therm_debye_dict['Svib'])
target_Svib = np.array(therm_dict['S_vib (J/mol/K)']['target'])

bvK_Svib = bvK_Svib[~np.isnan(debye_Svib)]
target_Svib = target_Svib[~np.isnan(debye_Svib)]
debye_Svib = debye_Svib[~np.isnan(debye_Svib)]

mae_debye_Svib = mean_absolute_error(target_Svib, debye_Svib)
mae_bvK_Svib = mean_absolute_error(target_Svib, bvK_Svib)

r2_debye_Svib = r2_score(target_Svib, debye_Svib)
r2_bvK_Svib = r2_score(target_Svib, bvK_Svib)

mad_debye_Svib = mean_absolute_deviation(target_Svib)
mad_bvK_Svib = mean_absolute_deviation(target_Svib)

'''
Plots of Debye and Born von Karman versus target distribution
'''

from jarvis.db.figshare import data as jdata

edos_pdos = jdata("edos_pdos")
dft_3d = jdata("dft_3d")

jid = 'JVASP-32'

match = next(i for i in dft_3d if i["jid"] == jid)
match_pdos = next(i for i in edos_pdos if i["jid"] == jid)

dos = np.array(match_pdos['pdos_elast'])
freq = np.linspace(0, 1000, len(dos))


int_target = pint.integrate_dos(freq,dos)
form_unit = pint.get_natoms_form_unit(match)
scale = int_target / form_unit / 3
dos = dos / scale

debye_dos = pint.debye_DOS(match, freq)
debye_k, debye_omega = pint.debye_dispersion(match)

bvk_dos = pint.BvK_DOS_2(match,freq)
bvk_k, bvk_omega = pint.BvK_dispersion(match)

fig, ax = plt.subplots(1, 3, figsize = (8, 3))
fig.tight_layout(w_pad = 2)

plt.subplot(1, 3, 1)
plt.plot(debye_k, debye_omega, 'xkcd:medium blue')
plt.plot(bvk_k, bvk_omega, 'xkcd:blood red')
fig.text(0.17, 0.7, 'Debye', va='center', fontsize = 14, color = 'xkcd:medium blue')
fig.text(0.22, 0.46, 'BvK', va='center', fontsize = 14, color = 'xkcd:blood red')
plt.xticks([0, 1], ['0', r'k$_{\mathrm{max}}$'])
plt.xlim([0, 1])
plt.ylim([0, 25])
plt.ylabel('Frequency (THz)')

plt.subplot(1, 3, 2)
plt.plot(freq, dos, 'xkcd:black')
plt.plot(freq, debye_dos, 'xkcd:medium blue')
plt.ylim(0, 0.07)
plt.xlabel(r'Frequency (cm$^{-1}$)')
plt.ylabel('Density of States')

plt.subplot(1, 3, 3)
plt.plot(freq, dos, 'xkcd:black')
plt.plot(freq, bvk_dos, 'xkcd:blood red')
plt.xlabel(r'Frequency (cm$^{-1}$)')
plt.ylim(0, 0.1)
plt.savefig('bvk_debye_comp.pdf', bbox_inches = 'tight')

'''
Transparent DOS Schematics
'''
plt.figure()
plt.plot(freq, dos, 'xkcd:black', linewidth = 4)
plt.plot(freq, debye_dos, 'xkcd:blood red', linewidth = 4)
plt.yticks([])
plt.xticks([])
plt.savefig('debye_disp_schem.pdf', bbox_inches = 'tight', transparent = True)

plt.figure()
plt.plot(freq, dos, 'xkcd:black', linewidth = 4)
plt.yticks([])
plt.xticks([])
plt.savefig('disp_schem.pdf', bbox_inches = 'tight', transparent = True)


#Longer frequency range
long_dos = np.pad(dos, [40,20])
long_freq = np.linspace(-800, 1400, len(long_dos))

fig = plt.figure(frameon = False)
ax = fig.add_axes([0, 0, 1, 1])
ax.axis('off')
plt.plot(long_freq, long_dos, 'xkcd:black', linewidth = 4)
plt.yticks([])
plt.xticks([])
plt.savefig('long_disp_schem.pdf', bbox_inches = 'tight', transparent = True)

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

plt.subplot(2, 2, 3)    

histo_svib = np.histogram2d(therm_dict[label]['target'], therm_dict[label]['prediction'], bins = 100)

cx = np.digitize(therm_dict[label]['target'], histo_svib[1])
cy = np.digitize(therm_dict[label]['prediction'], histo_svib[2])

cx = np.where(cx==101, 100, cx) - 1
cy = np.where(cy==101, 100, cy) - 1
pairs = [(x,y) for x,y in zip(cx, cy)]

c_svib = []

for p in pairs:
    c_svib.append(histo_svib[0][p])
    

plt.scatter(therm_dict[label]['target'], therm_dict[label]['prediction'], s = 3, c= c_svib, alpha=0.5, cmap = 'Spectral_r', vmax = 80)

x = np.linspace(0,4000,10)
y = np.linspace(0,4000,10)


plt.plot(x, y, linewidth=1, color='xkcd:black')

plt.plot(x, y - (np.quantile(therm_dict[label]['target'], 0.5) - np.quantile(therm_dict[label]['target'], 0.25)), linewidth=1, linestyle='--', color='xkcd:black')

plt.plot(x, y + (np.quantile(therm_dict[label]['target'], 0.75) - np.quantile(therm_dict[label]['target'], 0.5)), linewidth=1, linestyle='--', color='xkcd:black')


plt.xlim([0, 4000])
plt.ylim([0, 4000])

plt.ylabel(r'Predicted RT S$_{\mathrm{vib}}$ (J/mol/K)', fontsize = 9)
plt.xlabel(r'Target RT S$_{\mathrm{vib}}$ (J/mol/K)', fontsize = 9)
plt.colorbar(pad=0.01)


'''
HEAT CAPACITY -- Room Temperature
'''

plt.subplot(2,2,1)

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

plt.scatter(therm_dict[label]['target'], therm_dict[label]['prediction'], s = 3, c= c_cv, alpha=0.5, cmap = 'Spectral_r', vmax = 80)

x = np.linspace(0,300,10)
y = np.linspace(0,300,10)

plt.plot(x, y, linewidth=1, color='xkcd:black')

plt.plot(x, y - (np.quantile(therm_dict[label]['target'], 0.5) - np.quantile(therm_dict[label]['target'], 0.25)), linewidth=1, linestyle='--', color='xkcd:black')

plt.plot(x, y + (np.quantile(therm_dict[label]['target'], 0.75) - np.quantile(therm_dict[label]['target'], 0.5)), linewidth=1, linestyle='--', color='xkcd:black')

plt.ylim([0, 300])
plt.xlim([0,300])


plt.ylabel(r'Predicted RT C$_{\mathrm{V}}$ (J/mol/K)', fontsize = 9)
plt.xlabel(r'Target RT C$_{\mathrm{V}}$ (J/mol/K)', fontsize = 9)

plt.colorbar(pad=0.01)


'''
HEAT CAPACITY -- Half Debye Temperature
'''

plt.subplot(2,2,2)

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


plt.scatter(target, prediction, s = 3, c= c_cv_dT, alpha=0.5, cmap = 'Spectral_r', vmax = 80)

x = np.linspace(0,300,10)
y = np.linspace(0,300,10)

plt.plot(x, y, linewidth=1, color='xkcd:black')

plt.plot(x, y - (np.quantile(target, 0.5) - np.quantile(target, 0.25)), linestyle='--', linewidth=1, color='xkcd:black')

plt.plot(x, y + (np.quantile(target, 0.75) - np.quantile(target, 0.5)), linestyle='--', linewidth=1, color='xkcd:black')

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

plt.subplot(2,2,4)

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


plt.scatter(target, prediction, s = 3, c= c_iso, alpha=0.5, cmap = 'Spectral_r', vmax = 80)

x = np.linspace(0,50,10)
y = np.linspace(0,50,10)

plt.plot(x, y, linewidth=1, color='xkcd:black')

plt.plot(x, y - (np.quantile(target, 0.5) - np.quantile(target, 0.25)), linestyle='--', linewidth=1, color='xkcd:black')

plt.plot(x, y + (np.quantile(target, 0.75) - np.quantile(target, 0.5)), linestyle='--', linewidth=1, color='xkcd:black')

plt.ylim([0, 50])
plt.xlim([0,50])


plt.ylabel(r'Predicted $\tau^{-1}_{\mathrm{i}}$ (GHz)', fontsize = 9)
plt.xlabel(r'Target $\tau^{-1}_{\mathrm{i}}$ (GHz)', fontsize = 9)

plt.colorbar(pad=0.01)
plt.savefig('figures/{}_therm_props_multipanel_mockup.pdf'.format(run), bbox_inches = 'tight')

'''
Histograms
'''
plt.figure()

fig, ax = plt.subplots(1, 4, constrained_layout = True, figsize = (8, 2))
fig.tight_layout()

plt.subplot(1,4,3)
#fig.text(0.21, 0.25, '(e)', va='center', fontsize = 18)

plt.hist(therm_dict['S_vib (J/mol/K)']['target'], 50, color = 'xkcd:bluey green')
plt.xlabel(r'RT S$_{\mathrm{vib}}$ (J/mol/K)', fontsize = 12)

plt.subplot(1,4,4)
plt.hist(np.array(therm_dict['iso_scattering (Hz)']['target']) / 1e9, 200, color = 'xkcd:bluey green')
#plt.ylabel('Counts', fontsize = 12)
plt.xlabel(r'$\tau^{-1}_{\mathrm{i}}$ (GHz)', fontsize = 12)
plt.xlim([0,50])

plt.subplot(1,4,1)
plt.hist(therm_dict['Cp (J/mol/K)']['target'], 50, color = 'xkcd:bluey green')
plt.ylabel('Counts', fontsize = 12)
plt.xlabel(r'RT C$_{\mathrm{V}}$ (J/mol/K)', fontsize = 12)
#fig.text(0.45, 0.25, '(f)', va='center', fontsize = 18)

#Label peaks at 3NR



plt.subplot(1,4,2)
plt.hist(therm_dict['Cp (J/mol/K) DebT']['target'], 50, color = 'xkcd:bluey green')
plt.xlabel(r'0.5$\theta_{\mathrm{D}}$ C$_{\mathrm{V}}$ (J/mol/K)', fontsize = 12)


plt.savefig('figures/{}_therm_props_histograms.pdf'.format(run), bbox_inches = 'tight')


'''
Fractions of Debye Temperature
'''

mae_cv_list = np.load('output_files/{}_mae_cv_debyeT_list.npy'.format(run))
mae_svib_list = np.load('output_files/{}_mae_svib_debyeT_list.npy'.format(run))
debT_list = np.load('output_files/{}_debyeT_list.npy'.format(run))

frac = np.arange(0.1, 1.1, 0.1)

plt.figure()

fig, ax = plt.subplots(1, 2, gridspec_kw={'width_ratios' : [3,1]}, figsize = (8, 2))
fig.tight_layout(w_pad = 7.5)

plt.subplot(1,2,1)

ax2 = ax[0].twinx()

ax[0].scatter(frac, mae_cv_list, color = 'xkcd:black', s = 10, clip_on = False)
ax[0].plot(frac, mae_cv_list, color = 'xkcd:black')
plt.xlim([0,1])
ax[0].set_ylim([0, 3])

plt.xlabel(r'Fraction of $\theta_{\mathrm{D}}$')
ax[0].set_ylabel(r'MAE in C$_{\mathrm{V}}$ (J/mol/K)')

ax2.scatter(frac, mae_svib_list, color = 'xkcd:blood red', s = 10, clip_on = False)
ax2.plot(frac, mae_svib_list, color = 'xkcd:blood red')
ax2.set_ylim([0, 30])


ax[0].set_xlabel(r'Fraction of $\theta_{\mathrm{D}}$')
ax2.set_ylabel(r'MAE in S$_{\mathrm{vib}}$ (J/mol/K)', color = 'xkcd:blood red')
ax2.tick_params(axis='y', labelcolor='xkcd:blood red')

plt.subplot(1,2,2)

plt.hist(debT_list, 50, color = 'xkcd:bluey green')
plt.xlabel(r'$\theta_{\mathrm{D}}$ (K)')
plt.ylabel('Counts')

plt.savefig('figures/{}_Cv_MAE_T.pdf'.format(run), bbox_inches = 'tight')




