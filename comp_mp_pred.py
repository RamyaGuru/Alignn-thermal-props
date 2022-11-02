#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 06:57:46 2022

@author: rlg3

Compare ALIGNN and MP Results
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import json
from sklearn.metrics import mean_absolute_error, r2_score



def mean_absolute_deviation(data, axis=None):
    return np.mean(np.absolute(data - np.mean(data, axis)), axis)


'''
Load Data from CSV files
'''

mp_json = 'output_files/run21_thermal_props_mp_bins_mod.json'

pred_json = 'output_files/run21_thermal_props_dft3d_mod.json'

with open(mp_json) as infile:
    mp_dict = json.load(infile)
    
with open(pred_json) as jfile:
    pred_dict = json.load(jfile)

mp_df = pd.DataFrame.from_dict(mp_dict)

pred_df = pd.DataFrame.from_dict(pred_dict)

mp_mpid = np.array(mp_df['mpid'])

'''
Sort and Subset Prediction Dataset to Match MP Dataset
'''
pred_trim = pred_df.loc[pred_df['mpid'].isin(mp_mpid)]
pred_mpid = np.array(pred_trim['mpid'])
pred_trim_avg = pred_trim.groupby(['mpid']).mean()

mp_trim = mp_df.loc[mp_df['mpid'].isin(pred_mpid)]
mp_trim.set_index('mpid', inplace = True)

#Sort Dataframes first
mp_trim = mp_trim.sort_values(by = 'mpid')
pred_trim_avg = pred_trim_avg.sort_values(by = 'mpid')


'''
Drop scattering rates of zero from either dataset
'''
bool_1 = ~np.isnan(mp_trim['isotope_scatt'])
bool_2 = mp_trim['isotope_scatt'] != 0
bool_3 = ~np.isnan(pred_trim_avg['isotope_scatt'])
bool_4 = pred_trim_avg['isotope_scatt'] !=0

bool_full = bool_1 & bool_2 & bool_3 & bool_4

mp_iso = mp_trim.loc[bool_full]
pred_iso = pred_trim_avg.loc[bool_full] 

'''
Plot Histograms of Values
'''
plt.figure()

fig, ax = plt.subplots(1, 3, constrained_layout = True, figsize = (9, 2.5))
plt.tight_layout(h_pad = 1)

plt.subplot(2,3,1)
plt.title('C$_{\mathrm{V}}$ (J/mol/K)')
plt.ylabel('Counts')
plt.hist(mp_trim['Cv_mol'], 50, color = 'xkcd:black')
#plt.savefig('figures/mp_molar_Cv.pdf', bbox_inches = 'tight')
plt.hist(pred_trim_avg['Cv_mol'], 50, alpha = 0.6, color = 'xkcd:bluey green')
#plt.savefig('figures/pred_molar_Cv.pdf', bbox_inches = 'tight')

plt.subplot(2,3,2)
plt.title('S$_{\mathrm{vib}}$ (J/mol/K)')
plt.hist(mp_trim['Svib_mol'], 50, color = 'xkcd:black')
#plt.savefig('figures/mp_molar_Svib.pdf', bbox_inches = 'tight')

# plt.subplot(2,3,5)
# plt.title('ALIGNN Molar Svib')
plt.hist(pred_trim_avg['Svib_mol'], 50, alpha = 0.6, color = 'xkcd:bluey green')
#plt.savefig('figures/pred_molar_Svib.pdf', bbox_inches = 'tight')

plt.subplot(2,3,3)
plt.title(r'$\tau^{-1}_{\mathrm{i}}$ (GHz)')
plt.hist(mp_iso['isotope_scatt'] / 1E9, 50, color = 'xkcd:black', label = 'MP')
plt.xlim(0, 500)
#plt.savefig('figures/mp_molar_iso.pdf', bbox_inches = 'tight')

# plt.subplot(2,3,6)
# plt.title('ALIGNN Isotope Scattering')
plt.hist(pred_iso['isotope_scatt'] / 1E9, 50, alpha = 0.6, color = 'xkcd:bluey green', label = 'ALIGNN')
plt.xlim(0,500)
plt.legend()
plt.savefig('figures/pred_mp_prop_histograms_binned.pdf', bbox_inches = 'tight')



'''
Comparison between MP and ALIGNN Data
'''



plt.figure()

'''
Molar Heat Capacity
'''

fig, ax = plt.subplots(1, 3, constrained_layout = True, figsize = (9, 2.5))
plt.tight_layout(w_pad = 2)

plt.subplot(1,3,1)
plt.scatter(mp_trim['Cv_mol'], pred_trim_avg['Cv_mol'], alpha = 0.5, s = 3, c = 'xkcd:warm blue')

x = np.linspace(0,250,10)
y = np.linspace(0,250,10)


plt.plot(x, y, linewidth=1, color='xkcd:black')

plt.plot(x, y - (np.quantile(mp_trim['Cv_mol'], 0.5) - np.quantile(mp_trim['Cv_mol'], 0.25)), linestyle='--', linewidth=1, color='xkcd:black')

plt.plot(x, y + (np.quantile(mp_trim['Cv_mol'], 0.75) - np.quantile(mp_trim['Cv_mol'], 0.5)), linestyle='--', linewidth=1, color='xkcd:black')

plt.xlim([0, 250])
plt.ylim([0, 250])

plt.xlabel(r'MP C$_{\mathrm{V}}$ (J/mol/K)')
plt.ylabel(r'ALIGNN C$_{\mathrm{V}}$ (J/mol/K)')

'''
Molar Vibraitonal Entropy
'''

plt.subplot(1,3,2)
plt.scatter(mp_trim['Svib_mol'], pred_trim_avg['Svib_mol'], alpha = 0.5, s = 3, c='xkcd:warm blue')

x = np.linspace(0,600,10)
y = np.linspace(0,600,10)


plt.plot(x, y, linewidth=1, color='xkcd:black')
plt.plot(x, y - (np.quantile(mp_trim['Svib_mol'], 0.5) - np.quantile(mp_trim['Svib_mol'], 0.25)), linestyle='--', linewidth=1, color='xkcd:black')

plt.plot(x, y + (np.quantile(mp_trim['Svib_mol'], 0.75) - np.quantile(mp_trim['Svib_mol'], 0.5)), linestyle='--', linewidth=1, color='xkcd:black')

plt.xlim([0, 600])
plt.ylim([0, 600])

plt.xlabel(r'MP S$_{\mathrm{vib}}$ (J/mol/K)')
plt.ylabel(r'ALIGNN S$_{\mathrm{vib}}$ (J/mol/K)')

'''
Isotope Scattering
'''
plt.subplot(1,3,3)
plt.scatter(mp_iso['isotope_scatt'] / 1e9, pred_iso['isotope_scatt'] / 1e9, alpha = 0.5, s = 3, c = 'xkcd:warm blue')

x = np.linspace(0,500,10)
y = np.linspace(0,500,10)

plt.plot(x, y, linewidth=1, color='xkcd:black')
plt.plot(x, y - (np.quantile(mp_iso['isotope_scatt'] / 1e9, 0.5) - np.quantile(mp_iso['isotope_scatt'] / 1e9, 0.25)), linestyle='--', linewidth=1, color='xkcd:black')

plt.plot(x, y + (np.quantile(mp_iso['isotope_scatt'] / 1e9, 0.75) - np.quantile(mp_iso['isotope_scatt'] / 1e9, 0.5)), linestyle='--', linewidth=1, color='xkcd:black')

plt.xlim([0, 500])
plt.ylim([0, 500])

plt.xlabel(r'MP $\tau^{-1}_{\mathrm{i}}$ (GHz)')
plt.ylabel(r'ALIGNN $\tau^{-1}_{\mathrm{i}}$ (GHz)')

plt.savefig('figures/pred_mp_comp_binned.pdf', bbox_inches = 'tight')


'''
Highest Residual Vibrational Entropy: ALIGNN versus MP Spectra
'''

fig, ax = plt.subplots(2, 3, constrained_layout = True, figsize = (6, 3))
fig.tight_layout()

run = 'run21'

resid = np.abs(np.array(mp_trim['Svib_mol']) - np.array(pred_trim_avg['Svib_mol']))

resid_max = np.argsort(resid)[-6:]

'''
Error Statistics
'''
#MAE
Svib_mae = mean_absolute_error(mp_trim['Svib_mol'], pred_trim_avg['Svib_mol'])
Cv_mae = mean_absolute_error(mp_trim['Cv_mol'], pred_trim_avg['Cv_mol'])
tau_mae = mean_absolute_error(mp_iso['isotope_scatt'], pred_iso['isotope_scatt'])


#MAD
Svib_mad = mean_absolute_deviation(mp_trim['Svib_mol'])
Cv_mad = mean_absolute_deviation(mp_trim['Cv_mol'])
tau_mad = mean_absolute_deviation(mp_iso['isotope_scatt'])


#Ratio

mae_list = np.array([Svib_mae, Cv_mae, tau_mae])
mad_list = np.array([Svib_mad, Cv_mad, tau_mad])

ratio_list = mae_list / mad_list 

#R2 Score
Svib_r2 = r2_score(mp_trim['Svib_mol'], pred_trim_avg['Svib_mol'])
Cv_r2 = r2_score(mp_trim['Cv_mol'], pred_trim_avg['Cv_mol'])
tau_r2 = r2_score(mp_iso['isotope_scatt'], pred_iso['isotope_scatt'])


# Load the DOS spectra

# pred_file = '../../{}/pred_data_augmented.json'.format(run)

# freq_bin = np.linspace(0, 1000, 51)
# n = 1
# for indx in resid_max:
#     plt.subplot(2, 3, n)
#     plt.plot(freq_bin, np.array(mp_trim['scaled_dos'])[indx], color = 'xkcd:black')
#     mpid_temp = np.array(mp_trim['mpid'])[indx]
#     pred_dos = pred_trim.loc[pred_trim['mpid'] == mpid_temp]
#     plt.plot(freq_bin, list(pred_dos['scaled_dos'])[0], color = 'xkcd:purpley blue', alpha = 0.6)
#     n = n+1
    
# fig.text(-0.02, 0.5, 'Scaled DOS (a.u.)', va='center', rotation='vertical', fontsize = 12)
# fig.text(0.5, -0.02, r'Frequency (cm$^{-1}$)', ha='center', fontsize = 12)

# plt.savefig('figures/{}_Svib_high_residual_mp.pdf'.format(run), bbox_inches = 'tight')







