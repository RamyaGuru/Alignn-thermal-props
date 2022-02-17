#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 00:10:32 2022

@author: rlg3
"""

import pandas as pd
import matplotlib as mpl
mpl.rcdefaults()
import matplotlib.pyplot as plt
import numpy as np
import pdos_integration as pint
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from math import sqrt



def mean_absolute_deviation(data, axis=None):
    return np.mean(np.absolute(data - np.mean(data, axis)), axis)

'''
Load target and predicted dataframe. Plot against one another
'''
mpl.rcdefaults()
mpl.rcParams['font.size'] = 12

label = 'integrated_DOS'

pred_df = pd.read_csv('output_files/run11predictions_thermal_props_scale_2.csv')

target_df = pd.read_csv('output_files/run11target_thermal_props_scale_2.csv')

plt.figure()
plt.scatter(pred_df[label], target_df[label], s = 5, color= 'xkcd:medium blue', alpha=0.4)

x = np.linspace(0,50,10)
y = np.linspace(0,50,10)

plt.plot(x, y, linestyle='--', linewidth=2, color='xkcd:black')

plt.plot(x, y - np.std(target_df[label]), linestyle=':', linewidth=2, color='xkcd:black')

plt.plot(x, y + np.std(target_df[label]), linestyle=':', linewidth=2, color='xkcd:black')

plt.ylim([0, 50])
plt.xlim([0,50])

plt.xlabel(r'ML Integrated DOS', size = 14)
plt.ylabel(r'Target Integrated DOS', size = 14)

plt.savefig('true_vs_pred_intDOS.pdf', bbox_inches='tight')



'''
Plots of heat capacity and vibrational entropy scaling
'''


from jarvis.db.figshare import data as jdata

dft_3d = jdata("edos_pdos")

jid = 'JVASP-32'

match = next(i for i in dft_3d if i["jid"] == jid)


target = np.array(match['pdos_elast'])
freq = np.linspace(0, 1000, len(target))

plt.figure()
plt.plot(freq, target, label='DOS')
plt.plot(freq, pint.heat_capacity_scaling(freq, 300) * 100, color='xkcd:blood red', label=r'C$_{\mathrm{V}}$ scaling')
#plt.plot(freq, np.zeros(len(target)), linestyle=':', color = 'xkcd:grey')

ax = plt.gca()
ax.axes.get_yaxis().set_visible(False)


plt.xlim(0, 1000)
plt.ylabel('DOS (a.u.)')
plt.xlabel(r'Frequency cm$^{-1}$')
plt.legend()

plt.savefig('cv_scaling.pdf', bbox_inches='tight')


'''
Get other error statistics
'''

#R^2 values

r2_intDOS = r2_score(target_df['integrated_DOS'], pred_df['integrated_DOS'])

r2_Cp = r2_score(target_df['Cp (J/mol/K)'], pred_df['Cp (J/mol/K)'])

r2_Svib = r2_score(target_df['S_vib (J/mol/K)'], pred_df['S_vib (J/mol/K)'])


#RMSE values

rmse_intDOS = sqrt(mean_squared_error(target_df['integrated_DOS'], pred_df['integrated_DOS']))

rmse_Cp = sqrt(mean_squared_error(target_df['Cp (J/mol/K)'], pred_df['Cp (J/mol/K)']))

rmse_Svib = sqrt(mean_squared_error(target_df['S_vib (J/mol/K)'], pred_df['S_vib (J/mol/K)']))


#MAD values (Report for target distribution)

mad_intDOS = mean_absolute_deviation(target_df['integrated_DOS'])
mad_Cp = mean_absolute_deviation(target_df['Cp (J/mol/K)'])
mad_Svib = mean_absolute_deviation(target_df['S_vib (J/mol/K)'])


'''
Cv at same fraction of Debye temperature
'''

CV_at_debT = pd.read_csv('run11Cp_at_half_debyeT.csv')

plt.figure()
plt.scatter(CV_at_debT['Cp Target'],CV_at_debT['Cp Prediction'], s = 5, color= 'xkcd:medium blue', alpha=0.4)

x = np.linspace(0,250,10)
y = np.linspace(0,250,10)

plt.plot(x, y, linestyle='--', linewidth=2, color='xkcd:black')

plt.plot(x, y - np.std(CV_at_debT['Cp Target']), linestyle=':', linewidth=2, color='xkcd:black')

plt.plot(x, y + np.std(CV_at_debT['Cp Target']), linestyle=':', linewidth=2, color='xkcd:black')

plt.ylim([0, 250])
plt.xlim([0,250])

plt.xlabel(r'ML C$_{\mathrm{V}}$', size = 14)
plt.ylabel(r'Target C$_{\mathrm{V}}$', size = 14)

plt.savefig('Cp_at_debT.pdf', bbox_inches = 'tight')


Cp_target = CV_at_debT['Cp Target']
Cp_pred = CV_at_debT['Cp Prediction']

Cp_target = Cp_target[~np.isnan(Cp_target)]
Cp_pred = Cp_pred[~np.isnan(Cp_pred)]

mae_Cp_dT = mean_absolute_error(Cp_target, Cp_pred)
r2_Cp_dT = r2_score(Cp_target, Cp_pred)
rmse_Cp_dT = sqrt(mean_squared_error(Cp_target, Cp_pred))
mad_Cp_dT = mean_absolute_deviation(Cp_target)


'''
Plot Debye Temperature Histogram
'''

plt.figure()
plt.hist(CV_at_debT['DebyeT'], 20)
plt.xlabel('Debye Temperature')
plt.savefig('debyeT_histo.pdf', bbox_inches = 'tight')