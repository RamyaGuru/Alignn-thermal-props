#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 00:10:32 2022

@author: rlg3
"""

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pdos_integration as pint
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from math import sqrt

mpl.rcdefaults()
mpl.rcParams['font.sans-serif'] = 'Apple Symbols'
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['axes.xmargin'] = 0.1
mpl.rcParams['axes.formatter.useoffset'] = False


def mean_absolute_deviation(data, axis=None):
    return np.mean(np.absolute(data - np.mean(data, axis)), axis)

'''
Load target and predicted dataframe. Plot against one another
'''
mpl.rcdefaults()
#mpl.rcParams['font.size'] = 16
mpl.rcParams.update({'font.monospace': 'Ubuntu Mono', 'font.sans-serif': 'Open Sans', 'font.family': ['PT Sans'], 'font.size': 16.0})

label = 'integrated_DOS'

pred_df = pd.read_csv('output_files/run11predictions_thermal_props_orig.csv')

target_df = pd.read_csv('output_files/run11target_thermal_props_orig.csv')

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

#plt.savefig('true_vs_pred_intDOS.pdf', bbox_inches='tight')



'''
Plots of heat capacity scaling
'''


from jarvis.db.figshare import data as jdata

dft_3d = jdata("edos_pdos")

jid = 'JVASP-32'

match = next(i for i in dft_3d if i["jid"] == jid)


target = np.array(match['pdos_elast'])
freq = np.linspace(0, 1000, len(target))

fig = plt.figure()
plt.plot(freq, target, label='DOS', color = 'xkcd:black')
plt.plot(freq, pint.heat_capacity_scaling(freq, 300) * 10 + 1e-3, color='xkcd:blood red', label=r'Prefactor')
plt.fill_between(freq, pint.heat_capacity_scaling_trig(freq, 300) * 10 + 1e-3, color = 'xkcd:blood red', alpha = 0.3)

ax = plt.gca()
ax.axes.get_yaxis().set_visible(False)
fig.text(0.80, 0.7, r'$C_{\mathrm{V}}$', va='center', fontsize = 28)


plt.xlim(0, 1000)
plt.ylim(0, 0.009)
plt.ylabel('DOS (a.u.)')
plt.xlabel(r'Frequency (cm$^{-1}$)')
plt.legend(frameon = False)

plt.savefig('cv_scaling_trig.pdf', bbox_inches='tight')


'''
Plot of vibrational entropy scaling
'''
fig = plt.figure()
plt.plot(freq, target, label='DOS', color = 'xkcd:black')
plt.plot(freq, pint.vibrational_entropy_scaling(freq, 300) * 10 + 1e-3, color='xkcd:blood red', label=r'S$_{\mathrm{vib}}$ scaling')
plt.fill_between(freq, pint.vibrational_entropy_scaling_trig(freq, 300) * 10 + 1e-3, color = 'xkcd:blood red', alpha = 0.3)

ax = plt.gca()
ax.axes.get_yaxis().set_visible(False)

fig.text(0.78, 0.7, r'$S_{\mathrm{vib}}$', va='center', fontsize = 28)

plt.xlim(0, 1000)
plt.ylim(0, 0.009)
plt.ylabel('DOS (a.u.)')
plt.xlabel(r'Frequency (cm$^{-1}$)')
#plt.legend(frameon = False)

plt.savefig('svib_scaling.pdf', bbox_inches='tight')


'''
Plot of isotope scattering rate
'''
fig = plt.figure()
plt.plot(freq, target, label='DOS', color = 'xkcd:black')
plt.plot(freq, pint.isotopic_tau_scaling(match, freq) * 1e29, color='xkcd:blood red', label=r'S$_{\mathrm{vib}}$ scaling')
plt.fill_between(freq, pint.isotopic_tau_scaling(match, freq) * 1e29, color = 'xkcd:blood red', alpha = 0.3)

ax = plt.gca()
ax.axes.get_yaxis().set_visible(False)

fig.text(0.78, 0.7, r'$\tau^{-1}_{\mathrm{i}}$', va='center', fontsize = 28)

plt.xlim(0, 1000)
plt.ylim(0, 0.009)
plt.xlabel(r'Frequency (cm$^{-1}$)')
plt.savefig('iso_tau_scaling.pdf', bbox_inches='tight')


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

CV_at_debT = pd.read_csv('output_files/run11Cp_at_half_debyeT.csv')

fig = plt.figure(figsize = (8,4), constrained_layout = True)
plt.tight_layout()

plt.subplot(1,2,1)
plt.scatter(CV_at_debT['Cp Target'],CV_at_debT['Cp Prediction'], s = 5, color= 'xkcd:medium blue', alpha=0.4)

x = np.linspace(0,250,10)
y = np.linspace(0,250,10)

plt.plot(x, y, linewidth=2, color='xkcd:black')

plt.plot(x, y -(np.quantile(CV_at_debT['Cp Target'], 0.5) - np.quantile(CV_at_debT['Cp Target'], 0.25)),\
         linestyle='--', linewidth=2, color='xkcd:black')

plt.plot(x, y + (np.quantile(CV_at_debT['Cp Target'], 0.5) - np.quantile(CV_at_debT['Cp Target'], 0.25)),\
         linestyle='--', linewidth=2, color='xkcd:black')

fig.text(0.12, 0.89, r'R$^2$ = 0.759', va='center', fontsize = 18)
fig.text(0.12, 0.81, r'MAE = 13', va='center', fontsize = 18)
fig.text(0.12, 0.73, r'MAD = 32', va='center', fontsize = 18)

fig.text(0.44, 0.22, '(a)', va='center', fontsize = 22)
fig.text(0.94, 0.22, '(b)', va='center', fontsize = 22)

plt.ylim([0, 250])
plt.xlim([0,250])

plt.xlabel(r'Predicted DOS C$_{\mathrm{V}}$ (J/mol/K)')
plt.ylabel(r'Target DOS C$_{\mathrm{V}}$ (J/mol/K)')

plt.subplot(1,2,2)

plt.hist(CV_at_debT['DebyeT'], 20, color = 'xkcd:bluey green')
plt.xlabel('Debye Temperature')
plt.ylabel('Counts')


#plt.savefig('Cp_at_debT.pdf', bbox_inches = 'tight')


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

# plt.figure()
# plt.hist(CV_at_debT['DebyeT'], 20, color = 'xkcd:bluey green')
# plt.xlabel('Debye Temperature')
# plt.ylabel('Counts')
# plt.savefig('debyeT_histo.pdf', bbox_inches = 'tight')


'''
Plot the Isotope Scattering Rate
'''
gamma_target = np.array(CV_at_debT['Isotope Gamma Target']) / 1e9
gamma_pred = np.array(CV_at_debT['Isotope Gamma Prediction']) / 1e9

gamma_target = gamma_target[np.nonzero(gamma_target)]
gamma_pred = gamma_pred[np.nonzero(gamma_pred)]

debyeT_gamma = np.array(CV_at_debT['DebyeT'])[np.nonzero(gamma_target)]

debyeT_gen = np.array(CV_at_debT['DebyeT'])

mdiff = np.array(CV_at_debT['Mass Difference'])

mdiff_gamma = mdiff[np.nonzero(gamma_target)]

mae_gamma_dT = mean_absolute_error(gamma_target, gamma_pred)
r2_gamma_dT = r2_score(gamma_target, gamma_pred)
rmse_gamma_dT = sqrt(mean_squared_error(gamma_target, gamma_pred))
mad_gamma_dT = mean_absolute_deviation(gamma_target)

plt.figure()
plt.scatter(gamma_target, gamma_pred, s = 5, color= 'xkcd:medium blue', alpha=0.4)

x = np.linspace(0,100,10)
y = np.linspace(0,100,10)

plt.plot(x, y, linestyle='--', linewidth=2, color='xkcd:black')

plt.plot(x, y - np.std(gamma_target), linestyle=':', linewidth=2, color='xkcd:black')

plt.plot(x, y + np.std(gamma_target), linestyle=':', linewidth=2, color='xkcd:black')

plt.ylim([0, 100])
plt.xlim([0,100])

#plt.savefig('iso_tau_comparison.pdf', bbox_inches = 'tight')


'''
Heat Capacity versus Temperature
'''

jid = 'JVASP-317'

match2 = next(i for i in dft_3d if i["jid"] == jid)

target = np.array(match2['pdos_elast'])
freq = np.linspace(0, 1000, len(target))
scale = pint.get_natoms_from_db_entry(match2)
target = target * scale

debyeT = 372.2377

DP_limit = 3 * 8.3145 * 3

Cv_T = []

for temp in np.linspace(0, 700, 500):
    Cv_T.append(pint.heat_capacity(freq, target, T = temp))

fig = plt.figure()
plt.plot(np.linspace(0, 700, 500), Cv_T, color = 'xkcd:black', label = r'TiS$_2$')
plt.plot(np.ones(10) * debyeT, np.linspace(0, 85, 10), linestyle = ':', color = 'xkcd:black', linewidth = 3, alpha = 0.5)
fig.text(0.3, 0.4, r'$T_{\mathrm{D}}$ = 372 K', va='center', fontsize = 18, color = 'xkcd:black')
plt.plot(np.linspace(0, 700, 100), np.ones(100) * DP_limit, linestyle = ':', color = 'xkcd:blood red', linewidth = 3)
fig.text(0.13, 0.825, 'D-P Limit: 74.8 J/mol/K', va='center', fontsize = 16, color = 'xkcd:blood red')

plt.ylim([0, 85])
plt.xlim([0, 700])


jid = 'JVASP-1076'

match3 = next(i for i in dft_3d if i["jid"] == jid)

target = np.array(match3['pdos_elast'])
freq = np.linspace(0, 1000, len(target))
target = target * scale

debyeT = 457.8118

Cv_T = []
for temp in np.linspace(0, 700, 500):
    Cv_T.append(pint.heat_capacity(freq, target, T = temp))


plt.plot(np.linspace(0, 700, 500), Cv_T, color = 'xkcd:medium blue', label = r'MgSi$_2$')
fig.text(0.65, 0.6, r'$T_{\mathrm{D}}$ = 457 K', va='center', fontsize = 18, color = 'xkcd:medium blue')

plt.plot(np.ones(10) * debyeT, np.linspace(0, 85, 10), linestyle = ':', color = 'xkcd:medium blue', linewidth = 3, alpha = 0.5)

plt.ylabel(r'$C_{\mathrm{V}}$ (J/mol/K)')
plt.xlabel('Temperature (K)')
plt.legend(frameon = False)
#plt.plot(np.linspace(0, 600, 100), np.ones(100) * DP_limit, linestyle = ':', color = 'xkcd:blood red', linewidth = 3)

plt.savefig('Cv_vs_T.pdf', bbox_inches = 'tight')


'''
2x2 subplots with trends: Color-coded by idfferent attributes
'''
mpl.rcdefaults()

fig = plt.figure(constrained_layout=True)
plt.tight_layout()
#fig.subplots_adjust(hspace=2)
#Integrated DOS

label = 'integrated_DOS'

plt.subplot(2,2,1)

plt.scatter(target_df[label], pred_df[label], s = 3, c= mdiff, alpha=0.5, cmap = 'viridis', vmin=0, vmax=200)

x = np.linspace(0,50,10)
y = np.linspace(0,50,10)

plt.plot(x, y, linestyle='--', linewidth=2, color='xkcd:black')

plt.plot(x, y - (np.quantile(target_df[label], 0.5) - np.quantile(target_df[label], 0.25)), linestyle=':', linewidth=2, color='xkcd:black')

plt.plot(x, y + (np.quantile(target_df[label], 0.75) - np.quantile(target_df[label], 0.5)), linestyle=':', linewidth=2, color='xkcd:black')

plt.ylim([0, 50])
plt.xlim([0,50])

fig.text(0.125, 0.95, r'R$^2$ = 0.767', va='center', fontsize = 12)
plt.ylabel(r'Predicted Integrated DOS')
plt.xlabel(r'Target Integrated DOS')


#Heat Caoacity

label = 'Cp (J/mol/K)'

plt.subplot(2,2,2)

plt.scatter(target_df[label], pred_df[label], s = 3, c= mdiff, alpha=0.5, cmap = 'viridis', vmin=0, vmax=200)

x = np.linspace(0,300,10)
y = np.linspace(0,300,10)

plt.plot(x, y, linestyle='--', linewidth=2, color='xkcd:black')

plt.plot(x, y - (np.quantile(target_df[label], 0.5) - np.quantile(target_df[label], 0.25)), linestyle=':', linewidth=2, color='xkcd:black')

plt.plot(x, y + (np.quantile(target_df[label], 0.75) - np.quantile(target_df[label], 0.5)), linestyle=':', linewidth=2, color='xkcd:black')

plt.ylim([0, 300])
plt.xlim([0,300])

fig.text(0.565, 0.95, r'R$^2$ = 0.749', va='center', fontsize = 12)
plt.ylabel(r'Predicted DOS C$_{\mathrm{V}}$ (J/mol/K)')
plt.xlabel(r'Target  DOS C$_{\mathrm{V}}$ (J/mol/K)')


# Vibrational Entropy

label = 'S_vib (J/mol/K)'

plt.subplot(2,2,3)

plt.scatter(target_df[label], pred_df[label], s = 3, c= mdiff, alpha=0.5, cmap = 'viridis', vmin=0, vmax=200)

x = np.linspace(0,2000,10)
y = np.linspace(0,2000,10)

plt.plot(x, y, linestyle='--', linewidth=2, color='xkcd:black')

plt.plot(x, y - (np.quantile(target_df[label], 0.5) - np.quantile(target_df[label], 0.25)), linestyle=':', linewidth=2, color='xkcd:black')

plt.plot(x, y + (np.quantile(target_df[label], 0.75) - np.quantile(target_df[label], 0.5)), linestyle=':', linewidth=2, color='xkcd:black')

plt.ylim([0, 2000])
plt.xlim([0,2000])

fig.text(0.125, 0.45, r'R$^2$ = 0.716', va='center', fontsize = 12)
plt.ylabel(r'Predicted DOS S$_{\mathrm{vib}}$ (J/mol/K)')
plt.xlabel(r'Target  DOS S$_{\mathrm{vib}}$ (J/mol/K)')


# Isotopic Scattering Rate

plt.subplot(2,2,4)

plt.scatter(gamma_target, gamma_pred, s = 3, c= mdiff_gamma, alpha=0.5, cmap = 'viridis', vmin=0, vmax=200)

x = np.linspace(0,100,10)
y = np.linspace(0,100,10)

plt.plot(x, y, linestyle='--', linewidth=2, color='xkcd:black')

plt.plot(x, y - (np.quantile(gamma_target, 0.5) - np.quantile(gamma_target, 0.25)), linestyle=':', linewidth=2, color='xkcd:black')

plt.plot(x, y + (np.quantile(gamma_target, 0.75) - np.quantile(gamma_target, 0.5)), linestyle=':', linewidth=2, color='xkcd:black')

plt.ylim([0, 100])
plt.xlim([0,100])

fig.text(0.565, 0.45, r'R$^2$ = 0.860', va='center', fontsize = 12)
plt.ylabel(r'Predicted DOS $\tau^{-1}_{\mathrm{i}}$ (GHz)')
plt.xlabel(r'Target  DOS $\tau^{-1}_{\mathrm{i}}$ (GHz)')

plt.colorbar(label = 'Max Mass Difference')
#plt.savefig('target_vs_pred_plots_mdiff.pdf', bbox_inches = 'tight')



'''
Color by abundance of points in region
'''
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

#Integrated DOS

mpl.rcdefaults()


#fig = plt.figure(figsize = (8,8), constrained_layout=True)
fig, ax = plt.subplots(3,4, gridspec_kw={'height_ratios': [2, 2, 1.2]}, constrained_layout = True, figsize = (8,9))
plt.tight_layout()

label = 'integrated_DOS'

histo_idos = np.histogram2d(target_df[label], pred_df[label], bins = 100)

cx = np.digitize(target_df[label], histo_idos[1])
cy = np.digitize(pred_df[label], histo_idos[2])

cx = np.where(cx==101, 100, cx) - 1
cy = np.where(cy==101, 100, cy) - 1
pairs = [(x,y) for x,y in zip(cx, cy)]

c_idos = []

for p in pairs:
    c_idos.append(histo_idos[0][p])

plt.subplot(3,2,1)
x = np.linspace(0,50,10)
y = np.linspace(0,50,10)

plt.plot(x, y, linewidth=2, color='xkcd:black')

plt.plot(x, y - (np.quantile(target_df[label], 0.5) - np.quantile(target_df[label], 0.25)), linestyle='--', linewidth=2, color='xkcd:black')

plt.plot(x, y + (np.quantile(target_df[label], 0.75) - np.quantile(target_df[label], 0.5)), linestyle='--', linewidth=2, color='xkcd:black')

plt.scatter(target_df[label], pred_df[label], c = c_idos, s = 3, alpha = 0.5, cmap = 'Spectral_r', vmax = 100)

#fig.text(0.375, 0.725, '(a)', va='center', fontsize = 18)
fig.text(0.1, 0.95, r'R$^2$ = 0.525', va='center', fontsize = 16)
#plt.savefig('target_vs_pred_plots_conc.pdf', bbox_inches = 'tight')

plt.ylim([0, 50])
plt.xlim([0, 50])

plt.ylabel(r'Predicted Integrated DOS', fontsize = 12)
plt.xlabel(r'Target Integrated DOS', fontsize = 12)

#Heat Capacity

label = 'Cp (J/mol/K)'

histo_cv = np.histogram2d(target_df[label], pred_df[label], bins = 100)

cx = np.digitize(target_df[label], histo_cv[1])
cy = np.digitize(pred_df[label], histo_cv[2])

cx = np.where(cx==101, 100, cx) - 1
cy = np.where(cy==101, 100, cy) - 1
pairs = [(x,y) for x,y in zip(cx, cy)]

c_cv = []

for p in pairs:
    c_cv.append(histo_cv[0][p])


plt.subplot(3,2,2)


x = np.linspace(0,300,10)
y = np.linspace(0,300,10)

plt.plot(x, y, linewidth=2, color='xkcd:black')

plt.plot(x, y - (np.quantile(target_df[label], 0.5) - np.quantile(target_df[label], 0.25)), linestyle='--', linewidth=2, color='xkcd:black')

plt.plot(x, y + (np.quantile(target_df[label], 0.75) - np.quantile(target_df[label], 0.5)), linestyle='--', linewidth=2, color='xkcd:black')

#fig.text(0.815, 0.725, '(b)', va='center', fontsize = 18)
fig.text(0.6, 0.95, r'R$^2$ = 0.502', va='center', fontsize = 16)
plt.ylim([0, 300])
plt.xlim([0,300])

plt.scatter(target_df[label], pred_df[label], c = c_cv, s = 3, alpha = 0.5, cmap = 'Spectral_r', vmax = 100)

plt.ylabel(r'Predicted C$_{\mathrm{V}}$ (J/mol/K)', fontsize = 12)
plt.xlabel(r'Target C$_{\mathrm{V}}$ (J/mol/K)', fontsize = 12)



#Vibrational Entropy

label = 'S_vib (J/mol/K)'

histo_svib = np.histogram2d(target_df[label], pred_df[label], bins = 100)

cx = np.digitize(target_df[label], histo_svib[1])
cy = np.digitize(pred_df[label], histo_svib[2])

cx = np.where(cx==101, 100, cx) - 1
cy = np.where(cy==101, 100, cy) - 1
pairs = [(x,y) for x,y in zip(cx, cy)]

c_svib = []

for p in pairs:
    c_svib.append(histo_svib[0][p])
    
plt.subplot(3,2,3)


plt.scatter(target_df[label], pred_df[label], s = 3, c= c_svib, alpha=0.5, cmap = 'Spectral_r', vmax = 100)

x = np.linspace(0,2000,10)
y = np.linspace(0,2000,10)


plt.plot(x, y, linewidth=2, color='xkcd:black')

plt.plot(x, y - (np.quantile(target_df[label], 0.5) - np.quantile(target_df[label], 0.25)), linestyle='--', linewidth=2, color='xkcd:black')

plt.plot(x, y + (np.quantile(target_df[label], 0.75) - np.quantile(target_df[label], 0.5)), linestyle='--', linewidth=2, color='xkcd:black')

#fig.text(0.375, 0.37, '(c)', va='center', fontsize = 18)
fig.text(0.1, 0.61, r'R$^2$ = 0.560', va='center', fontsize = 16)
plt.xlim([0,2000])
plt.ylim([0,2000])

plt.ylabel(r'Predicted S$_{\mathrm{vib}}$ (J/mol/K)', fontsize = 12)
plt.xlabel(r'Target S$_{\mathrm{vib}}$ (J/mol/K)', fontsize = 12)

# Isotope Scattering

histo_scatt = np.histogram2d(gamma_target, gamma_pred, bins = 100)

cx = np.digitize(gamma_target, histo_scatt[1])
cy = np.digitize(gamma_pred, histo_scatt[2])

cx = np.where(cx==101, 100, cx) - 1
cy = np.where(cy==101, 100, cy) - 1
pairs = [(x,y) for x,y in zip(cx, cy)]

c_scatt = []

for p in pairs:
    c_scatt.append(histo_scatt[0][p])

plt.subplot(3,2,4)
print(max(c_scatt))

# plt.scatter(gamma_target, gamma_pred, s = 3, c= c_scatt, alpha=0.5, cmap = 'Spectral_r', vmax = 100)

# x = np.linspace(0,100,10)
# y = np.linspace(0,100,10)

# plt.plot(x, y, linewidth=2, color='xkcd:black')

# plt.plot(x, y - (np.quantile(gamma_target, 0.5) - np.quantile(gamma_target, 0.25)), linestyle='--', linewidth=2, color='xkcd:black')

# plt.plot(x, y + (np.quantile(gamma_target, 0.75) - np.quantile(gamma_target, 0.5)), linestyle='--', linewidth=2, color='xkcd:black')

# #fig.text(0.815, 0.37, '(d)', va='center', fontsize = 18)
# fig.text(0.555, 0.61, r'R$^2$ = 0.860', va='center', fontsize = 16)
# plt.ylim([0, 100])
# plt.xlim([0,100])


# plt.ylabel(r'Predicted DOS $\tau^{-1}_{\mathrm{i}}$ (GHz)', fontsize = 12)
# plt.xlabel(r'Target  DOS $\tau^{-1}_{\mathrm{i}}$ (GHz)', fontsize = 12)

plt.colorbar(label = 'Sample Count')

'''
Property histograms
'''
plt.subplot(3,4,9)
#fig.text(0.21, 0.25, '(e)', va='center', fontsize = 18)

plt.hist(target_df['integrated_DOS'], 20, color = 'xkcd:bluey green')
plt.ylabel('Counts', fontsize = 12)
plt.xlabel('Integrated DOS', fontsize = 12)

plt.subplot(3,4,10)
plt.hist(target_df['Cp (J/mol/K)'], 20, color = 'xkcd:bluey green')
plt.xlabel(r'C$_{\mathrm{V}}$ (J/mol/K)', fontsize = 12)
#fig.text(0.45, 0.25, '(f)', va='center', fontsize = 18)

plt.subplot(3,4,11)
plt.hist(target_df['S_vib (J/mol/K)'], 20, color = 'xkcd:bluey green')
plt.xlabel('S$_{\mathrm{vib}}$ (J/mol/K)', fontsize = 12)
#fig.text(0.675, 0.25, '(g)', va='center', fontsize = 18)

# plt.subplot(3,4,12)
# plt.hist(gamma_target, 20, color = 'xkcd:bluey green')
# plt.xlabel(r'$\tau^{-1}_{\mathrm{i}}$ (GHz)', fontsize = 12)
# #fig.text(0.932, 0.25, '(h)', va='center', fontsize = 18)
# plt.xlim([0, 100])

#plt.savefig('target_vs_pred_plots_conc.pdf', bbox_inches = 'tight')




'''
NEW FIGURE. Just showing vibraitonal entropy, heat capacity, and isotope scattering
'''

