#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 22:45:01 2022

@author: rlg3

Potential new main text figure to show the large sclae prediction of
new compound properties
"""

import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
from pyvalem.formula import Formula
import json


'''
Get the histograms for each predicted property
'''

pd.reset_option('display.float_format')#, '{:.2E}'.format)

run = 'run21'

json_file = 'output_files/{}_thermal_props_dft3d_mod.json'.format(run)

with open(json_file) as infile:
    dft3d_dict = json.load(infile)
    
dft3d_df = pd.DataFrame.from_dict(dft3d_dict)

'''
Histograms
'''

fig, ax = plt.subplots(1, 3, figsize = (9, 2))
fig.tight_layout(w_pad = 2)

plt.subplot(1, 3, 1)

plt.hist(dft3d_df['Cv_norm'], 50, color = 'xkcd:bluey green')
plt.xlabel(r'C$_{\mathrm{V}}$ (Jkg$^{-1}$K$^{-1}$)', fontsize = 12)
plt.ylabel('Counts')
plt.title('Heat Capacity', fontsize = 14)


plt.subplot(1, 3, 2)
plt.hist(dft3d_df['Svib_norm'], 50, color = 'xkcd:bluey green')
plt.xlabel(r'S$_{\mathrm{vib}}$ (Jkg$^{-1}$K$^{-1}$)' , fontsize = 12)
plt.title('Vibrational Entropy', fontsize = 14)

plt.subplot(1, 3, 3)
plt.hist(np.array(dft3d_df['isotope_scatt']) / 1e9, 250, color = 'xkcd:bluey green')
plt.xlabel(r'$\tau^{-1}$ (GHz)', fontsize = 12)
plt.xlim(0, 100)
plt.title('Isotopic Scattering Rate', fontsize = 14)


Svib_sorted = dft3d_df.sort_values(by = ['Svib_norm'])
Cv_sorted = dft3d_df.sort_values(by = ['Cv_norm'])
iso_scatt = dft3d_df.sort_values(by = ['isotope_scatt'])

iso_scatt = iso_scatt[iso_scatt['isotope_scatt'] > 0]

plt.savefig('pred_property_hist.pdf', bbox_inches = 'tight')


'''
Example DOS for high and low values
'''
freq = np.linspace(0, 1000, 51)

fig, ax = plt.subplots(2, 3, figsize = (9, 4))
fig.tight_layout(w_pad = 2, h_pad = 5)

plt.subplot(2, 3, 1)
plt.plot(freq, np.array(Cv_sorted['scaled_dos'][0:1])[0],\
         color = 'xkcd:black')
cf = Formula(np.array(Cv_sorted['form_unit'][0:1])[0])
jid = np.array(Cv_sorted['jid'][0:1])[0]
title_str = '$' + cf.latex + '$' + '; JID:' + jid
plt.title(title_str)

plt.subplot(2, 3, 4)

plt.plot(freq, np.array(Cv_sorted['scaled_dos'][-2:-1])[0],\
         color = 'xkcd:black')
cf = Formula(np.array(Cv_sorted['form_unit'][-2:-1])[0])
jid = np.array(Cv_sorted['jid'][-2:-1])[0]
title_str = '$' + cf.latex + '$' + '; JID:' + jid
plt.title(title_str)

plt.subplot(2, 3, 2)

plt.plot(freq, np.array(Svib_sorted['scaled_dos'][0:1])[0],\
         color = 'xkcd:black')
cf = Formula(np.array(Svib_sorted['form_unit'][0:1])[0])
jid = np.array(Svib_sorted['jid'][0:1])[0]
title_str = '$' + cf.latex + '$' + '; JID:' + jid
plt.title(title_str)

plt.subplot(2, 3, 5)

plt.plot(freq, np.array(Svib_sorted['scaled_dos'][-2:-1])[0],\
         color = 'xkcd:black')
cf = Formula(np.array(Svib_sorted['form_unit'][-2:-1])[0])
jid = np.array(Svib_sorted['jid'][-2:-1])[0]
title_str = '$' + cf.latex + '$' + '; JID:' + jid
plt.title(title_str)

plt.subplot(2, 3, 3)

plt.plot(freq, np.array(iso_scatt['scaled_dos'][0:1])[0],\
         color = 'xkcd:black')
cf = Formula(np.array(iso_scatt['form_unit'][0:1])[0])
jid = np.array(iso_scatt['jid'][0:1])[0]
title_str = '$' + cf.latex + '$' + '; JID:' + jid
plt.title(title_str)

plt.subplot(2, 3, 6)

plt.plot(freq, np.array(iso_scatt['scaled_dos'][-2:-1])[0],\
         color = 'xkcd:black')
cf = Formula(np.array(iso_scatt['form_unit'][-2:-1])[0])
jid = np.array(iso_scatt['jid'][-2:-1])[0]
title_str = '$' + cf.latex + '$' + '; JID:' + jid
plt.title(title_str)


fig.text(0.5, 1.05, 'DOS Corresponding to Lowest Property Value', ha='center', fontsize = 14)
fig.text(0.5, 0.5, 'DOS Corresponding to Highest Property Value', ha='center', fontsize = 14)
fig.text(-0.03, 0.5, 'Scaled DOS (a.u.)', va='center', rotation='vertical', fontsize = 12)
fig.text(0.5, -0.02, r'Frequency (cm$^{-1}$)', ha='center', fontsize = 12)

plt.savefig('pred_property_high_low.pdf', bbox_inches = 'tight')