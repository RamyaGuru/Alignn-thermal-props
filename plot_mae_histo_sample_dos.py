#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 22:24:28 2022

@author: rlg3

MAE histograms
"""

from sklearn.metrics import mean_absolute_error
import json
import matplotlib.pyplot as plt
import matplotlib as mpl
from jarvis.db.figshare import data as jdata
from jarvis.core.atoms import Atoms
import numpy as np
from pyvalem.formula import Formula

dft_3d = jdata("dft_3d")

mpl.rcdefaults()
#Calculate MAE for each spectrum pair

run = 'run21'
#folder = '{}/run_20_output'.format(run)

f_in ='../../{}/multi_out_predictions.json'.format(run)


MAE_list = [] # Stores the MAE for each sample DOS
atmV_list = []
jid_list = []

with open(f_in) as jfile:
    dos_dict = json.load(jfile)
    
    
# for d in dos_dict:
#     mae = mean_absolute_error(d['target'], d['predictions'])
#     d['MAE'] = mae
#     MAE_list.append(mae)
#     jid = d["id"]
#     start = jid.find('JVASP')
#     end = jid.find('.vasp')
#     jid = jid[start:end]
#     jid_list.append(jid)
#     match = next(d for d in dft_3d if d["jid"] == jid)
#     atoms = Atoms.from_dict(match['atoms'])
#     atmV = (atoms.volume / atoms.num_atoms)
#     atmV_list.append(atmV)
#     d['atmV'] = atmV
#     d['composition'] = atoms.composition.reduced_formula
#     d['nelem'] = len(set(atoms.elements))
#     spg = atoms.spacegroup()
#     n = int(spg[spg.find('(')+1:spg.find(')')])
#     if 0 < n < 3:
#         xtal = 'triclinic'
#     elif n < 16:
#         xtal =  'monoclinic'
#     elif n < 75:
#         xtal = 'orthorhombic'
#     elif n < 143:
#         xtal = 'tetragonal'
#     elif n < 168:
#         xtal = 'trigonal'
#     elif n < 195:
#         xtal = 'hexagonal' 
#     else:
#         xtal = 'cubic'
#     d['xtal_system'] =  xtal
    
# f_out = '../../{}/predictions_augmented.json'.format(run)  
 
# with open(f_out, 'w') as j_out:
#     json.dump(dos_dict, j_out)


'''
Load the dictionary from the json file
'''

in_file = '../../{}/predictions_augmented.json'.format(run)

with open(in_file) as json_file:
    dos_dict = json.load(json_file)
    
MAE_list = []
atmV_list = []
comp_list = []

nelem_mae = {1:[], 2:[], 3:[], 4:[], 5:[], 6:[]}


xtal_mae = {'cubic' : [],
             'hexagonal' : [],
             'trigonal' : [],
             'tetragonal' : [],
             'orthorhombic' : [],
             'monoclinic' : [],
             'triclinic' : []}

for d in dos_dict:
    comp_list.append(d['composition'])
    MAE_list.append(d['MAE'])
    atmV_list.append(d['atmV'])
    nelem_mae[d['nelem']].append(d['MAE'])
    xtal_mae[d['xtal_system']].append(d['MAE'])

#Plot histogram versus MAE
norm = mpl.colors.Normalize(vmin=0, vmax=0.25)
scale_map = plt.cm.ScalarMappable(norm, 'plasma')
scale_map.set_clim(0, 0.25)
cm = scale_map.get_cmap()

plt.figure(figsize = (7,3))
n, bins, patches = plt.hist(MAE_list, bins = 15, edgecolor = 'black')

bin_centers = 0.5 * (bins[:-1] + bins[1:])

# scale values to interval [0,1]
col = bin_centers - min(bin_centers)
col /= max(col)

for c, p in zip(col, patches):
    plt.setp(p, 'facecolor', cm(c))
    
    
plt.xlabel('Mean Absolute Error', fontsize = 14)
plt.ylabel('Number of Samples', fontsize = 14)
plt.savefig('figures/{}_mae_histogram.pdf'.format(run), bbox_inches = 'tight')
plt.show()


'''
MAE Trends
'''

histo_mae = np.histogram2d(atmV_list, MAE_list, bins = 100)

cx = np.digitize(atmV_list, histo_mae[1])
cy = np.digitize(MAE_list, histo_mae[2])

cx = np.where(cx==101, 100, cx) - 1
cy = np.where(cy==101, 100, cy) - 1

pairs = [(x,y) for x,y in zip(cx, cy)]

c_mae = []

for p in pairs:
    c_mae.append(histo_mae[0][p])

fig, ax = plt.subplots(3, 1, gridspec_kw={'height_ratios': [1.5, 1, 1]}, figsize = (5, 9))
fig.tight_layout(h_pad = 2)

plt.subplot(3, 1, 1)
plt.scatter(atmV_list, MAE_list, s=4, c=c_mae, cmap = 'viridis_r')
plt.ylabel('Mean Absolute Error', fontsize = 12)
plt.xlabel('Average Atomic Volume ($\AA ^3$)', fontsize = 12)
plt.xlim([5, 50])
plt.ylim([0, 0.25])

plt.colorbar().set_label(label = 'Sample Count', fontsize = 12)
#plt.savefig('mae_vs_atmV.pdf', bbox_inches = 'tight')


'''
Get Average MAE per number of elements
'''
#mpl.rcParams['font.size'] = 12

plt.subplot(3, 1, 2)
avg_mae = [sum(v) / len(v) for v in nelem_mae.values()]
count_list = [len(v) for v in nelem_mae.values()]
plt.bar(list(nelem_mae.keys()), avg_mae, color='xkcd:parchment', edgecolor='black', linewidth=1.2)
plt.xticks(fontsize = 14)
plt.xlabel('Number of Unique Elements', fontsize = 12)
ax2 = ax[1].twinx() 
plt.bar(list(nelem_mae.keys()), count_list, color = 'xkcd:bluish', edgecolor='black', linewidth=1.2, width = 0.4)
ax[1].set_ylabel('Mean Absolute Error', fontsize = 12)
ax2.set_ylabel('Counts in Test Set', fontsize = 12, color = 'xkcd:bluish')
ax2.tick_params(axis='y', colors='xkcd:bluish')


'''
Get Average MAE for each crystal system
'''
plt.subplot(3, 1, 3)
avg_mae = [sum(v) / len(v) for v in xtal_mae.values() if len(v) != 0]
count_list = [len(v) for v in xtal_mae.values() if len(v) != 0]
plt.bar(list(xtal_mae.keys()), avg_mae, color='xkcd:parchment', edgecolor='black', linewidth=1.2)
plt.xticks(rotation=30, ha='right', fontsize = 12)
#plt.xlabel('Crystal System', fontsize = 12)
ax3 = ax[2].twinx() 
plt.bar(list(xtal_mae.keys()), count_list, color = 'xkcd:bluish', edgecolor='black', linewidth=1.2, width = 0.4)
ax[2].set_ylabel('Mean Absolute Error', fontsize = 12)
ax3.set_ylabel('Counts in Test Set', fontsize = 12, color = 'xkcd:bluish')
ax3.tick_params(axis='y', colors='xkcd:bluish')

plt.savefig('figures/{}_mae_dependency.pdf'.format(run), bbox_inches = 'tight')

'''
Get list of compounds in each MAE bin.
Choose composition from the first 8 bins and plot example spectra
'''

MAE_bin = np.digitize(MAE_list, bins)

comp_bins = []
for bindx in range(1, 17):
    indx = np.where(MAE_bin == bindx)
    comp_bins.append(np.array(comp_list)[indx[0].astype(int)])
    # print('Compositions in MAE bin %d', bindx)
    # print(comp_bin)

'''
Generate grid of example DOS in first 8 bins
'''

fig, ax = plt.subplots(2, 4, figsize = (10, 5))
fig.tight_layout(h_pad = 2)
freq = np.linspace(-300, 1000, len(dos_dict[0]['target']))
for b in range(1,9):
    indx = np.where(MAE_bin == b)[0][10]
    dict_item = dos_dict[indx]
    plt.subplot(2, 4, b)
    plt.plot(freq, dict_item['target'],color = 'xkcd:black', linewidth = 2)
    plt.plot(freq, dict_item['predictions'], color = cm(col[b-1]), alpha = 0.8, linewidth = 2)
    f = Formula(dict_item['composition'])
    title_str = '$' + f.latex + '$' + '; MAE:' + str(round(dict_item['MAE'], 2))
    plt.title(title_str)

fig.text(-0.02, 0.5, 'Scaled DOS (a.u.)', va='center', rotation='vertical', fontsize = 14)
fig.text(0.5, -0.02, r'Frequency (cm$^{-1}$)', ha='center', fontsize = 14)

plt.savefig('figures/{}_sample_DOS_grid.pdf'.format(run), bbox_inches = 'tight')
# fig, ax = plt.subplots(2, 4, figsize = (10,5))
# fig.tight_layout()

# plt.subplot(2, 4, 1)
