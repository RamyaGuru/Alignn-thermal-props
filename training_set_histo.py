#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 15:25:45 2022

@author: rlg3

Training dataset hitogram
"""

import json
from jarvis.core.atoms import Atoms
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

run = "run21"

train_json='../../{}/PhononData-Freq-300_1000_20.json'.format(run)

with open(train_json) as f:
    dos_dict = json.load(f)
    
freq = np.linspace(-300, 1000, len(dos_dict[0]['pdos_elast']))
        
def label_dynamic_stability(freq, dos_dict):
    stable_list = []
    unstable_list = []
    zero_indx = np.where(freq > 0)[0][0] - 1
    neg_freq = np.linspace(-300, 0, zero_indx)
    for d in dos_dict:
        intdos_neg_target = np.trapz(neg_freq, d['pdos_elast'][:zero_indx])
        intdos_target = np.trapz(freq, d['pdos_elast'])
        if intdos_neg_target / intdos_target > 0.1:
            unstable_list.append(d)
        else:
            stable_list.append(d)
    return stable_list, unstable_list  
        

stable_list, unstable_list = label_dynamic_stability(freq, dos_dict)

#Histograms of average volume per atom
#List of number of unique elements in the compound
atmV_list = []
nary_list = []
xtal_list = {'cubic' : 0,#4290,
              'hexagonal' : 0,#9941,
              'trigonal' : 0,#8685,
              'tetragonal' : 0,#7216,
              'orthorhombic' : 0,#3510,
              'monoclinic' : 0,#1954,
              'triclinic' : 0}#525}
for d in stable_list:
    atoms = Atoms.from_dict(d['atoms'])
    atmV = (atoms.volume / atoms.num_atoms)
    atmV_list.append(atmV)
    nary_list.append(len(set(atoms.elements)))
    spg = atoms.spacegroup()
    n = int(spg[spg.find('(')+1:spg.find(')')])
    if 0 < n < 3:
        xtal_list['triclinic'] += 1
    elif n < 16:
        xtal_list['monoclinic'] += 1
    elif n < 75:
        xtal_list['orthorhombic'] += 1
    elif n < 143:
        xtal_list['tetragonal'] += 1
    elif n < 168:
        xtal_list['trigonal'] += 1
    elif n < 195:
        xtal_list['hexagonal'] += 1
    else:
        xtal_list['cubic'] += 1


with open('output_files/xtal_list.json', 'w') as out_f:
    json.dump(xtal_list, out_f)
    
mpl.rcdefaults()
mpl.rcParams['font.size'] = 16   
fig, ax = plt.subplots(3, 1, figsize = (5, 8))
fig.tight_layout(h_pad = 2)

plt.subplot(3, 1, 2)
plt.hist(atmV_list, bins = 20, color='xkcd:parchment', edgecolor='black', linewidth=1.2)
plt.ylabel('Sample Count', fontsize = 20)
plt.xlabel(r'Average Volume per Atom ($\AA^3$)')
plt.xlim([0, 80])
#plt.savefig('train_volume_histo.pdf', bbox_inches = 'tight')

nary_counts = {}
for i in np.array(range(max(nary_list))) + 1:
    nary_counts[i] = nary_list.count(i)
    
plt.subplot(3, 1, 1)
plt.bar(list(nary_counts.keys()), list(nary_counts.values()), color='xkcd:parchment', edgecolor='black', linewidth=1.2)
#plt.ylabel('Sample Count')
plt.xlabel('Number of Unique Elements')
plt.xticks([1,2,3,4,5,6])
#plt.savefig('train_nelem_histo.pdf', bbox_inches = 'tight')

plt.subplot(3, 1, 3)

plt.bar(list(xtal_list.keys()), list(xtal_list.values()), color='xkcd:parchment', edgecolor='black', linewidth=1.2)
plt.xticks(rotation=30, ha='right')
#plt.ylabel('Sample Count')

plt.savefig('training_set_histo_horiz_stable.pdf', bbox_inches = 'tight')


#Plot the relative abundance of the top 10 elements for each set

stable_elem = np.load("stable_elem_freq.npz")


unstable_elem = np.load("unstable_elem_freq.npz")



