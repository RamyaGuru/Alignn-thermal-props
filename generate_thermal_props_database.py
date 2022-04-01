#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 14:35:29 2022

@author: rlg3

Run Cv and Svib calculation for the entire database

Save as id_prop.csv and generate the POSCAR files to perform the ALIGNN run

Save binned DOS converted to inverse THz as well

Notes:
    1. Don't need to worry about any DOS normalization.
    2. However, should do the histogram transformation of the DOS before
    calculating heat capacity
    3. Finally, still need to do the scale factor on this data
"""

from jarvis.db.figshare import data as jdata
from jarvis.core.atoms import Atoms
import pdos_integration as pint
import numpy as np
from jarvis.core.spectrum import Spectrum

import json


icm_to_thz = 2.99792458e-2
datafile = '../alignn_props_run11/'
f = open(datafile + 'ids_train_val_test.json',)

ids_split = json.load(f)

jvasp_list = sum([v for v in ids_split.values()], [])

dft_3d = jdata("edos_pdos")

T = 300

jid_list = []
Svib_list = []
Cv_list = []
dos_ithz = []
dos_orig = []
dos_norm_db_max = []

def max_intensity_db(db):
    max_each = []
    for d in db:
        if d['pdos_elast'] != 'na':
            max_each.append(max(d['pdos_elast']))
    max_total = max(max_each)
    return max_total

#max_total = max_intensity_db(dft_3d)

for jvasp in jvasp_list:
    jid = jvasp
    start = jid.find('JVASP')
    end = jid.find('.vasp')
    jid = jid[start:end]
    jid_list.append(jid)  
    match = next(i for i in dft_3d if i["jid"] == jid)
    atoms = Atoms.from_dict(match["atoms"])
    atoms.write_poscar(datafile + jvasp)
    scale = pint.get_natoms_from_db_entry(match)
    freq = np.arange(0, 1000, 5)
    s = Spectrum(x= freq, y=np.array(match['pdos_elast']) * scale)
    Svib_list.append(pint.vibrational_entropy(s.x, s.y, T))
    Cv_list.append(pint.heat_capacity(s.x, s.y, T))
    dos_orig_list = np.array(match['pdos_elast'])
    dos_orig_str = ",".join(map(str, dos_orig_list))
    dos_orig.append(dos_orig_str)
    dos_ithz_list = np.array(match['pdos_elast']) / icm_to_thz
    dos_ithz_str = ",".join(map(str, dos_ithz_list))
    dos_ithz.append(dos_ithz_str)
    
f1 = open(datafile + 'id_prop_Cv.csv', 'w+')
f2 = open(datafile + 'id_prop_Svib.csv', 'w+')
f3 = open(datafile + 'id_prop_dos_ithz.csv', 'w+')

for i in range(len(jvasp_list)):
    f1.write("%s,%s\n" % (jvasp_list[i], Cv_list[i]))
    f2.write("%s,%s\n" % (jvasp_list[i], Svib_list[i]))
    f3.write("%s,%s\n" % (jvasp_list[i], dos_ithz[i]))

f1.close()
f2.close()
f3.close()
    
    
    
    


