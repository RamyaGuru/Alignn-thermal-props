#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 14:35:29 2022

@author: rlg3

Run Cv and Svib calculation for the entire database

Save as id_prop.csv and generate the POSCAR files to perform the ALIGNN run

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

f = open('../../run11/temp/ids_train_val_test.json',)

ids_split = json.load(f)

jvasp_list = sum([v for v in ids_split.values()], [])

dft_3d = jdata("edos_pdos")

T = 300

jid_list = []
Svib_list = []
Cv_list = []
for jid in jvasp_list:
    start = jid.find('JVASP')
    end = jid.find('.vasp')
    jid = jid[start:end]
    jid_list.append(jid)  
    match = next(i for i in dft_3d if i["jid"] == jid)
    scale = pint.get_natoms_from_db_entry(match)
    freq = np.arange(0, 1000, 5)
    s = Spectrum(x= freq, y=np.array(match['pdos_elast']) * scale)
    Svib_list.append(pint.vibrational_entropy(s.x, s.y, T))
    Cv_list.append(pint.heat_capacity(s.x, s.y, T))

f1 = open('id_prop_Cv.csv')
f2 = open('id_prop_Svib.csv')

for i in range(len(jvasp_list)):
    f1.write("%s,%s\n" % (jvasp_list[i], Cv_list[i]))
    f2.write("%s,%s\n" % (jvasp_list[i], Svib_list[i]))

f1.close()
f2.close()
    
    
    
    
    

