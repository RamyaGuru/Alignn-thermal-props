#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 16:44:32 2022

@author: rlg3

assemble the id_prop.csv from the JSON file
"""

import json
from jarvis.core.atoms import Atoms

minf = -300
maxf = 1000
step = 20

datafile = 'PhononData-Freq{}_{}_{}.json'.format(minf, maxf, step)


with open(datafile) as json_file:
    pdos_dict = json.load(json_file)
    

f1 = open('id_prop{}_{}_{}.csv'.format(minf, maxf, step), 'w+')

for d in pdos_dict:
    pos = 'POSCAR-{}.vasp'.format(d['jid'])
    atoms = Atoms.from_dict(d["atoms"])
    atoms.write_poscar(pos)    
    dos_str = ",".join(map(str, d['pdos_elast']))
    f1.write("%s,%s\n" % (pos, dos_str))
