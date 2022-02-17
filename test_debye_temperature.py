#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 16:22:44 2022

@author: rlg3


Test out the Debye temperature
"""

import numpy as np
from jarvis.db.figshare import data as jdata
import pandas as pd
from jarvis.core.atoms import Atoms
from jarvis.analysis.structure.spacegroup import Spacegroup3D
from jarvis.core.spectrum import Spectrum
from math import pi as pi
from jarvis.analysis.elastic.tensor import ElasticTensor


'''
Constants
'''

hbar = 1.0545718E-34
kB = 1.38064852E-23

def debye_temperature(atoms, vs):
    V0 = atoms.volume / sum([v for v in atoms.composition.to_dict().values()])
    print(atoms.volume)
    print(sum([v for v in atoms.composition.to_dict().values()]))
    print(V0)
    return (6 * pi**2 / (V0 * 1e-30)) ** (1/3) * vs * (hbar / kB)



dft_3d = jdata('dft_3d')
max_samples = 10


prop = 'elastic_tensor'

jid_list = ['JVASP-1002', 'JVASP-1103', 'JVASP-8185', 'JVASP-23'] #Si, PbTe, GaAs, CdTe

db_items=[]

for jid in jid_list:
    for i in dft_3d:
        if i['jid']==jid:
          db_items.append(i)


#Compare Debye temperatures

debyeT_jarvis = []
debyeT_toberer = []
vs_list = []

for x in db_items:
    et = ElasticTensor(x['elastic_tensor'])
    atoms = Atoms.from_dict(x['atoms'])
    vs = et.velocity_average(atoms)
    vs_list.append(vs)
    dt_jarvis = et.debye_temperature(atoms)
    dt_toberer = debye_temperature(atoms, vs)
    debyeT_jarvis.append(dt_jarvis)
    debyeT_toberer.append(dt_toberer)
    

