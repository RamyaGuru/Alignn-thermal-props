#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 18 22:27:10 2022

@author: rlg3
"""

from jarvis.io.phonopy.inputs import PhonopyInputs
from jarvis.core.atoms import Atoms
from jarvis.db.figshare import data as jdata
import os

#Change to directory with the phonon predictions
os.chdir('../../vasp_phonon_runs')

jid = ['JVASP-113634', 'JVASP-116263',\
       'JVASP-116267', 'JVASP-116273',\
           'JVASP-117476', 'JVASP-117483']


dft_3d = jdata("dft_3d")    
    
for j in jid:
    os.chdir('{}_PBEBO_FD_ELAST/MAIN-ELASTIC-bulk@{}'.format(j,j))    
    match = next(d for d in dft_3d if d["jid"] == j)
    
    atoms = Atoms.from_dict(match['atoms'])
    
    PI = PhonopyInputs(atoms)
    
    PI.mesh_dos()
    
    os.chdir('../../')