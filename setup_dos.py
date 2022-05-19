#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 18 22:27:10 2022

@author: rlg3
"""

from jarvis.io.phonopy.inputs import PhonopyInputs
from jarvis.core.atoms import Atoms
from jarvis.db.figshare import data as jdata

jid = 'JVASP-112289'

dft_3d = jdata("dft_3d")

match = next(d for d in dft_3d if d["jid"] == jid)

atoms = Atoms.from_dict(match['atoms'])

PI = PhonopyInputs(atoms)

PI.mesh_dos()