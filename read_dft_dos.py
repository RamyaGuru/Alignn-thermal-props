#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 19 07:53:26 2022

@author: rlg3


DFT DOS Thermal Proeprties
"""

import numpy as np
import pdos_integration as pint
from jarvis.core.atoms import Atoms
from jarvis.db.figshare import data as jdata
from jarvis.core.spectrum import Spectrum
import  matplotlib.pyplot as plt
from jarvis.core.specie import Specie

jid = 'JVASP-112289'

dft_3d = jdata("dft_3d")

match = next(d for d in dft_3d if d["jid"] == jid)

atoms = Atoms.from_dict(match['atoms'])

dos_file = 'output_files/total_dos.dat'

icm_to_thz = 2.99792458e-2


dos_array = np.genfromtxt(dos_file, skip_header = 1)


def transform_dos(match, dos_array):
    #Convert from THz to inv cm just so it works with methods
    dos_array[:,0] = dos_array[:,0] #/ icm_to_thz
    dos_array[:,1] = dos_array[:,1] #* icm_to_thz
    
    #Scale to 3N
    int_dos = pint.integrate_dos(dos_array[:,0], dos_array[:,1])
    form_unit = pint.get_natoms_form_unit(match)
    scale = (int_dos / form_unit) / 3.0
    dos_array[:,1] = dos_array[:,1] / scale 
    return dos_array


def get_molar_mass(atoms):
    mm = 0
    for k,v in atoms.composition.reduce()[0].items():
        el_sub = Specie(k)
        mm_sub = el_sub.atomic_mass * v
        mm = mm + mm_sub
    return mm


dos_array = transform_dos(match, dos_array)

zero_indx = np.where(dos_array[:,0] > 0)[0][0] - 1

min_indx = np.where(dos_array[:,0] > -300)[0][0] - 1

max_indx = np.where(dos_array[:,0] > 1000)[0][0] - 1

Cv = pint.heat_capacity(dos_array[zero_indx:,0], dos_array[zero_indx:,1], T = 300)

mm = get_molar_mass(atoms)

Cv_kg = Cv / mm * 1e3

plt.plot(dos_array[min_indx:max_indx, 0], dos_array[min_indx:max_indx, 1])
plt.xlabel(r'Frequency (cm$^{-1}$)')
plt.ylabel('DOS')
    
plt.savefig('test_dos_abridged.pdf', bbox_inches = 'tight')