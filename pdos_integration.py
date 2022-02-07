#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 14:45:10 2022

@author: rlg3


Analyze phonon DOS


Source: Course Notes by B. Fultz
http://www.its.caltech.edu/~matsci/btfgrp/Review_Main.pdf



DOS normalized by 3N, where N is the number of atoms in conv. unit cell

atoms_from_dict
"""

import numpy as np
from jarvis.db.figshare import data as jdata
import pandas as pd
from jarvis.core.atoms import Atoms
from jarvis.analysis.structure.spacegroup import Spacegroup3D
from jarvis.core.spectrum import Spectrum
#from jarvis.core.composition

'''
Constants
'''
kB = 8.617333262145e-5 # eV/K
hbar = 6.582119569e-16 # eV*s
e = 1.60217662e-19
Na = 6.0221409e23

icm_to_eV = 1.23981e-4


def histo_making(max_freq,min_freq,step,s):
    rrange=(max_freq-min_freq)
    nbins=int(rrange/step)
    hist0,bins=np.histogram(s.x,bins=nbins,range=(min_freq,max_freq))
    hist=[]
    for ii in range(nbins+1):
        hist.append(0.0)
    for ii in range(1,nbins-1):
        if (hist0[ii] > 0.000001):
           intensita=0.0
           for kk in range(len(s.y)):
               if (s.x[kk] >= bins[ii]) and (s.x[kk] < bins[ii+1]):
                  intensita=intensita+s.y[kk]
           hist[ii]=intensita
    return(hist,bins)

def transform_normalized_dos(dft_3d, norm_dos : dict, dos_label = 'target'):
    for n in range(len(norm_dos)):
        jid = norm_dos[n]["id"]
        start = jid.find('JVASP')
        end = jid.find('.vasp')
        jid = jid[start:end]
        match = next(i for i in dft_3d if i["jid"] == jid)
        freq = np.arange(0, 1000, 5)
        s = Spectrum(x= freq, y=match['pdos_elast'])
        hist, bins = histo_making(1000, 0, 5, s)
        max_intensity = np.max(hist)
        #max_intensity = np.max(match['pdos_elast'])
        norm_dos[n][dos_label] = np.array(norm_dos[n][dos_label]) * max_intensity
    return norm_dos


def get_nmol_from_db_entry(p):
    atoms = Atoms(lattice_mat = p['atoms']['lattice_mat'],\
                  coords = p['atoms']['coords'], elements = p['atoms']['elements'])
    form_unit = atoms.composition.reduce()
    red_atoms = sum([v for v in form_unit[0].values()])  
    ratio = atoms.num_atoms / red_atoms
    return red_atoms
    
    

def integrate_dos(omega, dos):
    omega = omega * icm_to_eV
    dos = dos / icm_to_eV
    return np.trapz(dos, omega)



def vibrational_entropy(omega, dos, natoms, T = 300):
    #Convert inv. cm to eV
    omega = omega * icm_to_eV
    dos = dos / icm_to_eV
    x = (omega) / (kB * T)
    n = 1 / (np.exp(x[1:]) - 1)
    S_vib = kB  * ((n + 1) * np.log(n + 1) + n * np.log(n)) * dos[1:]
    S_vib = np.insert(S_vib, 0, 0)
    return np.trapz(S_vib, omega) * e * Na * natoms


def heat_capacity(omega, dos, natoms, T = 300):
    #Convert inv. cm to eV
    omega = omega * icm_to_eV / natoms
    dos = dos / icm_to_eV * natoms
    x= (omega) / (kB * T)
    #Drop omega = 0 term?
    Cp = kB * x[1:]**2 * (np.exp(x[1:]) / (np.exp(x[1:]) - 1)**2) * dos[1:] #removed factor of 3
    Cp = np.insert(Cp, 0, 0)
    return np.trapz(Cp, omega) * e * Na * natoms #multiply by #atoms in formula unit??
    
    
    
if __name__ == '__main__':
    dft_3d = jdata("edos_pdos")
    max_samples = 10
    
    jid = 'JVASP-19985'
    
    x = []
    for i in dft_3d:
        if i['jid'] == jid:
            x.append(i)
    
    jid_list = []
    int_DOS = []
    S_vib = []
    Cp = []
    
    p=x[0]
    
    jid = p["jid"]
    
    n_mols = get_nmol_from_db_entry(p)
    jid_list.append(jid)
    target = np.array(p['pdos_elast'])
    freq = np.linspace(0, 1000, len(target))
    int_DOS.append(integrate_dos(freq, target))
    S_vib.append(vibrational_entropy(freq, target, n_mols, T = 1000))
    Cp.append(heat_capacity(freq, target, n_mols, T = 1000))
    
    output = {'JID' : jid_list,
              'integrated_DOS' : int_DOS,
              'S_vib' : S_vib,
              'Cp' : Cp}
    
    
    df = pd.DataFrame(output)
    
    df.to_csv('thermal_props.csv')
    
    
    def vib_scaling(omega, T = 300):
        omega = omega * icm_to_eV
        x = (omega) / (kB * T)
        n = 1 / (np.exp(x[1:]) - 1)
        S_vib = - kB  * (np.log(n))  
        S_vib = np.insert(S_vib, 0, 0)
        return S_vib
    
    import matplotlib.pyplot as plt
    Svib_scaling = vib_scaling(freq)
    plt.plot(freq * icm_to_eV, Svib_scaling)
        
    