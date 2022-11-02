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
from phonopy.structure.atoms import isotope_data
from math import pi, cos, sin, asin, acos
from jarvis.analysis.elastic.tensor import ElasticTensor
#from jarvis.core.composition

'''
Constants
'''
kB = 8.617333262145e-5 # eV/K
hbar = 6.582119569e-16 # eV*s
e = 1.60217662e-19
Na = 6.0221409e23

icm_to_eV = 1.23981e-4

icm_to_thz = 2.99792458e-2


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
        scale = get_natoms_from_db_entry(match)
        freq = np.arange(0, 1000, 5)
        s = Spectrum(x= freq, y=match['pdos_elast'])
        hist, bins = histo_making(1000, 0, 5, s)
        max_intensity = np.max(hist)
        #max_intensity = np.max(match['pdos_elast'])
        norm_dos[n][dos_label] = np.array(norm_dos[n][dos_label]) * max_intensity * scale
    return norm_dos


def transform_dos(dft_3d, norm_dos : dict, dos_label = 'target', unit_conv = 1):
    for n in range(len(norm_dos)):
        jid = norm_dos[n]["id"]
        start = jid.find('JVASP')
        end = jid.find('.vasp')
        jid = jid[start:end]
        match = next(i for i in dft_3d if i["jid"] == jid)
        scale = get_natoms_from_db_entry(match)
        #max_intensity = np.max(match['pdos_elast'])
        norm_dos[n][dos_label] = np.array(norm_dos[n][dos_label]) * unit_conv * scale
    return norm_dos    

def get_natoms_from_db_entry(p):
    atoms = Atoms.from_dict(p['atoms'])
    num_atoms = atoms.num_atoms
    try:
        spg = Spacegroup3D(atoms)
        cvn_atoms = spg.conventional_standard_structure.num_atoms
    except:
        cvn_atoms = num_atoms
    formula = atoms.composition.reduce()
    form_atoms = sum([v for v in formula[0].values()])
    #scale = formula[1]
    scale = (num_atoms / cvn_atoms) * form_atoms
    return scale
    
def get_natoms_form_unit(p):
    atoms = Atoms.from_dict(p['atoms'])
    num_atoms = atoms.num_atoms
    try:
        spg = Spacegroup3D(atoms)
        cvn_atoms = spg.conventional_standard_structure.num_atoms
    except:
        cvn_atoms = num_atoms
    formula = atoms.composition.reduce()
    form_atoms = sum([v for v in formula[0].values()])
    #scale = formula[1]
    scale = (num_atoms / cvn_atoms) * form_atoms
    return form_atoms

def debye_DOS(p, omega):
    '''
    Run up to the Debye temperature. Remaining values should be zeros.
    '''
    #Debye DOS
    #p is an element of the dft_3d database
    atoms = Atoms.from_dict(p['atoms'])
    et = ElasticTensor(p['elastic_tensor'])
    vs = et.velocity_average(atoms)
    omegaD = et.debye_temperature_toberer(atoms) * kB
    omegaD = omegaD / icm_to_eV
    print(omegaD)
    omega_new = omega[omega < omegaD]
    dos = 3 * omega_new ** 2 / (2 * pi ** 2 * vs ** 3)
    pad0 = len(omega) - len(dos)
    dos = np.pad(dos, (0,pad0))
    int_dos = integrate_dos(omega, dos)
    form_unit = get_natoms_form_unit(p)
    scale = (int_dos / form_unit) / 3.0
    dos = dos / scale
    return dos


def debye_dispersion(p):
    atoms = Atoms.from_dict(p['atoms'])
    et = ElasticTensor(p['elastic_tensor'])
    vs = et.velocity_average(atoms)  
    atmV = (atoms.volume / atoms.num_atoms) * 1e-30
    kmax = (6 * pi ** 2 / (atmV))**(1 / 3)
    k_list = np.linspace(0, kmax, 50)
    omega_list = vs * k_list
    return k_list / kmax, omega_list / 1e12 / (2 * pi)


def BvK_DOS(p, omega):
    '''
    Source: Supplemental Equation 
    '''
    #BvK DOS
    atoms = Atoms.from_dict(p['atoms'])
    et = ElasticTensor(p['elastic_tensor'])
    vs = et.velocity_average(atoms)
    atmV = (atoms.volume / atoms.num_atoms) * 1e-30
    form_unit = get_natoms_form_unit(p)
    kmax = (6 * pi ** 2 / (atmV))**(1 / 3)
    omega0 = (2 / pi) * vs * kmax * hbar
    omega_max = omega0 / icm_to_eV
    vs = et.velocity_average(atoms)
    # omegaD = et.debye_temperature_toberer(atoms) * kB
    # omegaD = omegaD / icm_to_eV
    omega_new = omega[omega < omega_max]
    omega_new[0] = 1e-10 
    k_omega = (2 * kmax / pi) * np.arcsin(omega_new / (omega0 / icm_to_eV))
    vg_bvk = (omega0 / hbar) * (pi / (2 * kmax)) *\
            np.cos(pi * k_omega / (2 * kmax))
    vp_bvk = (omega0 / hbar /  k_omega) * np.sin(pi * k_omega/ (2 * kmax))
    dos = 3 * (omega_new)** 2 / (2 * pi ** 2 * vg_bvk * vp_bvk ** 2)
    pad0 = len(omega) - len(dos)
    dos = np.pad(dos, (0,pad0))
    int_dos = integrate_dos(omega, dos)
    scale = (int_dos / form_unit) / 3.0
    dos = dos / scale 
    return dos


def BvK_DOS_2(p, omega):
    atoms = Atoms.from_dict(p['atoms'])
    et = ElasticTensor(p['elastic_tensor'])
    vs = et.velocity_average(atoms)
    atmV = (atoms.volume / atoms.num_atoms) * 1e-30
    form_unit = get_natoms_form_unit(p)
    kmax = (6 * pi ** 2 / (atmV))**(1 / 3)
    omega0 = (2 / pi) * vs * kmax * hbar
    omega0 = omega0 / icm_to_eV
    print(omega0)
    omega_new = omega[omega < omega0]
    dos = (1 / (2 * pi ** 2)) * ((2 / pi) *\
        np.arcsin(omega_new / (omega0)))**2 * (2 / pi) *\
        (1 / omega0) * (1 / (1 - (omega_new / omega0) ** 2) ** (1/2))  
    pad0 = len(omega) - len(dos)
    dos = np.pad(dos, (0,pad0))
    int_dos = integrate_dos(omega, dos)
    scale = (int_dos / form_unit) / 3.0
    dos = dos / scale 
    return dos

def BvK_dispersion(p):
    atoms = Atoms.from_dict(p['atoms'])
    et = ElasticTensor(p['elastic_tensor'])
    vs = et.velocity_average(atoms)  
    atmV = (atoms.volume / atoms.num_atoms) * 1e-30
    kmax = (6 * pi ** 2 / (atmV))**(1 / 3)
    k_list = np.linspace(0, kmax, 50)
    omega_list = (2 / pi) * vs * kmax * np.sin(pi * k_list / (2 * kmax))
    return k_list / kmax, omega_list / 1e12 / (2 * pi)   

def integrate_dos(omega, dos):
    omega = omega * icm_to_eV
    dos = dos / icm_to_eV
    return np.trapz(dos, omega)



def vibrational_entropy(omega, dos, T = 300):
    #Convert inv. cm to eV
    omega = omega * icm_to_eV
    dos = dos / icm_to_eV
    x = (omega) / (kB * T)
    n = 1 / (np.exp(x[1:]) - 1)
    S_vib = kB  * ((n + 1) * np.log(n + 1) - n * np.log(n)) * dos[1:]
    S_vib = np.insert(S_vib, 0, 0)
    return np.trapz(S_vib, omega) * e * Na

def vibrational_entropy_scaling(omega, T = 300):
    omega = omega * icm_to_eV
    x = (omega) / (kB * T)
    n = 1 / (np.exp(x[1:]) - 1)
    S_vib = kB  * ((n + 1) * np.log(n + 1) - n * np.log(n))
    S_vib = np.insert(S_vib, 0, S_vib[0])
    return S_vib  

def vibrational_entropy_scaling_trig(omega, T = 300):
    omega = omega * icm_to_eV
    x = (omega) / (kB * T)
    print(x)
    Svib = kB * ((x / 2) * (1 / np.tanh(x / 2)) - np.log(2 * np.sinh(x / 2)))
    return Svib

#def vibrational_entropy_trig(omega, dos, T = 300):
    

def heat_capacity(omega, dos, T = 300):
    #Convert inv. cm to eV
    omega = omega * icm_to_eV
    dos = dos / icm_to_eV
    x= (omega) / (kB * T)
    #Drop omega = 0 term?
    Cp = kB * x[1:]**2 * (np.exp(x[1:]) / (np.exp(x[1:]) - 1)**2) * dos[1:] #removed factor of 3
    Cp = np.insert(Cp, 0,0)
    return np.trapz(Cp, omega) * e * Na #multiply by #atoms in formula unit??

def heat_capacity_scaling(omega, T = 300):
    omega = omega * icm_to_eV
    x= (omega) / (kB * T)
    #Drop omega = 0 term?
    Cp = kB * x[1:]**2 * (np.exp(x[1:]) / (np.exp(x[1:]) - 1)**2)#removed factor of 3
    Cp = np.insert(Cp, 0, Cp[0])
    return Cp   

    
def heat_capacity_scaling_trig(omega, T = 300):
    omega = omega * icm_to_eV
    x = (omega) / (kB * T)
    Cp = kB * (x / 2)**2 * (1 / np.sinh(x / 2))**2
    #Cp = np.insert(Cp, 0, Cp[0])
    return Cp
    
def isotopic_gamma(p):
    atoms = Atoms.from_dict(p['atoms'])
    formula = atoms.composition.reduce()
    natoms = sum([v for v in formula[0].values()])
    ave_m = 0
    gamma = 0
    for k,v in formula[0].items():
        iso_list = isotope_data[k]
        ave_m_n = sum([iso[2] * iso[1] for iso in iso_list])
        gamma_n = sum([iso[2] * (iso[1] - ave_m_n)**2 for iso in iso_list])
        ave_m += ave_m_n * (v / natoms)
        gamma += gamma_n * (v / natoms)
    return gamma / (ave_m ** 2)

def isotopic_tau(p, omega, dos):
    gamma = isotopic_gamma(p)
    atoms = Atoms.from_dict(p['atoms'])
    atmV = (atoms.volume / atoms.num_atoms) * 1e-30
    omega = omega * icm_to_thz #* 1e12
    dos = dos / icm_to_thz / (atmV * atoms.num_atoms) #normalize by unit cell? Already has factor of 3?
    tau = (pi / 6) * (atmV * gamma * omega ** 2) * dos
    return np.trapz(tau, omega) * 1e12 # / 1e12) #Integrate over omega THz?

   

def isotopic_tau_scaling(p, omega):
    gamma = isotopic_gamma(p)
    atoms = Atoms.from_dict(p['atoms'])
    atmV = (atoms.volume / atoms.num_atoms) * 1e-30
    omega = omega * icm_to_thz #* 1e12
    tau = (pi / 6) * (atmV * gamma * omega ** 2)
    return tau # / 1e12) #Integrate over omega THz?   


def check_dynamical_stability(omega, dos):
    zero_indx = np.where(omega > 0)[0][0] - 1
    neg_freq = np.linspace(omega[0], 0, zero_indx)
    intdos_neg = np.trapz(neg_freq, dos[:zero_indx])
    intdos_full = np.trapz(omega, dos)
    stable = 1
    if intdos_neg / intdos_full > 0.1:
        stable = 0
    return stable
    
    
    
if __name__ == '__main__':
    dft_3d = jdata("dft_3d")
    max_samples = 10
    
    jid = 'JVASP-1076'
    
    x = []
    for i in dft_3d:
        if i['jid'] == jid:
            x.append(i)
    
    jid_list = []
    int_DOS = []
    S_vib = []
    Cp = []
    tau = []
    
    p=x[0]
    
    jid = p["jid"]
    
    #scale = get_scale_quantity_from_db_entry(p)
    # jid_list.append(jid)
    # target = np.array(p['pdos_elast'])
    freq = np.linspace(0, 1000, 200)
    # int_DOS.append(integrate_dos(freq, target))
    # S_vib.append(vibrational_entropy(freq, target, T = 1000))
    # Cp.append(heat_capacity(freq, target, T = 1000))
    # tau.append(isotopic_tau(p, freq, target))
    
    debye_dos = debye_DOS(p, freq)
    bvk_dos = BvK_DOS_2(p,freq)
    
    output = {'JID' : jid_list,
              'integrated_DOS' : int_DOS,
              'S_vib' : S_vib,
              'Cp' : Cp}
    
    
    #df = pd.DataFrame(output)
    
    #df.to_csv('thermal_props.csv')
    
    
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
        
    gamma = isotopic_gamma(x[0])
