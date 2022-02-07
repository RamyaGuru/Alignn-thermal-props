#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 13:00:16 2022

@author: rlg3

Thermal Properties on Francesca dataset
"""

import pdos_integration as pint
import numpy as np
import json
import pandas as pd
from jarvis.db.figshare import data as jdata

from sklearn.metrics import mean_absolute_error


'''
Compute Integrals on a Training Dataset
'''

run = 'run11'
label = 'target'

with open('../' + run + '/temp/multi_out_predictions.json') as json_file:
    dos_dict = json.load(json_file)

dft_3d = jdata("edos_pdos")
dos_dict = pint.transform_normalized_dos(dft_3d, dos_dict, dos_label = 'target')

# x = []
# for i in dos_dict:
#     if i['target'] != 'na':
#         x.append(i)

jid_list = []
int_DOS = []
S_vib = []
Cp = []
for i in dos_dict:
    jid = i["id"]
    start = jid.find('JVASP')
    end = jid.find('.vasp')
    jid = jid[start:end]
    jid_list.append(jid)
    target = np.array(i['target'])
    freq = np.linspace(0, 1000, len(target))
    int_DOS.append(pint.integrate_dos(freq, target))
    S_vib.append(pint.vibrational_entropy(freq, target, 1))
    Cp.append(pint.heat_capacity(freq, target, 1))

true_output = {'JID' : jid_list,
          'integrated_DOS' : int_DOS,
          'S_vib (J/mol/K)' : S_vib,
          'Cp (J/mol/K)' : Cp}
df = pd.DataFrame(true_output)

df.to_csv(run + label + '_thermal_props.csv')



'''
Compute Integrals on an ALIGNN Phonon DOS
'''

run = 'run11'

label = 'predictions'

with open('../' + run + '/temp/multi_out_predictions.json') as json_file:
    dos_dict = json.load(json_file)


dos_dict = pint.transform_normalized_dos(dft_3d, dos_dict, dos_label = 'predictions')
# x = []
# for i in dos_dict:
#     if i['predictions'] != 'na':
#         x.append(i)

jid_list = []
int_DOS = []
S_vib = []
Cp = []
for i in dos_dict:
    jid = i["id"]
    start = jid.find('JVASP')
    end = jid.find('.vasp')
    jid = jid[start:end]
    jid_list.append(jid)
    target = np.array(i['predictions'])
    freq = np.linspace(0, 1000, len(target))
    int_DOS.append(pint.integrate_dos(freq, target))
    S_vib.append(pint.vibrational_entropy(freq, target, 1))
    Cp.append(pint.heat_capacity(freq, target, 1))

pred_output = {'JID' : jid_list,
          'integrated_DOS' : int_DOS,
          'S_vib (J/mol/K)' : S_vib,
          'Cp (J/mol/K)' : Cp}
df = pd.DataFrame(pred_output)

df.to_csv(run + label + '_thermal_props.csv')


'''
Get subset of "ground truth" values corresponding to test set
'''
sort_true_output = {'JID' : [],
                    'integrated_DOS' : [],
                    'S_vib (J/mol/K)' : [],
                    'Cp (J/mol/K)' : []}

for jid in pred_output['JID']:
    indx = true_output['JID'].index(jid)
    sort_true_output['JID'].append(true_output['JID'][indx])
    sort_true_output['integrated_DOS'].append(true_output['integrated_DOS'][indx])
    sort_true_output['S_vib (J/mol/K)'].append(true_output['S_vib (J/mol/K)'][indx])
    sort_true_output['Cp (J/mol/K)'].append(true_output['Cp (J/mol/K)'][indx])


'''
Calculate MAE
'''
MAE = {'integrated_DOS' : [], 'S_vib (J/mol/K)' : [], 'Cp (J/mol/K)' : []}

MAE_avg = []
for prop in ['integrated_DOS', 'S_vib (J/mol/K)', 'Cp (J/mol/K)']:
    MAE[prop] =  mean_absolute_error(sort_true_output[prop], pred_output[prop], multioutput = 'raw_values')

MAE['MAE_int'] = MAE.pop('integrated_DOS')
MAE['MAE_Svib'] = MAE.pop('S_vib (J/mol/K)')
MAE['MAE_Cp'] = MAE.pop('Cp (J/mol/K)')

df1 = pd.DataFrame(sort_true_output)
df2 = pd.DataFrame(MAE)

df_tot = pd.concat([df1, df2], axis=1)

df_tot.to_csv('mae_' + run + '_thermal_props.csv')