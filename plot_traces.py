#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 12:07:52 2022

@author: rlg3

Plot the training and validation traces
"""


import json
import matplotlib.pyplot as plt
import numpy as np

folder = '../../run20/run_20_output/'

fval = open(folder + 'history_val.json')


val = json.load(fval)

ftrain = open(folder + 'history_train.json')

train = json.load(ftrain)

epochs = np.arange(1, 601)

'''
Plot MAE trace
'''


plt.figure()

plt.plot(epochs, val['mae'], label = 'Validation', color='xkcd:black')

plt.plot(epochs, train['mae'], label = 'Training', color = 'xkcd:red')

plt.ylabel('MAE')
plt.xlabel('Epochs')

#plt.ylim([0,10])

plt.legend()

plt.savefig('mae_traces_dos_orig.pdf', bbox_inches = 'tight')
