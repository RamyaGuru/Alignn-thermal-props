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

run = '21'

folder = '../../run{}/run_{}_output/'.format(run, run)

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

plt.savefig('figures/{}_mae_traces.pdf'.format(run), bbox_inches = 'tight')
