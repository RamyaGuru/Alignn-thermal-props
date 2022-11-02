#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 21:53:03 2022

@author: rlg3
"""

import numpy as np

#Plot the relative abundance of the top 10 elements for each set

#STABLE

stable_elem = np.load("stable_elem_freq.npz")

max_indx = np.argsort(stable_elem["arr_1"])[-10:]

elem = stable_elem["arr_0"][max_indx]

counts = np.array(stable_elem["arr_1"][max_indx])

total_elem = sum(stable_elem["arr_1"])

percentage = counts/total_elem



#UNSTABLE

unstable_elem = np.load("unstable_elem_freq.npz")

max_indx_2 = np.argsort(unstable_elem["arr_1"])[-10:]

counts_2 = np.array(unstable_elem["arr_1"][max_indx_2])

total_elem_2 = sum(unstable_elem["arr_1"])

elem_2 = unstable_elem["arr_0"][max_indx_2]

percentage_2 = counts_2 / total_elem