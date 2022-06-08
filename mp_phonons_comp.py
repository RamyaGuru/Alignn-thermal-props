#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 12:26:50 2022

@author: rlg3


Compare predictions to the MP phonons database on matminer
"""

from matminer.datasets import load_dataset

phondb = load_dataset('phonon_dielectric_mp')