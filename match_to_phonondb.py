#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 13:15:12 2022

@author: rlg3


Compare the ALIGNN Prediction Results to the PhononDB database

This script will get matches first.

Then, the DOS will be downloaded form phonondb.
"""

from jarvis.db.figshare import data
import json
import re


run = 'run21'
pred_file = '../../{}/pred_data.json'.format(run)

pdb_file = '../../writeup/ALIGNN_Phonons/phonondb/phondb.json'

dft_3d = data("dft_3d")

#Load data from the prediction set
with open(pred_file) as json_file:
    pred_data = json.load(json_file)
    
    
#Load data from the phonon database
with open(pdb_file) as json_file_2:
    phonondb = json.load(json_file_2)

pdb_mpids = []

for p in phonondb:
    reference = p['ref']
    start = reference.rfind('/')
    start = start + 1
    end = reference.find('.html')
    mpid = reference[start:end]
    pdb_mpids.append(mpid)
    
    
pred_mpids = []

for pd in pred_data:
    jarvis_entry = dft_3d[int(pd["id"])]
    mpid_pd = jarvis_entry['reference']
    pred_mpids.append(mpid_pd)
    

overlap = set(pdb_mpids) & set(pred_mpids)


#Write the matched mpids to a file



    
