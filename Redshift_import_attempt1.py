#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 14:52:10 2020

@author: jin
"""

file = open("Redshift_data.txt")
# splits data into lines
data = file.readlines()

# empty lists to organise the data
GRB_data_z = []
redshift = []
chosen_GRB = "970508"


for i in data:
    
    # splits data into rows
    row = i.split()
    GRB_name = row[0]
    
    # choses the data for the chosen GRB
    if GRB_name == chosen_GRB:
        
        redshift.append(row[9])
        GRB_data_z.append([row[0], row[9]])
        
    if GRB_name != chosen_GRB:
        continue

file.close()