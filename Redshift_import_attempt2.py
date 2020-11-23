#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 14:52:10 2020

@author: jin
"""

# open redshift data
with open("Redshift_Data.txt", encoding="ISO-8859-1") as file_2:
    # read the raw data
    data_2_raw = file_2.read()

# empty array for the data
data_2 = []

# split data into lines
for line in data_2_raw.split("\n"):
    if not line.strip():
            continue
    data_2.append(line.lstrip())


# empty lists to organise the data
GRB_data_2 = []
redshift = []

# loop to get the names and redshifts of GRBs
for i in data_2:
    
    # splits data into rows
    rows = i.split()
    
    # if there is data for the redshift, adds it to list
    if len(rows) > 10:
        redshift.append(rows[10])
        GRB_data_2.append([rows[0], rows[10]])
    
    # if there is no data for the redshift, adds None to list
    else:
        redshift.append(None)
        GRB_data_2.append([rows[0], None])
        
for j in GRB_data_2[1]:
    
    if j == 'Ê':
        GRB_data_2.append([rows[0], None])
    
    else:
        continue
    
print(GRB_data_2)