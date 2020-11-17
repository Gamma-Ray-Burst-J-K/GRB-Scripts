#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 18:44:06 2020

@author: jin
"""

# import packages
import numpy as np
import matplotlib.pyplot as plt

# import data from file
file = open("Radio_data.txt")
# splits data into lines
data = file.readlines()

# chosen GRB for light curve
chosen_GRB = "970508"

#Parameters for this GRB
z = 0.835 # redshift
beta = 1/3 # spectral index
d_l = 4271.0*3.08568e24 # luminosity distance in cm

# empty arrays to organise the data
time_array = []
freq_array = []
flux_array = []
flux_err_array = []
GRB_data = []

# for loop to itterate over data
for i in data:
    
    # splits data into rows
    row = i.split()
    #print(row)
    GRB_name = row[0]
    #print(GRB_name)
    
    # choses the data for the chosen GRB
    if GRB_name == chosen_GRB:
        
        time_array.append(row[5])
        freq_array.append(row[6])
        flux_array.append(row[7])
        flux_err_array.append(row[8])
        GRB_data.append([row[5], row[6], row[7], row[8]])
        
    if GRB_name != chosen_GRB:
        continue

time = []
flux = []
flux_err = []

for sublist in GRB_data:
    if sublist[2] == "-99.":
        continue
    if sublist[1] == "4.86":
        time.append(sublist[0])
        flux.append(sublist[2])
        flux_err.append(sublist[3])

print(time)
print(flux)
print(flux_err)

plt.scatter(time, flux)
plt.errorbar(time, flux, yerr = flux_err, fmt = ' ')
plt.xscale("log")
plt.xlabel("Time/days")
plt.yscale("log")
plt.ylabel("Flux")
plt.show()

#for j in flux:
    
#    def L(j):
#        L = j*4*np.pi*(d_l**2)*(1+z)**(beta-1)
#        return L

#    lum = L(j)

#print(lum)

file.close()