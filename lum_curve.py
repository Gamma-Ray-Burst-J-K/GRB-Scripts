#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 14:33:18 2020

@author: jin
"""

# import packages
import numpy as np
import matplotlib.pyplot as plt



# luminosity function
def L(flux, d_l, z, beta):
    
    # empty array for the luminosity
    lum = []
    
    # define luminosity function
    for j in flux:
        
        L = j*4*np.pi*(d_l**2)*(1+z)**(beta-1)
        lum.append(L)
        
    return lum

# luminosity error function
def L_err(flux_err, d_l, z, beta):

    # empty array for the error
    lum_err = []
    
    # define error function
    for k in flux_err:
        
        err = k*4*np.pi*(d_l**2)*(1+z)**(beta-1)
        lum_err.append(err)
        
    return lum_err



# import data from file
file = open("Radio_data.txt")
# splits data into lines
data = file.readlines()

# chosen GRB for light curve with redshift
chosen_GRB = "970508"

#Parameters for this GRB
z = 0.835 # redshift
beta = 1/3 # spectral index
d_l = 4271.0*3.08568e24 # luminosity distance in cm



# empty lists to organise the data
time_array = []
freq_array = []
flux_array = []
flux_err_array = []
GRB_data = []

# for loop to itterate over data
for i in data:
    
    # splits data into rows
    row = i.split()
    GRB_name = row[0]
    
    # choses the data for the chosen GRB
    if GRB_name == chosen_GRB:
        
        time_array.append(row[5])
        freq_array.append(row[6])
        flux_array.append(row[7])
        flux_err_array.append(row[8])
        GRB_data.append([row[5], row[6], row[7], row[8]])
        
    if GRB_name != chosen_GRB:
        continue



# new empty lists for the data
time_l = []
flux_l = []
flux_err_l = []

for sublist in GRB_data:
    
    # when there is no data recorded for this time
    if sublist[2] == "-99.":
        
        continue
    
    # chosen frequency
    if sublist[1] == "8.46":
        time_l.append(sublist[0])
        flux_l.append(sublist[2])
        flux_err_l.append(sublist[3])

# convert to list of floats
time = [float(x) for x in time_l]
flux = [float(x) for x in flux_l]
flux_err = [float(x) for x in flux_err_l]



# call on functions to convert to luminosity
lum = L(flux, d_l, z, beta)
lum_err = L_err(flux_err, d_l, z, beta)

# plot scatter graph of the flux
plt.scatter(time, flux)
plt.errorbar(time, flux, yerr = flux_err, fmt = ' ')
plt.xscale("log")
plt.xlabel("Time/days")
plt.yscale("log")
plt.ylabel("Flux")
plt.show()

# plot scatter graph of the luminosity
plt.scatter(time, lum)
plt.errorbar(time, lum, yerr = lum_err, fmt = ' ')
plt.xscale("log")
plt.xlabel("Time/days")
plt.yscale("log")
plt.ylabel("Luminosity")
plt.show()

# Close file
file.close()