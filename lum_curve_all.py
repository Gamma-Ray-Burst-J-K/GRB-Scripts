#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 14:33:18 2020

@author: jin
"""

"""
IMPORTS
"""
import numpy as np
import re
import matplotlib.pyplot as plt



"""
CONSTANTS
"""
beta = 1/3      # spectral index
H0 = 70         # Hubble parameter
WM = 0.3        # Omega(matter)
WV = 0.7        # Omega(Vaccuum)
WR = 0.         # Omega(radiation)
WK = 0.         # Omega curvaturve = 1-Omega(total)
c = 299792.458  # velocity of light in km/sec
Tyr = 977.8     # coefficent for converting 1/H into Gyr
DTT = 0.5       # time from z to now in units of 1/H0
age = 0.5       # age of Universe in units of 1/H0
zage = 0.1      # age of Universe at redshift z in units of 1/H0
a = 1.0         # 1/(1+z), the scale factor of the Universe
az = 0.5        # 1/(1+z(object))



"""
FUNCTIONS
"""
# luminosity function
def L(flux, d_l, z, beta):
    
    # empty array for the luminosity
    lum = []
    
    # define luminosity function
    for i in flux:
        
        F = i*1e-29 # convert units to erg s^-1 cm^-2 Hz^-1
        L = F*4*np.pi*(d_l**2)*(1+z)**(beta-1)
        lum.append(L)
        
    return lum

# luminosity error function
def L_err(flux_err, d_l, z, beta):

    # empty array for the error
    lum_err = []
    
    # define error function
    for j in flux_err:
        
        # convert units to erg s^-1 cm^-2 Hz^-1
        F_err = j*1e-29
        
        err = F_err*4*np.pi*(d_l**2)*(1+z)**(beta-1)
        lum_err.append(err)
        
    return lum_err



"""
REDSHIFTS
"""

# open afterglow data
with open("Radio_data.txt") as file_1:
    # splits data into lines
    data_1 = file_1.readlines()

names = []
GRB_sample = []         # names of all GRBs in sample

# generate list of GRBs in sample
for i in data_1:
    
    # splits data into rows
    row = i.split()
    names.append(row[0])

for GRB in names:
    if GRB not in GRB_sample:
        GRB_sample.append(GRB)
        
    
# open redshift data

GRB_redshift=[]
Redshift= []        # all redshifts
sample= []
z_sample = []       # redshifts in sample

with open('Redshifts.txt') as fp:

    for line in fp.readlines():
        row=line.split()
        
        try:
            if float(row[-1]) <10:
                
                z = float(row[-1])
                GRB_redshift.append([row[0], z])     # all GRBs with redshift
                Redshift.append(z)
                
        except:
            pass   


GRB_z = [item for item in GRB_redshift if item[0] in GRB_sample]


z_all = []          # redshifts of GRBs in sample
GRB_all = []        # GRBs in sample with redshifts

for data in GRB_z:
    z_all.append(data[1])
    GRB_all.append(data[0])



"""
LUMINOSITY DISTANCE
"""

# Define constants
H0 = 70         # Hubble parameter
WM = 0.3        # Omega(matter)
WV = 0.7        # Omega(Vaccuum)
WR = 0.         # Omega(radiation)
WK = 0.         # Omega curvaturve = 1-Omega(total)
c = 299792.458  # velocity of light in km/sec
Tyr = 977.8     # coefficent for converting 1/H into Gyr
DTT = 0.5       # time from z to now in units of 1/H0
age = 0.5       # age of Universe in units of 1/H0
zage = 0.1      # age of Universe at redshift z in units of 1/H0
a = 1.0         # 1/(1+z), the scale factor of the Universe
az = 0.5        # 1/(1+z(object))
h = H0/100.
WR = 4.165E-5/(h*h)   # includes 3 massless neutrino species, T0 = 2.72528
WK = 1-WM-WR-WV
age = 0.
n=1000         # number of points in integrals

# empty list for data
d_l_all = []        # all luminosity distances

for z in z_all:

    az = 1.0/(1+1.0*z)
    
    for i in range(n):
        a = az*(i+0.5)/n
        adot = np.sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
        age = age + 1./adot

    zage = az*age/n
    zage_Gyr = (Tyr/H0)*zage
    DTT = 0.0
    DCMR = 0.0

    # do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
    for i in range(n):
        a = az+(1-az)*(i+0.5)/n
        adot = np.sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
        DTT = DTT + 1./adot
        DCMR = DCMR + 1./(a*adot)

    DTT = (1.-az)*DTT/n
    DCMR = (1.-az)*DCMR/n
    age = DTT+zage
    age_Gyr = age*(Tyr/H0)
    DTT_Gyr = (Tyr/H0)*DTT
    DCMR_Gyr = (Tyr/H0)*DCMR
    DCMR_Mpc = (c/H0)*DCMR

    # tangential comoving distance

    ratio = 1.00
    x = np.sqrt(abs(WK))*DCMR
    if x > 0.1:
        if WK > 0:
            ratio =  0.5*(np.exp(x)-np.exp(-x))/x 
        else:
            ratio = np.sin(x)/x
    else:
        y = x*x
        if WK < 0: y = -y
        ratio = 1. + y/6. + y*y/120.
        DCMT = ratio*DCMR
        DA = az*DCMT
        DA_Mpc = (c/H0)*DA
        kpc_DA = DA_Mpc/206.264806
        DA_Gyr = (Tyr/H0)*DA
        DL = DA/(az*az)
        DL_Mpc = (c/H0)*DL

    # luminosuity distance in cm
    d_l = DL_Mpc * 3.08568e24
    
    d_l_all.append(d_l)

GRB_z_d_l = []      # GRB name with corresponding redshift and lum distance

for name, z, d_l in zip(GRB_all, z_all, d_l_all):
    GRB_z_d_l.append([name, z, d_l])
    



"""
ALL FLUX CURVES
"""

# empty array to hold data for 8.46GHz
GRB_name = []
time_array = []
flux_array = []
flux_err_array = []
data_freq = []
  
# empty array to hold data for 4.86GHz
GRB_name_4 = []
time_array_4 = []
flux_array_4 = []
flux_err_array_4 = []
data_freq_4 = []



# split data by GRB
for i in data_1:
    
    # splits data into rows
    row = i.split()
    freq = row[6]
    
    # Gets only data points for 8.46GHz
    if freq == '8.46':
        
        # removes no recorded data entries
        if row[7] == "-99.":
            continue
        
        # adds remaining data points to array
        else:
            GRB_name.append(row[0])
            time_array.append(float(row[5]))
            flux_array.append(float(row[7]))
            flux_err_array.append(float(row[8]))
            data_freq.append([row[0], float(row[5]), float(row[7]), float(row[8])])
    
    # Gets only points for 4.86GHz
    if freq == '4.86':
        
        # removes no recorded data entries
        if row[7] == "-99.":
            continue
        
        # adds remaining data points to array
        else:
            GRB_name_4.append(row[0])
            time_array_4.append(float(row[5]))
            flux_array_4.append(float(row[7]))
            flux_err_array_4.append(float(row[8]))
            data_freq_4.append([row[0], float(row[5]), float(row[7]), float(row[8])])
        
    else:
        continue



# first GRB in the radio data with a redshift
first_GRB =  GRB_name[0] 

# empty arrays to store data
time = []
flux = []
flux_err = []
GRB_plot = []

# loop through rows in data for 8.46GHz
for row in data_freq:
    
    # get name of GRB for that row
    GRB = row[0]
    
    # if GRB matches previous, add it to the array
    if re.match(GRB, first_GRB):
        
        # Adds data to array for that GRB
        time.append(row[1])
        flux.append(row[2])
        flux_err.append(row[3])

    
    # when it stops matching, plot the data and reset the arrays  
    else:
        
        if len(flux) >= 3:
            
            """
            # flux curve plot
            fig = plt.subplot()
            plt.title(f'Flux in the 8.56GHz band for GRB{GRB}')
            plt.scatter(time, flux)
            plt.errorbar(time, flux, yerr = flux_err, fmt = ' ')
            plt.xscale("log")
            plt.xlabel("Time [days]")
            plt.yscale("log")
            plt.ylabel(r'Flux in 8.5GHz band [$\mu$Jy]')
            plt.show()
            """
            
            # identify redshift and luminosity distance for that GRB
            for j in GRB_z_d_l:
                
                x = j[0]
                
                if re.match(GRB, x):
                    
                    z = j[1]
                    d_l = j[2]
                    
                    # call on luminosity function to plot lum curves
                    lum = L(flux, d_l, z, beta)
                    lum_err = L_err(flux_err, d_l, z, beta)
                    
                    # plot scatter graph of the luminosity
                    plt.title('Luminosity in the 8.56GHz band for GRB{GRB}')
                    plt.scatter(time, lum)
                    plt.errorbar(time, lum, yerr = lum_err, fmt = ' ')
                    plt.xscale("log")
                    plt.xlabel("Time [days]")
                    plt.yscale("log")
                    plt.ylabel("Luminosity in 8.5GHz band [erg s^-1]")
                    plt.show()
                    
                    # Adds GRB to list if it gets plot
                    GRB_plot.append(GRB)
            
        
        # reset the arrays
        time = []
        flux = []
        flux_err = []
        
        # move to next GRB
        first_GRB = GRB
        continue
   
print(len(GRB_plot))
