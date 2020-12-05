#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 14:52:10 2020

@author: jin
"""

# imports
import numpy as np



"""
REDSHIFTS
"""
# open redshift data

GRB_redshift=[]
Redshift=[]

with open('Redshifts.txt') as fp:

    for line in fp.readlines():
        row=line.split()
        try:
            if float(row[-1]) <10:  #Redshift if present is last entry of row and is a number, but last value may be a GCN/IAU circular number, redshift will always be a value less than 10 so only keep value if <10.
                GRB_redshift.append(row[0])
                Redshift.append(float(row[-1]))
        except:
            pass   



    
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
d_l_data = []

for z in Redshift:

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
    
    d_l_data.append(d_l)