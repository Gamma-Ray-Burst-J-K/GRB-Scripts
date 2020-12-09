#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 11:19:51 2020

@author: jin
"""

"""
IMPORTS
"""
import matplotlib.pyplot as plt



"""
OLD DATA
"""

# open afterglow data
with open("Radio_data.txt") as file_1:
    # splits data into lines
    data_1 = file_1.readlines()

# empty lists for data
freq_all = []
freq_flux_all = []


# split data by line
for i in data_1:
    
    # splits data into rows
    row = i.split()

    freq_all.append(row[6])  # array of frequency values
    freq_flux_all.append([row[6], row[7]])

# empty list for data
freq = []
    
for i in freq_flux_all[1:]:
    
    # removes count if no frequency recorded
    if i[1] == '-99':
        continue
    
    else:
        freq.append(float(i[0]))
        


"""
NEW DATA
"""

# open afterglow data
with open("Radio_data_new.txt") as file_2:
    # splits data into lines
    data_2 = file_2.readlines()
    
# empty lists for data
freq_all_new = []
freq_flux_all_new = []


# split data by line
for i in data_2:
    
    # splits data into rows
    row = i.split()

    freq_all_new.append(row[4])  # array of frequency values
    freq_flux_all_new.append([row[4], row[5]])



# empty list for data
freq_new = []
    
for i in freq_flux_all_new[1:]:
    
    # removes count if no frequency recorded
    if i[1] == '-99':
        continue
    
    else:
        freq_new.append(float(i[0]))



"""
PLOTS
"""

# plot histogram of frequencies for new data
plt.title('Histogram of frequencies for old data')
fig = plt.subplot
plt.hist(freq, bins = [0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024])
plt.xlabel("Frequency [GHz]")
plt.xscale("log")
plt.ylabel("No. of counts")
plt.show()      

# plot histogram of frequencies for old data
plt.title('Histogram of frequencies for updated data')
fig = plt.subplot
plt.hist(freq_new, bins = [0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024])
plt.xlabel("Frequency [GHz]")
plt.xscale("log")
plt.ylabel("No. of counts")
plt.show()