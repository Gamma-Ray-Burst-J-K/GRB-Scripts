#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import pandas as pd
import astropy.constants as c
import scipy
from scipy import stats
import math
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('font', family='serif')


grb = fits.open('gll_2flgc.fits')

table=grb[1].data
grb.close





# In[4]:



n=0
#data1 = [] #data with errors
#data2 = [] #data without errors

for i in table:
    if i['LUMINOSITY_DISTANCE'] !=0:
        n=n+1
        name = i['GCNNAME']
        z = i['REDSHIFT']
        dl = i['LUMINOSITY_DISTANCE']
        
        for k in i['LC_ENE_FLUX_ERR']:
            mask = np.where(i['LC_ENE_FLUX_ERR']!=0) #finding indexs of values where flux error is 0 and masking
            t_end = i['LC_END'][mask] 
            ene_flux = i['LC_ENE_FLUX'][mask] 
            ene_flux_err = i['LC_ENE_FLUX_ERR'][mask] 
            fluence = i['LC_FLUENCE'][mask] 
            flux = i['LC_FLUX'][mask] 
            flux_err = i['LC_FLUX_ERR'][mask] 
            indec = i['LC_INDEX'][mask]   # this is photon index, not spectral index;     
            #photon index = beta+1 e.g photon index 2 = beta +1, where beta = 1
            index_err = i['LC_INDEX_ERR'][mask] 
            median = i['LC_MEDIAN'][mask] 
            t_start = i['LC_START'][mask] 
            ts = i['LC_TS'][mask] 

        B = 1.9999999999 
        B=-1*indec   
        T_err = ((t_end-t_start)/2)
        T = t_start+T_err
        #print flux,ene_flux, dl, z, B,fluence
        lum = ene_flux*4*(np.pi)*(dl**2)*((1+z)**((B-1)-1))
        lum_err = ene_flux_err*4*(np.pi)*(dl**2)*((1+z)**((B-1)-1))

#        print lum,indec

        plt.figure(figsize = (12,8 ))
        plt.scatter(T/(1+z), lum, label=i['GCNNAME'] )
        plt.yscale('log')
        plt.xscale('log')
        plt.xlabel('log(T/(1+z)) [s]',fontsize = 18) 
        plt.ylabel('log(L) [erg/s]',fontsize = 18) 
        plt.title(f'Light Curves of LAT GRB{name}',fontsize = 18) 
        plt.tick_params(axis='both', which='major', labelsize=18)
        mpl.pyplot.errorbar(T/(1+z), lum, yerr=lum_err, xerr=T_err/(1+z), linestyle='',elinewidth=0.5)
        plt.legend()
        plt.show()
        exit()
            


# In[2]:


import matplotlib.pyplot as plt
plt.errorbar(1, 0.25, yerr=0.1, uplims=True, lolims=True, fmt='o')
plt.show()


# In[12]:


print(len(table['LC_INDEX'][0]))
print(len(table['LC_START'][0]))


# In[2]:


x = ['p','y','t','h','o','n']
print(x.index('o'))


# In[8]:


for i in table:
    if i['LUMINOSITY_DISTANCE'] !=0:
        #n=n+1
        name = i['GCNNAME']
        z = i['REDSHIFT']
        dl = i['LUMINOSITY_DISTANCE']
        #if np.all(i['LC_ENE_FLUX_ERR']>0) ==1 :
        for k in i['LC_ENE_FLUX_ERR']:
            mask = np.where(i['LC_ENE_FLUX_ERR']==0)
            t_end = i['LC_END'] #mutiple val
            


# In[11]:


mask


# In[9]:


t_end[mask]


# In[10]:


t_end


# In[13]:


indec = i['LC_INDEX']


# In[14]:


indec


# In[ ]:




