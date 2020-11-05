#!/usr/bin/env python
# coding: utf-8

# In[42]:


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
import corner


grb = fits.open('gll_2flgc.fits')

table=grb[1].data
grb.close


# In[59]:


Kc = []      #empty array for k corrections
Kc_std = []
n=0
redshift = []
diff_mean = []    #empty array for mean difference between corrected and uncorrected luminosities
diff_std = []
m=0
flux_diff = []

num1 = []
num2 = []
num3 = []



for i in table:
    if i['LUMINOSITY_DISTANCE'] !=0:
        n=n+1
        z = i['REDSHIFT']
        name = i['GCNNAME']
        dl = i['LUMINOSITY_DISTANCE']
        redshift.append(z)
        diff = []
        
        
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
            

        g = 1.6e-6
        B=-1*indec 
        k_correction = (1+z)**((B-1)-1)

        Kc.append(np.mean(k_correction))
        Kc_std.append(np.std(k_correction))
            
        lum = g*ene_flux*4*(np.pi)*(dl**2)*(k_correction)
        lum_err = g*ene_flux_err*4*(np.pi)*(dl**2)*(k_correction)
        T_err = ((t_end-t_start)/2)/(1+z)
        T = (t_start+T_err)/(1+z)
        
        t_1 = 5
        t_2 = 10
        t_3 = 20
        
        for j in lum:
            diff.append(abs((j/(4*(np.pi)*(dl**2)) - (g*ene_flux))))
            
        diff_mean.append(np.mean(diff))
        diff_std.append(np.std(diff))
        
        for i in T:
                if i >= (t_1-4) and i <=(t_1+4):
                    num1.append(i)
                    
                if i >= (t_2-5) and i <=(t_2+5):
                    num2.append(i)
                    
                if i >=(t_3-10) and i <=(t_3+10):
                    num3.append(i)
        
        if len(ene_flux)<=2:
            pass
        else:
            #print(len(ene_flux>1))
            plt.figure(figsize = (12,8 ))
            #print(n) 
            #plt.subplot(1,2,1)
            plt.scatter(T, lum, label=name )
            plt.yscale('log')
            plt.xscale('log')
            plt.xlabel('log(T-T0/(1+z)) [s]',fontsize = 14) 
            plt.ylabel('log(L) [erg/s]',fontsize = 14) 
            plt.title(f'Corrected Light Curves of LAT GRB{name}',fontsize = 14) 
            plt.tick_params(axis='both', which='major', labelsize=14)
            mpl.pyplot.errorbar(T, lum, yerr=lum_err, xerr=T_err, linestyle='',elinewidth=0.5)
            plt.axvline(t_1, label=f'{t_1}', color = 'r')
            plt.axvline(t_2, label=f'{t_2}', color = 'g')
            plt.axvline(t_3, label=f'{t_3}', color = 'b')
            plt.legend()
            plt.show()
        
            
            
            m=m+1
            
            
                    
            
                
        

        

print(len(num1))
print(len(num2))
print(len(num3))
print(m)


# In[54]:


plt.figure(figsize = (12,8 ))

plt.scatter(Kc, diff_mean)
plt.yscale('log')
plt.xscale('log')


plt.ylim(10e-10,10e-4)
#plt.xlim(0,26)
plt.xlabel('Mean k-correction',fontsize = 18) 
plt.ylabel('log(L/(4*(np.pi)*(dl**2)) - log(ene_flux) [erg/s]',fontsize = 18) 
plt.title('K-correction effect 1',fontsize = 18) 
#plt.tick_params(axis='both', which='major', labelsize=18)
#mpl.pyplot.errorbar(Kc, diff_mean, xerr=Kc_std, yerr=diff_std, linestyle='',elinewidth=0.5)
plt.show()


# In[55]:


plt.figure(figsize = (12,8 ))


plt.scatter(Kc, redshift)
plt.xlim(0,30)

plt.xlabel('Redshift',fontsize = 18) 
plt.ylabel('Mean k-correction',fontsize = 18) 
plt.title('K-correction effect 2',fontsize = 18) 


# In[52]:


plt.figure(figsize = (12,8 ))

z = [len(Kc), len(diff_mean)]



plt.contour([Kc, diff_mean],levels=[len(Kc)] )
#plt.yscale('log')
#plt.xscale('log')


plt.ylim(-10,1)
#plt.xlim(0,26)
plt.xlabel('Mean k-correction',fontsize = 18) 
plt.ylabel('Corrected - Not Corrected',fontsize = 18) 
plt.title('K-correction effect 1',fontsize = 18) 
#plt.tick_params(axis='both', which='major', labelsize=18)
#mpl.pyplot.errorbar(Kc, diff_mean, xerr=Kc_std, yerr=diff_std, linestyle='',elinewidth=0.5)
plt.show()


# In[46]:





# In[50]:





# In[ ]:





# In[ ]:





# In[13]:





# In[ ]:





# In[ ]:




