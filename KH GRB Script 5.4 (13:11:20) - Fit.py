#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import pandas as pd
import astropy.constants as c
import scipy
from scipy import stats, optimize
import math
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib as mpl

from pylab import *
from scipy import *
from scipy.optimize import curve_fit

import xlsxwriter



mpl.rc('font', family='serif')


grb = fits.open('gll_2flgc.fits')

table=grb[1].data
grb.close


# In[3]:


def power_law1(x, N, a): 
    '''normal power law 
    L = Nx^-a'''
    return N*np.power(x, -a)


# In[4]:


def power_law2(x, N, a, b):
    
    '''broken power law (a piecewise function)
    L = Nx^-a1 for x<xth
    L = Nx^-(a1-a2) for x>xth'''
    xth = x.max()
    
    for k in x:
        if k <=xth:
            return N*np.power(x, -a)
        else:
            return N*np.power(x, -(a-b))


# In[ ]:





# In[24]:


n=0
redshift = []
m=0
Ts = 10
parameters1 = []
parameters2 = []


#def luminosity(table):
Kc = []      #empty array for k corrections
Kc_std = []
n=0
redshift = []
m=0
Ts = 10
j = 0


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
        B_err = -1*index_err
        
        #calculated weighted mean of spectral index and error
        B_sum = []
        B_err_sum = []
        
        if len(B)==1:
            mean_B = B
            mean_B_err = abs(B_err)
        else:
            
            mean_B = np.average(B, weights=B_err)
            mean_B_err = abs(np.average(B_err))
                
                
        
        k_correction = (1+z)**((B-1)-1)

        lum = g*ene_flux*4*(np.pi)*(dl**2)*(k_correction)
        lum_err = g*ene_flux_err*4*(np.pi)*(dl**2)*(k_correction)

        T_err = ((t_end-t_start)/2)
        T = (t_start+T_err)/(1+z)
        T_err =  T_err/(1+z)


        #calculating the index of T either side of Ts
        diff = abs(T-Ts)
        idx1 = np.argmin(diff)
            #if T[idx1]
        if T[idx1]<Ts:
            idx2 = idx1+1
        if T[idx1]>Ts:
            idx2 = idx1 -1




        if T[idx1]> (Ts+30) or T[idx2]>30: #only choosing GRBs with points that fall within + 30 of Ts
            pass
        else:

            if len(ene_flux)<=2:
                pass
            else:
                plt.figure(figsize = (12,8 ))
                plt.scatter(T, lum, label=name )
                plt.yscale('log')
                plt.xscale('log')
                plt.xlabel('log(T-T0/(1+z)) [s]',fontsize = 14) 
                plt.ylabel('log(L) [erg/s]',fontsize = 14) 
                plt.title(f'Corrected Light Curves of LAT GRB{name}',fontsize = 14) 
                plt.tick_params(axis='both', which='major', labelsize=14)
                mpl.pyplot.errorbar(T, lum, yerr=lum_err, xerr=T_err, linestyle='',elinewidth=0.5)
                plt.axvline(Ts, color = 'red')#, label=f'{Ts}')   
                plt.legend()
                
                if idx1>Ts:
                    T_new = T[idx2:]
                    lum_new = lum[idx2:]

                else:
                    T_new = T[idx1:]
                    lum_new = lum[idx1:]



                p1= [1e55,1]
                p2= [1e57,10]
                #p3 = [1e55,1,1000]

                #power_law2_vec = np.vectorize(power_law2)

                pars1, cov1 = optimize.curve_fit(power_law1, xdata=T, ydata=lum, p0=p1, method='lm')
                pars2, cov2 = optimize.curve_fit(power_law1, xdata=T_new, ydata=lum_new, p0=p1, method='lm',maxfev = 1000,bounds=(-np.inf, np.inf) )
                #pars3, cov3 = optimize.curve_fit(power_law2, xdata=T_new, ydata=lum_new, p0=p3, method='lm',maxfev = 1000,bounds=(-np.inf, np.inf) )


                # Get the standard deviations of the parameters (square roots of the # diagonal of the covariance)
                stdevs1 = np.sqrt(np.diag(cov1))
                stdevs2 = np.sqrt(np.diag(cov2))
                
                
                # Calculate the residuals
                res1 = []
                res2 = []
                for i in range(len(lum)):
                    res1.append(abs(lum[i] - power_law1(T[i], *pars1)))
                    res2.append(abs(lum[i] - power_law1(T[i], *pars2)))


                #parameters1.append([pars1[0], pars1[1],stdevs1[0],stdevs1[1]])
                parameters1.append([name, pars1[0], stdevs1[0], pars1[1], stdevs1[1], mean_B, mean_B_err ])
                parameters2.append([name, pars2[0], stdevs2[0], pars2[1], stdevs2[1], mean_B, mean_B_err ])

                print(f'Weighted spectral index average b={np.round(mean_B,4)} ± {np.round(mean_B_err,4)}')
                if pars1[1]==pars2[1]:
                    plt.plot(T, power_law1(T, pars1[0],pars1[1]), linestyle='--', linewidth=2, label=f'Fit a={np.round(pars1[1],4)}', color='orange')
                    print(f'Whole and Ts fit returned the same temporal index a={np.round(pars1[1],4)} ± {np.round(stdevs1[1],4)}')
                    plt.show()
                    
                else:    
                    plt.plot(T, power_law1(T, pars1[0],pars1[1]), linestyle='--', linewidth=2, label=f'Whole Fit a={np.round(pars1[1],4)}',color='orange')
                    plt.plot(T, power_law1(T, pars2[0],pars2[1]), linestyle='--', linewidth=2, label=f'Ts Fit a={np.round(pars2[1],4)}',color='green')
                    print(f'Whole fit temporal index a={np.round(pars1[1],4)} ± {np.round(stdevs1[1],4)}')
                    print(f'Ts fit temporal index a={np.round(pars2[1],4)} ± {np.round(stdevs2[1],4)}')
                    plt.show()

                
                #plt.figure(figsize = (12,8 ))
                #plt.subplot(1,2,1)
                #plt.scatter(T,res1, color='orange')
                #print((res2))
                #plt.scatter(T,res2, color='green')
                #plt.ylim(0,0.2e54)
                #plt.yscale('log')
                #plt.xscale('log')
                
                #plt.subplot(1,2,2)
                #plt.scatter(T,res2)
                #plt.yscale('log')
                #plt.xscale('log')

                #print(f'Ts fit broken params are {pars3}')
                #plt.plot(T, power_law2(T, pars3[0],pars3[1], pars3[2]), linestyle='--', linewidth=2, label=f'Ts Fit a={np.round(pars3[1],4)}')

                

         


#saving to excel sheet (1-whole fit, 2-Ts fit)
#format [GRB Name, N, N error, a, a error, b, b error]
workbook1 = xlsxwriter.Workbook('GRB_params1.xlsx')
worksheet1 = workbook1.add_worksheet()
workbook2 = xlsxwriter.Workbook('GRB_params2.xlsx')
worksheet2 = workbook2.add_worksheet()
row = 0

#for col, data in enumerate(parameters1):
 #   worksheet1.write_row(col,row,data)
#workbook1.close()

#for col, data in enumerate(parameters2):
#    worksheet2.write_row(col,row,data)
#workbook2.close()

#for i in parameters1:
    #print(i[1])
  
#luminosity(table)        


# In[ ]:





# In[ ]:





# In[110]:


plt.figure(figsize = (12,8 ))
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
            
            B=-1*indec 
            B_err = -1*index_err
            plt.xlabel('Spectral Index Error',fontsize = 14) 
            plt.ylabel('Spectral Index',fontsize = 14) 
            plt.title('Spectral Index vs Error on Spectral Index',fontsize = 14) 
            plt.tick_params(axis='both', which='major', labelsize=14)
            #mpl.pyplot.errorbar(T, lum, yerr=lum_err, xerr=T_err, linestyle='',elinewidth=0.5)
            #plt.yscale('log')
            #plt.xscale('log')
            plt.scatter(B_err, B)


# In[115]:


plt.figure(figsize = (12,8 ))
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
            
            T_err = ((t_end-t_start)/2)
            T = (t_start+T_err)/(1+z)
            T_err =  T_err/(1+z)
            
            B=-1*indec 
            B_err = -1*index_err
            
            plt.scatter(T, abs(B))
            plt.xlabel('log(T-T0/(1+z)) [s]',fontsize = 14) 
            plt.ylabel('Spectral Index',fontsize = 14) 
            plt.title('Spectral Index vs Time',fontsize = 14) 
            plt.tick_params(axis='both', which='major', labelsize=14)
            #mpl.pyplot.errorbar(T, B, xerr=T_err, yerr=B_err, linestyle='',elinewidth=0.5)

            plt.yscale('log')
            plt.xscale('log')


# In[152]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




