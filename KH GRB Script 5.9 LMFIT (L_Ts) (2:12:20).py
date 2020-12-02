#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
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

import lmfit
from lmfit import Parameters, fit_report, minimize
from lmfit.models import ExponentialModel, GaussianModel, PowerLawModel, ExpressionModel, LinearModel

mpl.rc('font', family='serif')


grb = fits.open('gll_2flgc.fits')

table=grb[1].data
grb.close


# In[2]:


def power_law1(x, N, a): 
    '''simple power law 
    L = Nx^-a'''
    return N*np.power(x, a)


# In[3]:


def power_law2(x, N, a, b):
    
    '''broken power law (a piecewise function)
    L = Nx^-a1 for x<xth
    L = Nxth^-(a1-a2)x^-a2 for x>xth'''
    xth = np.where(x==x.max())
    
    #for k in x:
    if k <=x[xth]:
        return N*np.power(x, -a)
    else:
        return N*np.power(xth, -(a-b))*np.power(x,-b)


# In[4]:


def power_law3(x,N,a):
    '''simple power law in log space
    log(L) = log(N) - alogx'''
    return (a)*np.log10(x) + np.log10(N) 


# In[9]:


#def luminosity(table):
Kc = []      #empty array for k corrections
Kc_std = []

parameters1 = []
parameters2 = []
parameters3 = []
parameters4 = []

lum_ts = []
n=0
redshift = []
m=0
Ts = 10
j = 0
diff_mean = []
diff_std = []



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


        g = 1.6e-6               #extra factor missing in the LC_ENE_FLUX values in the .fits
        B=-1*indec               #calculating spectral index from measured photon indicies
        B_err = -1*index_err
        
        #calculating weighted mean of spectral index and error
        B_sum = []
        B_err_sum = []

        if len(B)==1:
            mean_B = B
            mean_B_err = abs(B_err)
        else:
            mean_B = np.average(B, weights=abs(B_err))
            mean_B_err = np.mean(abs(B_err))
           
                
            
        #calculating k correction
        k_correction = (1+z)**((mean_B-1)-1)

        #calculating corrected luminosities and their errors in normal and log space
        lum = g*ene_flux*4*(np.pi)*(dl**2)*(k_correction)
        lum_err = g*ene_flux_err*4*(np.pi)*(dl**2)*(k_correction)
        
        #calculating the asymmetric errors in log space
        lum_err_plus = np.log10(lum + lum_err) - np.log10(lum)        
        lum_err_minus= np.log10(lum) - np.log10(lum - lum_err)
        
        for j in lum:
            diff.append(abs((j/(4*(np.pi)*(dl**2)) - (g*ene_flux))))
            
        diff_mean.append(np.mean(diff))
        diff_std.append(np.std(diff))

        T_err = ((t_end-t_start)/2)
        T = (t_start+T_err)/(1+z)
        T_err =  T_err/(1+z)
        
        T_err_plus = np.log10(T + T_err) - np.log10(T)
        T_err_minus = np.log10(T) - np.log10(T - T_err)


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
                
                #plotting lightcurve from data
                plt.figure(figsize = (12,8 ))
                plt.scatter(T, lum, label=name )
                plt.yscale('log')
                plt.xscale('log')
                plt.xlabel('T-T0/(1+z) [s]',fontsize = 14) 
                plt.ylabel('L [erg/s]',fontsize = 14) 
                plt.title(f'Corrected Light Curves of LAT GRB{name}',fontsize = 14) 
                plt.tick_params(axis='both', which='major', labelsize=14)
                mpl.pyplot.errorbar(T, lum, yerr=lum_err, xerr=T_err, linestyle='',elinewidth=0.5)
                plt.axvline(Ts, color = 'red', linewidth=0.5)#, label=f'{Ts}')   
                plt.legend()

                if idx1>Ts:
                    T_new = T[idx2:]
                    lum_new = lum[idx2:]
                    lum_err_p_new = lum_err_plus[idx2:]
                    lum_err_m_new = lum_err_minus[idx2:]

                else:
                    T_new = T[idx1:]
                    lum_new = lum[idx1:]
                    lum_err_p_new = lum_err_plus[idx1:]
                    lum_err_m_new = lum_err_minus[idx1:]
                    
                #calculating weights
                weight_lum = 1/lum_err
                weight_lum_new = 1/lum_err_p_new

                #calculating luminosity at Ts
                if idx1 > idx2:
                    lum1 = lum[idx1]       #finding luminosity data points either side of Ts
                    lum2 = lum[idx2]
                    T1 = T[idx1]
                    T2 = T[idx2]
                    w1 = 1/lum_err[idx1]
                    w2 = 1/lum_err[idx2]
                else:
                    lum1 = lum[idx2]
                    lum2 = lum[idx1]
                    T1 = T[idx2]
                    T2 = T[idx1]
                    w1 = 1/lum_err[idx2]
                    w2 = 1/lum_err[idx1]

                log_T = np.log10(T)
                log_lum = np.log10(lum)
                log_Ts = np.log10(Ts)
                
                
                #calculating best fit parameters and covariances for the data lmfit
                #making powerlaw model for linear fits
                model1 = PowerLawModel(prefix='pow_')
                #making powerlaw model for log fits
                model2 = LinearModel(independent_vars=['x'])
                
                # make parameters with starting values:
                par1 = model1.make_params(pow_amplitude=1e55, pow_exponent=-1.0)
                par2 = model2.make_params(m=1,c=55)

                
                np.nan_to_num(weight_lum_new,copy=False)

                #running both sets of fits
                result1 = model1.fit(lum, par1, x=T, weights=weight_lum) #linear whole with weights
                result2 = model1.fit(lum_new, par1, x=T_new, weights=weight_lum_new) #linear ts with weight
                #result0 = model0.fit(lum, par1, x=T) #linear whole without weights
                result3 = model2.fit(np.log10(lum),par2,x=np.log10(T),weights=weight_lum) #log whole with weights
                result4 = model2.fit(np.log10(lum_new),par2,x=np.log10(T_new), weights=weight_lum_new) #log ts with weights
                result5 = model1.fit([lum1,lum2], par1, x=[T1,T2], weights=[w1,w1])
               
                a1, N1 = result1.best_values['pow_exponent'],result1.best_values['pow_amplitude'] #linear whole
                a2, N2 = result2.best_values['pow_exponent'],result2.best_values['pow_amplitude'] #linear ts
                a3, N3 = result3.best_values['slope'], result3.best_values['intercept'] #log whole
                a4, N4 = result4.best_values['slope'], result4.best_values['intercept'] #log ts
                a5, N5 = result5.best_values['pow_exponent'],result5.best_values['pow_amplitude'] #best fit for L(Ts)
                
                cova1, covN1, cova2, covN2  = [], [], [], []
                cova3, covN3, cova4, covN4, cova5, covN5 = [], [], [], [], [], []
                
                for i in result1.covar[1]:
                    cova1.append(i)
                    
                for i in result1.covar[0]:
                    covN1.append(i)

                for i in result2.covar[1]:
                    cova2.append(i)
                    
                for i in result2.covar[0]:
                    covN2.append(i)
                    
                for i in result3.covar[1]:
                    cova3.append(i)
                    
                for i in result3.covar[0]:
                    covN3.append(i)

                for i in result4.covar[1]:
                    cova4.append(i)
                    
                for i in result4.covar[0]:
                    covN4.append(i)
                    
                for i in result5.covar[1]:
                    cova5.append(i)
                    
                for i in result5.covar[0]:
                    covN5.append(i)
                    
                #calculating std
                stda1, stdN1 = np.sqrt(cova1[1]), np.sqrt(covN1[0])
                stda2, stdN2 = np.sqrt(cova2[1]), np.sqrt(covN2[0])
                stda3, stdN3 = np.sqrt(cova3[1]), np.sqrt(covN3[0])
                stda4, stdN4 = np.sqrt(cova4[1]), np.sqrt(covN4[0])
                stda5, stdN5 = np.sqrt(cova5[1]), np.sqrt(covN5[0])


                #if (result1.redchi <1.9):
                #    ci = lmfit.conf_interval(result1, result1, sigmas=[0.68])
                #    print( "Confidence intervals data 1")
                #    lmfit.report_ci(ci)

                #lmfit fit
                lum_fit1, log_lum_fit1 = power_law1(T, N1, a1), np.log10(power_law1(T, N1, a1))
                lum_fit2, log_lum_fit2 = power_law1(T, N2, a2), np.log10(power_law1(T, N2, a2))
                lum_fit3, log_lum_fit3 = power_law1(T, 10**N3, a3), np.log10(power_law1(T, 10**N3, a3))
                lum_fit4, log_lum_fit4 = power_law1(T, 10**N4, a4), np.log10(power_law1(T, 10**N4, a4))
                lum_fit5 = power_law1(T, N5, a5)
                
                #highlighting each point either side of ts
                plt.scatter([T1,T2],[lum1,lum2], color='red', marker='s')
                
                #creating two arrays with the parameters and their std 
                parameters1.append([name, N1, stdN1, a1, stda1, mean_B, mean_B_err ]) #linear whole
                parameters2.append([name, N2, stdN2, a2, stda2, mean_B, mean_B_err ]) #linear ts
                parameters3.append([name, N3, stdN3, a3, stda3, mean_B, mean_B_err ]) #log whole
                parameters4.append([name, N4, stdN4, a4, stda4, mean_B, mean_B_err ]) #log ts
                #parameters5.append([name, N5, stdN5, a5, stda5, mean_B, mean_B_err ])
                
                print(f'Weighted spectral index average b={np.round(mean_B,4)} ± {np.round(mean_B_err,4)}')
                print()

                #plotting best fits over the orignal data
                #plt.plot(T, lum_fit1, linestyle='--', linewidth=2, 
                #         label=f'PL Whole Fit a={np.round(a1,3)} ± {np.round(stda1,3)}',color='orange')
                plt.plot(T, lum_fit2, linestyle='--', linewidth=2, 
                         label=f'PL Ts_10 a={np.around(a2,3)} ± {np.around(stda2,3)}',color='green')
                #plt.plot(T, lum_fit3, linestyle='--', linewidth=2, 
                #         label=f'PL-SL Whole Fit a={np.round(a3,3)} ± {np.round(stda3,3)}', color='purple')
                #plt.plot(T, lum_fit4, linestyle='--', linewidth=2, 
                #         label=f'PL-SL Ts Fit a={np.round(a4,3)} ± {np.round(stda4,3)}', color='black')
                #plt.plot(T, (power_law1(T, N5,a5)), linestyle='--', linewidth=2, 
                #         label=f'L(Ts) Fit a={np.round(a5,3)} ± {np.round(stda5,3)}', color='blue')



                plt.tick_params(axis='both', which='major', labelsize=18)
                plt.legend()
                
                print(f'SPL Whole fit temporal index a={np.round(a1,3)} ± {np.round(stda1,3)}')
                print(f'SPL Ts fit temporal index a={np.round(a2,3)} ± {np.round(stda2,3)}')
                print()
                print(f'SPL-SL Whole fit temporal index a={np.round(a3,3)} ± {np.round(stda3,3)}')
                print(f'SPL-SL Ts fit temporal index a={np.round(a4,3)} ± {np.round(stda4,3)}')
                print()
                print(f'Luminsoity at {Ts}s is {power_law1(Ts, N5,a5)}')
                plt.axhline(power_law1(Ts, N5,a5),linewidth=2)
                plt.show()
                #plt.figure(figsize = (12,8 ))
                #plt.scatter(log_T, log_lum, label=name )
                #plt.xlabel('log(T-T0/(1+z)) [s]',fontsize = 14) 
                #plt.ylabel('log(L) [erg/s]',fontsize = 14) 
                #plt.title(f'Corrected Light Curves Log Space of LAT GRB{name}',fontsize = 14) 
                #plt.tick_params(axis='both', which='major', labelsize=18)
                #plt.errorbar(log_T, log_lum, yerr=[lum_err_minus,lum_err_plus],
                #             xerr=[T_err_minus,T_err_plus] , linestyle='',elinewidth=0.5)
                #plt.axvline(np.log10(Ts), color = 'red', linewidth=0.5)#, label=f'{Ts}')   
                #plt.plot(log_T, log_lum_fit1, linestyle='--', linewidth=2, 
                #         label=f'SPL Whole Fit a={np.round(a1,2)} ± {np.round(stda1,2)}', color='orange')
                #plt.plot(log_T, log_lum_fit2, linestyle='--', linewidth=2, 
                #         label=f'SPL Ts Fit a={np.round(a2,2)} ± {np.round(stda2,2)}', color='green')
                #plt.plot(log_T, log_lum_fit3, linestyle='--', linewidth=2, 
                #         label=f'SPL-SL Whole Fit a={np.round(a3,2)} ± {np.round(stda3,2)}', color='purple')
                #plt.plot(log_T, log_lum_fit4, linestyle='--', linewidth=2, 
                #         label=f'SPL-SL Ts Fit a={np.round(a4,2)} ± {np.round(stda4,2)}', color='black')
                #plt.plot(np.log(T), np.log(power_law1(T, pars5[0],pars5[1])), linestyle='--', linewidth=2, 
                #         label=f'SPL-SL Ts Fit a={np.round(pars5[1],4)} ± {np.round(stdevs5[1],4)}', color='blue')

                #plt.scatter(np.log10([T1,T2]),np.log10([lum1,lum2]), color='r', marker='s')

                #plt.axhline(np.log10(power_law1(Ts, N5,a5)),linewidth=2)
                #plt.legend()
                #plt.show()

                #np.savetxt(f'LAT_GRB{name}_data',np.column_stack((T,T_err,lum, lum_err)),delimiter=' ')

                #np.savetxt(f'LAT_GRB{name}_log_data',np.column_stack((log_T,T_err_plus,T_err_minus,
                                                                    #  log_lum, lum_err_plus, lum_err_minus)),delimiter=' ')

                    
                    




# In[49]:


np.sqrt(1.0330397e106)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[9]:





# In[ ]:





# In[26]:





# In[27]:





# 

# In[31]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




