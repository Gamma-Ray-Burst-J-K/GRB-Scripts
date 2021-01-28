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
import astropy
from astropy.modeling.powerlaws import BrokenPowerLaw1D

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


# In[5]:


def retrieve_grb_data(table, grb, z):
    '''
    table - fits files containing all data
    grb - name of grb
    
    want to retrieve all the data associated with a given grb'''
    for i in table:
        name = i['GCNNAME']
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
        #for n in i['GCNNAME']:
            #print(n)
        #n = np.where(i['GCNNAME'] == grb)
        #print(name)
        if name == grb:
            print(f'Match for GRB {name} found.')
            dl = i['LUMINOSITY_DISTANCE']
            return [ name,z, dl, t_end, ene_flux, ene_flux_err, fluence, flux, flux_err, indec, index_err, median, t_start, ts]
            
            
            

            
    
    
    
    


# In[6]:


#print(retrieve_grb_data(table,'080916C',5 ))
#new grb 140102A
#retrieve_grb_data(table,'140102A',5 )


# In[ ]:





# In[8]:


#def luminosity(table):


n=0
redshift = []
m=0
Ts = 20
j = 0
diff_mean = []
diff_std = []

#arrays for collecting terms
alpha = []
alpha_err = []

norm_alpha_tp, norm_lum_ts_tp = [],[]
norm_alpha_tp_err = []

norm_alpha_U, norm_lum_ts_U = [], []
norm_alpha_U_err = []

norm_alpha_per, norm_lum_ts_per = [], []
norm_alpha_per_err = []

Kc = []      #empty array for k corrections
Kc_std = []

all_parameters = []
all_parameterslog = []

lum_params = []



lum_ts = []
lum_ts_err = []
grb_name = []

for i in table:
    if i['LUMINOSITY_DISTANCE'] >0:
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


        if T[idx1]> (Ts+30) or T[idx2]>(Ts+30): #only choosing GRBs with points that fall within ± 30 of Ts
            pass
        else:

            if len(ene_flux)<=2:
                pass
            else:

                

                #calculating weights
                weight_lum = 1/lum_err
                
                if T[idx1]>Ts:                 #makign sure the weights are the correct way round
                    T_new = T[idx2:]
                    lum_new = lum[idx2:]
                    lum_err_p_new = lum_err_plus[idx2:]
                    lum_err_m_new = lum_err_minus[idx2:]
                    weight_ts = weight_lum[idx2:]
                    weight_lum_new = 1/lum_err_p_new
                else:
                    T_new = T[idx1:]
                    lum_new = lum[idx1:]
                    lum_err_p_new = lum_err_plus[idx1:]
                    lum_err_m_new = lum_err_minus[idx1:]
                    weight_ts = weight_lum[idx1:]
                    weight_lum_new = 1/lum_err_p_new
                
                #fitting from 10s-40s
                U = 10
                diff_minus1 = abs(T-(Ts-U))                   #lower T threshold
                diff_plus1 = abs(T-(Ts+(2*U)))                #upper T threshold
                idx_minus1 = np.argmin(diff_minus1)           #finding lower threshold index
                idx_plus1 = np.argmin(diff_plus1) +1          #finding upper threshold index
                T_min_plus1 = T[idx_minus1:idx_plus1]         #making new T array
                lum_min_plus1 = lum[idx_minus1:idx_plus1]     #making new L array
                weight_mp1 = weight_lum[idx_minus1:idx_plus1] #making new weights array
                
                #fitting from ±% Ts
                p = 0.3                   #percentage value
                d_minus2 = 10**((1-p)*np.log10(Ts))
                d_plus2 = 10**((1+p)*np.log10(Ts))
                diff_minus2 = abs(T-d_minus2)
                diff_plus2 = abs(T-d_plus2)
                idx_minus2 = np.argmin(diff_minus2)
                idx_plus2 = np.argmin(diff_plus2)+1
                T_min_plus2 = T[idx_minus2:idx_plus2]
                lum_min_plus2 = lum[idx_minus2:idx_plus2]
                weight_mp2 = weight_lum[idx_minus2:idx_plus2]

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
                par1 = model1.make_params(pow_amplitude=1e55, pow_exponent=-1.0) #linear powerlaw
                par2 = model2.make_params(m=1,c=55) #log 
                par3 = model1.make_params(pow_amplitude=1e5, pow_exponent=-1.0)
                par4 = model1.make_params(pow_amplitude=1e52, pow_exponent=-1.0)
                
                np.nan_to_num(weight_lum_new,copy=False)
                
                #running sets of fits
                result1 = model1.fit(lum, par1, x=T, weights=weight_lum) #linear whole with weights
                result2 = model1.fit(lum_new, par1, x=T_new, weights=weight_ts) #linear ts with weight
                result3 = model2.fit(np.log10(lum),par2,x=np.log10(T),weights=weight_lum) #log whole with weights
                result4 = model2.fit(np.log10(lum_new),par2,x=np.log10(T_new), weights=weight_lum_new) #log ts with weights
                
                result5 = model1.fit([lum1,lum2], par1, x=[T1,T2], weights=[w1,w1]) #two-point lum
                
                #running two-point fits
                resultn5 = model1.fit([lum1/1e50,lum2/1e50], par3, x=[T1,T2], weights=[w1,w1]) #norm two-point
                
                a1, N1 = result1.best_values['pow_exponent'],result1.best_values['pow_amplitude'] #linear whole params
                a2, N2 = result2.best_values['pow_exponent'],result2.best_values['pow_amplitude'] #linear ts params
                a3, N3 = result3.best_values['slope'], result3.best_values['intercept'] #log whole params
                a4, N4 = result4.best_values['slope'], result4.best_values['intercept'] #log ts params
                a5, N5 = result5.best_values['pow_exponent'],result5.best_values['pow_amplitude'] #two-point lum params
                
                
                #two-point luminosity parameters
                norma5, normN5 = resultn5.best_values['pow_exponent'],resultn5.best_values['pow_amplitude'] #norm two-point params

                #running ± U Ts fits and  ± 50 percent & parameters
                resultn3 = model1.fit(lum_min_plus1/1e50, par3, x=T_min_plus1, weights=weight_mp1) #±U Ts luminosity 
                resultn4 = model1.fit(lum_min_plus2/1e50, par3, x=T_min_plus2, weights=weight_mp2) #±30 percent lum 

                norma3, normN3 = resultn3.best_values['pow_exponent'],resultn3.best_values['pow_amplitude'] #±U Ts luminosity params
                norma4, normN4 = resultn4.best_values['pow_exponent'],resultn4.best_values['pow_amplitude'] #±30 percent lum params
                
                
                
                #finding the errors on the best fit
                ci1 = lmfit.conf_interval(result1, result1, sigmas=[0.68])
                ci2 = lmfit.conf_interval(result2, result2, sigmas=[0.68])
                ci3 = lmfit.conf_interval(result3, result3, sigmas=[0.68])
                ci4 = lmfit.conf_interval(result4, result4, sigmas=[0.68])
                
                normci5 = lmfit.conf_interval(resultn5, resultn5, sigmas=[0.68])
                normci3 = lmfit.conf_interval(resultn3, resultn3, sigmas=[0.68])
                normci4 = lmfit.conf_interval(resultn4, resultn4, sigmas=[0.68])
                print(normci4)
                
                #print(lmfit.fit_report(result5.params))
                #ci5 = lmfit.conf_interval(result5, result5, sigmas=[0.68], maxiter=500)
                
                #extracting sigma errors
                con1, con2, con3, con4, con5 = [],[],[],[],[]
                normcon5, normcon3, normcon4 = [],[],[]
                
                for key, value in ci1.items():
                    con1.append(value)
                    
                for key, value in ci2.items():
                    con2.append(value)
                
                for key, value in ci3.items():
                    con3.append(value)
                
                for key, value in ci4.items():
                    con4.append(value)

                for key, value in normci5.items():
                    normcon5.append(value)
                
                for key, value in normci3.items():
                    normcon3.append(value)
                
                for key, value in normci4.items():
                    normcon4.append(value)
                
                 
                #errors are an array with [sigma minus, sigma plus]
                #linear whole with weights errors
                a1_sig = [con1[1][0][1], con1[1][2][1]]
                N1_sig = [con1[0][0][1], con1[0][2][1]]
                
                #linear ts with weights errors
                a2_sig = [con2[1][0][1], con2[1][2][1]]
                N2_sig = [con2[0][0][1], con2[0][2][1]]
                
                #log whole with weights errors
                a3_sig = [abs(a3 - con3[1][0][1]), con3[1][2][1]]
                N3_sig = [abs(N3 - con3[0][0][1]), con3[0][2][1]]
                
                ##log ts with weights errors
                a4_sig = [con4[1][0][1], con4[1][2][1]]
                N4_sig = [con4[0][0][1], con4[0][2][1]]
                
                #normalised 2-point lum errors
                norma5_sig = [normcon5[1][0][1], normcon5[1][2][1]]
                normN5_sig = [normcon5[0][0][1], normcon5[0][2][1]]
                
                #normalised ±U Ts lum
                norma3_sig = [normcon3[1][0][1], normcon3[1][2][1]]
                normN3_sig = [normcon3[0][0][1], normcon3[0][2][1]]
                
                #normalised ±30% Ts lum
                norma4_sig = [normcon4[1][0][1], normcon4[1][2][1]]
                normN4_sig = [normcon4[0][0][1], normcon4[0][2][1]]
                
                #all_parameters.append([name, mean_B, mean_B_err, a1,abs(a1_sig[0]-a1),abs(a1_sig[1]-a1),N1,abs(N1_sig[0]-N1), abs(N1_sig[1]-N1),
                 #                     a2,abs(np.log10(a2_sig[0])-a2), abs(a2_sig[1]-a2),N2, abs(N2_sig[0]-N2), abs(N2_sig[1]-N2)])
                
                #all_parameterslog.append([name, mean_B, mean_B_err, a3,abs(np.log10(a3_sig[0])-a3), abs(np.log10(a3_sig[1])-a3),N3, abs(N3_sig[0]-N3), abs(N3_sig[1]-N3),
                #                      a4,abs(np.log10(a4_sig[0])-a4), abs(np.log10(a4_sig[1])-a4),N4, abs(N4_sig[0]-N4), abs(N4_sig[1]-N4)])
   
                
                print('normalised luminosity from two points')
                #print(lmfit.fit_report(resultn5.params))
                pln5 = power_law1(Ts, normN5,norma5)*1e50
                print(f'luminosity = {np.round(power_law1(Ts, normN5,norma5),5)*1e50}')# + {np.round(power_law1(Ts, normN5_sig_plus, a5),5)} - {np.round(power_law1(Ts, normN5_sig_minus, a5),5)} ')
                print(f'-{abs(power_law1(Ts,normN5_sig[0],norma5)*1e50 - pln5)}')
                print(f'+{abs(pln5 - power_law1(Ts,normN5_sig[1],norma5)*1e50)}')
                print()
                
                print('normalised luminosity Ts±10')
                #print(lmfit.fit_report(resultn3.params))
                pln3 = power_law1(Ts, normN3,norma3)*1e50
                print(f'luminosity = {power_law1(Ts, normN3, norma3)*1e50}')
                print(f'-{abs(power_law1(Ts,normN3_sig[0],norma3)*1e50 - pln3)}')
                print(f'+{abs(pln3 - power_law1(Ts,normN3_sig[1],norma3)*1e50)}')
                print()
                
                print('normalised luminosity Ts±30%')
                #print(lmfit.fit_report(resultn4.params))
                pln4 = power_law1(Ts, normN4,norma4)*1e50
                print(f'luminosity = {power_law1(Ts, normN4, norma4)*1e50}')
                print(f'-{abs(power_law1(Ts,normN4_sig[0], norma4)*1e50 - pln4)}')
                print(f'+{abs(pln4 - power_law1(Ts,normN4_sig[1], norma4)*1e50)}')
                print()
                
                #lum_params.append([name,power_law1(Ts, normN5,norma5)*1e50, abs(power_law1(Ts,normN5_sig[0],norma5)*1e50 - pln5), abs(pln5 - power_law1(Ts,normN5_sig[1],norma5)*1e50),
                #                   norma5, norma5_sig[0], norma5_sig[1], 
                #                  power_law1(Ts, normN3, norma3)*1e50,abs(power_law1(Ts,normN3_sig[0],norma3)*1e50 - pln3),
                #                   abs(pln3 - power_law1(Ts,normN3_sig[1],norma3)*1e50),norma3, norma3_sig[0], norma3_sig[1],
                #                  power_law1(Ts, normN4, norma4)*1e50,abs(power_law1(Ts,normN4_sig[0], norma4)*1e50 - pln4),
                 #                  abs(pln4 - power_law1(Ts,normN4_sig[1], norma4)*1e50),norma4, norma4_sig[0], norma4_sig[1],])

                #lmfit fit
                lum_fit1, log_lum_fit1 = power_law1(T, N1, a1), np.log10(power_law1(T, N1, a1))
                lum_fit2, log_lum_fit2 = power_law1(T, N2, a2), np.log10(power_law1(T, N2, a2))
                lum_fit3, log_lum_fit3 = power_law1(T, 10**N3, a3), np.log10(power_law1(T, 10**N3, a3))
                lum_fit4, log_lum_fit4 = power_law1(T, 10**N4, a4), np.log10(power_law1(T, 10**N4, a4))
                
                lum_fit5 = power_law1(T, N5, a5)
                
                alpha.append(a5)
                lum_ts.append(power_law1(Ts, N5,a5))

                
                #alpha_err.append(stda2)
                norm_lum_ts_per.append(pln4)
                norm_alpha_per.append(norma4)
                
                norm_lum_ts_tp.append(pln5)
                norm_alpha_tp.append(norma5)
                
                norm_lum_ts_U.append(pln3)
                norm_alpha_U.append(norma3)
                
                grb_name.append(name)
                
                
                #plotting different luminosity methods
                fig, axs = plt.subplots(3,figsize=(12,20), sharex=False, sharey=True)
                
                #highlighting each point either side of ts
                axs[0].scatter(T, lum, label=name )
                axs[0].set_yscale('log')
                axs[0].set_xscale('log')
                axs[0].set_xlabel(xlabel='T-T0/(1+z) [s]',fontsize = 14) 
                axs[0].set_ylabel(ylabel='L [erg/s]',fontsize = 14) 
                axs[0].set_title(f'GRB{name} for Ts={Ts}s Two-Point Luminosity Fit',fontsize = 14) 
                axs[0].tick_params(axis='both', which='major', labelsize=14)
                axs[0].errorbar(T, lum, yerr=lum_err, xerr=T_err, linestyle='',elinewidth=0.5)
                axs[0].axvline(Ts, color = 'red', linewidth=0.5)#, label=f'{Ts}')   
                axs[0].scatter([T1,T2],[lum1,lum2], color='red', marker='s')
                axs[0].axhline(power_law1(Ts, N5,a5),color='red', linewidth=1.5, label='2-Point Lum')
                axs[0].legend()
                
                #highlighting points for ±U Ts luminosity
                axs[1].scatter(T, lum, label=name )
                axs[1].set_yscale('log')
                axs[1].set_xscale('log')
                axs[1].set_xlabel(xlabel='T-T0/(1+z) [s]',fontsize = 14) 
                axs[1].set_ylabel(ylabel='L [erg/s]',fontsize = 14)
                axs[1].set_title(f'GRB{name} for Ts={Ts}s ±{U} Luminosity fit',fontsize = 14) 
                axs[1].tick_params(axis='both', which='major', labelsize=14)
                axs[1].errorbar(T, lum, yerr=lum_err, xerr=T_err, linestyle='',elinewidth=0.5)
                axs[1].axvline(Ts, color = 'red', linewidth=0.5)#, label=f'{Ts}')   
                axs[1].scatter(T_min_plus1, lum_min_plus1, color='black')
                axs[1].axhline(power_law1(Ts, normN3, norma3)*1e50, color='black', linewidth=1.5, label=f'±{U} Ts Lum')
                axs[1].legend()
                
                #highlighting points for ±50% Ts luminsoity
                axs[2].scatter(T, lum, label=name )
                axs[2].set_yscale('log')
                axs[2].set_xscale('log')
                axs[2].set_xlabel(xlabel='T-T0/(1+z) [s]',fontsize = 14) 
                axs[2].set_ylabel(ylabel='L [erg/s]',fontsize = 14) 
                axs[2].set_title(f'GRB{name} for Ts={Ts}s ±{100*p}% Luminosity fit',fontsize = 14) 
                plt.tick_params(axis='both', which='major', labelsize=14)
                axs[2].errorbar(T, lum, yerr=lum_err, xerr=T_err, linestyle='',elinewidth=0.5)
                plt.axvline(Ts, color = 'red', linewidth=0.5)#, label=f'{Ts}')   
                axs[2].legend()
                axs[2].scatter(T_min_plus2, lum_min_plus2, color='orange')
                axs[2].axhline(power_law1(Ts, normN4, norma4)*1e50,color='orange', linewidth=1.5, label='±50% Ts Lum')

                plt.show()
                
                
                #np.savetxt(f'LAT_GRB{name}_data',np.column_stack((lum, lum_err, T,T_err)),delimiter=' ')

                #np.savetxt(f'LAT_GRB{name}_log_data',np.column_stack((log_T,T_err_plus,T_err_minus,
                                                                    #  log_lum, lum_err_plus, lum_err_minus)),delimiter=' ')

                    
                  

                


# In[68]:


#np.savetxt('Multi-Fit-Params-Log.csv',all_parameterslog,fmt='%s', delimiter=',')
#np.savetxt('Multi-Fit-Params-Lineaar.csv',all_parameters,fmt='%s', delimiter=',')
#np.savetxt('Luminosities.csv',lum_params,fmt='%s', delimiter=',')


# In[16]:


#Spearman Ranks
#mask1 = np.where()
srank_tp = stats.spearmanr(norm_lum_ts_tp, norm_alpha_tp)

srank_U = stats.spearmanr(norm_lum_ts_U, norm_alpha_U)

srank_per = stats.spearmanr(norm_lum_ts_per, norm_alpha_per)


print(f'two point fit {srank_tp}') 
print()
print(f'10-40s fit {srank_U}')
print()
print(f'±30% fit {srank_per}')


# In[ ]:





# In[107]:


#Un-normalised fits
plt.figure(figsize = (12,8 ))

        
for i in range(len(grb_name)):
    plt.scatter(alpha[i], lum_ts[i],label=grb_name[i])
    

#print(np.log10(lum_ts_err)
plt.xlabel('Alpha',fontsize = 14) 
plt.ylabel('L_Ts=10 [erg/s]',fontsize = 14) 
plt.title(f'Luminosity at Ts={Ts} vs Decay Index From Un-normalised Fits',fontsize = 14) 
plt.tick_params(axis='both', which='major', labelsize=18)
#plt.errorbar(alpha, lum_ts, yerr=lum_ts_err, xerr=alpha_err, linestyle='',elinewidth=0.5)
plt.legend()

srank = stats.spearmanr(alpha, lum_ts)
print(srank)


# In[114]:


#Normalised Two-Point fits
plt.figure(figsize = (12,8 ))

        
for i in range(len(grb_name)):
    plt.scatter(norm_alpha_tp[i], norm_lum_ts_tp[i],label=grb_name[i])
    

#print(np.log10(lum_ts_err)
plt.xlabel('Alpha',fontsize = 14) 
plt.ylabel('L_Ts=10 [erg/s]',fontsize = 14) 
plt.title(f'Luminosity at Ts={Ts} vs Decay Index From Normalised Two-Point Fits',fontsize = 14) 
plt.tick_params(axis='both', which='major', labelsize=18)
#plt.errorbar(alpha, lum_ts, yerr=lum_ts_err, xerr=alpha_err, linestyle='',elinewidth=0.5)
plt.legend()



# In[17]:


#Normalised ±U fits
plt.figure(figsize = (12,8 ))

        
for i in range(len(grb_name)):
    plt.scatter(norm_alpha_U[i], norm_lum_ts_U[i],label=grb_name[i])
    

#print(np.log10(lum_ts_err)
plt.xlabel('Alpha',fontsize = 14) 
plt.ylabel('L_Ts=10 [erg/s]',fontsize = 14) 
plt.title(f'Luminosity at Ts={Ts} vs Decay Index From Normalised 10-40s Fits',fontsize = 14) 
plt.tick_params(axis='both', which='major', labelsize=18)
#plt.errorbar(alpha, lum_ts, yerr=lum_ts_err, xerr=alpha_err, linestyle='',elinewidth=0.5)

plt.legend()



# In[112]:


#Normalised ± 30% fits
plt.figure(figsize = (12,8 ))

        
for i in range(len(grb_name)):
    plt.scatter(norm_alpha_per[i], norm_lum_ts_per[i],label=grb_name[i])
    

#print(np.log10(lum_ts_err)
plt.xlabel('Alpha',fontsize = 14) 
plt.ylabel('L_Ts=10 [erg/s]',fontsize = 14) 
plt.title(f'Luminosity at Ts={Ts} vs Decay Index From Normalised ± 30% Fits',fontsize = 14) 
plt.tick_params(axis='both', which='major', labelsize=18)
#plt.xlim(-4,4)
#plt.ylim(-3,3)
#plt.errorbar(alpha, lum_ts, yerr=lum_ts_err, xerr=alpha_err, linestyle='',elinewidth=0.5)
plt.legend()



# In[111]:



fig, (ax1, ax2,ax3) = plt.subplots(1, 3,figsize=(25,10))#, sharex=True, sharey=True)
fig.suptitle('Lum vs Alpha', fontsize=19)

for i in range(len(grb_name)):
    if norm_alpha_tp[i]>-10:
        ax1.scatter(norm_alpha_tp[i], norm_lum_ts_tp[i])#,label=grb_name[i])
        ax2.scatter(norm_alpha_U[i], norm_lum_ts_U[i])#,label=grb_name[i])
        ax3.scatter(norm_alpha_per[i], norm_lum_ts_per[i])#,label=grb_name[i])
    else:
        ax2.scatter(norm_alpha_U[i], norm_lum_ts_U[i])#,label=grb_name[i])
        ax3.scatter(norm_alpha_per[i], norm_lum_ts_per[i])#,label=grb_name[i])
    

#print(np.log10(lum_ts_err)
ax1.set_xlabel('Alpha',fontsize = 14) 
ax2.set_xlabel('Alpha',fontsize = 14) 
ax3.set_xlabel('Alpha',fontsize = 14) 
ax1.set_ylabel('L_Ts=10 [erg/s]',fontsize = 14) 
ax1.set_title('Two-Point Lum')
ax2.set_title('10-40s Lum')
ax3.set_title('±30% Ts Lum')
#plt.title(f'Luminosity at Ts={Ts} vs Decay Index From Normalised ± 30% Fits',fontsize = 14) 
#plt.tick_params(axis='both', which='major', labelsize=18)


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




