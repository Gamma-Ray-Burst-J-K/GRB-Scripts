#!/usr/bin/env python
# coding: utf-8

# In[6]:


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

import sklearn

mpl.rc('font', family='serif')


grb = fits.open('gll_2flgc.fits')

table=grb[1].data


# In[7]:


def power_law1(x, N, a): 
    '''simple power law 
    L = Nx^-a'''
    return N*np.power(x, a)


# In[ ]:





# In[8]:


def power_law3(x,N,a):
    '''simple power law in log space
    log(L) = log(N) - alogx'''
    return (a)*np.log10(x) + np.log10(N) 


# In[9]:


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
            
            
            

            
    
    
    
    


# In[10]:


#print(retrieve_grb_data(table,'080916C',5 ))
#new grb 140102A
#retrieve_grb_data(table,'140102A',5 )


# In[86]:


#def luminosity(table):

v=0
n=0
redshift = []
m=0

Ts = 20
upper_bound = 40
lower_bound = 5
        
j = 0
diff_mean = []
diff_std = []

#arrays for collecting terms
alpha = []
alpha_err = []

norm_alpha_tp, norm_lum_tp, norm_amp_tp = [],[], []
norm_alpha_tp_err = []

norm_alpha_ts, norm_lum_ts, norm_amp_ts = [],[], []
norm_alpha_ts_err = []

norm_alpha_U, norm_lum_ts_U, norm_amp_U = [], [], []
norm_alpha_U_err = []

norm_alpha_per, norm_lum_ts_per, norm_amp_per = [], [],[]
norm_alpha_per_err = []

norm_alpha_peak, norm_lum_ts_peak, norm_amp_peak = [], [], []
norm_alpha_peak_err = []

all_T, all_lum, all_T_err, all_lum_err = [], [], [], []
all_a_ts, all_a_U, all_a_per, all_a_peak, all_a_tp = [], [], [], [], []

brightness_ts, brightness_U, brightness_per, brightness_peak, brightness_tp = [], [], [], [], []


all_parameters = []

lum_params = []

hist_lum_1, hist_lum_10, hist_lum_100 = [], [], []     #for plotting histogram, inputting luminosities at 1s,10s, 20s up to 50s

all_names = []

lum_ts = []
lum_ts_err = []
grb_name = []

for i in table:
    if i['LUMINOSITY_DISTANCE'] >0:
        z = i['REDSHIFT']
        name = i['GCNNAME']
        dl = i['LUMINOSITY_DISTANCE']
        diff = []
        all_names.append(name)
        
        n=n+1
        #print(n)


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
        B = -1*indec               #calculating spectral index from measured photon indicies
        B_err = index_err
        
        #calculating weighted mean of spectral index and error
        B_sum = []
        B_err_sum = []

        if len(B)==1:
            mean_B = B
            #mean_B_err = abs(B_err)
        else:
            mean_B = np.average(B, weights=abs(1/B_err))
            #mean_B_err = np.mean(abs(B_err))
           
                
            
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

        if idx1+1 == len(T):
            if T[idx1]<Ts:
                idx2 = idx1 
                idx1 = idx1 - 1
            if T[idx1]>Ts:
                idx2 = idx1 - 1

        else:
            if T[idx1]<Ts:
                idx2 = idx1+1
            if T[idx1]>Ts:
                idx2 = idx1 -1

            
        #if T[idx1]> (Ts+50) or T[idx2]>(Ts+50): #only choosing GRBs with points that fall within ± 30 of Ts
        #    pass
        #else:

        #if T[idx1] > len(T) or T[idx2] > len(T) or T[idx1]>(Ts+50) or T[idx2]<(Ts+50):

        if len(ene_flux)<=2:
            pass
        else:
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


            #fitting from ±U
            diff_minus1 = abs(T-lower_bound)                   #lower T threshold
            diff_plus1 = abs(upper_bound-T)                #upper T threshold
            idx_minus1 = np.argmin(diff_minus1)+1          #finding lower threshold index
            idx_plus1 = np.argmin(diff_plus1)       #finding upper threshold index
            T_min_plus1 = T[idx_minus1:idx_plus1]         #making new T array
            lum_min_plus1 = lum[idx_minus1:idx_plus1]     #making new L array
            weight_mp1 = weight_lum[idx_minus1:idx_plus1]    #making new weights array
            

            #fitting from ±% Ts
            p = 0.3                 #percentage value
            d_minus2 = 10**((1-p)*np.log10(Ts)) 
            d_plus2 =  10**((1+p)*np.log10(Ts))
            diff_minus2 = abs(T-d_minus2)
            diff_plus2 = abs(T-d_plus2)
            idx_minus2 = np.argmin(diff_minus2)
            idx_plus2 = np.argmin(diff_plus2)+1
            T_min_plus2 = T[idx_minus2:idx_plus2]
            lum_min_plus2 = lum[idx_minus2:idx_plus2]
            weight_mp2 = weight_lum[idx_minus2:idx_plus2]
            
            #fitting from peak luminosity onwards
            peak_idx = np.argmax(lum)
            peak_lum = lum[peak_idx:]
            peak_T = T[peak_idx:]
            peak_weight = weight_lum[peak_idx:]

            
            #not going ahead with fits where the indexes are too small to fit
            if len(T_min_plus1)<=0 or len(T_min_plus2)<=0:
                pass
            else:

                print(name)
                
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

                # make parameters with starting values:
                par3 = model1.make_params(pow_amplitude=1e6, pow_exponent=-1.0)


                np.nan_to_num(weight_lum_new,copy=False)
                
                #running ts onwards fits & params
#                 resultn2 = model1.fit(lum_new/1e50, par3, x=T_new, weights=weight_ts) #linear ts with weight
#                 norma2, normN2 = resultn2.best_values['pow_exponent'],resultn2.best_values['pow_amplitude'] #norm two-point params
#                 normci2 = lmfit.conf_interval(resultn2, resultn2, sigmas=[0.68])


                
                #running ± U Ts fits & params
#                 resultn3 = model1.fit(lum_min_plus1/1e50, par3, x=T_min_plus1, weights=weight_mp1) #±U Ts luminosity 
#                 norma3, normN3 = resultn3.best_values['pow_exponent'],resultn3.best_values['pow_amplitude'] #±U Ts luminosity params
#                 normci3 = lmfit.conf_interval(resultn3, resultn3, sigmas=[0.68])
                
                
                #running ± 30 percent fits & parameters
                resultn4 = model1.fit(lum_min_plus2/1e50, par3, x=T_min_plus2, weights=weight_mp2) #±30 percent lum 
                norma4, normN4 = resultn4.best_values['pow_exponent'],resultn4.best_values['pow_amplitude'] #±30 percent lum params
                normci4 = lmfit.conf_interval(resultn4, resultn4, sigmas=[0.68])
                                

                #running two-point fits & params
                #resultn5 = model1.fit([lum1/1e50,lum2/1e50], par3, x=[T1,T2], weights=[w1,w1]) #norm two-point
                #norma5, normN5 = resultn5.best_values['pow_exponent'],resultn5.best_values['pow_amplitude'] #norm two-point params
                #normci5 = lmfit.conf_interval(resultn5, resultn5, sigmas=[0.68])
                
                
                #running peak fits & params
                resultn6 = model1.fit(peak_lum/1e50, par3, x=peak_T, weights=peak_weight)
                norma6, normN6 = resultn6.best_values['pow_exponent'],resultn6.best_values['pow_amplitude']
                normci6 = lmfit.conf_interval(resultn6, resultn6, sigmas=[0.68])


                #extracting sigma errors
                normcon2, normcon5, normcon3, normcon4, normcon6 = [],[],[],[], []

                for key, value in normci2.items():
                    normcon2.append(value)
                    
                #for key, value in normci5.items():
                #    normcon5.append(value)

                for key, value in normci3.items():
                    normcon3.append(value)

                for key, value in normci4.items():
                    normcon4.append(value)
                    
                for key, value in normci6.items():
                    normcon6.append(value)


                #errors are an array with [sigma minus, sigma plus]
                #normalised ts onwards lum errors
                norma2_sig = [normcon2[1][0][1], normcon2[1][2][1]]
                normN2_sig = [normcon2[0][0][1], normcon2[0][2][1]]
                
                #normalised 2-point lum errors
                #norma5_sig = [normcon5[1][0][1], normcon5[1][2][1]]
                #normN5_sig = [normcon5[0][0][1], normcon5[0][2][1]]

#                 #normalised ±U Ts lum
#                 norma3_sig = [normcon3[1][0][1], normcon3[1][2][1]]
#                 normN3_sig = [normcon3[0][0][1], normcon3[0][2][1]]

                #normalised ±30% Ts lum
                norma4_sig = [normcon4[1][0][1], normcon4[1][2][1]]
                normN4_sig = [normcon4[0][0][1], normcon4[0][2][1]]
                
                #normalised peak-Ts lum
                norma6_sig = [normcon6[1][0][1], normcon6[1][2][1]]
                normN6_sig = [normcon6[0][0][1], normcon6[0][2][1]]

                #all_parameters.append([name, mean_B, mean_B_err, a1,abs(a1_sig[0]-a1),abs(a1_sig[1]-a1),N1,abs(N1_sig[0]-N1), abs(N1_sig[1]-N1),
                 #                     a2,abs(np.log10(a2_sig[0])-a2), abs(a2_sig[1]-a2),N2, abs(N2_sig[0]-N2), abs(N2_sig[1]-N2)])

                #all_parameterslog.append([name, mean_B, mean_B_err, a3,abs(np.log10(a3_sig[0])-a3), abs(np.log10(a3_sig[1])-a3),N3, abs(N3_sig[0]-N3), abs(N3_sig[1]-N3),
                #                      a4,abs(np.log10(a4_sig[0])-a4), abs(np.log10(a4_sig[1])-a4),N4, abs(N4_sig[0]-N4), abs(N4_sig[1]-N4)])


                #pln5 = power_law1(Ts, normN5,norma5)*1e50
                
#                 print('normalised luminosity from ts onwards')
#                 #print(lmfit.fit_report(resultn5.params))
#                 pln2 = power_law1(Ts, normN2,norma2)*1e50
#                 print(f'luminosity = {np.round(power_law1(Ts, normN2,norma2),5)*1e50}')
#                 print(f'-{abs(power_law1(Ts,normN2_sig[0],norma2)*1e50 - pln2)}')
#                 print(f'+{abs(pln2 - power_law1(Ts,normN2_sig[1],norma2)*1e50)}')
#                 print()

#                 print('normalised luminosity Ts±10')
#                 #print(lmfit.fit_report(resultn3.params))
#                 pln3 = power_law1(Ts, normN3,norma3)*1e50
#                 print(f'luminosity = {power_law1(Ts, normN3, norma3)*1e50}')
#                 print(f'-{abs(power_law1(Ts,normN3_sig[0],norma3)*1e50 - pln3)}')
#                 print(f'+{abs(pln3 - power_law1(Ts,normN3_sig[1],norma3)*1e50)}')
#                 print()

                print('normalised luminosity Ts±30%')
                #print(lmfit.fit_report(resultn4.params))
                pln4 = power_law1(Ts, normN4,norma4)*1e50
                print(f'luminosity = {power_law1(Ts, normN4, norma4)*1e50}')
                print(f'-{abs(power_law1(Ts,normN4_sig[0], norma4)*1e50 - pln4)}')
                print(f'+{abs(pln4 - power_law1(Ts,normN4_sig[1], norma4)*1e50)}')
                print()
                
                print('normalised luminosity peak-Ts')
                #print(lmfit.fit_report(resultn4.params))
                pln6 = power_law1(Ts, normN6,norma6)*1e50
                print(f'luminosity = {power_law1(Ts, normN6, norma6)*1e50}')
                print(f'-{abs(power_law1(Ts,normN6_sig[0], norma6)*1e50 - pln6)}')
                print(f'+{abs(pln6 - power_law1(Ts,normN6_sig[1], norma6)*1e50)}')
                print()

                
                
                lum_fitn2 = power_law1(T,normN2, norma2) *1e50
                #lum_fitn5 = power_law1(T,normN5, norma5) *1e50
                #lum_fitn3 = power_law1(T,normN3, norma3) *1e50
                lum_fitn4 = power_law1(T,normN4, norma4) *1e50
                lum_fitn6 = power_law1(T,normN6, norma6) *1e50


                #alpha_err.append(stda2)
                norm_lum_ts_per.append(pln4)
                norm_alpha_per.append(norma4)
                norm_amp_per.append(normN4)
                norm_alpha_per_err.append(norma4_sig[1])

                ##norm_lum_ts.append(pln2)
                #norm_alpha_ts.append(norma2)
                #norm_amp_ts.append(normN2)

#                 norm_lum_ts_U.append(pln3)
#                 norm_alpha_U.append(norma3)
#                 norm_amp_U.append(normN3)
#                 norm_alpha_U_err.append(norma3_sig[1])
                
                #norm_lum_tp.append(pln5)
                #norm_alpha_tp.append(norma5)
                #norm_amp_tp.append(normN5)
                #norm_alpha_tp_err.append(norma5_sig[1])
                
                norm_lum_ts_peak.append([name,pln6])
                norm_alpha_peak.append(norma6)
                norm_amp_peak.append(normN6)
                norm_alpha_peak_err.append(norma6_sig[1])

                grb_name.append(name)


                v=v+1
                #print(v)
                #plotting different luminosity methods
                #fig, axs = plt.subplots(4,figsize=(12,30), sharex=False, sharey=True)

                #Ts Luminosity
                #plt.figure(figsize = (12,8 ))
                #plt.scatter(T, lum, label=name )
                #plt.yscale('log')
                #plt.xscale('log')
                #plt.xlabel(xlabel='T-T0/(1+z) [s]',fontsize = 14) 
                #plt.ylabel(ylabel='L [erg/s]',fontsize = 14)
                #plt.title(f'GRB{name} for Ts={Ts}s onwards Luminosity fit',fontsize = 14) 
                #plt.tick_params(axis='both', which='major', labelsize=14)
                #plt.errorbar(T, lum, yerr=lum_err, xerr=T_err, linestyle='',elinewidth=0.5)
                #plt.axvline(Ts, color = 'red', linewidth=0.5)#, label=f'{Ts}')   
                #plt.scatter(T_new, lum_new, color='green')
                #plt.axhline(power_law1(Ts, normN2, norma2)*1e50, color='green', linewidth=1.5)
                #axs[0].plot(T, lum_fitn2, linestyle='--', linewidth=2,color='green')
                #plt.legend()
                
                '''±U Ts Luminosity'''
#                 plt.figure(figsize = (12,8 ))
#                 plt.scatter(T, lum, label=name )
#                 plt.yscale('log')
#                 plt.xscale('log')
#                 plt.xlabel(xlabel='T-T0/(1+z) [s]',fontsize = 14) 
#                 plt.ylabel(ylabel='L [erg/s]',fontsize = 14)
#                 plt.title(f'GRB{name} for Ts={Ts}s {lower_bound}-{upper_bound}s Luminosity fit',fontsize = 14) 
#                 plt.tick_params(axis='both', which='major', labelsize=14)
#                 plt.errorbar(T, lum, yerr=lum_err, xerr=T_err, linestyle='',elinewidth=0.5)
#                 plt.axvline(Ts, color = 'red', linewidth=0.5)#, label=f'{Ts}')   
#                 plt.scatter(T_min_plus1, lum_min_plus1, color='black')
#                 plt.axhline(power_law1(Ts, normN3, norma3)*1e50, color='black', linewidth=1.5)
#                 plt.plot(T, lum_fitn3, linestyle='--', linewidth=2,color='black')
#                 plt.legend()

                '''±30% Ts Luminosity'''
                plt.figure(figsize = (12,8 ))
                plt.scatter(T, lum, label=name )
                plt.yscale('log')
                plt.xscale('log')
                plt.xlabel(xlabel='T-T0/(1+z) [s]',fontsize = 14) 
                plt.ylabel(ylabel='L [erg/s]',fontsize = 14) 
                plt.title(f'GRB{name} for Ts={Ts}s ±{100*p}% Luminosity fit',fontsize = 14) 
                plt.tick_params(axis='both', which='major', labelsize=14)
                plt.errorbar(T, lum, yerr=lum_err, xerr=T_err, linestyle='',elinewidth=0.5)
                plt.axvline(Ts, color = 'red', linewidth=0.5)#, label=f'{Ts}')   
                plt.legend()
                plt.scatter(T_min_plus2, lum_min_plus2, color='orange')
                plt.axhline(power_law1(Ts, normN4, norma4)*1e50,color='orange', linewidth=1.5)
                plt.plot(T, lum_fitn4, linestyle='--', linewidth=2,color='orange')
                plt.show()
                
                '''Peak Luminosity'''
                plt.figure(figsize = (12,8 ))
                plt.scatter(T, lum, label=name )
                plt.yscale('log')
                plt.xscale('log')
                plt.xlabel(xlabel='T-T0/(1+z) [s]',fontsize = 14) 
                plt.ylabel(ylabel='L [erg/s]',fontsize = 14) 
                plt.title(f'GRB{name} for Ts={Ts}s Peak Luminosity fit',fontsize = 14) 
                plt.tick_params(axis='both', which='major', labelsize=14)
                plt.errorbar(T, lum, yerr=lum_err, xerr=T_err, linestyle='',elinewidth=0.5)
                plt.axvline(Ts, color = 'red', linewidth=0.5)#, label=f'{Ts}')   
                plt.legend()
                plt.scatter(peak_T, peak_lum, color='red')
                plt.axhline(power_law1(Ts, normN6, norma6)*1e50,color='red', linewidth=1.5)
                plt.plot(T, lum_fitn6, linestyle='--', linewidth=2,color='red')
                plt.show()
                
                

                for i in T:
                    all_T.append(i)

                for i in T_err:
                    all_T_err.append(i)

                for i in lum:
                    all_lum.append(i)
                    #all_a_ts.append(norma2)
                    #all_a_tp.append(norma5)
                    #all_a_U.append(norma3)
                    all_a_per.append(norma4)
                    all_a_peak.append(norma6)
                    #brightness_ts.append(power_law1(Ts, normN2, norma2))
                    #brightness_U.append(power_law1(Ts, normN3, norma3))
                    brightness_per.append(power_law1(Ts, normN4, norma4))
                    brightness_peak.append(power_law1(Ts, normN6, norma6))

                for i in lum_err:
                    all_lum_err.append(i)



                #np.savetxt(f'LAT_GRB{name}_data',np.column_stack((lum, lum_err, T,T_err)),delimiter=' ')

                #np.savetxt(f'LAT_GRB{name}_log_data',np.column_stack((log_T,T_err_plus,T_err_minus,
                                                                    #  log_lum, lum_err_plus, lum_err_minus)),delimiter=' ')




                


# In[87]:


len(grb_name)


# In[105]:


#def luminosity(table):

v=0
n=0
redshift = []
m=0

Ts = 20
lower_bound = 5
upper_bound = 40

all_a_whole, all_lum_whole = [], []
j = 0
diff_mean = []
diff_std = []

plt.figure(figsize = (12,8 ))

lum_ts = []
lum_ts_err = []
#grb_name = []

for i in table:
    if i['LUMINOSITY_DISTANCE'] >0:
        z = i['REDSHIFT']
        name = i['GCNNAME']
        dl = i['LUMINOSITY_DISTANCE']
        diff = []
        
        n=n+1
        #print(n)


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
        B = -1*indec               #calculating spectral index from measured photon indicies
        B_err = index_err
        
        #calculating weighted mean of spectral index and error
        B_sum = []
        B_err_sum = []

        if len(B)==1:
            mean_B = B
            #mean_B_err = abs(B_err)
        else:
            mean_B = np.average(B, weights=abs(1/B_err))
            #mean_B_err = np.mean(abs(B_err))
           
                
            
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

        if idx1+1 == len(T):
            if T[idx1]<Ts:
                idx2 = idx1 
                idx1 = idx1 - 1
            if T[idx1]>Ts:
                idx2 = idx1 - 1

        else:
            if T[idx1]<Ts:
                idx2 = idx1+1
            if T[idx1]>Ts:
                idx2 = idx1 -1

            
        #if T[idx1]> (Ts+50) or T[idx2]>(Ts+50): #only choosing GRBs with points that fall within ± 30 of Ts
        #    pass
        #else:

        #if T[idx1] > len(T) or T[idx2] > len(T) or T[idx1]>(Ts+50) or T[idx2]<(Ts+50):

        if len(ene_flux)<=2:
            pass
        else:
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


            #fitting from ±U
#             diff_minus1 = abs(T-lower_bound)                   #lower T threshold
#             diff_plus1 = abs(upper_bound-T)                #upper T threshold
#             idx_minus1 = np.argmin(diff_minus1)+1          #finding lower threshold index
#             idx_plus1 = np.argmin(diff_plus1)       #finding upper threshold index
#             T_min_plus1 = T[idx_minus1:idx_plus1]         #making new T array
#             lum_min_plus1 = lum[idx_minus1:idx_plus1]     #making new L array
#             weight_mp1 = weight_lum[idx_minus1:idx_plus1]    #making new weights array
            

            #fitting from ±% Ts
            p = 0.3                 #percentage value
            d_minus2 = 10**((1-p)*np.log10(Ts)) 
            d_plus2 =  10**((1+p)*np.log10(Ts))
            diff_minus2 = abs(T-d_minus2)
            diff_plus2 = abs(T-d_plus2)
            idx_minus2 = np.argmin(diff_minus2)
            idx_plus2 = np.argmin(diff_plus2)+1
            T_min_plus2 = T[idx_minus2:idx_plus2]
            lum_min_plus2 = lum[idx_minus2:idx_plus2]
            weight_mp2 = weight_lum[idx_minus2:idx_plus2]
            
            #fitting from peak luminosity onwards
            peak_idx = np.argmax(lum)
            peak_lum = lum[peak_idx:]
            peak_T = T[peak_idx:]
            peak_weight = weight_lum[peak_idx:]

            
            #not going ahead with fits where the indexes are too small to fit
            if len(T_min_plus1)<=1 or len(T_min_plus2)<=1:
                pass
            else:

                v=v+1
                print(name)
                print(v)
                
                
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

                # make parameters with starting values:
                par3 = model1.make_params(pow_amplitude=1e6, pow_exponent=-1.0)


                np.nan_to_num(weight_lum_new,copy=False)
                
                #running whole fits & params
#                 resultn1 = model1.fit(lum/1e50, par3, x=T, weights=weight_lum) 
#                 norma1, normN1 = resultn1.best_values['pow_exponent'],resultn1.best_values['pow_amplitude'] #norm whole params
#                 normci1 = lmfit.conf_interval(resultn1, resultn1, sigmas=[0.68])
                
                
                #running ts onwards fits & params
#                 resultn2 = model1.fit(lum_new/1e50, par3, x=T_new, weights=weight_ts) #linear ts with weight
#                 norma2, normN2 = resultn2.best_values['pow_exponent'],resultn2.best_values['pow_amplitude'] #norm two-point params
#                 normci2 = lmfit.conf_interval(resultn2, resultn2, sigmas=[0.68])


                
                #running ± U Ts fits & params
#                 resultn3 = model1.fit(lum_min_plus1/1e50, par3, x=T_min_plus1, weights=weight_mp1) #±U Ts luminosity 
#                 norma3, normN3 = resultn3.best_values['pow_exponent'],resultn3.best_values['pow_amplitude'] #±U Ts luminosity params
#                 normci3 = lmfit.conf_interval(resultn3, resultn3, sigmas=[0.68])
                
                
                #running ± 30 percent fits & parameters
                resultn4 = model1.fit(lum_min_plus2/1e50, par3, x=T_min_plus2, weights=weight_mp2) #±30 percent lum 
                norma4, normN4 = resultn4.best_values['pow_exponent'],resultn4.best_values['pow_amplitude'] #±30 percent lum params
                normci4 = lmfit.conf_interval(resultn4, resultn4, sigmas=[0.68])
                                

                #running two-point fits & params
                #resultn5 = model1.fit([lum1/1e50,lum2/1e50], par3, x=[T1,T2], weights=[w1,w1]) #norm two-point
                #norma5, normN5 = resultn5.best_values['pow_exponent'],resultn5.best_values['pow_amplitude'] #norm two-point params
                #normci5 = lmfit.conf_interval(resultn5, resultn5, sigmas=[0.68])
                
                
                #running peak fits & params
                resultn6 = model1.fit(peak_lum/1e50, par3, x=peak_T, weights=peak_weight)
                norma6, normN6 = resultn6.best_values['pow_exponent'],resultn6.best_values['pow_amplitude']
                normci6 = lmfit.conf_interval(resultn6, resultn6, sigmas=[0.68])


                #lum_fitn1 = power_law1(T,normN1, norma1) *1e50
                #lum_fitn2 = power_law1(T,normN2, norma2) *1e50
                #lum_fitn5 = power_law1(T,normN5, norma5) *1e50
                #lum_fitn3 = power_law1(T,normN3, norma3) *1e50
                lum_fitn4 = power_law1(T,normN4, norma4) *1e50
                lum_fitn6 = power_law1(T,normN6, norma6) *1e50
                
                #all_a_whole.append(norma1)
                #all_lum_whole.append(power_law1(Ts, normN1, norma1)*1e50)
                
                if norma6 < -1.1:         #3 for ±U, 4 for ±30%, 5 for tp, 6 for peak

                    plt.plot(T, lum_fitn6, linestyle='--', linewidth=2,color='red')
                else:
                    plt.plot(T, lum_fitn6, linestyle='--', linewidth=2,color='blue')
                    
                if norma6 > 0:
                    plt.plot(T,lum_fitn6, linestyle='--', linewidth=2, color='g', 
                             label=f'{name} - alpha={np.round(norma6,3)}' )

                plt.yscale('log')
                plt.xscale('log')
                plt.scatter(T, lum, color='lightgray')
                #plt.ylim(1e39,1e60)
                plt.xlabel(xlabel='T-T0/(1+z) [s]',fontsize = 14) 
                plt.ylabel(ylabel='L [erg/s]',fontsize = 14) 
                #plt.title(f'All GRBs from {lower_bound}-{upper_bound}s Fit Luminosity Fit',fontsize = 14) 
                #plt.title(f'All GRBs from Ts={Ts}s ±30% Luminosity Fit',fontsize = 14)
                plt.title(f'All GRBs from Ts={Ts}s Peak Onward Luminosity Fit',fontsize = 14)
                #plt.title(f'All GRBs from Ts={Ts}s Two-Point Luminosity Fit',fontsize = 14)
                plt.tick_params(axis='both', which='major', labelsize=14)

                #plt.legend()
                plt.axvline(Ts, color = 'red', linewidth=0.5)


# In[ ]:





# In[95]:




print(len(grb_name))


# In[107]:


'''Plotting Lum at Ts vs Alpha from Racusin et al'''

ras_alpha = [0.96, 1.59, 1.11, 1.21, 1.28]

grbs = ['080916C', '90510.0', '090926A', '110731A', '130427A']

norm_lum_ts_peak[0]

data = []
n=0

#matching the GRBs based on name to then plot their lum_ts vs Rascun alpha's
for i in table:
    for k in grbs:
        for j in range(len(norm_lum_ts_peak)):
            if k == i['GCNNAME'] and k == norm_lum_ts_peak[j][0]:
                data.append([i['GCNNAME'], norm_lum_ts_peak[j][1], ras_alpha[n]])
                n=n+1

    

    
plt.figure(figsize = (8,5 ))
luminosity = []
for l in data:
    plt.scatter(l[2], l[1], label = l[0])
    luminosity.append(l[1])



plt.xlabel('Alpha',fontsize = 14) 
plt.ylabel('Luminosity',fontsize = 14) 
plt.yscale('log')
plt.title(f'Ts={Ts}s for Rascun et. al. Alphas Luminosity vs Decay ',fontsize = 14) 
plt.tick_params(axis='both', which='major', labelsize=14)
plt.legend()
plt.show()
    
stats.spearmanr(ras_alpha, np.log10(luminosity))


# In[ ]:





# In[22]:


'''Finding Mean and Standard Deviations of Alpha'''
alpha_per_weights, alpha_U_weights, alpha_tp_weights, alpha_peak_weights = [], [], [], []

for i in range(len(norm_alpha_per)):
    alpha_per_weights.append(1/abs(norm_alpha_per_err[i]-norm_alpha_per[i]))
    
mean_per_alpha = np.average(norm_alpha_per, weights=alpha_per_weights)

# print(f'Per no weights mean = {np.mean(norm_alpha_per)}')
# print(f'Per with weights mean = {mean_per_alpha}')
# print(f'Per std = {np.std(norm_alpha_per)}')
# print()


for i in range(len(norm_alpha_U)):
    alpha_U_weights.append(1/abs(norm_alpha_U_err[i]-norm_alpha_U[i]))
    
mean_U_alpha = np.average(norm_alpha_U, weights=alpha_U_weights)

# print(f'U no weights mean = {np.mean(norm_alpha_U)}')
# print(f'U with weights mean = {mean_U_alpha}')
# print(f'U std = {np.std(norm_alpha_U)}')
# print()


for i in range(len(norm_alpha_peak)):
    alpha_peak_weights.append(1/abs(norm_alpha_peak_err[i]-norm_alpha_peak[i]))
    
mean_peak_alpha = np.average(norm_alpha_peak, weights=alpha_peak_weights)

print(f'Peak no weights mean = {np.mean(norm_alpha_peak)}')
print(f'Peak with weights mean = {mean_peak_alpha}')
print(f'Peak std = {np.std(norm_alpha_peak)}')


# In[9]:


plt.figure(figsize = (12,8 ))
plt.title(f'Ts={Ts}s Alpha Distribution')
#plt.ylabel('Number of GRBs')
plt.xlabel('Alpha')
plt.hist(norm_alpha_U, color='black', label=f'{lower_bound}-{upper_bound}s Fit')
plt.hist(norm_alpha_per, color='orange', label='Ts±30% Fit ')
plt.hist(norm_alpha_peak, color='red', label='Peak Fit')
plt.legend()
plt.show()


# In[10]:


plt.figure(figsize = (8,5 ))
plt.title(f'Alpha Distribution')
#plt.ylabel('Number of GRBs')
plt.xlabel('Alpha')
plt.hist(norm_alpha_U, color='black')

plt.show()


# In[11]:


plt.figure(figsize = (8,5 ))
plt.title(f'Alpha Distribution')
#plt.ylabel('Number of GRBs')
plt.xlabel('Alpha')

plt.hist(norm_alpha_per, color='orange')

plt.show()


# In[12]:


plt.figure(figsize = (8,5 ))
plt.title(f'Alpha Distribution')
#plt.ylabel('Number of GRBs')
plt.xlabel('Alpha')

plt.hist(norm_alpha_peak, color='red')
plt.show()


# In[ ]:





# In[ ]:





# In[108]:


'''Flux Version/ Not Corrected to rest frame'''

#def luminosity(table):

v=0
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



#norm_alpha_peak, norm_lum_ts_peak, norm_amp_peak = [], [], []



brightness_ts, brightness_U, brightness_per, brightness_peak, brightness_tp = [], [], [], [], []


all_parameters = []

lum_params = []

hist_lum_1, hist_lum_10, hist_lum_100 = [], [], []     #for plotting histogram, inputting luminosities at 1s,10s, 20s up to 50s

plt.figure(figsize = (12,8 ))

lum_ts = []
lum_ts_err = []
grb_name1 = []

for i in table:
    #if i['LUMINOSITY_DISTANCE'] >0:
    z = i['REDSHIFT']
    name = i['GCNNAME']
    dl = i['LUMINOSITY_DISTANCE']
    diff = []

    n=n+1
    #print(n)


    for k in i['LC_FLUX_ERR']:
        mask = np.where(i['LC_FLUX_ERR']!=0) #finding indexs of values where flux error is 0 and masking
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


    g = 1.6e-6                     #extra factor missing in the LC_ENE_FLUX values in the .fits
    B = -1*indec               #calculating spectral index from measured photon indicies
    B_err = index_err

    #calculating weighted mean of spectral index and error
    B_sum = []
    B_err_sum = []

    if len(B)==1:
        mean_B = B
        #mean_B_err = abs(B_err)
    else:
        #mean_B = np.average(B, weights=abs(1/B_err))
        mean_B = np.mean(B)
        #mean_B_err = np.mean(abs(B_err))



    #calculating k correction
    k_correction = (1+z)**((mean_B-1)-1)

    #calculating corrected luminosities and their errors in normal and log space
    lum = g*ene_flux*4*(np.pi)*(dl**2)*(k_correction)
    lum_err = g*ene_flux_err*4*(np.pi)*(dl**2)*(k_correction)

    #calculating the asymmetric errors in log space
    lum_err_plus = np.log10(lum + lum_err) - np.log10(lum)        
    lum_err_minus= np.log10(lum) - np.log10(lum - lum_err)
    
    flux_err_plus = np.log10(flux + flux_err) - np.log10(flux)        
    flux_err_minus = np.log10(flux) - np.log10(flux - flux_err)


    T_err = ((t_end-t_start)/2)
    t = t_start+T_err
    T_err =  T_err

    T_err_plus = np.log10(t + T_err) - np.log10(t)
    T_err_minus = np.log10(t) - np.log10(t - T_err)
    
        #calculating the index of T either side of Ts
    diff = abs(t-Ts)
    #print(diff)
    
    
    if len(diff)==0:
        pass
    else:
        idx1 = np.argmin(diff)

        if idx1+1 == len(t):
            if t[idx1]<Ts:
                idx2 = idx1 
                idx1 = idx1 - 1
            if t[idx1]>Ts:
                idx2 = idx1 - 1

        else:
            if t[idx1]<Ts:
                idx2 = idx1+1
            if t[idx1]>Ts:
                    idx2 = idx1 -1


        if t[idx1] > (Ts+100) or t[idx2] > (Ts+100):
            pass
        else:
            
            if len(flux)<=2:
                pass
            else:
                #weight_lum = 1/lum_err

                weight_flux = 1/flux_err

                if t[idx1]>Ts:                 #makign sure the weights are the correct way round
                    T_new = t[idx2:]
                    flux_new = flux[idx2:]
                    flux_err_p_new = flux_err_plus[idx2:]
                    flux_err_m_new = flux_err_minus[idx2:]
                    weight_ts = weight_flux[idx2:]
                    weight_flux_new = 1/flux_err_p_new
                else:
                    T_new = t[idx1:]
                    flux_new = flux[idx1:]
                    flux_err_p_new = flux_err_plus[idx1:]
                    flux_err_m_new = flux_err_minus[idx1:]
                    weight_ts = weight_flux[idx1:]
                    weight_lum_new = 1/flux_err_p_new


                #fitting from ±U
                upper_bound = 60
                lower_bound = 10
                diff_minus1 = abs(t-lower_bound)                   #lower T threshold
                diff_plus1 = abs(upper_bound-t)                #upper T threshold
                idx_minus1 = np.argmin(diff_minus1)+1          #finding lower threshold index
                idx_plus1 = np.argmin(diff_plus1)       #finding upper threshold index
                T_min_plus1 = t[idx_minus1:idx_plus1]         #making new T array
                flux_min_plus1 = flux[idx_minus1:idx_plus1]     #making new L array
                weight_mp1 = weight_flux[idx_minus1:idx_plus1]    #making new weights array


                #fitting from ±% Ts
                p = 0.3                 #percentage value
                d_minus2 = 10**((1-p)*np.log10(Ts)) 
                d_plus2 =  10**((1+p)*np.log10(Ts))
                diff_minus2 = abs(t-d_minus2)
                diff_plus2 = abs(t-d_plus2)
                idx_minus2 = np.argmin(diff_minus2)
                idx_plus2 = np.argmin(diff_plus2)+1
                T_min_plus2 = t[idx_minus2:idx_plus2]
                flux_min_plus2 = flux[idx_minus2:idx_plus2]
                weight_mp2 = weight_flux[idx_minus2:idx_plus2]

                #fitting from peak luminosity onwards
                peak_idx = np.argmax(flux)
                peak_flux = flux[peak_idx:]
                peak_T = t[peak_idx:]
                peak_weight = weight_flux[peak_idx:]


                #calculating luminosity at Ts
                if idx1 > idx2:
                    lum1 = flux[idx1]       #finding luminosity data points either side of Ts
                    lum2 = flux[idx2]
                    T1 = t[idx1]
                    T2 = t[idx2]
                    w1 = 1/flux_err[idx1]
                    w2 = 1/flux_err[idx2]
                else:
                    lum1 = flux[idx2]
                    lum2 = flux[idx1]
                    T1 = t[idx2]
                    T2 = t[idx1]
                    w1 = 1/flux_err[idx2]
                    w2 = 1/flux_err[idx1]

                #calculating best fit parameters and covariances for the data lmfit
                #making powerlaw model for linear fits
                model1 = PowerLawModel(prefix='pow_')

                # make parameters with starting values:
                par3 = model1.make_params(pow_amplitude=1e6, pow_exponent=-1.0)

                # make parameters with starting values:
                par4 = model1.make_params(pow_amplitude=1e-3, pow_exponent=-1.0)


                np.nan_to_num(weight_flux_new,copy=False)

                #running ts onwards fits & params
                result1 = model1.fit(flux, par4, x=t, weights=(1/flux_err)) #linear ts with weight
                a1, N1 = result1.best_values['pow_exponent'],result1.best_values['pow_amplitude'] #norm two-point params
                ci1 = lmfit.conf_interval(result1, result1, sigmas=[0.68])




                #running peak fits & params
                resultn6 = model1.fit(peak_flux, par4, x=peak_T, weights=peak_weight)
                norma6, normN6 = resultn6.best_values['pow_exponent'],resultn6.best_values['pow_amplitude']
                normci6 = lmfit.conf_interval(resultn6, resultn6, sigmas=[0.68])


                #extracting sigma errors
                normcon2, normcon5, normcon3, normcon4, normcon6 = [],[],[],[], []

         
           
                pln6 = power_law1(Ts, normN6,norma6)
                

                lum_fit1 = power_law1(t, N1, a1)

                lum_fitn6 = power_law1(t,normN6, norma6)


             

                #norm_lum_ts_peak.append(pln6)
                #norm_alpha_peak.append(norma6)
                #norm_amp_peak.append(normN6)

                grb_name1.append(name)


                v=v+1
                #print(v)
                

                '''±U Ts Luminosity'''
                #plt.figure(figsize = (12,8 ))
                #plt.scatter(t, flux, label=name )
                #plt.yscale('log')
                #plt.xscale('log')
                #plt.xlabel(xlabel='T-T0/(1+z) [s]',fontsize = 14) 
                #plt.ylabel(ylabel='L [erg/s]',fontsize = 14)
                #plt.title(f'GRB{name} ',fontsize = 14) 
                #plt.tick_params(axis='both', which='major', labelsize=14)
                ##plt.errorbar(T, lum, yerr=lum_err, xerr=T_err, linestyle='',elinewidth=0.5)
                #plt.axvline(Ts, color = 'red', linewidth=0.5)#, label=f'{Ts}')   
                #plt.scatter(T_min_plus1, lum_min_plus1, color='black')
                #plt.axhline(power_law1(Ts, normN3, norma3)*1e50, color='black', linewidth=1.5)
                #plt.plot(t, lum_fit1, linestyle='--', linewidth=2,color='black')
                #plt.plot(t, lum_fitn6, linestyle='--', linewidth=2,color='red')
                #plt.legend()
                #plt.xlim(1e0, 5e1)
                #plt.show()



                if norma6 < -1.2:         #3 for ±U, 4 for ±30%, 5 for tp, 6 for peak

                    plt.plot(t, lum_fitn6, linestyle='--', linewidth=2,color='red', label = v)
                else:
                    plt.plot(t, lum_fitn6, linestyle='--', linewidth=2,color='blue', label = v)
                    
                if norma6 > 0:
                    plt.plot(t,lum_fitn6, linestyle='--', linewidth=2, color='g', 
                             label=f'{name} - alpha={np.round(norma6,3)}' )
                
                plt.yscale('log')
                plt.xscale('log')
                plt.scatter(t, flux, color='lightgray')
                #plt.ylim(1e39,1e60)
                plt.xlabel(xlabel='T -T0 [s]',fontsize = 14) 
                plt.ylabel(ylabel='Flux [ph/cm^2s]',fontsize = 14) 
                #plt.title(f'All GRBs from {lower_bound}-{upper_bound}s Luminosity Fit',fontsize = 14) 
                #plt.title(f'All GRBs from Ts ±30% Luminosity Fit',fontsize = 14)
                plt.title(f'All GRBs from Peak Onward Flux Fit',fontsize = 14)
                #plt.title(f'All GRBs from Two-Point Luminosity Fit',fontsize = 14)
                plt.tick_params(axis='both', which='major', labelsize=14)

                #plt.legend()
                plt.axvline(Ts, color = 'red', linewidth=0.5)
                
                



                


# In[ ]:





# In[ ]:





# In[90]:


'''Fast GRBs'''
alpha_tp_fast = []
lum_tp_fast = []

alpha_U_fast = []
lum_U_fast = []

alpha_per_fast = []
lum_per_fast = []

alpha_peak_fast = []
lum_peak_fast = []

threshold = -1.1


#for i in range(len(norm_alpha_tp)):
#    if norm_alpha_tp[i]< threshold:
#        alpha_tp_fast.append(norm_alpha_tp[i])
#        lum_tp_fast.append(norm_lum_tp[i])
        
#plt.figure(figsize = (8,5 ))
#plt.scatter(alpha_tp_fast, lum_tp_fast)
#plt.xlabel('Alpha',fontsize = 14) 
#plt.ylabel('Luminosity',fontsize = 14) 
#plt.title(f'Two-Point Luminosity vs Decay',fontsize = 14) 
#plt.tick_params(axis='both', which='major', labelsize=14)
#plt.show()
#print(f' Two Point Fit {stats.spearmanr(alpha_tp_fast, lum_tp_fast)}')
#print()




# for i in range(len(norm_alpha_U)):
#     if norm_alpha_U[i] < threshold:
#         alpha_U_fast.append(norm_alpha_U[i])
#         lum_U_fast.append(norm_lum_ts_U[i])
        
# plt.figure(figsize = (8,5 ))
# plt.scatter(alpha_U_fast, lum_U_fast, label=f'Threshold={threshold}')
# plt.xlabel('Alpha',fontsize = 14) 
# plt.ylabel('Luminosity',fontsize = 14) 
# plt.yscale('log')
# plt.title(f'Ts={Ts}s for {lower_bound}-{upper_bound}s Fit Luminosity vs Decay ',fontsize = 14) 
# plt.tick_params(axis='both', which='major', labelsize=14)
# plt.legend()
# plt.show()
# print(f'{lower_bound}-{upper_bound}s Fit {stats.spearmanr(alpha_U_fast, np.log10(lum_U_fast))}')
# print(len(alpha_U_fast))    
    
    
    
for i in range(len(norm_alpha_per)):
    if norm_alpha_per[i] < threshold:
        alpha_per_fast.append(norm_alpha_per[i])
        lum_per_fast.append(norm_lum_ts_per[i])
        
plt.figure(figsize = (8,5 ))
plt.scatter(alpha_per_fast, lum_per_fast, label=f'Threshold={threshold}')
plt.xlabel('Alpha',fontsize = 14) 
plt.ylabel('Luminosity',fontsize = 14)
plt.yscale('log')
plt.title(f'Ts={Ts}s for Ts ± 30%  Fit Luminosity vs Decay',fontsize = 14) 
plt.tick_params(axis='both', which='major', labelsize=14)
plt.legend()
plt.show()

print(f' 30% Fit {stats.spearmanr(alpha_per_fast, np.log10(lum_per_fast))}')
print(len(alpha_per_fast))
        
    
    
for i in range(len(norm_alpha_peak)):
    if norm_alpha_peak[i] < threshold:
        alpha_peak_fast.append(norm_alpha_peak[i])
        lum_peak_fast.append(norm_lum_ts_peak[i][1])

plt.figure(figsize = (8,5 ))
plt.scatter(alpha_peak_fast, lum_peak_fast, label=f'Threshold={threshold}')
plt.xlabel('Alpha',fontsize = 14) 
plt.ylabel('Luminosity',fontsize = 14) 
plt.yscale('log')
plt.title(f'Ts={Ts}s Peak Onwards Fit Luminosity vs Decay',fontsize = 14) 
plt.tick_params(axis='both', which='major', labelsize=14)
plt.legend()
plt.show()

print(f' Peak Fit {stats.spearmanr(alpha_peak_fast, np.log10(lum_peak_fast))}')
print(len(alpha_peak_fast))




#print(lum_U_fast)
#print(alpha_U_fast)

        


# In[ ]:





# In[92]:


'''Slow GRBs'''
alpha_ts_slow = []
lum_ts_slow = []

alpha_U_slow = []
lum_U_slow = []

alpha_per_slow = []
lum_per_slow = []

alpha_peak_slow = []
lum_peak_slow = []



#for i in range(len(norm_alpha_ts)):
#    if norm_alpha_ts[i] > threshold:
#        alpha_ts_slow.append(norm_alpha_ts[i])
#        lum_ts_slow.append(norm_lum_ts[i])
        
#plt.figure(figsize = (8,5 ))
#plt.scatter(alpha_ts_slow, lum_ts_slow)
#plt.xlabel('Alpha',fontsize = 14) 
#plt.ylabel('Luminosity',fontsize = 14) 
#plt.title(f'Ts Onwards',fontsize = 14) 
#plt.tick_params(axis='both', which='major', labelsize=14)
#plt.show()
#print(f' Peak Fit {stats.spearmanr(alpha_ts_slow, lum_ts_slow)}')
#print()


# for i in range(len(norm_alpha_U)):
#     if norm_alpha_U[i] > threshold:
#         alpha_U_slow.append(norm_alpha_U[i])
#         lum_U_slow.append(norm_lum_ts_U[i])
        
# plt.figure(figsize = (8,5 ))
# plt.scatter(alpha_U_slow, lum_U_slow)
# plt.xlabel('Alpha',fontsize = 14) 
# plt.yscale('log')
# plt.ylabel('Luminosity',fontsize = 14) 
# plt.title(f'Ts ±U ',fontsize = 14) 
# plt.tick_params(axis='both', which='major', labelsize=14)
# plt.show()
# print(f' 5-40s Fit {stats.spearmanr(alpha_U_slow, np.log10(lum_U_slow))}')
# print(len(alpha_U_slow))    
    
    
    
for i in range(len(norm_alpha_per)):
    if norm_alpha_per[i] > threshold:
        alpha_per_slow.append(norm_alpha_per[i])
        lum_per_slow.append(norm_lum_ts_per[i])
        
plt.figure(figsize = (8,5 ))
plt.scatter(alpha_per_slow, lum_per_slow)
plt.xlabel('Alpha',fontsize = 14) 
plt.ylabel('Luminosity',fontsize = 14) 
plt.yscale('log')
plt.title(f'±30% Onwards',fontsize = 14) 
plt.tick_params(axis='both', which='major', labelsize=14)
plt.show()

print(f' Per Fit {stats.spearmanr(alpha_per_slow, np.log10(lum_per_slow))}')
print(len(alpha_per_slow))  
    
    
for i in range(len(norm_alpha_peak)):
    if norm_alpha_peak[i] > threshold:
        alpha_peak_slow.append(norm_alpha_peak[i])
        lum_peak_slow.append(norm_lum_ts_peak[i][1])

plt.figure(figsize = (8,5 ))
plt.scatter(alpha_peak_slow, lum_peak_slow)
plt.xlabel('Alpha',fontsize = 14) 
plt.yscale('log')
plt.ylabel('Luminosity',fontsize = 14) 
plt.title(f'Peak Onwards',fontsize = 14) 
plt.tick_params(axis='both', which='major', labelsize=14)
plt.show()

print(f' Peak Fit {stats.spearmanr(alpha_peak_slow, np.log10(lum_peak_slow))}')
print(len(alpha_peak_slow))  




# In[40]:


#Spearman Ranks not separated 
srank_ts = stats.spearmanr(norm_lum_ts, norm_alpha_ts)

srank_U = stats.spearmanr(norm_lum_ts_U, norm_alpha_U)

srank_per = stats.spearmanr(norm_lum_ts_per, norm_alpha_per)

srank_peak = stats.spearmanr(norm_lum_ts_peak, norm_alpha_peak)


#print(f'ts fit {srank_ts}')
#print()
print(f'{lower_bound}-{upper_bound}s fit {srank_U}')
print()
print(f'±30% fit {srank_per}')
print()
print(f'peak fit {srank_peak}') 
print()


# In[ ]:





# In[66]:


'''± U Fit Luminosity at 1s, 10s, 100s'''

U_lum_1,U_lum_10, U_lum_100 = [], [], []

for i in range(len(norm_amp_U)):
    U_lum_1.append(power_law1(1, norm_amp_U[i],norm_alpha_U[i])*1e50)
    U_lum_10.append(power_law1(10, norm_amp_U[i],norm_alpha_U[i])*1e50)
    U_lum_100.append(power_law1(100, norm_amp_U[i],norm_alpha_U[i])*1e50)
    

plt.figure(figsize = (8,5 ))
plt.title(f'{lower_bound}-{upper_bound}s Luminosity Fit Distribution at 1s')
plt.ylabel('Number of GRBs')
plt.xlabel('log(L) at 1s')
plt.hist(np.log10(U_lum_1))
plt.show()

plt.figure(figsize = (8,5 ))
plt.title(f'{lower_bound}-{upper_bound}s Luminosity Fit Distribution at 10s')
plt.ylabel('Number of GRBs')
plt.xlabel('log(L) at 10s')
plt.hist(np.log10(U_lum_10))
plt.show()

plt.figure(figsize = (8,5 ))
plt.title(f'{lower_bound}-{upper_bound}s Luminosity Fit Distribution at 100s')
plt.ylabel('Number of GRBs')
plt.xlabel('log(L) at 100s')
plt.hist(np.log10(U_lum_100))
plt.show()


# In[65]:


'''±30% Fit Luminosity at 1s, 10s, 100s'''

per_lum_1,per_lum_10, per_lum_100 = [], [], []

for i in range(len(norm_amp_per)):
    per_lum_1.append(power_law1(1, norm_amp_per[i],norm_alpha_per[i])*1e50)
    per_lum_10.append(power_law1(10, norm_amp_per[i],norm_alpha_per[i])*1e50)
    per_lum_100.append(power_law1(100, norm_amp_per[i],norm_alpha_per[i])*1e50)
    

plt.figure(figsize = (8,5 ))
plt.title('Ts ±30% Luminosity Distribtuion at 1s')
plt.ylabel('Number of GRBs')
plt.xlabel('log(L) at 1s')
plt.hist(np.log10(per_lum_1))
plt.show()

plt.figure(figsize = (8,5 ))
plt.title('Ts ±30% Luminosity Distribution at 10s')
plt.ylabel('Number of GRBs')
plt.xlabel('log(L) at 10s')
plt.hist(np.log10(per_lum_10))
plt.show()

plt.figure(figsize = (8,5 ))
plt.title('Ts ±30% Luminosity Distribution at 100s')
plt.ylabel('Number of GRBs')
plt.xlabel('log(L) at 100s')
plt.hist(np.log10(per_lum_100))
plt.show()


# In[64]:


'''Peak Fit Luminosity at 1s, 10s, 100s'''

peak_lum_1,peak_lum_10, peak_lum_100 = [], [], []

for i in range(len(norm_amp_peak)):
    peak_lum_1.append(power_law1(1, norm_amp_peak[i],norm_alpha_peak[i])*1e50)
    peak_lum_10.append(power_law1(10, norm_amp_peak[i],norm_alpha_peak[i])*1e50)
    peak_lum_100.append(power_law1(100, norm_amp_peak[i],norm_alpha_peak[i])*1e50)
    

plt.figure(figsize = (8,5 ))
plt.title('Peak Fit Luminosity Distribution at 1s')
plt.ylabel('Number of GRBs')
plt.xlabel('log(L) at 1s')
plt.hist(np.log10(peak_lum_1))
plt.show()

plt.figure(figsize = (8,5 ))
plt.title('Peak Fit Luminosity Distribution at 10s')
plt.ylabel('Number of GRBs')
plt.xlabel('log(L) at 10s')
plt.hist(np.log10(peak_lum_10))
plt.show()

plt.figure(figsize = (8,5 ))
plt.title('Peak Fit Luminosity Distribution at 100s')
plt.ylabel('Number of GRBs')
plt.xlabel('log(L) at 100s')
plt.hist(np.log10(peak_lum_100))
plt.show()


# In[ ]:





# In[ ]:


'''TS Fit Luminosity at 1s, 10s, 100s'''

ts_lum_1,ts_lum_10, ts_lum_100 = [], [], []

for i in range(len(norm_amp_ts)):
    ts_lum_1.append(power_law1(1, norm_amp_ts[i],norm_alpha_ts[i])*1e50)
    ts_lum_10.append(power_law1(10, norm_amp_ts[i],norm_alpha_ts[i])*1e50)
    ts_lum_100.append(power_law1(100, norm_amp_ts[i],norm_alpha_ts[i])*1e50)
    

plt.figure(figsize = (8,5 ))
plt.title('Luminsoity Distribution at 1s')
plt.ylabel('Number of GRBs')
plt.xlabel('log(L) at 1s')
plt.hist(np.log10(ts_lum_1))
plt.show()

plt.figure(figsize = (8,5 ))
plt.title('Luminosity Distrubiton at 10s ')
plt.ylabel('Number of GRBs')
plt.xlabel('log(L) at 10s')
plt.hist(np.log10(ts_lum_10))
plt.show()

plt.figure(figsize = (8,5 ))
plt.title('Luminosity at 100s')
plt.ylabel('Number of GRBs')
plt.xlabel('log(L) at 100s')
plt.hist(np.log10(ts_lum_100))
plt.show()


# In[ ]:





# In[ ]:





# In[22]:


#np.savetxt('Multi-Fit-Params-Log.csv',all_parameterslog,fmt='%s', delimiter=',')
#np.savetxt('Multi-Fit-Params-Lineaar.csv',all_parameters,fmt='%s', delimiter=',')
#np.savetxt('Luminosities.csv',lum_params,fmt='%s', delimiter=',')


# In[72]:


'''Ts onwards'''
#plt.figure(figsize = (12,8 ))
#plt.scatter(all_T, all_lum, c=all_a_ts, s=50, cmap='Greys')
#plt.yscale('log')
#plt.xscale('log')
#plt.xlabel(xlabel='T-T0/(1+z) [s]',fontsize = 14) 
#plt.ylabel(ylabel='L [erg/s]',fontsize = 14) 
#plt.title(f'All GRBs Ts Luminosity Fit',fontsize = 14) 
#plt.tick_params(axis='both', which='major', labelsize=14)
#axs[0].errorbar(all_T, all_lum, yerr=all_lum_err, xerr=all_T_err, linestyle='',elinewidth=0.5)
#plt.axvline(Ts, color = 'red', linewidth=0.5)
#plt.colorbar()
#plt.show()

'''10-40s'''
plt.figure(figsize = (12,8 ))
plt.scatter(all_T, all_lum, c=all_a_U, s=50, cmap='coolwarm')
plt.yscale('log')
plt.xscale('log')
plt.xlabel(xlabel='T-T0/(1+z) [s]',fontsize = 14) 
plt.ylabel(ylabel='L [erg/s]',fontsize = 14) 
plt.title(f'Colour Map for All GRBs from {lower_bound}-{upper_bound}s Luminosity Fit',fontsize = 14) 
plt.tick_params(axis='both', which='major', labelsize=14)
#axs[0].errorbar(all_T, all_lum, yerr=all_lum_err, xerr=all_T_err, linestyle='',elinewidth=0.5)
plt.axvline(Ts, color = 'red', linewidth=0.5)
plt.colorbar()
plt.show()

'''±30%'''
plt.figure(figsize = (12,8 ))
plt.scatter(all_T, all_lum, c=all_a_per, s=50, cmap='coolwarm')
plt.yscale('log')
plt.xscale('log')
plt.xlabel(xlabel='T-T0/(1+z) [s]',fontsize = 14) 
plt.ylabel(ylabel='L [erg/s]',fontsize = 14) 
plt.title(f'Colour Map for All GRBs from ±30% Luminosity Fit',fontsize = 14) 
plt.tick_params(axis='both', which='major', labelsize=14)
#axs[0].errorbar(all_T, all_lum, yerr=all_lum_err, xerr=all_T_err, linestyle='',elinewidth=0.5)
plt.axvline(Ts, color = 'red', linewidth=0.5)
plt.colorbar()
plt.show()

'''Peak'''
plt.figure(figsize = (12,8 ))
#highlighting each point either side of ts
plt.scatter(all_T, all_lum, c=all_a_peak, s=50, cmap='coolwarm')
plt.yscale('log')
plt.xscale('log')
plt.xlabel(xlabel='T-T0/(1+z) [s]',fontsize = 14) 
plt.ylabel(ylabel='L [erg/s]',fontsize = 14) 
plt.title(f'Colour Map for All GRBs From Peak Luminosity Fit',fontsize = 14) 
plt.tick_params(axis='both', which='major', labelsize=14)
#axs[0].errorbar(all_T, all_lum, yerr=all_lum_err, xerr=all_T_err, linestyle='',elinewidth=0.5)
plt.axvline(Ts, color = 'red', linewidth=0.5) 
plt.colorbar()
plt.show()


# In[98]:


'''Ts onwardss'''
#plt.figure(figsize = (12,8 ))
#plt.scatter(all_T, all_lum, c=brightness_ts, s=50, cmap='Blues')
#plt.yscale('log')
#plt.xscale('log')
#plt.xlabel(xlabel='T-T0/(1+z) [s]',fontsize = 14) 
##plt.ylabel(ylabel='L [erg/s]',fontsize = 14) 
#plt.title(f'All GRBs Ts Luminosity Fit',fontsize = 14) 
#plt.tick_params(axis='both', which='major', labelsize=14)
#axs[0].errorbar(all_T, all_lum, yerr=all_lum_err, xerr=all_T_err, linestyle='',elinewidth=0.5)
#plt.axvline(Ts, color = 'red', linewidth=0.5)
#plt.colorbar()
#plt.show()

# '''10-40s'''
# plt.figure(figsize = (12,8 ))
# plt.scatter(all_T, all_lum, c=brightness_U, s=35, cmap='coolwarm')
# plt.yscale('log')
# plt.xscale('log')
# plt.xlabel(xlabel='T-T0/(1+z) [s]',fontsize = 14) 
# plt.ylabel(ylabel='L [erg/s]',fontsize = 14) 
# plt.title(f'Colour Map for All GRBs from {lower_bound}-{upper_bound}s Luminosity Fit',fontsize = 14) 
# plt.tick_params(axis='both', which='major', labelsize=14)
# #axs[0].errorbar(all_T, all_lum, yerr=all_lum_err, xerr=all_T_err, linestyle='',elinewidth=0.5)
# plt.axvline(Ts, color = 'red', linewidth=0.5)
# plt.colorbar()
# plt.show()

'''±30%'''
plt.figure(figsize = (12,8 ))
plt.scatter(all_T, all_lum, c=brightness_per, s=35, cmap='rainbow')
plt.yscale('log')
plt.xscale('log')
plt.xlabel(xlabel='T-T0/(1+z) [s]',fontsize = 14) 
plt.ylabel(ylabel='L [erg/s]',fontsize = 14) 
plt.title(f'Colour Map for All GRBs from ±30% Luminosity Fit',fontsize = 14) 
plt.tick_params(axis='both', which='major', labelsize=14)
#axs[0].errorbar(all_T, all_lum, yerr=all_lum_err, xerr=all_T_err, linestyle='',elinewidth=0.5)
plt.axvline(Ts, color = 'red', linewidth=0.5)
plt.colorbar()
plt.show()

'''Peak'''
plt.figure(figsize = (12,8 ))
#highlighting each point either side of ts
plt.scatter(all_T, all_lum, c=brightness_peak, s=35, cmap='rainbow')
plt.yscale('log')
plt.xscale('log')
plt.xlabel(xlabel='T-T0/(1+z) [s]',fontsize = 14) 
plt.ylabel(ylabel='L [erg/s]',fontsize = 14) 
plt.title(f'Colour Map for All GRBs From Peak Luminosity Fit',fontsize = 14) 
plt.tick_params(axis='both', which='major', labelsize=14)
#axs[0].errorbar(all_T, all_lum, yerr=all_lum_err, xerr=all_T_err, linestyle='',elinewidth=0.5)
plt.axvline(Ts, color = 'red', linewidth=0.5) 
plt.colorbar()
plt.show()

#print(len(brightness_peak))


# In[ ]:





# In[ ]:


fig, axs = plt.subplots(4,figsize=(12,22), sharex=False, sharey=True)
#fig.suptitle(f'Lum vs Alpha Ts={Ts}s', fontsize=19)

print(len(norm_alpha_per))
for i in range(len(grb_name)):
    axs[0].scatter(norm_alpha_ts[i], norm_lum_ts[i])
    axs[1].scatter(norm_alpha_U[i], norm_lum_ts_U[i])#,label=grb_name[i])
    axs[2].scatter(norm_alpha_per[i], norm_lum_ts_per[i])#,label=grb_name[i])
    axs[3].scatter(norm_alpha_peak[i], norm_lum_ts_peak[i])#,label=grb_name[i])


#print(np.log10(lum_ts_err)
axs[0].set_xlabel('Alpha',fontsize = 14) 
axs[1].set_xlabel('Alpha',fontsize = 14) 
axs[2].set_xlabel('Alpha',fontsize = 14)
axs[3].set_xlabel('Alpha',fontsize = 14) 
#ax3.set_xlabel('Alpha',fontsize = 14) 
axs[0].set_ylabel('L_Ts=10 [erg/s]',fontsize = 14) 
#ax1.set_title('Two-Point Lum')
axs[0].set_title('Ts Lum')
axs[1].set_title('5-50s Lum')
axs[2].set_title('±30% Ts Lum')
axs[3].set_title('Peak Ts Lum')
#plt.title(f'Luminosity at Ts={Ts} vs Decay Index From Normalised ± 30% Fits',fontsize = 14) 
#plt.tick_params(axis='both', which='major', labelsize=18)


# In[ ]:





# In[ ]:





# In[14]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


no_z_names = []

for i in table:
    if i['LUMINOSITY_DISTANCE'] == 0:
            no_z_names.append(i['GCNNAME'])
    
np.savetxt('Grb-names-no-z.txt', no_z_names, delimiter=',', fmt='%s')

print(len(no_z_names))

