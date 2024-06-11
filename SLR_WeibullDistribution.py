# -*- coding: utf-8 -*-
"""
Weibull Distribution implementation to generate probability distributions from the decadal SLR scenarios downloaded from the New Zealand SeaRise Project (https://www.searise.nz/)
Linear interpolation to generate a comprehensive projected uniform SLR distribution with 5000 values per year from 2020 to 2100.

@author: Charline Dalinghaus
"""

import os
import pickle

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from scipy.stats import weibull_min

os.chdir("C:/")

#%%
###############################################################################
######################## Part A: Weibull distribution #########################
###############################################################################

RCPs=pd.read_csv('file.csv') #File downloaded from the New Zealand SeaRise Project

SLRdict = {
    'RCP4.5+VLM': RCPs[(RCPs['scenario'] == 'SSP2-4.5 + VLM (medium confidence)') & (RCPs['year'] > 2005) & (RCPs['year'] <= 2100)].reset_index(drop=True),
    'RCP8.5+VLM': RCPs[(RCPs['scenario'] == 'SSP5-8.5 + VLM (medium confidence)') & (RCPs['year'] > 2005) & (RCPs['year'] <= 2100)].reset_index(drop=True) 
} #Dictionary for the SLR data

# Create array k
k_values = np.linspace(0.001, 5, num=10000) #Must be >0

final_values = []
colors = iter(cm.seismic(np.linspace(0, 1, 18)))

plt.figure (figsize=(10, 6))

for scenario, data in SLRdict.items():
    for _, row in data.iterrows():
        year = row['year']
        p17 = row['p17']
        p50 = row['p50']
        p83 = row['p83']
        
        # Calculate Least Squares to determine the parameters of the Weibull distribution
        E = (0.17 - 1+np.exp(-(p17/(p50/((np.log(2))**(1/k_values))))**k_values))**2 + \
            (0.83 - 1+np.exp(-(p83/(p50/((np.log(2))**(1/k_values))))**k_values))**2
        
        min_idx = np.nanargmin(E)
        shape_param = k_values[min_idx]
        scale_param = p50 / (np.log(2) ** (1 / shape_param))

        p17_check = 1 - np.exp(-(p17 / scale_param) ** shape_param)
        p50_check = 1 - np.exp(-(p50 / scale_param) ** shape_param)
        p83_check = 1 - np.exp(-(p83 / scale_param) ** shape_param)
        
        final_values.append([scenario, year, scale_param, shape_param, p17_check, p50_check, p83_check])
        
        # Generate the Weibull distribution function
        dist = weibull_min(shape_param, scale=scale_param)
        values = dist.rvs(size=5000)
        
        x = np.linspace(values.min(), values.max(), 100)
        y = dist.pdf(x)
     
        c = next(colors)  #Generate colors for each iteration
        plt.plot(x, y, '-', c=c, lw=2, label=f'{scenario} {year}')
        
plt.xlim(0, 2.1); plt.ylim(0, 12)
plt.xlabel('Relative MSL (m)', fontsize = 16, fontweight = 'bold'); plt.ylabel('Probability density function\n' + '(Weibull Distribution)', fontsize = 16, fontweight = 'bold')
plt.legend(loc='best', ncol=2, prop={'size': 14,'weight':'bold'})

plt.savefig("SLR_WeibullDist.png", dpi=300)

w_values = pd.DataFrame(final_values, columns =['RCP', 'year', 'l', 'k', 'p17', 'p50', 'p83'])
w_values.to_csv('weibull_values.csv')

#%%
###############################################################################
######################## Part B: Linear Interpolation #########################
###############################################################################

## Interpolating the Weibull distribution parameters (k and l)
WD = {'RCP4.5+VLM': {},
      'RCP8.5+VLM': {}} #Dictionary for the Weibull Distribution 

colors = iter(cm.seismic(np.linspace(0, 1, 2)))

for scenario in SLRdict:
    weibull_k = w_values.loc[w_values['RCP'] == scenario, 'k']
    weibull_l = w_values.loc[w_values['RCP'] == scenario, 'l']
    t = SLRdict[scenario]['year']
    
    t_interp = np.arange(2020, 2101) #Set of years at which to interpolate
    weibull_interp_k= np.interp(t_interp, t, weibull_k)  #Interpolate k and l values yearly
    weibull_interp_l= np.interp(t_interp, t, weibull_l)
        
    # Generate Weibull distribution for each year
    fig, ax = plt.subplots(figsize=(10, 6))
    c = next(colors)
    for k, l, t_int in zip(weibull_interp_k, weibull_interp_l, t_interp):
            W_dist = weibull_min(k, scale=l) #Generate the Weibull distribution function
            W_values = W_dist.rvs(size=5000).tolist()
            WD[scenario][t_int.tolist()] = W_values
            ax.plot(np.full((5000,), t_int), W_values, '.', color=c) 

    ax.set_xlabel('Date', fontsize = 16, fontweight = 'bold'); ax.set_ylabel('Relative MSL (m)', fontsize = 16, fontweight = 'bold')
    ax.set_xlim(2019, 2101); ax.set_ylim(0, 2.1)
    ax.legend([scenario], loc='upper left', prop={'size': 14, 'weight': 'bold'})

    plt.savefig("SLR_WeibullDist"+str(scenario)+".png", dpi=300)

pickle.dump(WD, open("weibull_distribution", "wb"))

