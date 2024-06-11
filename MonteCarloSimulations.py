# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 16:29:57 2023

Monte Carlo simulations for Near Future Data

@author: Charline Dalinghaus
"""

import os
import pickle
import random
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from tqdm.contrib.concurrent import process_map #parallelize using all CPUs and show a smart progress meter

os.chdir("C:/")

#%%
## ============================== Read TWL data ===============================

# Set the data directory
data_dir = Path.cwd()

# Read TWL data
filepathSLR = data_dir / "Data/weibull_distribution"
slr = pd.read_pickle(open(filepathSLR, "rb"))

filepathAT = data_dir / "Data/at_dailyAmp"
at = pd.read_pickle(open(filepathAT, "rb"))

filepathSS = data_dir / "Data/ss"
ss = pd.read_pickle(open(filepathSS, "rb"))

filepathWAVES = data_dir / "Data/waves"
waves = pd.read_pickle(open(filepathWAVES, "rb"))

filepathFSLOPE = data_dir / "Data/GeneralForeshoreSlope.txt"
foreshore_slope = np.loadtxt(filepathFSLOPE)

#%%
## =================== Monte Carlo - Future Data 2026 - 2045 ==================

del waves['Hist']

start = '2026-01-01'; end = '2045-12-31' #Data in common between all datasets
date_range = pd.date_range(start=start, end=end, freq='Y')
years = pd.Series(date_range.year)

mask_at = (at['time'] >= start) & (at['time'] < '2046-01-01') #Create a mask for the in-common time data between all datasets
dates_at = at.loc[mask_at] #Apply the mask

# Iterate over each day in the dataframe
def mc_fut1(year):
    twl_100000samples_fut1 = pd.DataFrame()
    for i in range(100000):
        
        # Select a random scenario and model
        scenario = random.choice(list(waves.keys())) # Get a random KEY (SCENARIO) from the outer dictionary
        model = random.choice(list(waves[scenario].keys())) # Get a random KEY (MODEL) from the outer dictionary
        
        # Select the WAVE DATA referred to the selected model
        mask_w = (waves[scenario][model]['time'] >= start) & (waves[scenario][model]['time'] < '2046-01-01') #Create a mask for the in-common time data between all datasets
        wave_data = waves[scenario][model].loc[mask_w]
        
        # Pick a random wave data from the same year of the "year" variable
        wave_data_same_year = wave_data[wave_data['time'].dt.year == year] #Filter the "wave_data" to get all data for the same month
        random_day = wave_data_same_year.sample(n=1)['time'].values[0] #Randomly select a day from this year  
        wave_day_random = wave_data_same_year[wave_data_same_year['time'] == random_day] #Get the corresponding wave data using the selected day
      
        # Select the SS DATA referred to the same wave random day
        ss_data = ss[scenario][model]
        random_day_date = str(random_day)[:10]  # Convert the numpy datetime64 to a string and extract the date part
        ss_day_random = ss_data[ss_data['time'].astype(str).str[:10] == random_day_date]
        
        # Select the AT DATA referred to the same day
        at_day_random = dates_at.loc[dates_at['time'].dt.date.isin(wave_day_random['time'].dt.date)]
            
        # Randomly select an SLR referred to the same year
        rand_slr = random.choice(slr[scenario+'+VLM'][wave_day_random['time'].dt.year.iloc[0]])
        
        # Calculate Runup using the selected Wave data
        # Variables
        H0 = wave_day_random['Hs'].iloc[0]
        Tp = wave_day_random['Tp'].iloc[0]
        L0 = 1.56*(Tp**2)
        βf = random.choice(foreshore_slope)
        ξ0 = βf/((H0/L0)**0.5)
        D50 = 0.16/1000 #In meters. D50 from Reinen-Hamill et al., 2006 & Schofield, 1970
    
        # Choose Wave Setup equation
        equationsWS = ['0.355*H0*(ξ0**0.5)', #Dalinghaus et al. (2022)
                      'H0/4.08*(ξ0/3.25 + ξ0/(ξ0 + 0.64) + ξ0/(1625*D50 + ξ0))', #Dalinghaus et al. (2022)
                      '0.220*(H0**0.629)*(L0**0.371)*(βf**0.538)'] #Ji et al. (2018)
        wave_setupEQ = random.choice(equationsWS)
        wave_setup = eval(wave_setupEQ)
    
        # Choose Swash equation #Total Swash=((Sig2 + Sinc2)**0.5)/2 #Stockdon et al. (2006)
        equationsS = ['((((0.19+0.008*Omega)*((H0*L0*βf)**0.5))**2+(Sinc)**2)**0.5)', #Gomes da Silva et al. (2018, 2019)
                      '146.737*(βf**2)+((Tp*(H0**3))/(5.8+10.595*(H0**3)))-4397.838*(βf**4)', #Passarella et al. (2017)
                      '(((0.06*((H0*L0)**0.5))**2)+((0.75*βf*((H0*L0)**0.5))**2))**0.5'] #Stockdon et al. (2006)
        swashEQ = random.choice(equationsS)
             
        # Parameters formula Gomes da Silva et al. (2019):
        if swashEQ == '((((0.19+0.008*Omega)*((H0*L0*βf)**0.5))**2+(Sinc)**2)**0.5)':
            Ws = 273*D50**1.1 #0.1 < D50 < 1 mm
            Omega = H0/(Ws*Tp)
            Sinc=[]
            if Omega >= 5.5:
                Sinc = ((2.83*(βf**2.12))*((H0/L0)**-0.82))*H0 #dissipative
            else:
                if Omega <= 1.5:
                    Sinc = ((0.50*(βf**-0.37))*((H0/L0)**-0.15))*H0 #reflective
                else:
                    Sinc = ((0.15*(βf**0.56))*((H0/L0)**-0.64))*H0 #intermediate
        
        swash = eval(swashEQ)
    
        # Runup
        ru2 = 1.1*(wave_setup + (swash/2)) #Stockdon et al. (2006)
            
        # Calculate Fut1 TWL
        twl = rand_slr + at_day_random.iloc[0,1] + ss_day_random.iloc[0,1] + ru2
        twl_sample_fut1 = pd.DataFrame({ 'Date': [wave_day_random.iloc[0,0]],
                                         'TWL': [twl],
                                         'SLR': [rand_slr],
                                         'AT': [at_day_random.iloc[0,1]],
                                         'SS': [ss_day_random.iloc[0,1]],
                                         'Runup': [ru2],
                                         'Setup Eq.': [wave_setupEQ],
                                         'Swash Eq.': [swashEQ],
                                         'Hs0': [H0],
                                         'Tp': [Tp],
                                         'βf': [βf],
                                         'Scenario': [scenario],
                                         'Model': [model]})
    
        twl_100000samples_fut1 = pd.concat([twl_100000samples_fut1, twl_sample_fut1,])
    return twl_100000samples_fut1

if __name__ == '__main__':
    twl_samples_fut1 = pd.concat(process_map(mc_fut1, years, chunksize = 1, max_workers = 5))
    
    pickle.dump(twl_samples_fut1, open("MonteCarlo/TWL_MCsimul_Fut1", "wb"))

