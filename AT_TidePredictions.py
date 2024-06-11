# -*- coding: utf-8 -*-
"""
Extraction of tidal constituents, reconstruction and prediction of Tides using the UTIDE implementation for Python (https://github.com/wesleybowman/UTide)

@author: Charline Dalinghaus
"""

import osxx
import pickle

import matplotlib.pyplot as plt
import pandas as pd
from utide import solve, reconstruct

os.chdir("C:/")

#%%
## ======================== Extract Tidal Constituents ========================

# Read AT data
at = pickle.load(open('at_NIWA', "rb"))

tide = at['tide']
time_present = at['time']
lat = -36.584896  # Latitude of the location where the data was collected

coef = solve(
    time_present,
    tide,
    lat=lat,
    method="ols",
    conf_int="MC",
    verbose=False,
)

# print(coef.keys())
# print(coef['name'])

constituents = pd.DataFrame({
    'Constituents': coef['name'],
    'Amplitude': coef['A'],
    'Amplitude CI': coef['A_ci'],
    'Phase': coef['g'],
    'Phase CI': coef['g_ci'],
    '% Energy': coef['PE'],
    'Signal-to-Noise power Ratio': coef['SNR']
})

#Plot Constituents
plt.figure(figsize=(15, 5))
plt.scatter(constituents['Constituents'], constituents['Amplitude'])
plt.xlabel('Tidal Constituent')
plt.ylabel('Amplitude')
plt.title('Amplitudes of Harmonic Constituents')
plt.xticks(rotation=90)
plt.ylim(0, 1.25)
plt.show()

constituents.to_csv('TidalConstituents.txt', index=False, sep='\t')

#%%
## ============================= Reconstruct Tide =============================

#Reconstruct (check if the reconstruction match with the data)
tide_rec = reconstruct(time_present, coef, verbose=False)
print(tide_rec.keys())

#Plot Tide
fig, (ax0, ax1, ax2) = plt.subplots(figsize=(17, 5), nrows=3, sharey=True, sharex=True)

ax0.plot(time_present, tide, label="NIWA Tide", color="C0")
ax1.plot(time_present, tide_rec.h, label="Reconstruction", color="C1")
ax2.plot(time_present, tide - tide_rec.h, label="Residual", color="C2")
fig.legend(ncol=3, loc="upper center");

#%%
## =============================== Predict Tide ===============================

#Predict
time_future = pd.date_range(start='2030-12-31 11:10:00', end='2100-12-31 23:50:00', freq='10min')
time_future = pd.Series(time_future)

tide_pred = reconstruct(time_future, coef, verbose=False)
print(tide_pred.keys())

#Plot
time_present = time_present.dt.tz_convert(None)
tide_present = pd.DataFrame({'time': time_present, 'tide': tide})
tide_future = pd.DataFrame({'time': time_future, 'tide': tide_pred['h']})
at_joined = pd.concat([tide_present, tide_future], axis=0).reset_index(drop=True)

time_joined = at_joined['time']
tide_joined = at_joined['tide']

fig, (ax0, ax1, ax2) = plt.subplots(figsize=(17, 5), nrows=3, sharey=True, sharex=True)

ax0.plot(time_present, tide, label="Reconstruction", color="C0")
ax1.plot(time_future, tide_pred.h, label="Prediction", color="C1")
ax2.plot(time_joined, tide_joined, label="Complete Tide Series", color="C2")
fig.legend(ncol=3, loc="upper center")

pickle.dump(at_joined, open("at", "wb"))
plt.savefig('.../AT.png', dpi=300)  
