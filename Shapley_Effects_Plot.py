# -*- coding: utf-8 -*-
"""

Plot Shapley Effects

@author: Charline Dalinghaus
"""

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

os.chdir("C:/")
maindir= "C:/"
path = 'TWL_MCsimul_'

#%% 

# Colormap and variable order
variables_fut = ['MSL', 'AT', 'SS', 'Runup','Setup Eq', 'Swash Eq', 'βf', 'Scenario', 'Model']

variable_order = ['Model', 'βf', 'Swash Eq', 'Setup Eq', 'Runup', 'SS', 'AT', 'MSL', 'Scenario']
variable_colors = {
    'Model': '#FFD089',
    'βf': '#FEA16B',
    'Swash Eq': '#FF4343',
    'Setup Eq': '#E98DBB',
    'Runup': '#7A4C93',
    'SS': '#003795',
    'AT': '#0B99C8',
    'MSL': '#72F2E9',
    'Scenario': '#4FC496',
}

#%% 
########################### Future Data 2026 - 2045 ###########################

# Load data for near future analysis
datadir_fut1 = os.path.join(maindir, path + 'Fut1')
common_name = "shapley_effects_"
all_files_fut1 = os.listdir(datadir_fut1)
csv_files_fut1 = [file for file in all_files_fut1 if common_name in file and file.endswith('.csv')]

shapley_effects_fut1 = {}
for file in csv_files_fut1:
    file_path = os.path.join(datadir_fut1, file)
    file_number = file.replace(common_name, '').split('.')[0]
    df = pd.read_csv(file_path)
    shapley_effects_fut1[file_number] = df

# Extract and process data
sumYEAR = []
data_fut1 = {
    'MSL': [],
    'AT': [],
    'SS': [],
    'Runup': [],
    'Setup Eq': [],
    'Swash Eq': [],
    'βf': [],
    'Scenario': [],
    'Model': []
}

for key in shapley_effects_fut1:
    for i, variable in enumerate(variables_fut):
        value = shapley_effects_fut1[key].iloc[i, 0]
        data_fut1[variable].append(value)
    sumYEAR.append(int(key))

data_fut1 = pd.DataFrame(data_fut1)
data_fut1.insert(0, 'Year', sumYEAR)


# PLOT!!
#Stacked bar chart with moving averages over a 3-year window for each variable
window_size = 3
data_smoothed = data_fut1.rolling(window=window_size, min_periods=1, center=True).mean()
data_smoothed['Year']=data_fut1['Year']

plt.figure(figsize=(10, 6))
for variable in variable_order:
    if variable in data_fut1:
        if variable == 'Model':
            bottom = np.zeros(len(sumYEAR))
        plt.bar(data_smoothed['Year'], data_smoothed[variable], label=variable, color=variable_colors[variable], bottom=bottom)
        bottom += data_smoothed[variable]

handles, labels = plt.gca().get_legend_handles_labels()
plt.xticks(np.arange(min(data_smoothed['Year']), max(data_smoothed['Year'])+1, 2))
plt.legend(reversed(handles), reversed(labels), loc='upper left')
plt.title('Stacked Bar Chart of Shapley Effects with 3 Moving Average (Future Data 2026-2045)')
plt.xlabel('Date', fontweight='bold', fontsize = 16)
plt.ylabel('Shapley Effects', fontweight='bold', fontsize = 16)
plt.savefig('Shapley/' + path + 'Fut1.jpg')         