# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 17:57:32 2023

@author: dmsel
"""

# Import 

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
data = pd.read_excel('C:/Users/dmsel/440_HW/440_Project/RNAseq_Dataframe_Input.xlsx')

#%% Analysis starts here

t0_data = data[data.columns[0:3]]


AL_young = t0_data[t0_data['Group'] == '06mo_AL']

AL_old = t0_data[t0_data['Group'] == '19mo_AL']
AL_old = AL_old.reset_index(drop=True)

#AL_comparison = AL_young[data.columns[0]]
AL_comparison = AL_young.copy()
AL_comparison = AL_comparison.drop(columns=['Group','00_A'])
AL_comparison['fold change'] = AL_old['00_A']/AL_young['00_A']
AL_comparison.replace(np.inf,np.nan,inplace=True)
AL_comparison.dropna(inplace=True)


#%%
from bioinfokit import analys, visuz
# load dataset as pandas dataframe

df = analys.get_data('volcano').data

visuz.GeneExpression.volcano(df=AL_comparison, lfc='fold change', show=True)


#%% Plotting circadian rhythms of BMAL1 gene for all feeding groups

# Isolate data that corresponds to BMAL1 gene
BMAL1_data = data.loc[data['GeneID'] == 'ENSMUSG00000055116.8_Arntl']
BMAL1_data = BMAL1_data.reset_index(drop=True)

#Set up new dataframe that we can process further
BMAL1_Avg = BMAL1_data.copy()

#Generate list of timepoints for the dataframe
column_names = list(range(0,48,4))

#Loop through the 12 timepoints, averaging together the two replicates from the dataset
for i in range(12):
    averages = (BMAL1_Avg.iloc[:,2*i+2]+BMAL1_Avg.iloc[:,2*i+3])/2 #Find average of the two replicates
    averages = pd.DataFrame(averages,columns = [column_names[i]]) #Convert averages to a dataframe with the timepoint as the column header
    BMAL1_Avg = pd.concat([BMAL1_Avg,averages], axis=1) #Append the new column to the existing dataframe

#Delete the columns that had the replicate data, leaving us with only the average data
BMAL1_Avg.drop(BMAL1_Avg.columns[[list(range(2,26))]], axis=1, inplace=True)


#%% Plot dataframe

#List of subplot titles that correspond to the feeding group
titles = ['AL','CR.day.12h','CR.day.2h','CR.night.12h','CR.night.2h','CR.spread']

#Set up 3x2 subplot
fig, axs = plt.subplots(3,2,figsize=(6, 8),sharex=True,sharey=True)
count = 0

#Loop through the dataframe, separating each feeding group and combining the two ages on each plot
for i in range(3):
    for k in range(2):
        axs[i,k].plot(BMAL1_Avg.iloc[count,2:]) #Plot for young mice
        axs[i,k].plot(BMAL1_Avg.iloc[count+6,2:]) #Plot for old mice
        axs[i,k].set_title(titles[count]) #Set the title to match the feeding group
        count +=1

#Set figure labels
axs[2,0].set_xlabel('Hours')
axs[2,1].set_xlabel('Hours')
fig.supylabel('mRNA level (RPKM)')
fig.legend(['Young mice','Old mice'],loc='lower center',ncol=2)
plt.show()

