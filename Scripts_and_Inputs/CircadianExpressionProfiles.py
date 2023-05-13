# Script to generate plots of gene expression profiles of any gene

# Import  libraries
import pandas as pd
import matplotlib.pyplot as plt

# Upload the RPKM-normalized gene counts
data = pd.read_excel('C:/Users/dmsel/440_HW/440_Project/Data/RNAseq_Dataframe_Input.xlsx') #Change this file path to match where you saved your data file


#%% Process dataframe to isolate data that corresponds to a user-selected gene and calculate the average across the two replicates at each timepoint

# Isolate data that corresponds to gene of interest
gene = data.loc[data['GeneID'] == 'ENSMUSG00000027533.10_Fabp5']
gene = gene.reset_index(drop=True)
gene_name = gene.iloc[0,0].split('_')[1]
    
#Set up new dataframe that we can process further
gene_avg = gene.copy()
    
#Generate list of timepoints for the dataframe
column_names = list(range(0,48,4))
    
#Loop through the 12 timepoints, averaging together the two replicates from the dataset
for i in range(12):
    averages = (gene_avg.iloc[:,2*i+2]+gene_avg.iloc[:,2*i+3])/2 #Find average of the two replicates
    averages = pd.DataFrame(averages,columns = [column_names[i]]) #Convert averages to a dataframe with the timepoint as the column header
    gene_avg = pd.concat([gene_avg,averages], axis=1) #Append the new column to the existing dataframe
    
#Delete the columns that had the replicate data, leaving us with only the average data
gene_avg.drop(gene_avg.columns[[list(range(2,26))]], axis=1, inplace=True)

plt.plot(gene_avg.iloc[7,2:]) #Plot for day mice
plt.plot(gene_avg.iloc[9,2:]) #Plot for night mice

#Set figure labels
plt.xlabel('Hours')
plt.ylabel('mRNA level (RPKM)')
plt.legend(['Day 12h','Night 12h'])
plt.title(gene_name)
plt.show()