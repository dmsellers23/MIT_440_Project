# MIT_440_Project

### Overview

This repository contains the data, code for analysis, and figures for my MIT course 20.440 project. The data is bulk RNA seq data published from Acosta-Rodriguez et al., (doi: 10.1126/science.abk0297) which has already been normalized to reads per kilobase per million (RPKM). The analysis includes Python (and eventually R) scripts for analyzing the published data. The data and code here will be used to look further into how eating time affects gene expression in old and young mice. 

### Data

The data provided is from doi: 10.1126/science.abk0297. In this paper, authors sought to break down the independent effects of calorie restriction, time restricted feeding, and circadian feeding. They used a mouse model with automatic feeding schedules, and broke the mice up into 6 groups:
1.	Ad libitum: mice had unlimited access to food
2.	Calorie restriction, spread: mice were provided 30% fewer calories than the ad libitum group, but were able to eat across a 24 hour period.
3.	Calorie restriction, 12 hours, daytime: same calorie restriction, but mice only had access to food for 12 hours during the day
4.	Calorie restriction, 2 hours, daytime: same calorie restriction, but mice only had access to food for 2 hours during the day
5.	Calorie restriction, 12 hours, nighttime: same calorie restriction, but mice only had access to food for 12 hours during the night
6.	Calorie restriction, 2 hours, nighttime: same calorie restriction, but mice only had access to food for 12 hours during the night
The authors harvested liver tissue from 2 mice in each group at young and old ages, and when livers were harvested, they were harvested from mice every 4 hours over a 48 hour period. Bulk RNA seq was performed to assess differences in gene expression across ages, time, and feeding groups.

### Folder Structure

The "Data" folder contains all of the raw data files. The data files provided are the raw counts, the normalized counts (normalized to RPKM), and the input data file for running the code. Use the RNAseq_Dataframe_Input file to run the code.

The "Code" folder contains the code for analyzing the data. 

The "Figures" folder contains the figures generated for this project. So far, only one figure has been generated which comes from the Project_analysis.py file in the  Code folder.


### Installation
In order to run the code:
1. Download the .py file
2. Make sure that you have pandas, numpy, and matplotlib packages installed
3. Download the input data file (RNAseq_Dataframe_Input) and make sure the file path in the code matches where the input data file is saved
4. Run the code!
