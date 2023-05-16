# MIT_440_Project

### Overview

This repository contains the data, code for analysis, and figures for my MIT course 20.440 project. The data is bulk RNA seq data published from Acosta-Rodriguez et al., (doi: 10.1126/science.abk0297). Both the raw and RPKM normalized counts were posted in the Scripts_and_Inputs folder. The analysis includes Python (and eventually R) scripts for analyzing the published data. The data and code here will be used to look further into how eating time affects gene expression in old and young mice. 

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

The "Scripts_and_Inputs" folder contains all of the raw data files. The data files provided are the raw counts (starts with GE_RawData), the normalized counts (normalized to RPKM), and the column headers (starts with ColData) for DESeq2. This folder also contains the scripts for running differential expression analysis (DiffExp_Analysis_R.R) and plotting the longitudinal gene expression profiles (CircadianExpressionProfiles.py). When running the differential gene expression analysis, be sure to pair the ColData to the correct raw counts file.

The "Figures" folder contains the figures generated for this project. This includes the PCA figure, all of the volcano plots for differential expression, the gene ontology figure, and the gene expression profiles.

The "Results" folder contains the differential expression results from each pairwise comparison run in DiffExp_Analysis_R.R. These all show the genes and their corresponding log2 fold change and p values. There is also a table that summarizes the results of all comparisons neatly in one file.

The "GSEA_inputs" folder contains the two input files needed to run GSEA. These two files are the raw counts for the night vs. day 12 hour feeding groups as well as the phenotype labels. 
