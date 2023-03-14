# MIT_440_Project
440 Course Project
READ ME file

Overview
This repository contains the data, analysis code, and figures for our project. The data is bulk RNA seq data published from Acosta-Rodriguez et al., (doi: 10.1126/science.abk0297) which has already been normalized to reads per kilobase per million (RPKM). The analysis includes Python (and eventually R) scripts for analyzing the published data. 

Data
The data provided is from doi: 10.1126/science.abk0297. In this paper, authors sought to break down the independent effects of calorie restriction, time restricted feeding, and circadian feeding. They used a mouse model with automatic feeding schedules, and broke the mice up into 6 groups:
1.	Ad libitum: mice had unlimited access to food
2.	Calorie restriction, spread: mice were provided 30% fewer calories than the ad libitum group, but were able to eat across a 24 hour period.
3.	Calorie restriction, 12 hours, daytime: same calorie restriction, but mice only had access to food for 12 hours during the day
4.	Calorie restriction, 2 hours, daytime: same calorie restriction, but mice only had access to food for 2 hours during the day
5.	Calorie restriction, 12 hours, nighttime: same calorie restriction, but mice only had access to food for 12 hours during the night
6.	Calorie restriction, 2 hours, nighttime: same calorie restriction, but mice only had access to food for 12 hours during the night
The authors harvested liver tissue from 2 mice in each group at young and old ages, and when livers were harvested, they were harvested from mice every 4 hours over a 48 hour period. Bulk RNA seq was performed to assess differences in gene expression across ages, time, and feeding groups.

Folder Structure
