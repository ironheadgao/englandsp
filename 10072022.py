import math
import pandas as pd
import numpy as np
from sympy import Symbol,cos,integrate
import itertools as it
from pandas import DataFrame
import datetime
from datetime import timedelta
import glob
from OSGridConverter import grid2latlong
from OSGridConverter import latlong2grid
import utm
names = locals()

region =['ANGLIAN']
#creat data_to
for i in region:#for every region
    data_to = pd.read_csv("\\04062022_Spear_"+i.title() +"_taxon_data_final_SPECIES.csv",encoding='gbk')
    sumoftotal = pd.read_csv("\\" + i.title() + "_taxon_data_final_SPECIES.csv",encoding='gbk')
    if i == 'ANGLIAN':
        data_to = data_to.rename(columns={'SumOfTOTAL_ABUNDANCE': 'Rank','Rank':'Sensitivity','Sensitivity':'Generation','Generation':'Refuge','Refuge':'Exposed','Exposed':'NONE'})
        data_to.insert(7,'SumOfTOTAL_ABUNDANCE',sumoftotal['SumOfTOTAL_ABUNDANCE'])
        data_to = data_to.drop('NONE', axis=1)
    elif i =='MIDLANDS':
        data_to= data_to.iloc[:,1:]
        data_to = data_to.drop({ 'Sensitivity', 'Generation', 'Refuge', 'Exposed'},axis=1)
        data_to = data_to.rename(
            columns={ 'TOTAL_ABUNDANCE':'Sensitivity', 'REGION':'Generation', 'SumOfTOTAL_ABUNDANCE':'Refuge', 'Rank':'Exposed'})
        data_to.insert(7, 'SumOfTOTAL_ABUNDANCE', sumoftotal['SumOfTOTAL_ABUNDANCE'])
        data_to.insert(8, 'SumOfTOTAL_ABUNDANCE_Copy', sumoftotal['SumOfTOTAL_ABUNDANCE'])
    else:
        data_to.insert(6, 'SumOfTOTAL_ABUNDANCE_Copy', sumoftotal['SumOfTOTAL_ABUNDANCE'])


    data_spear = data_to.drop_duplicates(subset=['ANALYSIS_ID'])
    data_spear['SPEAR'] = 0
    data_to['AtRisk'] = 0
    data_to['LOG_up'] = 0
    data_to['LOG_down'] = 0
    # the loop need to be optimized
    for j in range(0,data_to.shape[0]):
        data_to.iloc[j,-3] = 1 if (data_to.iloc[j,9] > -0.36) and (float(data_to.iloc[j,10]) > 0.5) else 0
        data_to.iloc[j, -2] = np.log10(data_to.iloc[j,7]+1) * data_to.iloc[j,-3]
        data_to.iloc[j, -1] = np.log10(data_to.iloc[j, 7] + 1)
        print(j,i,'AtRisk')
    # generate final spear dataframe
    for m in  range(0,data_spear.shape[0]):
        if m < data_spear.shape[0]-1:
            n=m+1
            data_spear.iloc[m, -1] = np.sum(data_to.iloc[data_spear.iloc[m, 0]:data_spear.iloc[n, 0], -2])/np.sum(data_to.iloc[data_spear.iloc[m, 0]:data_spear.iloc[n, 0], -1])
        else :
            n = data_to.shape[0]-1
            data_spear.iloc[m, -1] = np.sum(data_to.iloc[data_spear.iloc[m, 0]:n, -2])/np.sum(data_to.iloc[data_spear.iloc[m, 0]:n, -1])
        print(m, i, 'SPEAR')
    data_spear.to_csv("D:\\HydroClass\\ceh\\04072022\\06072022_Spear_" + i.title() + "_taxon_data_final_SPECIES.csv",encoding='gbk')
    print(i)

