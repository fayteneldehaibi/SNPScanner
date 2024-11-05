# -*- coding: utf-8 -*-
'''
#SNPControlCheck

@author: Fayten El-Dehaibi
Created on 27 Nov 2018
'''

import sys
import os
import xlsxwriter
import numpy
import pandas
import scipy
from scipy import stats, optimize
import time
import re
from itertools import compress
import random

print('Reading...')
startTime = time.time()
data_ = pandas.read_csv('C:/Users/Monster Workstation/Documents/Gargantuan Genome Dataset/Pitt Metabolomics+Sc2i.csv',
                        engine='python', header=0, iterator=True, chunksize=15000,skipinitialspace=True)
data = pandas.concat(data_,ignore_index=True)
countTable_ = pandas.read_csv('C:/Users/Monster Workstation/Documents/vodovotz-database-query/SNPScanner/All Data Fixed Table.csv',
                              engine='python', header=0, iterator=True, chunksize=15000,skipinitialspace=True)
countTable = pandas.concat(countTable_,ignore_index=True)
goodPats_ = pandas.read_csv('C:/Users/Monster Workstation/Documents/vodovotz-database-query/SNPScanner Original Pool.csv',
                            engine='python', header=0, iterator=True, chunksize=15000,skipinitialspace=True)
goodPats = pandas.concat(goodPats_,ignore_index=True)
print('Read!')

outputTitle = 'Even ControlSNP Randomized Clinical Analysis--correct.xlsx'

#pull index of 0 cyt snps from data.Names
lastIndex = countTable[countTable['Number Significant Mediators']==1].index[0]
snpNames = list(countTable.loc[:lastIndex,'SNP Name'])
controlIndex = [data[data.Name==n].index[0] for n in snpNames] #DO NOT SORT THIS LIST!
#KEEP LIST OF CONTROLSNPINDEX

#screen for proper patient columns
blanks = []
checkForBlanks = pandas.DataFrame(pandas.isnull(data.iloc[controlIndex,12:]))
checkForBlanksList = list(checkForBlanks.all())
for x in range(len(list(checkForBlanks.columns))):
            if checkForBlanksList[x]:
                blanks.append(list(checkForBlanks.columns)[x])
data = data.drop(blanks,axis='columns')
timeIndices = data[data['patient_id'].map(lambda x: x == 'time')].index
preCytokinesCols = list(set(data.loc[timeIndices[0]:controlIndex[0], list(data.columns)[0]].dropna(axis=0)))
cytokinesCols = [c for c in preCytokinesCols if(str(c) == 'time' or str(c).startswith('plasma_2_'))]
cyts = [cyt for cyt in cytokinesCols if cyt != 'time']
cytokineIndices = {}
for cyt in range(len(cytokinesCols)):
    cytokineIndices[cytokinesCols[cyt]] = list(data[data['patient_id'].map(lambda x: x == cytokinesCols[cyt])].index)
cytIndVals =[*cytokineIndices.values()]
deathRow = data[data.patient_id == 'discharge_discharged_to']
deathRow = deathRow[deathRow.isin(['Death'])]
deathRow = deathRow.dropna(axis='columns', how='all')
deadPats = list(deathRow.columns)
survivors = [c for c in list(data.columns)[12:] if c not in deadPats]
#survivors = [c for c in list(goodPats.columns)[4:]] #Used for original 380 patient cohort only
goodCols = ['patient_id','Name'] + survivors
data = data.loc[:,goodCols]

#don't need to run checkAA/checkBB, because these SNPs are already screened for equivalent distribution
#reduce total data index to ICULOS,total LOS, vent, and MODS + controlIndex
params = ['discharge_total_icu_days','discharge_total_hospital_days','discharge_total_ventilator_days',
         'MOD1','MOD2','MOD3','MOD4','MOD5','MOD6','MOD7']
paramIndex = [data[data.patient_id==p].index[0] for p in params]
dataIndex = paramIndex + controlIndex
data = data.iloc[dataIndex,:]
print('DataFrame reduced')
#KEEP LIST OF PARAMETERINDEX
paramCols = ['ICU','TotalLOS','Vent','MOD1','MOD2','MOD3','MOD4','MOD5','MOD6','MOD7']
patientsDF = pandas.DataFrame(columns=paramCols)
patientSeries_ = []
for s in survivors:
    patientSeries = pandas.Series({'ID':s,
                                   'ICU':data.loc[paramIndex[0],s],
                                   'TotalLOS':data.loc[paramIndex[1],s],
                                   'Vent':data.loc[paramIndex[2],s],
                                   'MOD1':data.loc[paramIndex[3],s],
                                   'MOD2':data.loc[paramIndex[4],s],
                                   'MOD3':data.loc[paramIndex[5],s],
                                   'MOD4':data.loc[paramIndex[6],s],
                                   'MOD5':data.loc[paramIndex[7],s],
                                   'MOD6':data.loc[paramIndex[8],s],
                                   'MOD7':data.loc[paramIndex[9],s]})
    patientSeries_.append(patientSeries)
patientsDF = pandas.concat(patientSeries_,axis=1)
patientsDF.columns = patientsDF.iloc[0,:]
patientsDF = patientsDF.drop(patientsDF.index[0],axis='index')
#NOW do checkAA/checkBB to get patientAA/BB {name:[<columnList>]}
checkAA = data.isin(['AA'])
checkBB = data.isin(['BB'])
SNPsAndPatientsAA = {}
SNPsAndPatientsBB = {}
for ctlIndex,SNPname in zip(controlIndex,snpNames):
    lineAA = checkAA.loc[ctlIndex,:][checkAA.loc[ctlIndex,:]]
    lineBB = checkBB.loc[ctlIndex,:][checkBB.loc[ctlIndex,:]]
    SNPsAndPatientsAA[SNPname] = list(lineAA.index)
    SNPsAndPatientsBB[SNPname] = list(lineBB.index)
overallResults = []
shortResults = []
startTime2 = time.time()
results = []
for snp in snpNames:
    dataDFs = [pandas.DataFrame(),pandas.DataFrame()]
    dataDicts = [SNPsAndPatientsAA,SNPsAndPatientsBB]
    for gtype in range(2):
        theseCols = dataDicts[gtype][snp] #list of columns
        dataDFs[gtype] = patientsDF.loc[:,theseCols].T
    dataDFAA = dataDFs[0]
    dataDFBB = dataDFs[1]
    mwu0,pICU = scipy.stats.mannwhitneyu(pandas.Series(dataDFAA['ICU'],dtype=float).dropna(),
                                         pandas.Series(dataDFBB['ICU'],dtype=float).dropna(),alternative='two-sided')
    mwu1,pLOS = scipy.stats.mannwhitneyu(pandas.Series(dataDFAA['TotalLOS'],dtype=float).dropna(),
                                         pandas.Series(dataDFBB['TotalLOS'],dtype=float).dropna(),alternative='two-sided')
    mwu2,pVent = scipy.stats.mannwhitneyu(pandas.Series(dataDFAA['Vent'],dtype=float).dropna(),
                                          pandas.Series(dataDFBB['Vent'],dtype=float).dropna(),alternative='two-sided')
    # v This fixes MODS
    modSeriesAA = pandas.Series(pandas.concat([dataDFAA['MOD1'],dataDFAA['MOD2'],dataDFAA['MOD3'],
                                 dataDFAA['MOD4'],dataDFAA['MOD5'],dataDFAA['MOD6'],
                                 dataDFAA['MOD7']],axis=0),dtype=float).dropna()
    modSeriesBB = pandas.Series(pandas.concat([dataDFBB['MOD1'],dataDFBB['MOD2'],dataDFBB['MOD3'],
                                 dataDFBB['MOD4'],dataDFBB['MOD5'],dataDFBB['MOD6'],
                                 dataDFBB['MOD7']],axis=0),dtype=float).dropna()
    mwu3,pMODS = scipy.stats.mannwhitneyu(modSeriesAA,modSeriesBB,alternative='two-sided')
    results.append(pandas.Series({'SNP':snp,
                                       'ICU pVal':pICU,
                                       'TotalLOS pVal':pLOS,
                                       'Vent pVal':pVent,
                                       'MODS pVal':pMODS}))
batchResultDF = pandas.concat(results,axis=1,ignore_index=True).T
batchResultDF['Significance Count'] = batchResultDF.iloc[:,1:].lt(0.05,axis=0).sum(axis=1)
overallResults.append(batchResultDF)
finalResults = pandas.concat(overallResults,ignore_index=True)
print(finalResults.loc[:,'Significance Count'].value_counts())

writer =  pandas.ExcelWriter(outputTitle, engine='xlsxwriter')
finalResults.to_excel(writer,'Control SNP Batch MWU',index=False)
writer.close()
print(time.time()-startTime)
print('Done!')