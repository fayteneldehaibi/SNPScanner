# -*- coding: utf-8 -*-
'''
HaplotypeCheck

@author: Fayten El-Dehaibi
Created on 10 Nov 2023
'''

import sys
import os
import PyQt5
from PyQt5 import QtCore, QtGui, QtWidgets
import xlsxwriter
import numpy
import pandas
import multiprocessing
import scipy
from scipy import stats, optimize
import re
import time

app = QtWidgets.QApplication(sys.argv)

def MannWhitney(inputList):
    param = inputList[0]
    patients = inputList[1]
    paramDict = inputList[2]
    timeCytPairs1 = dict()
    timeCytPairs2 = dict()
    pairsDicts = [timeCytPairs1, timeCytPairs2]
    for group in range(2):
        for patient in patients[group]:
            patientCytDF = paramDict[patient]
            for x in range(len(patientCytDF['time'])):
                val = patientCytDF.loc[x, param]
                timeVal = str(patientCytDF.loc[x, 'time'])
                if all([val is not None, val == val]):
                    pairsDicts[group][timeVal] = val
    if (len(timeCytPairs1) == 0 or len(timeCytPairs2) == 0):
        return None
    group1 = pandas.Series(data=timeCytPairs1, dtype=float)
    group1 = group1.dropna()
    group2 = pandas.Series(data=timeCytPairs2, dtype=float)
    group2 = group2.dropna()
    mannWhitneyValue, pVal = scipy.stats.mannwhitneyu(group1, group2, alternative='two-sided')
    return [param,pVal,mannWhitneyValue]

class Form(QtWidgets.QDialog):
    def __init__(self, parent = None):
        super(Form, self).__init__(parent)

        pathLine = QtWidgets.QLineEdit()
        pathLine.setReadOnly(True)
        self.getFileButton = QtWidgets.QPushButton('...')
        sheetLine = QtWidgets.QLineEdit('Sheet1')
        self.snpBox = QtWidgets.QTextEdit()
        self.snpBox.resize(80,80)
        nHaploSize = QtWidgets.QSpinBox()
        nHaploSize.setMinimum(2)
        nHaploSize.setMaximum(4)
        paramList = QtWidgets.QListWidget()
        paramList.setToolTip('Double-click a row to remove it.')
        newParam = QtWidgets.QLineEdit()
        paramBtn = QtWidgets.QPushButton('Add to List')
        self.outputLine = QtWidgets.QLineEdit('Output Name')
        runBtn = QtWidgets.QPushButton('Run!')

        fileLayout = QtWidgets.QHBoxLayout()
        fileLayout.addWidget(QtWidgets.QLabel('File:'))
        fileLayout.addWidget(pathLine)
        fileLayout.addWidget(self.getFileButton)
        pageLayout = QtWidgets.QHBoxLayout()
        pageLayout.addWidget(QtWidgets.QLabel('Sheet:'))
        pageLayout.addWidget(sheetLine)
        haploLayout = QtWidgets.QHBoxLayout()
        haploLayout.addWidget(QtWidgets.QLabel('Haplotype Size:'))
        haploLayout.addWidget(nHaploSize)
        paramLayout = QtWidgets.QHBoxLayout()
        paramLayout.addWidget(QtWidgets.QLabel('New Parameter:'))
        paramLayout.addWidget(newParam)
        paramLayout.addWidget(paramBtn)
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(QtWidgets.QLabel('Detect haplotypes from list of SNPs.'))
        layout.addLayout(fileLayout)
        layout.addLayout(pageLayout)
        layout.addWidget(QtWidgets.QLabel(''))
        layout.addWidget(QtWidgets.QLabel('List of SNP Names (space delimited)'))
        layout.addWidget(self.snpBox)
        layout.addLayout(haploLayout)
        layout.addWidget(QtWidgets.QLabel('Parameters to Check'))
        layout.addWidget(paramList)
        layout.addLayout(paramLayout)
        layout.addWidget(QtWidgets.QLabel(''))
        layout.addWidget(self.outputLine)
        layout.addWidget(runBtn)
        self.setLayout(layout)
        self.setWindowTitle('Haplotype Check')

        self.getFileButton.clicked.connect(lambda: self.getPath(pathLine))
        paramBtn.clicked.connect(lambda: paramList.addItem(newParam.text()) if len(newParam.text()) > 1 else None)
        paramList.itemDoubleClicked.connect(lambda: paramList.takeItem(paramList.currentRow()))
        runBtn.clicked.connect(lambda: self.errNoFile() if pathLine.text().isspace()
                                       else self.run(nHaploSize.value(),paramList,pathLine.text(),sheetLine.text()))

    def getPath(self,path):
        fileBox = QtWidgets.QFileDialog(self)
        inFile = QtCore.QFileInfo(fileBox.getOpenFileName(filter='*.xls* *.csv')[0])
        filepath = inFile.absoluteFilePath()
        if any([filepath.endswith('.xls'),filepath.endswith('.xlsx'),filepath.endswith('.csv')]):
            path.setText(filepath)

    def readFile(self,path,sheet):
        self.data = pandas.DataFrame()
        if(path.endswith('.xls') or path.endswith('.xlsx')):
            self.data = pandas.read_excel(path,sheet_name=sheet)
        if(path.endswith('.csv')):
            df = pandas.read_csv(path,engine='python',header=0,iterator=True,
                                 chunksize=15000,infer_datetime_format=True)
            self.data = pandas.DataFrame(pandas.concat(df,ignore_index=True))

    def errNoFile(self):
        message = QtWidgets.QMessageBox.warning(self,'Error: No File Selected',
                                            'Please select an input file.')

    def paramNotFound(self):
        message = QtWidgets.QMessageBox.warning(self, 'Error: No Parameters Found',
                                                'Parameters given not found in data set.')

    def processPatientsCyts(self,col,cytDict):
        badCharRegex = re.compile(r'>|<|%*')
        patDF = pandas.DataFrame()
        for cyt in cytDict.keys():
            values = []
            for x in range(len(cytDict['time'])):
                vals = cytDict[cyt]
                val = col.loc[vals[x]]
                if all([val is not None, val == val]):
                    if cyt == 'time':
                        if ('h' in val):
                            if val.startswith('h'):
                                val = val[1:] + 'h'
                            if ('r' in val):
                                val = ''.join(val.split('r'))
                            if (len(val) < 3):
                                val = '0' + val
                    else:
                        if any(['>' in val, '<' in val, '*' in val]):
                            val = ''.join(re.split(badCharRegex, val))
                elif all([val != val,cyt != 'time']):
                    val = None
                values.append(val)
            patDF[cyt] = values
        patDF = patDF.dropna(axis='index',how='all')
        return patDF

    def get2SNPCombos(self,snpList):
        # get combinations of SNPs
        comboList = []
        for x in range(len(snpList)-1):
            last = snpList.pop()
            subset = [[last, y] for y in snpList]
            comboList = comboList + subset
        return comboList

    def checkHaplotypesCytokine(self,patients,cyts,cytDict):
        sigCount = []
        to_Output = []
        with multiprocessing.Pool(4) as pool:
            chunkForSNP = pandas.DataFrame()
            toMWU = [[cyt, patients, cytDict] for cyt in cyts]
            results_list = pool.map(MannWhitney, iterable=toMWU)
            series_list = [pandas.DataFrame(data={'Mediator':r[0],'p Value':r[1],'Mann-Whitney U Score':r[2]},index=[0])
                           for r in results_list if r is not None]
            chunkForSNP = pandas.concat(series_list, axis=0, ignore_index=True)
            chunkForSNP = chunkForSNP.dropna(how='any')
            to_Output.append(chunkForSNP)
        output = pandas.concat(to_Output,ignore_index=True,axis=0)
        return output

    def run(self,nHaplo,paramList,path,sheet):
        params = [paramList.item(i).text() for i in range(paramList.count())]
        snps_ = re.split('\s',self.snpBox.toPlainText())
        snps_ = list(set(snps_))#Removes duplicates
        snps = [s for s in snps_ if len(s) > 0]#Removes null entries Python dectects
        if len(snps) < 2:
            message = QtWidgets.QMessageBox.warning(self, 'Not Enough SNPs',
                                                     'Enter at least 2 SNPs to be compared.')
            return 0
        print('Reading...')
        startTime = time.time()
        self.readFile(path, sheet)
        print('Import time: ' + str(time.time() - startTime))
        paramsPresent = [param for param in params if any(self.data.iloc[:, 0].str.contains(param))]
        if len(params) > len(paramsPresent):
            if len(paramsPresent) > 0:
                message = QtWidgets.QMessageBox.question(self, 'Not All Parameters Found',
                                                         'Some parameters given were not found. Analyze data with found parameters?',
                                                         QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                                                         QtWidgets.QMessageBox.Yes)
                if message == QtWidgets.QMessageBox.No:
                    return 0
            else:
                self.paramNotFound()
                return 0
        #Reduce database to necessary patients, SNPs, and parameters
        startTime = time.time()
        startSNPS = self.data[self.data.Chr.notnull()].index[0]
        snpIndices = [self.data[self.data.Name==n].index[0] for n in snps]
        '''
        blanks = []
        checkForBlanks = pandas.DataFrame(pandas.isnull(self.data.iloc[0:startSNPS, 12:]))
        checkForBlanksList = list(checkForBlanks.all())
        for x in range(len(list(checkForBlanks.columns))):
            if checkForBlanksList[x]:
                blanks.append(list(checkForBlanks.columns)[x])
        self.data = self.data.drop(blanks, axis='columns')
        deathRow = self.data[self.data.patient_id == 'discharge_discharged_to']
        deathRow = deathRow[deathRow.isin(['Death'])]
        deathRow = deathRow.dropna(axis='columns', how='all')
        deadPats = list(deathRow.columns)
        '''
        timeIndices = list(self.data[self.data['patient_id'].map(lambda x: x == 'time')].index)
        preCytokinesCols = list(set(self.data.iloc[timeIndices[0]:startSNPS, 0].dropna(axis=0)))
        cyts = [c for c in preCytokinesCols if any(re.findall(r'|'.join(paramsPresent), c, re.IGNORECASE))]
        if len(cyts) < 1:
            self.paramNotFound()
            return 0
        cytokineIndices = {}
        for cyt in cyts:
            cytokineIndices[cyt] = list(self.data[self.data['patient_id'].map(lambda x: x == cyt)].index)
        cytokineIndices['time'] = timeIndices
        cytIndVals = [*cytokineIndices.values()]
        allRows = sorted([val for cyt in cytIndVals for val in cyt])
        #v Remove later
        pittFile = pandas.read_excel('H:/SNP Scanner/rs10404939 Demographic Table.xlsx',
                                     sheet_name='Deriv Total',header=0)
        hr1pats_ = list(pittFile.loc[:,'patient_id'])
        hr1pats = [str(h) for h in hr1pats_ if h==h]
        #^ Remove later
        keepCols = [c for c in list(self.data.columns) if
                    #all([c not in deadPats,len(self.data.loc[allRows, c]) > 0,c in hr1pats])]  # survivors with cytokines
                    c in hr1pats]#remove later
        print('Patients: ' + str(len(keepCols)))
        allRows.extend(snpIndices)
        self.data = self.data.loc[allRows, keepCols]
        print('Database Reduction: ' + str(time.time() - startTime))
        patCytDict = {}
        time0 = time.time()
        for patient in keepCols:
            pColumn = pandas.Series(self.data.loc[:, patient], index=self.data.index)
            patCytDict[patient] = self.processPatientsCyts(pColumn, cytokineIndices)
        print('Cytokine Processing: ' + str(time.time() - time0))
        #Get list of patients for each genotype for each SNP
        time1 = time.time()
        SNPsAndPats = {}
        for snp,snpInd in zip(snps,snpIndices):
            line = self.data.loc[snpInd,:]
            try:
                SNPsAndPats[snp] = [list(line[line=='AA'].columns),
                                    list(line[line=='AB'].columns),
                                    list(line[line=='BB'].columns)]
            except(AttributeError):
                SNPsAndPats[snp] = [list(line[line == 'AA'].index),
                                    list(line[line == 'AB'].index),
                                    list(line[line == 'BB'].index)]
        allCombos = []
        if nHaplo==2:
            allCombos = self.get2SNPCombos(snps)
        print(str(len(allCombos)) + ' SNP Combinations Acquired: '+ str(time.time()-time1))
        time2 = time.time()
        resultDF = pandas.DataFrame(columns=['SNP 1', 'gtype 1', 'n 1', 'SNP 2', 'gtype 2', ' n 2',
                                             'n Intersection', '% 1 Represented', '% 2 Represented',
                                             'Intersection'])
        results = []
        gtypes = ['AA','AB','BB']
        for combo in allCombos:
            snp1 = combo[0]
            snp2 = combo[1]
            for gtype1 in range(3): #0 == AA, 1 == AB, 2==BB
                pats1 = SNPsAndPats[snp1][gtype1]
                if len(pats1) < 1:
                    continue
                for gtype2 in range(3):
                    pats2 = SNPsAndPats[snp2][gtype2]
                    if len(pats2) < 1:
                        continue
                    haplotype = sorted(list(set(pats1).intersection(set(pats2))))
                    percent1 = len(haplotype)/len(pats1)
                    percent2 = len(haplotype)/len(pats2)
                    result = pandas.Series({'SNP 1':snp1,
                                            'gtype 1':gtypes[gtype1],
                                            'n 1':len(pats1),
                                            'SNP 2':snp2,
                                            'gtype 2':gtypes[gtype2],
                                            'n 2':len(pats2),
                                            'n Intersection':len(haplotype),
                                            '% 1 Represented':percent1,
                                            '% 2 Represented':percent2,
                                            'Intersection':haplotype})
                    results.append(result)
        resultDF = pandas.concat(results,axis=1).T
        dfResult = resultDF[resultDF['n Intersection'] >= 20]
        resultDF = resultDF.sort_values('n Intersection',ascending=False)
        print(str(time.time()-time2) + 's Processing Time')
        resultDF.to_csv(self.outputLine.text() + '.csv',index=False,sep=',',mode='w',chunksize=15000)
        #compare inflammation of haplotypes with their 'opposites' (e.g. (SNP1 AA & SNP2 BB) vs. (SNP1 BB & SNP2 AA))
        dfResult = dfResult[(dfResult['gtype 1'] != 'AB') & (dfResult['gtype 2'] != 'AB')]
        results2 = []
        shortResults = []
        print('Beginning cyt comparison...')
        step = int((len(dfResult.index)) / 100)
        time1 = time.time()
        i = 0
        for part in list(dfResult.index):
            if (step > 0 and i % step == 0):
                print(str(i / step) + '% complete: ' + str(time.time() - time1), str(i))
            i += 1
            result2 = pandas.DataFrame()
            shortResult = pandas.DataFrame()
            try:
                SNP1 = dfResult.loc[part, 'SNP 1']
            except(KeyError):
                continue
            gtype1 = dfResult.loc[part, 'gtype 1']
            SNP2 = dfResult.loc[part, 'SNP 2']
            gtype2 = dfResult.loc[part, 'gtype 2']
            patientsPart = dfResult.loc[part, 'Intersection']
            try:
                if(gtype1=='AA' and gtype2=='AA') or (gtype1=='BB' and gtype2=='BB'):#For pairs where both are AA or BB
                    gtype2 = 'BB' if gtype1=='AA' else 'AA'
                    counterpart = dfResult[(((dfResult['SNP 1'] == SNP1) & (dfResult['gtype 1'] == gtype1))
                                            & ((dfResult['SNP 2'] == SNP2) & (dfResult['gtype 2'] == gtype1)))
                                           | (((dfResult['SNP 1'] == SNP2) & (dfResult['gtype 1'] == gtype2))
                                              & ((dfResult['SNP 2'] == SNP1) & (dfResult['gtype 2'] == gtype2)))].index[0]
                    name1 = SNP1 + ' ' + gtype1 + ' & ' + SNP2 + ' ' + gtype1
                    name2 = SNP1 + ' ' + gtype2 + ' & ' + SNP2 + ' ' + gtype2
                else:#For pairs where one is AA and the other is BB
                    counterpart = dfResult[(((dfResult['SNP 1'] == SNP1) & (dfResult['gtype 1'] == gtype2))
                                            & ((dfResult['SNP 2'] == SNP2) & (dfResult['gtype 2'] == gtype1)))
                                           | (((dfResult['SNP 1'] == SNP2) & (dfResult['gtype 1'] == gtype1))
                                              & ((dfResult['SNP 2'] == SNP1) & (dfResult['gtype 2'] == gtype2)))].index[0]
                    name1 = SNP1 + ' ' + gtype1 + ' & ' + SNP2 + ' ' + gtype2
                    name2 = SNP1 + ' ' + gtype2 + ' & ' + SNP2 + ' ' + gtype1
            except(IndexError):
                continue
            patientsCounter = dfResult.loc[counterpart, 'Intersection']
            # send patient lists to MWU
            result2 = self.checkHaplotypesCytokine([patientsPart,patientsCounter],cyts,patCytDict)
            sigParamCount = result2['p Value'].lt(0.05).sum()
            shortResult = pandas.Series(data = {'Group 1':name1,
                                                'n 1':dfResult.loc[part,'n Intersection'],
                                                'Group 2':name2,
                                                'n 2': dfResult.loc[counterpart, 'n Intersection'],
                                                '# Significant Cyts':sigParamCount})
            result2.insert(0,'HType 1',name1)
            result2.insert(1,'HType 2',name2)
            result2['Count 1'] = dfResult.loc[part,'n Intersection']
            result2['Count 2'] = dfResult.loc[counterpart,'n Intersection']
            results2.append(result2)
            shortResults.append(shortResult)
            dfResult = dfResult.drop(counterpart,axis='index')
        cytResult = pandas.concat(results2,axis=0,ignore_index=True)
        cytResult = cytResult.sort_values('p Value',ascending=False)
        shortDF = pandas.concat(shortResults,axis=1,ignore_index=True).T
        shortDF = shortDF.sort_values('# Significant Cyts')
        print('Cyt Comparison: ' + str((time.time()-time1)))
        cytResult.to_csv(self.outputLine.text() + ' HType Report.csv',index=False,sep=',',mode='w',chunksize=15000)
        shortDF.to_csv(self.outputLine.text() + ' HType Table.csv',index=False,sep=',',mode='w',chunksize=15000)
        self.close()

if __name__=="__main__":
    form = Form()
    form.show()
    app.exec_()