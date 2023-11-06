# -*- coding: utf-8 -*-
'''
SNP Scanner

@author: Fayten El-Dehaibi
Created on 27 Nov 2018
'''

import sys
import os
import PyQt5
from PyQt5 import QtCore, QtGui, QtWidgets
import xlsxwriter
import numpy
import pandas
import scipy
from scipy import stats, optimize
import time
import matplotlib.pyplot
import re

app = QtWidgets.QApplication(sys.argv)

class Form(QtWidgets.QDialog):
    def __init__(self, parent = None):
        super(Form, self).__init__(parent)
        
        pathLine = QtWidgets.QLineEdit()
        pathLine.setReadOnly(True)
        self.getFileButton = QtWidgets.QPushButton('...')
        sheetLine = QtWidgets.QLineEdit('Sheet1')
        cellLine = QtWidgets.QLineEdit('A1')
        #listbox of variables
        
        self.outputName = QtWidgets.QLineEdit()
        self.scanButton = QtWidgets.QPushButton('Find Control SNPs')

        fileLayout = QtWidgets.QHBoxLayout()
        fileLayout.addWidget(QtWidgets.QLabel('File:'))
        fileLayout.addWidget(pathLine)
        fileLayout.addWidget(self.getFileButton)
        pageLayout = QtWidgets.QHBoxLayout()
        pageLayout.addWidget(QtWidgets.QLabel('Sheet:'))
        pageLayout.addWidget(sheetLine)
        pageLayout.addWidget(QtWidgets.QLabel('Cell:'))
        pageLayout.addWidget(cellLine)
        outputLayout = QtWidgets.QHBoxLayout()
        outputLayout.addWidget(QtWidgets.QLabel('Output File Name:'))
        outputLayout.addWidget(self.outputName)
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(QtWidgets.QLabel('Find candidates for control SNPs in your data sheet.'))
        layout.addLayout(fileLayout)
        layout.addLayout(pageLayout)
        layout.addLayout(outputLayout)
        layout.addWidget(QtWidgets.QLabel(''))
        layout.addWidget(self.scanButton)
        self.setLayout(layout)
        self.setWindowTitle('SNP Scanner')

        self.getFileButton.clicked.connect(lambda: self.getPath(pathLine))
        self.scanButton.clicked.connect(lambda: self.errNoFile() if len(pathLine.text())<1 else self.scan(pathLine.text(),sheetLine.text(),cellLine.text()))

    def getPath(self,path):
        fileBox = QtWidgets.QFileDialog(self)
        inFile = QtCore.QFileInfo(fileBox.getOpenFileName(filter='*.xls* *.csv')[0])
        filepath = inFile.absoluteFilePath()
        if any([filepath.endswith('.xls'),filepath.endswith('.xlsx'),filepath.endswith('.csv')]):
            path.setText(filepath)

    def getXLSCoords(self,cell):
        if(cell.isspace()):
            cell.setText('A1')
        row,col = xlsxwriter.utility.xl_cell_to_rowcol(cell)
        a = None
        b = None
        if(row > 0):
            b = range(row)
        if(col > 0):
            a = range(col,16384)
        return a, b

    def readFile(self,path,sheet,cell):
        self.data = pandas.DataFrame()
        if(path.endswith('.xls') or path.endswith('.xlsx')):
            XLScoordR, XLScoordC = self.getXLSCoords(cell)
            self.data = pandas.read_excel(path,sheetname=sheet,
                                     skiprows=XLScoordR,parse_cols=XLScoordC)
        if(path.endswith('.csv')):
            df = pandas.read_csv(path,engine='python',header=0,iterator=True,
                                 chunksize=15000,infer_datetime_format=True)
            self.data = pandas.DataFrame(pandas.concat(df,ignore_index=True))

    def errNoFile(self):
        message = QtWidgets.QMessageBox.warning(self,'Error: No File Selected',
                                            'Please select a file.')

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

    def sortCytokineForSNP(self,snps,cyts,patientsAA,patientsBB,cytDict,startTime):
        print(str(len(snps))+' Candidate SNPs')
        step = int((len(snps))/100)
        significanceCount = pandas.DataFrame(columns=['SNP Name','Number Significant Mediators'],dtype='object')
        for snp in range(len(snps)):
            if(step > 0 and snp%step==0):
                print(str(snp/step)+'% complete: '+str(time.time()-startTime), str(snp))
            name = snps[snp]
            sigCytCount = 0
            for cyt in cyts:
                timeCytPairsAA = dict()
                timeCytPairsBB = dict()
                patients = [patientsAA[name],patientsBB[name]]
                pairsDicts = [timeCytPairsAA,timeCytPairsBB]
                ratio = str(len(patients[0])) + 'AA : ' + str(len(patients[1])) + 'BB'
                for gtype in range(2):
                    for p in patients[gtype]:
                        patientDF = cytDict[p]
                        for x in range(len(patientDF['time'])):
                            val = patientDF.loc[x,cyt]
                            timeVal = str(patientDF.loc[x,'time'])
                            if all([val is not None, val==val]):
                                pairsDicts[gtype][timeVal] = val
                if(len(timeCytPairsAA)==0 or len(timeCytPairsBB)==0):
                    continue
                AA = pandas.Series(data=timeCytPairsAA,dtype=float)
                AA = AA.dropna()
                BB = pandas.Series(data=timeCytPairsBB,dtype=float)
                BB = BB.dropna()
                mannWhitneyValue,p = scipy.stats.mannwhitneyu(AA,BB,alternative='two-sided')
                self.output = self.output.append(pandas.Series({'SNP':name,
                                                                'Mediator':cyt,
                                                                'p Value':p,
                                                                'Mann-Whitney U Score':mannWhitneyValue,
                                                                'ratio':ratio}),
                                                ignore_index=True)
                if p < 0.05:
                    sigCytCount += 1
            significanceCount = significanceCount.append(pandas.Series({'SNP Name':name,
                                                                        'Number Significant Mediators':sigCytCount}),
                                                                       ignore_index=True)
        if len(significanceCount.index) > 0:
            significanceCount = significanceCount.sort_values('Number Significant Mediators',axis=0,ascending=True)
            significanceCount.to_csv(self.outputName.text()+' Table.csv',index=False,sep=',',mode='w',chunksize=15000)
            graph = significanceCount.plot(kind='bar',title=self.outputName.text()+': Number Statistically Significant Inflammatory Mediators per SNP',fontsize=14)
            graph.set_xlabel('SNPs',fontsize=12)
            graph.set_ylabel('Frequency',fontsize=12)
            matplotlib.pyplot.show()
        else:
            print('No statistically significant SNP-mediator permutations found.')

    def scan(self,path,sheet,cell):
        print('Reading...')
        startTime = time.time()
        self.readFile(path,sheet,cell)
        print('Import time: ' + str(time.time()-startTime))
        startTime = time.time()
        startSNPS = self.data[self.data.Chr.notnull()].index[0]
        endSNPs = self.data[self.data.Normalization.notnull()].index[0]
        print(str(endSNPs - startSNPS) + ' SNPs')
        blanks = []
        checkForBlanks = pandas.DataFrame(pandas.isnull(self.data.iloc[0:startSNPS,12:]))
        checkForBlanksList = list(checkForBlanks.all())
        for x in range(len(list(checkForBlanks.columns))):
            if checkForBlanksList[x]:
                blanks.append(list(checkForBlanks.columns)[x])
        self.data = self.data.drop(blanks,axis='columns')
        timeIndices = self.data[self.data['patient_id'].map(lambda x: x == 'time')].index
        preCytokinesCols = list(set(self.data.iloc[timeIndices[0]:startSNPS, 0].dropna(axis=0)))
        cytokinesCols = [c for c in preCytokinesCols if(c == 'time' or c.startswith('plasma_2_'))]
        cyts = [cyt for cyt in cytokinesCols if cyt != 'time']
        cytokineIndices = {}
        for cyt in range(len(cytokinesCols)):
            cytokineIndices[cytokinesCols[cyt]] = list(self.data[self.data['patient_id'].map(lambda x: x == cytokinesCols[cyt])].index)
        cytIndVals =[*cytokineIndices.values()]
        allRows = [val for cyt in cytIndVals for val in cyt]
        allRows.extend([*range(startSNPS,endSNPs,1)])
        deathRow = self.data[self.data.patient_id == 'discharge_discharged_to']
        deathRow = deathRow[deathRow.isin(['Death'])]
        deathRow = deathRow.dropna(axis='columns', how='all')
        deadPats = list(deathRow.columns)
        deadPats_ = [c for c in list(self.dataTable.columns) if c not in deadPats]
        #deadPats_ = list(self.data.columns[:12]) + deadPats  # for nonsurvivor runs
        #deadPats_ = list(self.data.columns[:12]) + ['HR-056','HR-341','HR-423','HR-426','HR-454','HR-1005','HR-1139','HR-1194','HR-1283',
        #                                            'HR-198','HR-224','HR-260']  # for survivor matched run #HR-1261?
        self.data = self.data.loc[allRows,deadPats_]
        print('Checking genotypes...')
        checkAA = self.data[self.data.isin(['AA'])]
        checkBB = self.data[self.data.isin(['BB'])]
        patAA_ = checkAA.any(axis='index')
        patBB_ = checkBB.any(axis='index')
        patAA = list(patAA_[patAA_].index)
        patBB = list(patBB_[patBB_].index)
        patients = [p for p in list(self.data.columns) if (p in patAA or p in patBB)]
        print(str(time.time()-startTime) + ' for Database Reduction')
        print('Initializing patients...')
        #create 3d DF of patient -> cyt -> values
        patCytDict = {}
        time2 = time.time()
        for patient in range(len(patients)):
            p = patients[patient]
            pColumn = pandas.Series(self.data.loc[:,p],index=self.data.index)
            patCytDict[p] = self.processPatientsCyts(pColumn,cytokineIndices)
        print(str(time.time()-time2) + ' for Patient Preprocessing')
        SNPsAndPatientsAA = {}
        SNPsAndPatientsBB = {}
        print('Analyzing...')
        time2 = time.time()
        goodSNPs = []
        for i in range(startSNPS,endSNPs):
        #for i in range(startSNPS,startSNPS+5500):
            totalCount = self.data.loc[i,:].count()-12
            a = checkAA.loc[i,:].count()
            b = checkBB.loc[i,:].count()
            AA = a / totalCount
            BB = b / totalCount
            #if (a > 20 and b > 20) and ((((0.05 * AA) <= BB) and ((0.1 * AA) >= BB)) or (((0.05 * BB) <= AA) and ((0.1 * BB) >= AA))):  # rare SNPs
            #if (a > 0 and b > 0) and ((((0.05 * AA) <= BB) and ((0.1 * AA) >= BB)) or (((0.05 * BB) <= AA) and ((0.1 * BB) >= AA))):  # 'rare' SNP nonsurvivors
            if all([a > 20, b > 20, AA*0.9 <= BB, AA*1.1 >= BB]):#equal distribution SNPs
                SNPname = self.data.loc[i,'Name']
                goodSNPs.append(SNPname)
                SNPsAndPatientsAA[SNPname] = list(checkAA.loc[i,:].dropna(axis=0,how='all').index)
                SNPsAndPatientsBB[SNPname] = list(checkBB.loc[i,:].dropna(axis=0,how='all').index)
        self.output = pandas.DataFrame(columns=['SNP', 'Mediator', 'p Value', 'Mann-Whitney U Score','ratio'],dtype='object')
        print(str(time.time()-time2) + ' for SNP Screening')
        self.sortCytokineForSNP(goodSNPs, cyts, SNPsAndPatientsAA,SNPsAndPatientsBB,patCytDict,time.time())
        self.output = self.output.sort_values('p Value',axis=0,ascending=False)
        self.output.to_csv(self.outputName.text()+' Report.csv',index=False,sep=',',mode='w',chunksize=15000)
        print('Processing time: ' + str(time.time()-startTime))
        self.close()

form = Form()
form.show()
app.exec_()
