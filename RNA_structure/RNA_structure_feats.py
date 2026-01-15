import os
import oligo_melting as OligoMelt
import sys
from collections import deque
import RNA
import pandas as pd
import numpy as np
from time import time
scaffold = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
class Features:

    def __init__(self, scaffold, gRNA, dataset=None):
        self.scaffold = scaffold
        self.scaffoldStruct = RNA.fold_compound(scaffold).mfe()[0]
        self.guide = gRNA
        self.mergedStruct = RNA.fold_compound(gRNA + scaffold).mfe()[0]
        self.mergedEnergy = RNA.fold_compound(gRNA + scaffold).mfe()[1]
        self.lstIndexs = []
        self.initHeads()  # find the indexes of the three heads.
        self.firstHeadBegin = self.lstIndexs[0]
        self.firstHeadEnding = self.lstIndexs[1]
        self.secondHeadBegin = self.lstIndexs[2]
        self.secondHeadEnding = self.lstIndexs[3]
        self.thirdHeadBegin = self.lstIndexs[4]
        self.thirdHeadEnding = self.lstIndexs[5]
        self.dataset = dataset

    ''' 
     initHeads - 
     when creating instance of class Features, it will automatically 
     finds the "heads" of the scaffold structure.
    '''

    def initHeads(self):
        index = 0
        myStack = deque()
        while index < len(self.scaffoldStruct) - 1:  # go over scaffold clause struct

            while self.scaffoldStruct[index] != '(':  # locate the beginning index of the next head.
                index += 1
            self.lstIndexs.append(index)
            myStack.append(self.scaffoldStruct[index])
            while myStack:  # as long as the stack is not empty, we didnt locate the ending index of the current head
                index += 1
                if self.scaffoldStruct[index] == '(':
                    myStack.append(self.scaffoldStruct[index])
                if self.scaffoldStruct[index] == ')':
                    myStack.pop()
            self.lstIndexs.append(index)

    '''
    calculate guide energy
    '''

    def guideEnergy(self):
        ene = RNA.fold_compound(self.guide).mfe()[1]
        return ene

    '''
    calculate both guide & scaffold energy
    '''

    def guideAndScaffoldEnergy(self):
        return self.mergedEnergy, self.mergedStruct

    '''
    check if the connection went "well", in other words
    make sure the structure is achieving "base pairs"
    '''

    def basePairs(self):
        guide = self.mergedStruct[0:len(self.guide)]
        myStack = deque()
        index = 0
        while index < len(guide):
            if guide[index] == '(':
                myStack.append(guide[index])
            if guide[index] == ')':
                myStack.pop()
            index += 1

        if myStack:  # if the stack is not empty, means it connected well to the scaffold -> return true
            return True
        else:  # else = stack is empty, connection failed.
            return False

    '''
    confirm first head stays in place
    '''

    def isFirstHeadOk(self):
        originalHead = self.scaffoldStruct[self.lstIndexs[0]:self.lstIndexs[1]]
        newHead = self.mergedStruct[len(self.guide) + self.lstIndexs[0]: len(self.guide) + self.lstIndexs[1]]
        return originalHead == newHead

    '''
    confirm second head stays in place
    '''

    def isSecondHeadOk(self):
        originalHead = self.scaffoldStruct[self.lstIndexs[2]:self.lstIndexs[3]]
        newHead = self.mergedStruct[len(self.guide) + self.lstIndexs[2]: len(self.guide) + self.lstIndexs[3]]
        return originalHead == newHead

    '''
    confirm third head stays in place
    '''

    def isThirdHeadOk(self):
        originalHead = self.scaffoldStruct[self.lstIndexs[4]:self.lstIndexs[5]]
        newHead = self.mergedStruct[len(self.guide) + self.lstIndexs[4]: len(self.guide) + self.lstIndexs[5]]
        return originalHead == newHead

    '''
    confirm all heads stay in place
    '''

    def allHeadAreOk(self):
        return self.isFirstHeadOk() and self.isSecondHeadOk() and self.isThirdHeadOk()

def calculate(seq):
    try:
        seq_as_rna = seq.replace('T', 'U')
        (_, _, _, _, g_RNADNA, _) = OligoMelt.Duplex.calc_tm(seq_as_rna,celsius=True,tt_mode='RNA:DNA')
        (_, _, _, _, g_DNADNA, _) = OligoMelt.Duplex.calc_tm(seq,celsius=True,tt_mode='DNA:DNA')
        g_RNADNA = round(100*g_RNADNA)/100
        g_DNADNA = round(100*g_DNADNA)/100
        # g_DNADNA = Oligo_Melting.calculate_oligo_dna(seq)
        # g_RNADNA = Oligo_Melting.calculate_oligo_rna(seq_as_rna)
        f1 = Features(scaffold, seq)
        # print(f1)
        r1 = round(10*f1.guideEnergy())/10
        r2, r21 = f1.guideAndScaffoldEnergy()
        r2 = round(100*r2)/100
        r3 = f1.basePairs()
        r4 = np.bool(f1.isFirstHeadOk())
        r5 = np.bool(f1.isSecondHeadOk())
        r6 = np.bool(f1.isThirdHeadOk())
        r7 = np.bool(r5 and r6 and r4)
        # can either return results and str or as list.
        result = [g_DNADNA,   g_RNADNA,       r1,         r2,            r21, r3, r4, r5, r6, r7]
                #['g_DNADNA', 'g_RNADNA', 'guideEne', 'guide&scafEne', 'clause',
                   #r3,             r4,       r5,     r6,     r7
                 #'isBasePairs', 'Head1', 'Head2', 'Head3', '1&2&3']
        # result = [seq, r1, r2, r21, r3, r4, r5, r6, r7]
        return result
        # return str(seq) +" " + str(g_DNADNA) + " "+ str(g_RNADNA) + " " + str(r1) + " "+str(r2)+" "+str(r21)+" "+str(r3)+" "+ str(r4) +" "+str(r5)+" "+str(r6)+" "+str(r7)
    except Exception as e:
        print(e)
        return ""

isana_cols = ['g_DNADNA','g_RNADNA', 'guideEne', 'guide&scafEne', 'clause', 'isBasePairs', 'Head1', 'Head2', 'Head3', '1&2&3']
isana_cols_types = [np.nan,np.nan,np.nan,np.nan,'',np.bool,np.bool,np.bool,np.bool,np.bool]

# isana_cols = ['g_DNADNA','g_RNADNA', 'guideEne', 'guide&scafEne', 'isBasePairs', 'Head1', 'Head2', 'Head3', '1&2&3']
# isana_cols_types = [np.nan,np.nan,np.nan,np.nan,np.bool,np.bool,np.bool,np.bool,np.bool]
def make_RNA_struct_df(df):
    df = df.drop(columns=[col for col in df.columns if any(x in col for x in isana_cols)])
    new_cols_df = pd.DataFrame({col: col_type  for col,col_type in zip(isana_cols,isana_cols_types)}, index=df.index)
    df = pd.concat([df, new_cols_df], axis=1)
    return df

def calc_RNA_structure_feats(df):
    print('started calculating RNA structure feats!')
    df = make_RNA_struct_df(df)
    for i in df.index:
        if i%1000==0:
            print(f'{i/df.shape[0]:.3f}')
        seq = df.loc[i,'target_seq']
        seq = seq[:20]
        seq = seq.upper().replace('U','T')
        res = calculate(seq)
        df.loc[i,isana_cols] = res
    print('finished calculating RNA structure feats!')
    return df
