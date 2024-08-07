import pickle
import numpy as np
from datetime import datetime

import pandas as pd
import statsmodels.api as sm
import sys
import re
import warnings
from scipy.stats import ks_2samp
import glob
from pathlib import Path
import os
import json

#This script checks if U, L, y0,z0,y1 values reach equilibrium and writes summary file of dist


argErrCode=2
sameErrCode=3
if (len(sys.argv)!=3):
    print("wrong number of arguments")
    exit(argErrCode)


jsonFromSummaryLast=json.loads(sys.argv[1])
jsonDataFromConf=json.loads(sys.argv[2])

TDirRoot=jsonFromSummaryLast["TDirRoot"]
U_dist_dataDir=jsonFromSummaryLast["U_dist_dataDir"]
effective_data_num_required=int(jsonDataFromConf["effective_data_num_required"])


summary_U_distFile=TDirRoot+"/summary_U_dist.txt"


def sort_data_files_by_loopEnd(oneDir):
    dataFilesAll=[]
    loopEndAll=[]
    for oneDataFile in glob.glob(oneDir+"/*.csv"):
        # print(oneDataFile)
        dataFilesAll.append(oneDataFile)
        matchEnd=re.search(r"loopEnd(\d+)",oneDataFile)
        if matchEnd:
            indTmp=int(matchEnd.group(1))
            loopEndAll.append(indTmp)
    endInds=np.argsort(loopEndAll)
    sortedDataFiles=[dataFilesAll[i] for i in endInds]
    return sortedDataFiles


def parseSummaryU_Dist():
    startingFileInd=-1
    startingVecPosition=-1

    summaryFileExists=os.path.isfile(summary_U_distFile)
    if summaryFileExists==False:
        return startingFileInd,startingVecPosition

    with open(summary_U_distFile,"r") as fptr:
        lines=fptr.readlines()
    for oneLine in lines:
        #match startingFileInd
        matchStartingFileInd=re.search(r"startingFileInd=(\d+)",oneLine)
        if matchStartingFileInd:
            startingFileInd=int(matchStartingFileInd.group(1))

        #match startingVecPosition
        matchStartingVecPosition=re.search(r"startingVecPosition=(\d+)",oneLine)
        if matchStartingVecPosition:
            startingVecPosition=int(matchStartingVecPosition.group(1))

    return startingFileInd, startingVecPosition

def auto_corrForOneColumn(colVec):
    """

    :param colVec: a vector of data
    :return:
    """
    same=False
    eps=1e-3
    NLags=int(len(colVec)*3/4)
    with warnings.catch_warnings():
        warnings.filterwarnings("error")
    try:
        acfOfVec=sm.tsa.acf(colVec,nlags=NLags)
    except Warning as w:
        same=True
    acfOfVecAbs=np.abs(acfOfVec)
    minAutc=np.min(acfOfVecAbs)

    lagVal=-1
    if minAutc<=eps:
        lagVal=np.where(acfOfVecAbs<=eps)[0][0]
    return same,lagVal

def ksTestOneColumn(colVec,lag):
    """

    :param colVec: a vector of data
    :param lag: auto-correlation length
    :return:
    """
    colVecSelected=colVec[::lag]

    lengthTmp=len(colVecSelected)
    if lengthTmp%2==1:
        lengthTmp-=1
    lenPart=int(lengthTmp/2)

    colVecToCompute=colVecSelected[-lengthTmp:]

    #ks test
    selectedVecPart0=colVecToCompute[:lenPart]
    selectedVecPart1=colVecToCompute[lenPart:]
    result=ks_2samp(selectedVecPart0,selectedVecPart1)
    return result.pvalue, lenPart*2

def checkU_distDataFilesForOneT(U_dist_csv_dir):
    """

    :param dist_csv_dir:
    :return:
    """
    U_dist_sortedDataFilesToRead=sort_data_files_by_loopEnd(U_dist_csv_dir)
    if len(U_dist_sortedDataFilesToRead)==0:
        print("no data for U_dist.")
        exit(0)

    startingFileInd,startingVecPosition=parseSummaryU_Dist()

    if startingFileInd<0:
        #we guess that the equilibrium starts at this file
        startingFileInd=int(len(U_dist_sortedDataFilesToRead)*5/7)
    startingFileName=U_dist_sortedDataFilesToRead[startingFileInd]
    # print(startingFileName)
    #read the starting U_dist csv file
    in_dfStart=pd.read_csv(startingFileName)

    in_nRowStart,in_nColStart=in_dfStart.shape
    if startingVecPosition<0:
        #we guess equilibrium starts at this position
        startingVecPosition=int(in_nRowStart/2)



    distArray=np.array(in_dfStart)

    #read the rest of the csv files
    for csv_file in U_dist_sortedDataFilesToRead[(startingFileInd+1):]:
        in_df=pd.read_csv(csv_file)

        in_array=np.array(in_df)

        distArray=np.r_[distArray,in_array]
        # y0Vec=np.r_[y0Vec,in_y0]
        # z0Vec=np.r_[z0Vec,in_z0]
        # y1Vec=np.r_[y1Vec,in_y1]




    dist_same=[]
    dist_lags=[]



    _,nCol=distArray.shape
    # print("nrow="+str(nrow))
    # print("ncol="+str(ncol))
    for j in range(0,nCol):
        sameTmp,lagTmp=auto_corrForOneColumn(distArray[:,j])
        dist_same.append(sameTmp)
        dist_lags.append(lagTmp)


    same_exist=any(dist_same)

    #all values are the same, exit with err code
    if same_exist==True:
        with open(summary_U_distFile,"w+") as fptr:
            msg="error: same\n"
            fptr.writelines(msg)
            exit(sameErrCode)
    pThreshHold=0.01
    #if one lag==-1, then the auto-correlation is too large
    if np.min(dist_lags)>0:
        lagMax=np.max(dist_lags)

        pValsAll=[]
        lengthValAll=[]
        for j in range(0,nCol):
            pTmp,lengthTmp=ksTestOneColumn(distArray[:,j],lagMax)
            pValsAll.append(pTmp)
            lengthValAll.append(lengthTmp)
        numDataPoints=np.min(lengthValAll)
        if np.min(pValsAll)>=pThreshHold and numDataPoints>=200:
            if numDataPoints>=effective_data_num_required:
                newDataPointNum=0
            else:
                newDataPointNum=effective_data_num_required-numDataPoints

            msg="equilibrium\n" \
                +"lag="+str(lagMax)+"\n" \
                +"numDataPoints="+str(numDataPoints)+"\n" \
                +"startingFileInd="+str(startingFileInd)+"\n" \
                +"startingVecPosition="+str(startingVecPosition)+"\n" \
                +"newDataPointNum="+str(newDataPointNum)+"\n"

            with open(summary_U_distFile,"w+") as fptr:
                fptr.writelines(msg)
            exit(0)

        continueMsg="continue\n"
        if np.min(pValsAll)<pThreshHold:
            #not the same distribution
            continueMsg+="p value: "+str(np.min(pValsAll))+"\n"
        if numDataPoints<200:
            #not enough data number

            continueMsg+="numDataPoints="+str(numDataPoints)+" too low\n"
            continueMsg+="lag="+str(lagMax)+"\n"
        with open(summary_U_distFile,"w+") as fptr:
            fptr.writelines(continueMsg)
        exit(0)

    else:
        msg="high correlation"
        with open(summary_U_distFile,"w+") as fptr:
            fptr.writelines(msg)
        exit(0)










checkU_distDataFilesForOneT(U_dist_dataDir)