import sys
import glob
import re
import json
from decimal import Decimal
import pandas as pd
import numpy as np
import subprocess


#this script loads previous data
numArgErr=4
if (len(sys.argv)!=3):
    print("wrong number of arguments.")
    exit(numArgErr)


jsonDataFromConf =json.loads(sys.argv[1])
jsonFromSummary=json.loads(sys.argv[2])

potential_function_name=jsonDataFromConf["potential_function_name"]
U_dist_dataDir=jsonFromSummary["U_dist_dataDir"]
startingFileInd=jsonFromSummary["startingFileInd"]
startingVecPosition=jsonFromSummary["startingVecPosition"]


#search and read U_dist files
#give arbitrary values to L, y0, z0, y1 without reading data
UInit=6

y0Init=1
z0Init=1
y1Init=1
LInit=y0Init+z0Init+y1Init+0.7

loopLastFile=-1

#search csv files
csvFileList=[]
loopEndAll=[]
for file in glob.glob(U_dist_dataDir+"/*.csv"):
    csvFileList.append(file)
    matchEnd=re.search(r"loopEnd(\d+)",file)
    if matchEnd:
        loopEndAll.append(int(matchEnd.group(1)))


def create_loadedJsonData(UVal,LVal,y0Val,z0Val,y1Val,loopLastFileVal):
    """

    :param UVal:
    :param LVal:
    :param y0Val:
    :param z0Val:
    :param y1Val:
    :param loopLastFileVal:
    :return: loadedJsonData as string
    """
    initDataDict={
        "U":str(UVal),
        "L":str(LVal),
        "y0":str(y0Val),
        "z0":str(z0Val),
        "y1":str(y1Val),
        "loopLastFile":str(loopLastFileVal)
    }
    return json.dumps(initDataDict)

#if no data found, return the arbitrary values
if len(csvFileList)==0:
    loadedJsonDataStr=create_loadedJsonData(UInit,LInit,y0Init,z0Init,y1Init,loopLastFile)
    loadedJsonData_stdout="loadedJsonData="+loadedJsonDataStr
    print(loadedJsonData_stdout)
    exit(0)


#if found csv data
sortedEndInds=np.argsort(loopEndAll)
sortedLoopEnd=[loopEndAll[ind] for ind in sortedEndInds]
sortedCsvFileNames=[csvFileList[ind] for ind in sortedEndInds]
loopLastFile=sortedLoopEnd[-1]

lastFileName=sortedCsvFileNames[-1]

def get_last_row(csv_file):
    result = subprocess.run(['tail', '-n', '1', csv_file], stdout=subprocess.PIPE)
    last_row = result.stdout.decode('utf-8').strip()
    return last_row


csvLastRowStr=get_last_row(lastFileName)
valsInLastRow = [float(value) for value in csvLastRowStr.split(',')]

UInit=valsInLastRow[0]
LInit=valsInLastRow[1]
y0Init=valsInLastRow[2]
z0Init=valsInLastRow[3]
y1Init=valsInLastRow[4]


loadedJsonDataStr=create_loadedJsonData(UInit,LInit,y0Init,z0Init,y1Init,loopLastFile)
loadedJsonData_stdout="loadedJsonData="+loadedJsonDataStr
print(loadedJsonData_stdout)
exit(0)