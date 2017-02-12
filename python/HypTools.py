#This module serves the purpose of converting back and forth between file line number and hypothesis in log space

import numpy as np
import SterileSearchPy as ssp


def Make2DHypDicts(FileName):
    Data=np.loadtxt(FileName)
    PointToHypID={}
    HypIDToPoint={}
    MaxHypID=[]
    AllSin2s=[]
    AllDm2s=[]
    for i in range(0,len(Data)):
        HypID=i+1
        Sin2=np.round(Data[i][3],4)
        Dm2=np.round(Data[i][2],4)
        PointToHypID[(Sin2,Dm2)]=HypID
        HypIDToPoint[HypID]=(Sin2,Dm2)
        MaxHypID=HypID
        if(not Sin2 in AllSin2s):
            AllSin2s.append(Sin2)
        if(not Dm2 in AllDm2s):
            AllDm2s.append(Dm2)
    ReturnDict={}
    ReturnDict['ParamsToHypothesisNumber']=PointToHypID
    ReturnDict['HypothesisNumberToParams']=HypIDToPoint
    ReturnDict['MaxHypID']=MaxHypID
    ReturnDict['AllParamVals']=[AllSin2s,AllDm2s]
    return ReturnDict



def Make2DScanFile(FileName, FirstLogDM2=-2, LastLogDM2=2., LogDM2Step=0.05, FirstLogSin2=-3, LastLogSin2=0, LogSin2Step=0.05):

    File=open(FileName,'w')
    for LogDM2 in np.arange(FirstLogDM2,LastLogDM2,LogDM2Step):
        for Sin2 in np.arange(FirstLogSin2,LastLogSin2,LogSin2Step):
            Thet=np.arcsin(pow(pow(10,Sin2),0.5))/2.
            File.write(str(np.round(pow(10,LogDM2),6))+" " + str(np.round(Thet,6)) + " " + str(np.round( LogDM2,6)) + " " +str(np.round(Sin2,6)) +'\n')
    File.close()


def MakeSNP(ScanValsFile,LineNumber,filetype='2D'):
    snp=ssp.SterileNuParams()
    if filetype=='2D':
        # Key 0 is always the null hypothesis.
        if LineNumber==0:
            snp.del14=0
            snp.del24=0
            snp.dm41sq=0
            snp.modelId=0
            snp.th14=0
            snp.th24=0
            snp.th34=0
            return snp
        else:
            HypDict=Make2DHypDicts(ScanValsFile)
            [Log10Sin2, Log10Dm2]=HypDict['HypothesisNumberToParams'][LineNumber]
            snp.modelId=LineNumber
            snp.dm41sq=np.round(pow(10.,Log10Dm2),6)
            snp.th24=np.round(np.arcsin(pow(pow(10.,Log10Sin2),0.5))/2.,6)
            snp.del14=0
            snp.del24=0
            snp.th14=0
            snp.th34=0
            return snp
    else:
        #Here we will eventually extend to 2D file lists.
        print("Unsupported scan values file type", filetype)
    
