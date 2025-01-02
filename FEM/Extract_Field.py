# -*- coding: utf-8 -*-
from visualization import *
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from connectorBehavior import *
from abaqus import *
from abaqusConstants import *
import random
import time
import numpy as np
from odbAccess import *
from abaqusConstants import *
import string
import scipy.io as scio
import numpy as np

# scipy.io.savemat('file_name.mat', {'data':my_data})
jobname=getInput("input the inp file name without the .inp:")
jobname1='Label.mat'
jobnameS11='S11.mat'
jobnameS22='S22.mat'
jobnameS33='S33.mat'
jobnameS12='S12.mat'
jobnameE11='E11.mat'
jobnameE22='E22.mat'
jobnameE33='E33.mat'
jobnameE12='E12.mat'

jobname=jobname+'.odb'
odb=openOdb(jobname)
assembly = odb.rootAssembly
instance1 = odb.rootAssembly.instances['A-1']

EOUT = instance1
f=open(jobname1,'w')

frame = odb.steps['Step-1'].frames[-1]
allFields = frame.fieldOutputs
stressSet1 = allFields['S'] 
TstressSet = stressSet1.getSubset(region=EOUT)
stressS = TstressSet.values

strainSet1 = allFields['E'] 
TstrainSet = strainSet1.getSubset(region=EOUT)
strainE = TstrainSet.values

num_data = len(stressS)
All_S11 = []
All_S22 = []
All_S33 = []
All_S12 = []

All_E11 = []
All_E22 = []
All_E33 = []
All_E12 = []

All_L = []
for n_data in range(num_data):
	All_L.append(stressS[n_data].elementLabel)
	All_S11.append(stressS[n_data].data[0])
	All_S22.append(stressS[n_data].data[1])
	All_S33.append(stressS[n_data].data[2])
	All_S12.append(stressS[n_data].data[3])
	All_E11.append(strainE[n_data].data[0])
	All_E22.append(strainE[n_data].data[1])
	All_E33.append(strainE[n_data].data[2])
	All_E12.append(strainE[n_data].data[3])

All_S11=np.array(All_S11)
All_S22=np.array(All_S22)
All_S33=np.array(All_S33)
All_S12=np.array(All_S12)
All_E11=np.array(All_E11)
All_E22=np.array(All_E22)
All_E33=np.array(All_E33)
All_E12=np.array(All_E12)
All_L=np.array(All_L)

scio.savemat(jobname1, {'data':All_L})

scio.savemat(jobnameS11, {'data':All_S11})
scio.savemat(jobnameS22, {'data':All_S22})
scio.savemat(jobnameS33, {'data':All_S33})
scio.savemat(jobnameS12, {'data':All_S12})

scio.savemat(jobnameE11, {'data':All_E11})
scio.savemat(jobnameE22, {'data':All_E22})
scio.savemat(jobnameE33, {'data':All_E33})
scio.savemat(jobnameE12, {'data':All_E12})

# odb.close()