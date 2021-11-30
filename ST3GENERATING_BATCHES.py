import os
import numpy as np
import random
from time      import time
from DefGeneral       import *

cwd = os.getcwd()


ICDir         = "/n/holystore01/LABS/guenette_lab/Users/lrogers/IC/"
NEW_PATH      = "/n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/data/SmallerSphere3/"
Size          = "NEXT_NEW/"

Regions = ["MUONS"]


Njobs = int(10)
Nevents  = int(300)

MacrosDir  =  NEW_PATH 
ScriptDir  =  NEW_PATH 


region="MUONS"
for i in range(100,200,1):

    eventnum     =   (i+1)*Nevents

    TemplateCONF  = cwd + "/detsim.conf"
    ICType        = "detsim"
    InputType     = "MCmuonsDiscrete"
    WorkType      = "Muons"
    OutputDir  =  NEW_PATH + ICType+"/Output/"
    LogDir     =  NEW_PATH + ICType+"/Logs/"


    GEN_CONFIGURATION(i, TemplateCONF, WorkType,ICType,InputType, MacrosDir, region)  
    
    ScriptGen(i, Nevents, WorkType,ICType,InputType, MacrosDir, OutputDir,
                  ScriptDir, LogDir, ICDir,region)

    '''TemplateCONF  = cwd + "/hypathia.conf"
    ICType        = "hypathia"
    InputType     = "detsim"
    WorkType      = "detsim"
    OutputDir  =  NEW_PATH + ICType+"/Output/"
    LogDir     =  NEW_PATH + ICType+"/Logs/"


    GEN_CONFIGURATION(i, TemplateCONF, WorkType,ICType,InputType, MacrosDir, region)
    
    ScriptGen2(i, Nevents, WorkType,ICType,InputType, MacrosDir, OutputDir,
                  ScriptDir, LogDir, ICDir,region)'''
    
    
    TemplateCONF  = cwd + "/diomira_0.conf"                                                                                            
    ICType        = "diomira"                                                                                                             
    InputType     = "detsim"                                                                                                                 
    WorkType      = "detsim"                                                                                                                 
    OutputDir  =  NEW_PATH + ICType+"/Output/"                                                                                              
    LogDir     =  NEW_PATH + ICType+"/Logs/"                                                                                              
    

    GEN_CONFIGURATION(i, TemplateCONF, WorkType,ICType,InputType, MacrosDir, region)                                                       
    
    ScriptGen2(i, Nevents, WorkType,ICType,InputType, MacrosDir, OutputDir,                                                          
                  ScriptDir, LogDir, ICDir,region)


    TemplateCONF  = cwd + "/irene2.conf"
    ICType        = "irene"
    InputType     = "diomira"
    WorkType      = "diomira"
    OutputDir  =  NEW_PATH + ICType+"/Output/"
    LogDir     =  NEW_PATH + ICType+"/Logs/"


    GEN_CONFIGURATION(i, TemplateCONF, WorkType,ICType,InputType, MacrosDir, region)

    ScriptGen2(i, Nevents, WorkType,ICType,InputType, MacrosDir, OutputDir,
                  ScriptDir, LogDir, ICDir,region)
    
    '''
    TemplateCONF  = cwd + "/penthesilea.conf"
    ICType        = "penthesilea"
    InputType     = "hypathia"
    WorkType      = "hypathia"
    OutputDir  =  NEW_PATH + ICType+"/Output/"
    LogDir     =  NEW_PATH + ICType+"/Logs/"
    '''



    TemplateCONF  = cwd + "/penthesilea2.conf"
    ICType        = "penthesilea"
    InputType     = "irene"
    WorkType      = "irene"
    OutputDir  =  NEW_PATH + ICType+"/Output/"
    LogDir     =  NEW_PATH + ICType+"/Logs/"


    GEN_CONFIGURATION(i, TemplateCONF, WorkType,ICType,InputType, MacrosDir, region)

    ScriptGen2(i, Nevents, WorkType,ICType,InputType, MacrosDir, OutputDir,
                  ScriptDir, LogDir, ICDir,region)



    TemplateCONF  = cwd + "/esmeralda.conf"
    ICType        = "esmeralda"
    InputType     = "penthesilea"
    WorkType      = "penthesilea"
    OutputDir  =  NEW_PATH + ICType+"/Output/"
    LogDir     =  NEW_PATH + ICType+"/Logs/"


    GEN_CONFIGURATION(i, TemplateCONF, WorkType,ICType,InputType, MacrosDir, region)
    
    ScriptGen2(i, Nevents, WorkType,ICType,InputType, MacrosDir, OutputDir,
                  ScriptDir, LogDir, ICDir,region)


