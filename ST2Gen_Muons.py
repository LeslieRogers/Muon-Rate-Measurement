import os
import numpy as np
import random
from time      import time
from Def3       import *

cwd = os.getcwd()

TemplateCONF  = cwd + "/ST1NEW-Mu_0.config.mac"
TemplateINIT  = cwd + "/NEW-Muons_0_t2.init.mac"

WorkType      = "Muon.verticle"
WorkType      = "Muons"
NEW_PATH      = "/n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/data/SmallerSphere3/MCmuonsDiscrete"
Size          = "NEXT_NEW/"
#Sim           = "MC.nexus_p5_05_00/"

#NexusDir      = "/n/holystore01/LABS/guenette_lab/Lab/software/next/nexus/nexus-v6_01_01/build/source/nexus"
NexusDir      = "/n/home10/lrogers/packages/nexus/build/source/nexus"

## Put the config specific variables in a dictionary
## That way you can only put the ones you need
config_vars = {}
config_vars['active_diam']    = 260
config_vars['active_length']  = 260
config_vars['fcage_thickn']   = 1
config_vars['ics_thickn']     = 12
config_vars['vessel_thickn']  = 2
config_vars['gas_temperature'] = 300
config_vars['gas_pressure']    = 10.1
config_vars['Xe136DecayMode']  = 1 # Decay0 interface for BB decays ... (BB0nu: DecayMode 1), (BB2nu: DecayMode 4)
config_vars['threshold']      = 2.3
config_vars['min_eng']        = 100
config_vars['max_eng']        = 100
region="MUONS"
#config_vars['region'] = MUONS
grp_ran=4
SEED = random.sample(range(10, 10000000), k=100005)
seed_index = 0

Energies = [.1,1,10,50,100,200,400,500]
#for nrg in Energies:
    


Nevents  = int(300)
    
MacrosDir  =  NEW_PATH + "/Config/"
OutputDir  =  NEW_PATH + "/Output/"
ScriptDir  =  NEW_PATH + "/Scripts/"
LogDir     =  NEW_PATH + "/Logs/"

config_vars['min_eng']        = 1 
config_vars['max_eng']        = 1     
config_vars['region']         ="MUONS"
#print("Generating the files for " + str(nrg))
for i in range(100,200):
    seed_index += 1
    config_vars['seed']     = SEED[seed_index]
    eventnum     =   (i+1)*Nevents
    config_vars['event_id'] = eventnum


    GEN_CONFIGURATION(i, TemplateCONF, WorkType, MacrosDir, **config_vars)
    
    GEN_INITIALIZATION(i, TemplateINIT, WorkType, MacrosDir, region)
    
    ScriptGen(i, Nevents, WorkType, MacrosDir, OutputDir,
                  ScriptDir, LogDir, NexusDir, **config_vars)
    grp_ran=grp_ran+1
    #submitJobs(WorkType)
