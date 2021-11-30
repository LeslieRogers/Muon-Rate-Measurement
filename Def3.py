import numpy as np

def GEN_CONFIGURATION(ijob, TemplateCONF, WorkType, MacrosDir, **kwargs):
    try:
        region = kwargs['region']
    except KeyError:
        print('Give me a generation region, please')
        raise

    Seed         = np.random.randint(1000000)
    Num          = str(int(ijob)).zfill(4)

    ConfigPath   = MacrosDir   + WorkType+'-'+Num+'-'+region+".config.mac"
    InitPath     = MacrosDir   + WorkType+'-'+Num+'-'+region+".init.mac"
    outputFile   = "/scratch/" + WorkType+'-'+Num+'-'+region

    #ConfigPath   = MacrosDir   + WorkType+'-'+Num+'-'+region+".config.mac"
    #InitPath     = MacrosDir   + WorkType+'-'+Num+'-'+region+".init.mac"
    #outputFile   = "/scratch/" + WorkType+'-'+Num+'-'+region
    with open(TemplateCONF) as config, open(ConfigPath, 'w') as configN:
        for line in config:
            if 'DIAMETER' in line:
                try:
                    line = line.replace('DIAMETER',
                                        str(kwargs['active_diam']) + '.')
                except KeyError:
                    print('Template requires active_diam')
                    raise
            elif 'LENGTH' in line:
                try:
                    line = line.replace('LENGTH',
                                        str(kwargs['active_length']) + '.')
                except KeyError:
                    print('Template requires active_length')
                    raise
            elif 'FC_THICKNESS' in line:
                try:
                    line = line.replace('FC_THICKNESS',
                                        str(kwargs['fcage_thickn']) + '.')
                except KeyError:
                    print('Template requires fcage_thickn')
                    raise
            elif 'ICS_THICKNESS' in line:
                try:
                    line = line.replace('ICS_THICKNESS',
                                        str(kwargs['ics_thickn']) + '.')
                except KeyError:
                    print('Template requires ics_thickn')
                    raise
            elif 'VESSEL_THICKNESS' in line:
                try:
                    line = line.replace('VESSEL_THICKNESS',
                                        str(kwargs['vessel_thickn']) + '.')
                except KeyError:
                    print('Template requires vessel_thickn')
                    raise
            elif 'TEMPERATURE' in line:
                try:
                    line = line.replace('TEMPERATURE',
                                        str(kwargs['gas_temperature']) + '.')
                except KeyError:
                    print('Template requires gas_temperature')
                    raise
            
            elif 'PRESSURE' in line:
                try:
                    line = line.replace('PRESSURE',
                                        str(kwargs['gas_pressure']) + '.')
                except KeyError:
                    print('Template requires gas_pressure')
                    raise

            elif 'REGION' in line:
                line = line.replace('REGION', region)
            elif 'ATOMIC_NUMBER' in line:
                try:
                    line = line.replace('ATOMIC_NUMBER',
                                        str(kwargs['atomic_number']))
                except KeyError:
                    print('Template requires atomic_number')
                    raise
            elif 'ATOMIC_MASS' in line:
                try:
                    line = line.replace('ATOMIC_MASS',
                                        str(kwargs['mass_number']))
                except KeyError:
                    print('Template requires mass_number')
                    raise
            elif 'DECAY_MODE' in line:
                try:
                    line = line.replace('DECAY_MODE',
                                        str(kwargs['Xe136DecayMode']))
                except KeyError:
                    print('Template requires Xe136DecayMode')
                    raise
            elif 'SEED' in line:
                try:
                    line = line.replace('SEED',
                                        str(kwargs['seed']) )
                except KeyError:
                    print('Template requires seed')
                    raise
            elif "EVENT_ID" in line:
                try:
                    line = line.replace("EVENT_ID",
                                        str(kwargs['event_id']) )
                except KeyError:
                    print('Template requires event_id')
                    raise

                #line = line.replace('SEED', str(Seed))

            elif 'THRESHOLD' in line:
                try:
                    line = line.replace('THRESHOLD',
                                        str(kwargs['threshold']))
                except KeyError:
                    print('Template requires threshold')
                    raise
            elif 'OUTPUT' in line:
                line = line.replace('OUTPUT', outputFile)
            elif 'MIN_ENG' in line:
                try:
                    line = line.replace('MIN_ENG',
                                        str(kwargs['min_eng']))
                except KeyError:
                    print('Template requires minimum_energy')
                    raise
            elif 'MAX_ENG' in line:
                try:
                    line = line.replace('MAX_ENG',
                                        str(kwargs['max_eng']))
                except KeyError:
                    print('Template requires maximum_energy')
                    raise
            configN.write(line)


def GEN_INITIALIZATION(ijob, TemplateINIT, WorkType, MacrosDir,region):
    Num          = str(int(ijob)).zfill(4)
    ConfigPath   = MacrosDir + WorkType+'-'+Num+'-'+region+'.config.mac'
    InitPath     = MacrosDir + WorkType+'-'+Num+'-'+region+'.init.mac'
    with open(TemplateINIT) as config, open(InitPath, 'w') as configN:
        for line in config:
            if 'CONFIG_PATH' in line:
                line = line.replace('CONFIG_PATH', ConfigPath)
            elif 'BiMAC' in line:
                delayed = ""
                if ("Bi214" in WorkType):
                    delayed = "/nexus/RegisterDelayedMacro /n/holystore01/LABS/guenette_lab/Lab/software/next/nexus/nexus-v5_07_00/macros/physics/Bi214.mac"
                line = line.replace('BiMAC', delayed)
            configN.write(line)


def ScriptGen(ijob, Nevents, WorkType, MacrosDir, OutputDir,
              ScriptDir, LogDir, NexusDir, **kwargs):
    try:
        region = kwargs['region']
    except KeyError:
        print('Give me a generation region, please')
        raise

    #SOURCE     = "source /data4/NEXT/sw/Releases/NEXT_v1_05_02/setup.sh"
    NEXUS      = NexusDir+" -b -n "
    Num        = str(int(ijob)).zfill(4)
    ConfigPath = MacrosDir + WorkType+'-'+Num+'-'+region+".config.mac"
    InitPath   = MacrosDir + WorkType+'-'+Num+'-'+region+".init.mac"
    outputName = WorkType+'-'+Num+'-'+region+".h5"
    outputFile = OutputDir + WorkType+'-'+Num+'-'+region
    LogName    = WorkType+'-'+Num+'-'+region+'.log'
    LogPath    = LogDir    + WorkType+'-'+Num+'-'+region+'.log'
    LogScratch = "/scratch/"+LogName
    ScriptPath = ScriptDir + WorkType+'-'+Num+'-'+region+'.sh'

    Command    = NEXUS+str(Nevents)+' '+InitPath+' >& '+LogScratch
    MEMS = "#SBATCH --mem=8000          # Memory pool for all cores (see also --mem-per-cpu)\n"
    if WorkType == "Muon.verticle" :
        MEMS = "#SBATCH --mem=8000          # Memory pool for all cores (see also --mem-per-cpu)\n"
    with open(ScriptPath, 'w') as jobF:
        jobF.write("#!/bin/bash\n")
        jobF.write("#SBATCH -J nexus        # A single job name for the array \n")
        jobF.write("#SBATCH -c 2                # Number of cores \n")
        jobF.write("#SBATCH -N 1                # Ensure that all cores are on one machine\n")
        jobF.write("#SBATCH --mem 8000          #Memory request\n")
        jobF.write("#SBATCH -t 0-06:00          # Runtime in D-HH:MM, minimum of 10 minutes\n")
        jobF.write("#SBATCH -p guenette         # Partition to submit to \n")
        jobF.write(MEMS)
        #jobF.write("#SBATCH -o %A_%a.out    # Standard output \n")
        #jobF.write("#SBATCH -e %A_%a.err    # Standard error \n")
        # Setup Nexus
        jobF.write("source ~/packages/nexus/setup_nexus.sh\n")
        jobF.write("source /n/holystore01/LABS/guenette_lab/Lab/software/next/ups_products/setup \n")
        
        # Configure scisoft softward products
        #jobF.write(". /n/holystore01/LABS/guenette_lab/Lab/software/next/scisoft/setup \n")
        # Setup some UPS products: ROOT and HDF5
        #jobF.write("setup cmake  v3_14_3 \n")
        #jobF.write("setup geant4 v4_10_6_p01 -q e19:prof \n")
        #jobF.write("setup gsl v2_5 -q prof \n")
        #jobF.write("setup root v6_18_04 -q e19:prof \n")
        #jobF.write("setup hdf5 v1_10_5 -q e19 \n")
        # Setup GATE 
        jobF.write("export LD_LIBRARY_PATH=/n/holystore01/LABS/guenette_lab/Lab/software/next/GATE/2.0/lib:$LD_LIBRARY_PATH \n")
        #jobF.write(SOURCE+" \n")
        jobF.write("\n")
        jobF.write("cd /n/holystore01/LABS/guenette_lab/Users/lrogers/JOB_CONTROL/JUNK \n")
        jobF.write(Command)
        jobF.write("\n")
        jobF.write("mv "+LogScratch+" "+LogDir)
        jobF.write("\n")
        jobF.write("mv "+"/scratch/"+outputName+" "+OutputDir)
        jobF.write("\n")
