import numpy as np

def GEN_CONFIGURATION(ijob, TemplateCONF, WorkType,ICtype,InputType, MacrosDir, region):

    #ICtype       = 'diomira'
    #InputType    = 'detsim'
    Seed         = np.random.randint(1000000)
    Num          = str(int(ijob)).zfill(4)

    #ConfigPath   = MacrosDir +ICtype+'Discrete/Config/'+ ICtype+'-'+Num+'.conf'
    ConfigPath   = MacrosDir +ICtype+'/Config/'+ ICtype+'-'+Num+'.conf'
    inputFile    = MacrosDir +InputType+'/Output/'  + WorkType+ '-' +Num+'-'+region+'.h5'
    outputFile   = "/scratch/" + ICtype+'-'+Num+'-'+region+'.h5'

    with open(TemplateCONF) as config, open(ConfigPath, 'w') as configN:
        for line in config:
            if 'OUTPUT' in line:
                line = line.replace('OUTPUT', outputFile)
            elif 'INPUT' in line:
                line = line.replace('INPUT', inputFile)
            elif 'RUNNUM' in line:
                line = line.replace('RUNNUM', str(ijob))

            configN.write(line)
def ScriptGen(ijob, Nevents, WorkType,ICtype,InputType, MacrosDir, OutputDir,
              ScriptDir, LogDir, ICDir,region):

    #SOURCE     = "source /data4/NEXT/sw/Releases/NEXT_v1_05_02/setup.sh"
    IC         = ICDir
    Num        = str(int(ijob)).zfill(4)
    #ConfigPath   = MacrosDir +ICtype+'Discrete/Config/'+ ICtype+'-'+Num+'.conf'
    ConfigPath   = MacrosDir +ICtype+'/Config/'+ ICtype+'-'+Num+'.conf'
    outputName = ICtype+'-'+Num+'-'+region+'.h5'
    #outputFile = MacrosDir +ICtype+'Discrete/Output/'+ ICtype+'-'+Num
    outputFile = MacrosDir +ICtype+'/Output/'+ ICtype+'-'+Num
    LogName    = WorkType+'-'+Num+'.log'
    #LogPath    = MacrosDir + ICtype+'Discrete/Logs/'    + ICtype+'-'+Num+'.log'
    LogPath    = MacrosDir + ICtype+'/Logs/'    + ICtype+'-'+Num+'.log'
    LogScratch = "/scratch/"+LogName
    #ScriptPath = MacrosDir + ICtype+'Discrete/Scripts/'+ICtype+'-'+Num+'.sh'
    #ScriptPath = MacrosDir + ICtype+'/Scripts/'+ICtype+'-'+Num+'.sh'
    ScriptPath = MacrosDir + 'AllCities/Scripts/Allcities-'+Num+'.sh'
    Command    = 'city '+ICtype+' '+ConfigPath+' >& '+LogScratch
    MEMS = "#SBATCH --mem=10000          # Memory pool for all cores (see also --mem-per-cpu)\n"

    #print(ScriptPath)
    with open(ScriptPath, 'w') as jobF:
        jobF.write("#!/bin/bash\n")
        jobF.write("#SBATCH -J "+ ICtype+"       # A single job name for the array \n")
        jobF.write("#SBATCH -n 2                # Number of cores \n")
        jobF.write("#SBATCH -N 1                # Ensure that all cores are on one machine\n")
        jobF.write("#SBATCH -t 0-05:00          # Runtime in D-HH:MM, minimum of 10 minutes\n")
        jobF.write("#SBATCH -p guenette         # Partition to submit to \n")
        jobF.write(MEMS)
        #jobF.write("#SBATCH -o %A_%a.out    # Standard output \n")
        #jobF.write("#SBATCH -e %A_%a.err    # Standard error \n")



        # Setup UPS
        jobF.write("source /n/holystore01/LABS/guenette_lab/Lab/software/next/ups_products/setup \n")

        # Configure scisoft softward products
        #jobF.write(". /n/holystore01/LABS/guenette_lab/Lab/software/next/scisoft/setup \n")

        #Setup Conda
        #jobF.write(". /n/holystore01/LABS/guenette_lab/Lab/software/next/miniconda/etc/profile.d/conda.sh \n")
        #jobF.write("export ICTDIR=/n/holystore01/LABS/guenette_lab/Users/lrogers/IC \n")
        #jobF.write("export ICDIR=$ICTDIR/invisible_cities \n")
        #jobF.write("export PATH=$ICTDIR/bin:$PATH \n")
        #jobF.write("export PYTHONPATH=$ICTDIR:$PYTHONPATH \n")

        # Setup IC
        #jobF.write("conda activate IC-3.7-2020-06-16 \n")

        jobF.write("source /n/holystore01/LABS/guenette_lab/Users/lrogers/initstuff.sh \n")
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



def ScriptGen2(ijob, Nevents, WorkType,ICtype,InputType, MacrosDir, OutputDir,
              ScriptDir, LogDir, ICDir,region):


    Num        = str(int(ijob)).zfill(4)
    ConfigPath   = MacrosDir +ICtype+'/Config/'+ ICtype+'-'+Num+'.conf'
    outputName = ICtype+'-'+Num+'-'+region+'.h5'
    LogName    = WorkType+'-'+Num+'.log'
    LogPath    = MacrosDir + ICtype+'/Logs/'    + ICtype+'-'+Num+'.log'
    LogScratch = "/scratch/"+LogName
    Command    = 'city '+ICtype+' '+ConfigPath+' >& '+LogScratch

    ScriptPath = MacrosDir + 'AllCities/Scripts/Allcities-'+Num+'.sh'
    with open(ScriptPath, 'a') as jobF:
        jobF.write("\n")
        jobF.write(Command)
        jobF.write("\n")
        jobF.write("mv "+LogScratch+" "+LogDir)
        jobF.write("\n")
        jobF.write("mv "+"/scratch/"+outputName+" "+OutputDir)
        jobF.write("\n")
