#!/bin/bash
#SBATCH -J detsim       # A single job name for the array 
#SBATCH -n 2                # Number of cores 
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 0-05:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p guenette         # Partition to submit to 
#SBATCH --mem=18000          # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o %A_%a.out                                                                                                                                                  
#SBATCH -e %A_%a.err                                                                                                                                               

JOBNUMBER=${SLURM_ARRAY_TASK_ID}

source /n/holystore01/LABS/guenette_lab/Lab/software/next/ups_products/setup
source /n/holystore01/LABS/guenette_lab/Users/lrogers/initstuff.sh

#source /n/holystore01/LABS/guenette_lab/Lab/software/next/ups_products/setup 
#source /n/holystore01/LABS/guenette_lab/Users/lrogers/initstuff.sh 
#export LD_LIBRARY_PATH=/n/holystore01/LABS/guenette_lab/Lab/software/next/GATE/2.0/lib:$LD_LIBRARY_PATH 

#cd /n/holystore01/LABS/guenette_lab/Users/lrogers/JOB_CONTROL/JUNK 
city detsim /n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/data/SmallerSphere2/detsim/Config/detsim-${JOBNUMBER}.conf >& /scratch/Muons-${JOBNUMBER}.log
mv /scratch/Muons-${JOBNUMBER}.log /n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/data/SmallerSphere2/detsim/Logs/
mv /scratch/detsim-${JOBNUMBER}-MUONS.h5 /n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/data/SmallerSphere2/detsim/Output/

city diomira /n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/data/SmallerSphere2/diomira/Config/diomira-${JOBNUMBER}.conf >& /scratch/detsim-${JOBNUMBER}.log
mv /scratch/detsim-${JOBNUMBER}.log /n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/data/SmallerSphere2/diomira/Logs/
mv /scratch/diomira-${JOBNUMBER}-MUONS.h5 /n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/data/SmallerSphere2/diomira/Output/

city irene /n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/data/SmallerSphere2/irene/Config/irene-${JOBNUMBER}.conf >& /scratch/diomira-${JOBNUMBER}.log
mv /scratch/diomira-${JOBNUMBER}.log /n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/data/SmallerSphere2/irene/Logs/
mv /scratch/irene-${JOBNUMBER}-MUONS.h5 /n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/data/SmallerSphere2/irene/Output/

city penthesilea /n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/data/SmallerSphere2/penthesilea/Config/penthesilea-${JOBNUMBER}.conf >& /scratch/irene-${JOBNUMBER}.log
mv /scratch/irene-${JOBNUMBER}.log /n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/data/SmallerSphere2/penthesilea/Logs/
mv /scratch/penthesilea-${JOBNUMBER}-MUONS.h5 /n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/data/SmallerSphere2/penthesilea/Output/

city esmeralda /n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/data/SmallerSphere2/esmeralda/Config/esmeralda-${JOBNUMBER}.conf >& /scratch/penthesilea-${JOBNUMBER}.log
mv /scratch/penthesilea-${JOBNUMBER}.log /n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/data/SmallerSphere2/esmeralda/Logs/
mv /scratch/esmeralda-${JOBNUMBER}-MUONS.h5 /n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/data/SmallerSphere2/esmeralda/Output/
