#!/bin/bash                                                                                                                                                                   
#SBATCH -J MakingMu        # A single job name for the array                                                                                                                  
#SBATCH -c 1               # Number of cores                                                                                                                                  
#SBATCH -N 1                # Ensure that all cores are on one machine                                                                                                        
#SBATCH --mem 8000          #Memory request                                                                                                                                   
#SBATCH -t 0-04:00          # Runtime in D-HH:MM, minimum of 10 minutes                                                                                                       
#SBATCH -p guenette         # Partition to submit to                                                                                                                          
#SBATCH -o %A_%a.out                                                                                                                                                          
#SBATCH -e %A_%a.err                                                                                                                                                           

JOBNUMBER=${SLURM_ARRAY_TASK_ID}



#source ~/packages/nexus/setup_nexus.sh
#source /n/holystore01/LABS/guenette_lab/Lab/software/next/ups_products/setup 
#export LD_LIBRARY_PATH=/n/holystore01/LABS/guenette_lab/Lab/software/next/GATE/2.0/lib:$LD_LIBRARY_PATH 

#cd /n/holystore01/LABS/guenette_lab/Users/lrogers/JOB_CONTROL/JUNK 
/n/home10/lrogers/packages/nexus/build/source/nexus -b -n 300 /n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/data/SmallerSphere3/MCmuonsDiscrete/Config/Muons-0${JOBNUMBER}-MUONS.init.mac >& /scratch/Muons-0${JOBNUMBER}-MUONS.log
mv /scratch/Muons-0${JOBNUMBER}-MUONS.log /n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/data/SmallerSphere3/MCmuonsDiscrete/Logs/
mv /scratch/Muons-0${JOBNUMBER}-MUONS.h5 /n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/data/SmallerSphere3/MCmuonsDiscrete/Output/


