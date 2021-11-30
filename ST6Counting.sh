#!/bin/bash                                                                                  
#SBATCH -J Counting        # A single job name for the array
#SBATCH -c 1               # Number of cores                             
#SBATCH -N 1                # Ensure that all cores are on one machine                
#SBATCH --mem 10000          #Memory request                                                  
#SBATCH -t 0-01:00          # Runtime in D-HH:MM, minimum of 10 minutes                        
#SBATCH -p guenette         # Partition to submit to                                              
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

JOBNUMBER=${SLURM_ARRAY_TASK_ID}
#JOBNUMBER=1

#. /n/holystore01/LABS/guenette_lab/Lab/software/next/miniconda/etc/profile.d/conda.sh
#Setup IC                                                                        
#source /n/holystore01/LABS/guenette_lab/Users/lrogers/initstuff.sh

#python /n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/MuonSelctionForTruths.py 
python /n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/MuonSelctionForTruths2.py ${JOBNUMBER}
#python /n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/FindStartPointsetc.py  ${JOBNUMBER}
#python $script $p1 $p2 $min $max > ${p2}log_${min}_${max}.txt
