#!/bin/bash
#SBATCH -J Muons # A single job name for the array
#SBATCH -c 1 # Number of cores
#SBATCH -p shared # Partition
#SBATCH --mem 4000 # Memory request (6Gb)
#SBATCH -t 0-3:00 # Maximum execution time (D-HH:MM)
#SBATCH -o Muons_%A_%a.out # Standard output
#SBATCH -e Muons_%A_%a.err # Standard error

echo "Initialising NEXUS environment" 2>&1 | tee -a log_nexus_"${SLURM_ARRAY_TASK_ID}".txt
start=`date +%s`

# Create the directory
cd $SCRATCH/guenette_lab/Users/$USER/
mkdir -p Muons/jobid_"${SLURM_ARRAY_TASK_ID}"
cd Muons/jobid_"${SLURM_ARRAY_TASK_ID}"

# Copy the files over
cp ~/packages/MuonGenerator/macros_/NEXT100_muons_hallA.config.mac .
cp ~/packages/MuonGenerator/macros_/NEXT100_muons_hallA.init.mac .
cp ~/packages/nexus/macros/physics/Xe137.mac .
cp ~/packages/nexus/data/MuonAnaAllRuns.csv .

# Edit the file configs
sed -i "s#.*outputFile.*#/nexus/persistency/outputFile Next100Muons_hallA_example.next.h5#" NEXT100_muons_hallA.config.mac

# Setup nexus and run
echo "Setting Up NEXUS" 2>&1 | tee -a log_nexus_"${SLURM_ARRAY_TASK_ID}".txt
source ~/packages/nexus/setup_nexus.sh

# Also setup IC
source ~/packages/IC/setup_IC.sh

echo "Running NEXUS" 2>&1 | tee -a log_nexus_"${SLURM_ARRAY_TASK_ID}".txt
for i in {1..2}; do

	# Replace the seed in the file	
	echo "The seed number is: $((1111111*${SLURM_ARRAY_TASK_ID}+$i))" 2>&1 | tee -a log_nexus_"${SLURM_ARRAY_TASK_ID}".txt
	sed -i "s#.*random_seed.*#/nexus/random_seed $((1111111*${SLURM_ARRAY_TASK_ID}+$i))#" NEXT100_muons_hallA.config.mac
	
	nexus -n 5000 NEXT100_muons_hallA.init.mac 2>&1 | tee -a log_nexus_"${SLURM_ARRAY_TASK_ID}".txt

	# Rename the output file
	python ~/packages/MuonGenerator/slim_file.py "Next100Muons_hallA_example.next.h5" "$(basename Next100Muons_hallA_example.next.h5 .next.h5)_${SLURM_ARRAY_TASK_ID}_$i.next.h5"
	#mv Next100Muons_hallA_example.next.h5 "$(basename Next100Muons_hallA_example.next.h5 .next.h5)_${SLURM_ARRAY_TASK_ID}_$i.next.h5"
	echo; echo; echo;
done

# Cleaning up
rm -v Next100Muons_hallA_example.next.h5 2>&1 | tee -a log_nexus_"${SLURM_ARRAY_TASK_ID}".txt

echo "FINISHED....EXITING" 2>&1 | tee -a log_nexus_"${SLURM_ARRAY_TASK_ID}".txt
end=`date +%s`
runtime=$((end-start))
echo "$runtime s" 2>&1 | tee -a log_nexus_"${SLURM_ARRAY_TASK_ID}".txt
