#!/bin/bash
#SBATCH -J myjob           # Job name
#SBATCH -o myjob.out       # Name of stdout output file
#SBATCH -e myjob.err       # Name of stderr error file
#SBATCH -p normal         # Queue (partition) name
#SBATCH -N 1              # Total # of nodes (must be 1 for serial)
#SBATCH -n 1              # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 8:00:00       # Run time (hh:mm:ss)
#SBATCH --mail-user=dreambetter@gmail.com
#SBATCH --mail-type=all   # Send email at begin and end of job

export PATH=~/myTools/HiChIP/HiC-Pro/bin/utils/:$PATH
export PYTHONPATH=""
genomesize=/home1/04935/shaojf/myTools/HiChIP/HiC-Pro/annotation/chrom_hg19.sizes
# hicpro2higlass -i INPUT -r RESOLUTION -c CHROMSIZE [-n] [-h]
## Convert matrix file into .cool file
# HICPRO_PATH/bin/utils/hicpro2higlass.sh -i hic_results/matrix/dixon_2M/raw/1000000/dixon_2M_1000000.matrix -r 1000000 -c chrom_hg19.sizes -n

# mkdir higlass.mcool.files
# cd higlass.mcool.files
# ln -s ../HiCPro.output/*/hic_results/data/*/*_allValidPairs .
## Convert allValidPairs file into .cool file
# for f in HiCPro.output/*/hic_results/data/*/*_allValidPairs
source activate py36
for f in *_allValidPairs
do
	hicpro2higlass.sh -i $f -r 5000 -c $genomesize 1> $f.log 2>&1 &
	sleep 10m
done
source deactivate py36

# docker exec higlass-container python higlass-server/manage.py  ingest_tileset --filename /tmp/SRR5831497.mcool --datatype matrix --filetype cooler
for f in *_allValidPairs
do
	~/Documents/GitHub/HiC-Pro/bin/utils/hicpro2higlass.sh -i $f -r 5000 -c hg19.chrom.sizes
done
