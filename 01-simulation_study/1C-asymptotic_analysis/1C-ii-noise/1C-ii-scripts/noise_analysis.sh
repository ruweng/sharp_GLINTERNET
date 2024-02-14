#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=10:mem=10gb
#PBS -N noise_analysis
#PBS -J 0-1999

cd $PBS_O_WORKDIR

module load anaconda3/personal
source activate r413

n_cores=10
node=$PBS_ARRAY_INDEX

Rscript noise_analysis.R $n_cores $node