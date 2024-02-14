#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=10:mem=10gb
#PBS -N correlation_analysis
#PBS -J 0-1199

cd $PBS_O_WORKDIR

module load anaconda3/personal
source activate r413

n_cores=10
node=$PBS_ARRAY_INDEX

Rscript correlation_analysis.R $n_cores $node