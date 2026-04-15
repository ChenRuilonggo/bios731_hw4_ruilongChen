#!/bin/bash
#SBATCH --job-name=bios731_mixsim
#SBATCH --output=/projects/YangLabData/Ruilong/731/logs/mixsim_n${N}_%A_%a.out
#SBATCH --error=/projects/YangLabData/Ruilong/731/logs/mixsim_n${N}_%A_%a.err
#SBATCH --array=1-500
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --partition=day-long-cpu,month-long-cpu,yanglab
#SBATCH --mail-user=ruilong.chen@emory.edu
#SBATCH --mail-type=END,FAIL

mkdir -p /projects/YangLabData/Ruilong/731/logs
mkdir -p /projects/YangLabData/Ruilong/731/results/raw

module load R

if [ -z "${N}" ]; then
  echo "Error: N is not set."
  echo "Submit like: sbatch --export=ALL,N=100 /projects/YangLabData/Ruilong/731/slurm/run_simulation.sh"
  exit 1
fi

BATCH_SIZE=1
TASK_ID=${SLURM_ARRAY_TASK_ID}

echo "Running scenario: n=${N}, task_id=${TASK_ID}, batch_size=${BATCH_SIZE}"

Rscript /projects/YangLabData/Ruilong/731/run_simulation.R "${N}" "${TASK_ID}" "${BATCH_SIZE}"