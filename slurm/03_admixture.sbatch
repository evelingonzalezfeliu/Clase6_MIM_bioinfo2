#!/bin/bash
#SBATCH --job-name=03_admixture
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/03_admixture_%j.out
#SBATCH --error=logs/03_admixture_%j.err

# Rutas y número de componentes
PROYECTO="$SLURM_SUBMIT_DIR"
ADMIXTURE_SIF="$PROYECTO/../gwas-class-containers/admixture_1.3.0--0.sif"
RESULTADOS="$PROYECTO/results"
K=3

cd "$RESULTADOS"
cp chr22_pruned.fam "admixture_K${K}.fam"

# ADMIXTURE con K=3
singularity exec --bind "$PROYECTO:$PROYECTO" "$ADMIXTURE_SIF" \
  admixture --cv=10 -j"$SLURM_CPUS_PER_TASK" chr22_pruned.bed "$K"
