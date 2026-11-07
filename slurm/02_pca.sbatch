#!/bin/bash
#SBATCH --job-name=02_pca
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/02_pca_%j.out
#SBATCH --error=logs/02_pca_%j.err

# Rutas de entrada y salida
PROYECTO="$SLURM_SUBMIT_DIR"
PLINK_SIF="$PROYECTO/../gwas-class-containers/plink_1.90b7.7--h18e278d_1.sif"
ENTRADA="$PROYECTO/results/chr22_pruned"
SALIDA="$PROYECTO/results/pca"

# PCA con 10 componentes principales
singularity exec --bind "$PROYECTO:$PROYECTO" "$PLINK_SIF" \
  plink \
  --bfile "$ENTRADA" \
  --pca 10 header tabs \
  --out "$SALIDA"
