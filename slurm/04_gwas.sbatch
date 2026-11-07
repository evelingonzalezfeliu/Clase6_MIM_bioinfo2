#!/bin/bash
#SBATCH --job-name=04_gwas
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/04_gwas_%j.out
#SBATCH --error=logs/04_gwas_%j.err

# Rutas y semilla de la simulación
PROYECTO="$SLURM_SUBMIT_DIR"
PLINK_SIF="$PROYECTO/../gwas-class-containers/plink_1.90b7.7--h18e278d_1.sif"
GENOTIPOS="$PROYECTO/results/chr22_qc"
PCA="$PROYECTO/results/pca.eigenvec"
METADATOS="$PROYECTO/data/integrated_call_samples_v3.20130502.ALL.panel"
FENOTIPO="$PROYECTO/results/fenotipo.tsv"
SEMILLA=20261106

# Simulación reproducible del fenotipo
Rscript "$PROYECTO/R/00_simulate_phenotype.R" \
  "$METADATOS" "$GENOTIPOS.fam" "$FENOTIPO" "$SEMILLA"

# GWAS sin componentes principales
singularity exec --bind "$PROYECTO:$PROYECTO" "$PLINK_SIF" \
  plink \
  --bfile "$GENOTIPOS" \
  --pheno "$FENOTIPO" \
  --pheno-name PHENO \
  --logistic hide-covar \
  --allow-no-sex \
  --out "$PROYECTO/results/gwas_sin_pc"

# GWAS ajustado por PC1-PC5
singularity exec --bind "$PROYECTO:$PROYECTO" "$PLINK_SIF" \
  plink \
  --bfile "$GENOTIPOS" \
  --pheno "$FENOTIPO" \
  --pheno-name PHENO \
  --covar "$PCA" \
  --covar-number 1-5 \
  --logistic hide-covar \
  --allow-no-sex \
  --out "$PROYECTO/results/gwas_con_pc"
