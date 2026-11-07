#!/bin/bash
#SBATCH --job-name=01_plink
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/01_plink_%j.out
#SBATCH --error=logs/01_plink_%j.err

set -euo pipefail

# Rutas de entrada y salida
PROYECTO="$SLURM_SUBMIT_DIR"
VCF="$PROYECTO/data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
PANEL="$PROYECTO/data/integrated_call_samples_v3.20130502.ALL.panel"
PLINK_SIF="$PROYECTO/../gwas-class-containers/plink_1.90b7.7--h18e278d_1.sif"
SALIDA="$PROYECTO/results/chr22_qc"

mkdir -p "$PROYECTO/results"

# Selección de cuatro poblaciones
awk 'BEGIN {OFS="\t"} NR > 1 && ($2 == "CEU" || $2 == "YRI" || $2 == "CHB" || $2 == "PEL") {print $1, $1}' \
  "$PANEL" > "$PROYECTO/results/muestras.keep"

# Preparación y control de calidad con PLINK
singularity exec --bind "$PROYECTO:$PROYECTO" "$PLINK_SIF" \
  plink \
  --vcf "$VCF" \
  --double-id \
  --keep "$PROYECTO/results/muestras.keep" \
  --snps-only just-acgt \
  --biallelic-only strict \
  --set-missing-var-ids '@:#:$1:$2' \
  --geno 0.02 \
  --mind 0.02 \
  --maf 0.05 \
  --make-bed \
  --out "$SALIDA"

# LD pruning para PCA y ADMIXTURE
singularity exec --bind "$PROYECTO:$PROYECTO" "$PLINK_SIF" \
  plink \
  --bfile "$SALIDA" \
  --indep-pairwise 50 5 0.2 \
  --out "$PROYECTO/results/ld"

singularity exec --bind "$PROYECTO:$PROYECTO" "$PLINK_SIF" \
  plink \
  --bfile "$SALIDA" \
  --extract "$PROYECTO/results/ld.prune.in" \
  --make-bed \
  --out "$PROYECTO/results/chr22_pruned"

# Tablas para los gráficos de QC en R
singularity exec --bind "$PROYECTO:$PROYECTO" "$PLINK_SIF" \
  plink \
  --bfile "$SALIDA" \
  --freq \
  --out "$PROYECTO/results/qc"

singularity exec --bind "$PROYECTO:$PROYECTO" "$PLINK_SIF" \
  plink \
  --bfile "$SALIDA" \
  --missing \
  --out "$PROYECTO/results/qc"
