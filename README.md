# Práctico de GWAS y ancestría

La guía completa del práctico está en [ACTIVIDAD.md](ACTIVIDAD.md).

La descripción de los comandos y archivos generados está en
[DETALLE_ANALISIS.md](DETALLE_ANALISIS.md).

## Flujo del práctico

```text
VCF del cromosoma 22 + información de las muestras
                         ↓
             Preparación previa con PLINK
                         ↓
          ┌──────────────┼──────────────┐
          ↓              ↓              ↓
         PCA         ADMIXTURE      Mini-GWAS
          ↓              ↓              ↓
       Gráfico         Gráfico      Manhattan y QQ
          └──────────────┴──────────────┘
                         ↓
              Interpretación en RStudio
```

La guía incluye la descarga de archivos, los cuatro trabajos Slurm y tres
scripts de R para PCA, ADMIXTURE, QC y mini-GWAS.
