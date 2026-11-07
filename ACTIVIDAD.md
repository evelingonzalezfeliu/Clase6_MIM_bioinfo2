# Práctico de bioinformática: GWAS y ancestría

## Índice

1. [Descargar datos y contenedores](#parte-1-descargar-datos-y-contenedores)
2. [Ejecutar los cuatro trabajos Slurm](#parte-2-ejecutar-los-trabajos-slurm)
3. [Abrir RStudio y ejecutar los tres scripts de R](#parte-3-explorar-los-resultados-en-rstudio)
4. [Responder las preguntas a partir de los gráficos](#preguntas)

## Parte 1. Descargar datos y contenedores

En esta parte descargaremos los contenedores y los datos que utilizaremos
durante el práctico. Copie y pegue cada bloque de comandos en la terminal del
clúster.

Todavía no usaremos `sbatch`.

## Paso 1. Entrar al repositorio

Cambie `/RUTA/AL/REPOSITORIO` por la ruta donde clonó el repositorio:

```bash
cd /RUTA/AL/REPOSITORIO/Clase6_MIM_bioinfo2
```

## Paso 2. Descargar los tres contenedores

Utilizaremos imágenes de BioContainers:

- [PLINK 1.90b7.7](https://quay.io/repository/biocontainers/plink)
- [ADMIXTURE 1.3.0](https://quay.io/repository/biocontainers/admixture)
- [bcftools 1.22](https://quay.io/repository/biocontainers/bcftools)

Copie y pegue el bloque completo:

```bash
export CONTAINER_DIR="$PWD/../gwas-class-containers"
mkdir -p "$CONTAINER_DIR"

singularity pull \
  "$CONTAINER_DIR/plink_1.90b7.7--h18e278d_1.sif" \
  docker://quay.io/biocontainers/plink:1.90b7.7--h18e278d_1

singularity pull \
  "$CONTAINER_DIR/admixture_1.3.0--0.sif" \
  docker://quay.io/biocontainers/admixture:1.3.0--0

singularity pull \
  "$CONTAINER_DIR/bcftools_1.22--h3a4d415_1.sif" \
  docker://quay.io/biocontainers/bcftools:1.22--h3a4d415_1
```

Los contenedores quedarán fuera del repositorio, en la carpeta
`gwas-class-containers`.

## Paso 3. Descargar los datos

Usaremos el cromosoma 22 de 1000 Genomes Project, fase 3, en GRCh37. Los
archivos originales están disponibles en el
[sitio de 1000 Genomes](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/).

Copie y pegue el bloque completo:

```bash
mkdir -p data

curl --fail --location --retry 3 --continue-at - \
  --output data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
  https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

curl --fail --location --retry 3 --continue-at - \
  --output data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi \
  https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi

curl --fail --location --retry 3 \
  --output data/integrated_call_samples_v3.20130502.ALL.panel \
  https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

curl --fail --location --retry 3 \
  --output data/20140625_related_individuals.txt \
  https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/20140625_related_individuals.txt
```

Ya están descargados los archivos que se utilizarán en las siguientes partes
del práctico.

---

## Parte 2. Ejecutar los trabajos Slurm

En esta parte seguiremos este flujo:

```text
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

Ejecute los trabajos desde la carpeta principal del repositorio. Espere que
termine cada trabajo antes de enviar el siguiente.

### Paso 1. Preparación con PLINK

Selecciona muestras CEU, YRI, CHB y PEL, convierte el VCF a formato PLINK,
aplica filtros básicos y realiza LD pruning. Genera `chr22_qc` para el
mini-GWAS y `chr22_pruned` para PCA y ADMIXTURE.

```bash
sbatch slurm/01_plink.sbatch
```

### Paso 2. PCA

Calcula diez componentes principales usando los SNP después del LD pruning.

```bash
sbatch slurm/02_pca.sbatch
```

### Paso 3. ADMIXTURE

Estima tres componentes (`K=3`) usando los mismos SNP empleados en el PCA.

```bash
sbatch slurm/03_admixture.sbatch
```

### Paso 4. Mini-GWAS

Simula un fenotipo caso-control con una semilla fija y ejecuta dos modelos de
asociación: uno sin componentes principales y otro ajustado por PC1-PC5.

```bash
sbatch slurm/04_gwas.sbatch
```

Al finalizar, los archivos necesarios para hacer los gráficos estarán en la
carpeta `results/`. Los utilizaremos en la siguiente parte desde RStudio.

---

## Parte 3. Explorar los resultados en RStudio

Abra RStudio desde Open OnDemand. En la consola de R, indique la ruta del
repositorio:

```r
setwd("/RUTA/AL/REPOSITORIO/Clase6_MIM_bioinfo2")
```

### Paso 1. PCA

Este script muestra PC1 frente a PC2, colorea las muestras por población y
guarda la figura como `results/PCA.png`.

```r
source("R/01_pca.R")
```

### Paso 2. ADMIXTURE

Este script representa los tres componentes estimados para cada individuo y
guarda la figura como `results/ADMIXTURE_K3.png`.

```r
source("R/02_admixture.R")
```

Los componentes 1, 2 y 3 no deben interpretarse automáticamente como
ancestrías continentales específicas.

### Paso 3. QC y mini-GWAS

Este script muestra el QC de muestras y variantes. Después compara los
resultados de asociación sin PC y con PC mediante gráficos de asociación y QQ.
Guarda `results/QC.png` y `results/GWAS_Manhattan_QQ.png`.

```r
source("R/03_gwas_qc.R")
```

Los gráficos de asociación contienen solamente el cromosoma 22; por lo tanto,
no son gráficos Manhattan de todo el genoma.

---

## Preguntas

### PCA

Observe `results/PCA.png`:

1. ¿Las cuatro poblaciones forman grupos separados o se superponen?
2. ¿Por qué un PCA basado solamente en el cromosoma 22 puede ser menos estable
   que un PCA de todo el genoma?

### ADMIXTURE

Observe `results/ADMIXTURE_K3.png`:

1. ¿Qué diferencias observa entre las proporciones de las cuatro poblaciones?
2. ¿Por qué no debemos asignar automáticamente un continente a cada color?

### QC y mini-GWAS

Observe `results/QC.png` y `results/GWAS_Manhattan_QQ.png`:

1. ¿Qué muestran las distribuciones de MAF y missingness después del QC?
2. ¿Cómo cambia el gráfico QQ al incluir PC1-PC5 en el modelo?
3. ¿Por qué no podemos obtener conclusiones sobre todo el genoma usando
   solamente el cromosoma 22?
