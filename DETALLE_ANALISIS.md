# Guía de comandos y archivos generados

Esta guía explica qué hace cada trabajo Slurm, cuáles son sus entradas y qué
significa cada archivo creado en `results/`.

## Flujo de archivos

```text
VCF + panel de muestras
        │
        ▼
01_plink.sbatch
        ├── chr22_qc.* ────────────────► 04_gwas.sbatch
        │                                  ├── fenotipo.tsv
        │                                  ├── gwas_sin_pc.*
        │                                  └── gwas_con_pc.*
        │
        └── chr22_pruned.*
                 ├─────────────────────► 02_pca.sbatch
                 │                          ├── pca.eigenvec
                 │                          └── pca.eigenval
                 │
                 └─────────────────────► 03_admixture.sbatch
                                            ├── chr22_pruned.3.Q
                                            └── chr22_pruned.3.P
```

Los archivos con el mismo prefijo pertenecen al mismo conjunto de datos. Por
ejemplo, `chr22_qc.bed`, `chr22_qc.bim` y `chr22_qc.fam` deben permanecer
juntos.

## Árbol resumido de `results/`

```text
results/
├── muestras.keep                  # Muestras CEU, YRI, CHB y PEL
├── chr22_qc.bed/.bim/.fam         # Genotipos después del QC
├── ld.prune.in                    # SNP conservados después del LD pruning
├── ld.prune.out                   # SNP descartados por LD
├── chr22_pruned.bed/.bim/.fam     # Genotipos para PCA y ADMIXTURE
├── qc.frq                         # Frecuencias y MAF
├── qc.imiss                       # Missingness por individuo
├── qc.lmiss                       # Missingness por variante
├── pca.eigenvec                   # PC1-PC10 por individuo
├── pca.eigenval                   # Eigenvalores del PCA
├── admixture_K3.fam               # Orden de individuos del archivo Q
├── chr22_pruned.3.Q               # Proporciones de los tres componentes
├── chr22_pruned.3.P               # Frecuencias alélicas estimadas
├── fenotipo.tsv                   # Fenotipo caso-control simulado
├── gwas_sin_pc.assoc.logistic     # Asociación sin ajuste por PC
├── gwas_con_pc.assoc.logistic     # Asociación ajustada por PC1-PC5
├── *.log                          # Registros de PLINK
└── *.nosex                        # Muestras sin sexo definido en PLINK
```

Después de ejecutar los scripts de R también aparecerán `PCA.png`,
`ADMIXTURE_K3.png`, `QC.png` y `GWAS_Manhattan_QQ.png`.

---

## 1. Preparación y QC con PLINK

Archivo Slurm: `slurm/01_plink.sbatch`

### Entradas

- VCF del cromosoma 22: genotipos originales de 1000 Genomes.
- Panel de muestras: indica la población de cada individuo.
- Contenedor de PLINK 1.90b7.7.

### Selección de muestras

El comando `awk` lee el panel y selecciona individuos de cuatro poblaciones:

- `CEU`: población de Utah con ancestría del norte y oeste de Europa.
- `YRI`: Yoruba de Ibadan, Nigeria.
- `CHB`: Han de Beijing, China.
- `PEL`: población peruana de Lima.

La selección se guarda en `muestras.keep`. Sus dos columnas corresponden a
`FID` e `IID`, los identificadores que utiliza PLINK.

### Conversión y filtros

El comando principal es:

```bash
plink \
  --vcf VCF \
  --double-id \
  --keep muestras.keep \
  --snps-only just-acgt \
  --biallelic-only strict \
  --set-missing-var-ids '@:#:$1:$2' \
  --geno 0.02 \
  --mind 0.02 \
  --maf 0.05 \
  --make-bed \
  --out chr22_qc
```

Significado de los parámetros:

- `--vcf`: indica el archivo de genotipos de entrada.
- `--double-id`: usa el ID de la muestra como `FID` e `IID`.
- `--keep`: conserva solamente las muestras seleccionadas.
- `--snps-only just-acgt`: conserva SNP con alelos A, C, G o T.
- `--biallelic-only strict`: conserva variantes con exactamente dos alelos.
- `--set-missing-var-ids`: asigna un ID a las variantes que no lo tienen,
  usando cromosoma, posición y alelos.
- `--geno 0.02`: elimina variantes con más de 2 % de genotipos ausentes.
- `--mind 0.02`: elimina individuos con más de 2 % de genotipos ausentes.
- `--maf 0.05`: conserva variantes con frecuencia del alelo menor de al menos
  5 % en el conjunto analizado.
- `--make-bed`: genera el formato binario de PLINK.
- `--out`: define el prefijo de los archivos de salida.

### Salidas `chr22_qc`

- `chr22_qc.bed`: matriz binaria de genotipos. No se abre como texto.
- `chr22_qc.bim`: información de las variantes. Incluye cromosoma, ID,
  posición y alelos.
- `chr22_qc.fam`: información y orden de las muestras.
- `chr22_qc.log`: comando ejecutado, número de muestras, variantes y filtros.
- `chr22_qc.nosex`: muestras cuyo sexo no está definido en el archivo PLINK.
  Su presencia no significa que el trabajo haya fallado.

El conjunto `chr22_qc` conserva las variantes después del QC y será utilizado
en el mini-GWAS.

### LD pruning

El comando de poda por desequilibrio de ligamiento es:

```bash
plink --bfile chr22_qc --indep-pairwise 50 5 0.2 --out ld
```

Los tres números significan:

- `50`: ventana de 50 variantes.
- `5`: avance de 5 variantes entre ventanas.
- `0.2`: umbral de correlación `r²`; se retiran variantes redundantes cuando
  superan este valor.

Salidas:

- `ld.prune.in`: variantes relativamente independientes que se conservan.
- `ld.prune.out`: variantes descartadas por LD.
- `ld.log`: resumen de la ejecución.
- `ld.nosex`: archivo auxiliar de sexo desconocido.

Las variantes de `ld.prune.in` se extraen para crear:

- `chr22_pruned.bed`
- `chr22_pruned.bim`
- `chr22_pruned.fam`
- `chr22_pruned.log`
- `chr22_pruned.nosex`

Este conjunto reducido se utiliza en PCA y ADMIXTURE. La poda evita que grupos
de SNP muy correlacionados tengan un peso excesivo.

### Tablas para QC

PLINK también ejecuta `--freq` y `--missing`:

- `qc.frq`: frecuencia de los alelos y MAF de cada variante.
- `qc.imiss`: missingness por individuo; `F_MISS` es la fracción ausente.
- `qc.lmiss`: missingness por variante; `F_MISS` es la fracción ausente.
- `qc.log`: registro del comando de QC.
- `qc.nosex`: muestras con sexo no definido en PLINK.

Estas tablas son leídas por `R/03_gwas_qc.R` para crear `QC.png`.

### Qué deberíamos observar

- Los tres archivos `chr22_qc.bed/.bim/.fam` deben existir juntos.
- `chr22_pruned.bim` debe contener menos variantes que `chr22_qc.bim`.
- Los valores de `F_MISS` corresponden al conjunto que ya pasó los filtros.
- La MAF de las variantes conservadas debería ser igual o superior a 0.05,
  salvo pequeñas diferencias de representación o redondeo.

---

## 2. PCA

Archivo Slurm: `slurm/02_pca.sbatch`

### Entrada

El PCA utiliza `chr22_pruned.bed/.bim/.fam`, es decir, los SNP después del QC
y del LD pruning.

### Comando principal

```bash
plink --bfile chr22_pruned --pca 10 header tabs --out pca
```

Parámetros:

- `--bfile chr22_pruned`: lee el conjunto binario de PLINK.
- `--pca 10`: calcula diez componentes principales.
- `header`: agrega nombres de columnas.
- `tabs`: separa las columnas mediante tabulaciones.

### Salidas

- `pca.eigenvec`: una fila por muestra, con `FID`, `IID` y PC1-PC10.
- `pca.eigenval`: un eigenvalor para cada componente calculado.
- `pca.log`: resumen del análisis.
- `pca.nosex`: archivo auxiliar de sexo desconocido.

`R/01_pca.R` combina `pca.eigenvec` con el panel de muestras y genera
`PCA.png`.

### Qué deberíamos observar

En el gráfico PC1 versus PC2 pueden aparecer agrupamientos relacionados con la
estructura genética de CEU, YRI, CHB y PEL. La separación exacta debe
observarse en el resultado real, no asumirse de antemano.

Este PCA utiliza solamente el cromosoma 22. Por ello contiene menos marcadores
y puede ser menos estable que un PCA de todo el genoma.

---

## 3. ADMIXTURE

Archivo Slurm: `slurm/03_admixture.sbatch`

### Entrada

ADMIXTURE utiliza `chr22_pruned.bed/.bim/.fam`, el mismo conjunto podado por LD
que se utiliza para PCA.

Antes de ejecutar el programa se copia `chr22_pruned.fam` como
`admixture_K3.fam`. Esta copia conserva el orden exacto de las muestras que
corresponde a las filas del archivo `.Q`.

### Comando principal

```bash
admixture --cv=10 -j4 chr22_pruned.bed 3
```

Parámetros:

- `--cv=10`: calcula el error de validación cruzada usando diez particiones.
- `-j4`: utiliza cuatro CPU.
- `chr22_pruned.bed`: genotipos de entrada.
- `3`: ajusta un modelo con `K=3` componentes.

### Salidas

- `chr22_pruned.3.Q`: una fila por individuo y tres columnas. Cada valor es la
  proporción estimada para un componente; las tres proporciones de una fila
  deberían sumar aproximadamente 1.
- `chr22_pruned.3.P`: frecuencias alélicas estimadas para las variantes en los
  tres componentes.
- `admixture_K3.fam`: orden e identidad de las filas del archivo `.Q`.
- `logs/03_admixture_JOBID.out`: progreso del programa y error de validación
  cruzada para `K=3`.

`R/02_admixture.R` une `.Q`, `.fam` y el panel de muestras para generar
`ADMIXTURE_K3.png`.

### Qué deberíamos observar

El gráfico presenta una barra por individuo, dividida en tres colores. Se
pueden comparar patrones entre poblaciones y variación entre individuos.

Los colores son componentes estadísticos sin nombre. No deben etiquetarse
automáticamente como continentes o ancestrías específicas. Además, el análisis
usa solamente el cromosoma 22 y un único valor de `K`.

---

## 4. Fenotipo y mini-GWAS

Archivo Slurm: `slurm/04_gwas.sbatch`

### Entradas

- `chr22_qc.bed/.bim/.fam`: genotipos después del QC, sin LD pruning.
- `pca.eigenvec`: componentes principales calculados previamente.
- Panel de muestras: población y superpoblación de cada individuo.
- `R/00_simulate_phenotype.R`: simulación del fenotipo.

### Fenotipo simulado

La simulación usa la semilla fija `20261106`. En `fenotipo.tsv`:

- `FID` e `IID`: identificadores de PLINK.
- `pop`: población.
- `super_pop`: superpoblación.
- `PHENO`: `1` para controles y `2` para casos.

Las probabilidades de ser caso difieren entre superpoblaciones para crear un
ejemplo de confusión por estructura poblacional. Cada población contiene casos
y controles.

### GWAS sin PC

```bash
plink \
  --bfile chr22_qc \
  --pheno fenotipo.tsv \
  --pheno-name PHENO \
  --logistic hide-covar \
  --allow-no-sex \
  --out gwas_sin_pc
```

- `--pheno`: indica la tabla de fenotipos.
- `--pheno-name PHENO`: selecciona la columna que contiene caso/control.
- `--logistic`: ajusta una regresión logística para cada variante.
- `hide-covar`: no imprime pruebas adicionales para las covariables.
- `--allow-no-sex`: permite analizar muestras sin sexo registrado en `.fam`.

### GWAS con PC

El segundo modelo agrega:

```bash
--covar pca.eigenvec --covar-number 1-5
```

- `--covar`: indica el archivo de covariables.
- `--covar-number 1-5`: utiliza PC1, PC2, PC3, PC4 y PC5.

Los dos modelos usan exactamente las mismas muestras y variantes. La única
diferencia es el ajuste por los cinco componentes principales.

### Salidas

- `fenotipo.tsv`: fenotipo reproducible utilizado por ambos modelos.
- `gwas_sin_pc.assoc.logistic`: asociación sin ajustar por PC.
- `gwas_con_pc.assoc.logistic`: asociación ajustada por PC1-PC5.
- `gwas_sin_pc.log` y `gwas_con_pc.log`: resumen de cada modelo.
- `gwas_sin_pc.nosex` y `gwas_con_pc.nosex`: muestras sin sexo registrado.

En los archivos `.assoc.logistic`, las columnas principales son:

- `CHR`: cromosoma.
- `SNP`: identificador de la variante.
- `BP`: posición.
- `A1`: alelo evaluado.
- `TEST`: prueba; el script de R utiliza `ADD`.
- `NMISS`: número de muestras incluidas en esa prueba.
- `OR`: odds ratio estimado.
- `STAT`: estadístico de la regresión.
- `P`: valor p.

`R/03_gwas_qc.R` genera `GWAS_Manhattan_QQ.png` para comparar los modelos sin
y con PC.

### Qué deberíamos observar

El gráfico QQ permite observar si el modelo sin PC presenta mayor desviación
respecto de la diagonal y si el ajuste por PC reduce esa desviación. El cambio
real debe evaluarse en los gráficos generados.

El fenotipo es simulado y el análisis contiene únicamente el cromosoma 22. El
gráfico de asociación no es un Manhattan de todo el genoma y no debe utilizarse
para concluir que existen asociaciones biológicas reales.

---

## Archivos auxiliares comunes

- `.log`: registra el comando, la versión del programa, cantidades procesadas
  y mensajes del análisis. Es el primer archivo que se debe revisar ante un
  problema.
- `.nosex`: lista muestras cuyo sexo no está definido en PLINK. En este
  ejercicio es un archivo esperado y no representa por sí solo un error.
- `logs/*.out`: salida normal de los trabajos Slurm.
- `logs/*.err`: errores y advertencias de los trabajos Slurm.

Los archivos `.bed` son binarios. No deben abrirse ni modificarse con un editor
de texto. Los archivos `.bim`, `.fam`, `.Q`, `.P`, `.eigenvec`, `.eigenval`,
`.frq`, `.imiss`, `.lmiss`, `.tsv` y `.assoc.logistic` sí son tablas de texto.
