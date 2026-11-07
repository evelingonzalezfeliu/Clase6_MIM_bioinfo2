# Leer los resultados del PCA y la información de las muestras
pca <- read.table("results/pca.eigenvec", header = TRUE,
                  stringsAsFactors = FALSE, check.names = FALSE)
metadatos <- read.delim(
  "data/integrated_call_samples_v3.20130502.ALL.panel",
  stringsAsFactors = FALSE
)

# Agregar la población de cada muestra
names(pca)[1:2] <- c("FID", "IID")
pca$poblacion <- metadatos$pop[match(pca$IID, metadatos$sample)]

# Colores para las cuatro poblaciones
colores <- c(
  CEU = "#377EB8",
  YRI = "#E41A1C",
  CHB = "#4DAF4A",
  PEL = "#984EA3"
)

# Gráfico PC1 versus PC2
plot(
  pca$PC1, pca$PC2,
  col = colores[pca$poblacion], pch = 19,
  xlab = "PC1", ylab = "PC2",
  main = "PCA: cromosoma 22"
)
legend("topright", legend = names(colores), col = colores,
       pch = 19, bty = "n")

# Guardar el gráfico
dev.copy(png, "results/PCA.png", width = 1100, height = 850, res = 130)
dev.off()
