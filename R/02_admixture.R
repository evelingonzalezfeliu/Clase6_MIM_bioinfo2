# Leer ADMIXTURE, el orden de las muestras y sus poblaciones
K <- 3
q <- read.table("results/chr22_pruned.3.Q")
fam <- read.table("results/admixture_K3.fam", stringsAsFactors = FALSE)
metadatos <- read.delim(
  "data/integrated_call_samples_v3.20130502.ALL.panel",
  stringsAsFactors = FALSE
)

# Agregar la población y ordenar las muestras
names(fam)[1:2] <- c("FID", "IID")
fam$poblacion <- metadatos$pop[match(fam$IID, metadatos$sample)]
orden <- order(fam$poblacion, fam$IID)
fam <- fam[orden, ]
q <- q[orden, ]

# Gráfico de los tres componentes
colores <- c("#E69F00", "#56B4E9", "#009E73")
posiciones <- barplot(
  t(as.matrix(q)), col = colores, border = NA, space = 0,
  ylab = "Proporción estimada", xlab = "Individuos",
  main = "ADMIXTURE: K=3, cromosoma 22",
  names.arg = rep("", nrow(q))
)

# Escribir los nombres de las poblaciones bajo el gráfico
inicios <- which(!duplicated(fam$poblacion))
finales <- c(inicios[-1] - 1, nrow(fam))
centros <- (posiciones[inicios] + posiciones[finales]) / 2
axis(1, at = centros, labels = fam$poblacion[inicios], tick = FALSE)
legend("topright", legend = paste("Componente", 1:K),
       fill = colores, bty = "n")

# Guardar el gráfico
dev.copy(png, "results/ADMIXTURE_K3.png",
         width = 1500, height = 750, res = 130)
dev.off()
