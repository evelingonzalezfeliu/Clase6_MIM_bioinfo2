# Leer los resultados del control de calidad
frecuencias <- read.table("results/qc.frq", header = TRUE)
missing_snp <- read.table("results/qc.lmiss", header = TRUE)
missing_muestra <- read.table("results/qc.imiss", header = TRUE)
fam <- read.table("results/chr22_qc.fam")
bim <- read.table("results/chr22_qc.bim")

# Gráficos de QC
par(mfrow = c(2, 2))

hist(frecuencias$MAF, breaks = 30, col = "steelblue",
     main = "Frecuencia del alelo menor", xlab = "MAF")

hist(missing_snp$F_MISS, breaks = 30, col = "orange",
     main = "Missingness por variante", xlab = "F_MISS")

hist(missing_muestra$F_MISS, breaks = 30, col = "seagreen",
     main = "Missingness por muestra", xlab = "F_MISS")

plot.new()
text(0.5, 0.65, paste("Muestras después del QC:", nrow(fam)))
text(0.5, 0.45, paste("Variantes después del QC:", nrow(bim)))
text(0.5, 0.25, "Datos del cromosoma 22")

# Guardar los gráficos de QC
dev.copy(png, "results/QC.png", width = 1400, height = 1000, res = 130)
dev.off()

# Leer los resultados de asociación
sin_pc <- read.table("results/gwas_sin_pc.assoc.logistic", header = TRUE)
con_pc <- read.table("results/gwas_con_pc.assoc.logistic", header = TRUE)

# Conservar el modelo aditivo y los valores P válidos
sin_pc <- sin_pc[sin_pc$TEST == "ADD" & !is.na(sin_pc$P) & sin_pc$P > 0, ]
con_pc <- con_pc[con_pc$TEST == "ADD" & !is.na(con_pc$P) & con_pc$P > 0, ]

# Asociación y QQ sin PC y con PC
par(mfrow = c(2, 2))

plot(sin_pc$BP, -log10(sin_pc$P), pch = 20, col = "steelblue",
     xlab = "Posición en el cromosoma 22", ylab = "-log10(P)",
     main = "Asociación sin PC")
abline(h = -log10(5e-8), lty = 2, col = "red")

esperado <- -log10(ppoints(nrow(sin_pc)))
observado <- -log10(sort(sin_pc$P))
plot(esperado, observado, pch = 20,
     xlab = "Esperado", ylab = "Observado", main = "QQ sin PC")
abline(0, 1, lty = 2, col = "red")

plot(con_pc$BP, -log10(con_pc$P), pch = 20, col = "steelblue",
     xlab = "Posición en el cromosoma 22", ylab = "-log10(P)",
     main = "Asociación con PC1-PC5")
abline(h = -log10(5e-8), lty = 2, col = "red")

esperado <- -log10(ppoints(nrow(con_pc)))
observado <- -log10(sort(con_pc$P))
plot(esperado, observado, pch = 20,
     xlab = "Esperado", ylab = "Observado", main = "QQ con PC")
abline(0, 1, lty = 2, col = "red")

# Guardar los gráficos del mini-GWAS
dev.copy(png, "results/GWAS_Manhattan_QQ.png",
         width = 1500, height = 1100, res = 130)
dev.off()
