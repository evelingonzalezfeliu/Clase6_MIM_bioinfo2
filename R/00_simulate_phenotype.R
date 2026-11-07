args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("uso: Rscript 00_simulate_phenotype.R metadata fam salida semilla")
metadata_file <- args[1]; fam_file <- args[2]; output_file <- args[3]; seed <- as.integer(args[4])

meta <- read.delim(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
fam <- read.table(fam_file, stringsAsFactors = FALSE)
if (ncol(fam) < 2L) stop("FAM inválido")
names(fam)[1:2] <- c("FID", "IID")
idx <- match(fam$IID, meta$sample)
if (anyNA(idx)) stop("Hay IDs del FAM ausentes en metadatos")
d <- data.frame(FID=fam$FID, IID=fam$IID, pop=meta$pop[idx],
                super_pop=meta$super_pop[idx], stringsAsFactors=FALSE)
rates <- c(AFR=.70, AMR=.60, EAS=.40, EUR=.30, SAS=.50)
set.seed(seed)
d$PHENO <- 1L
for (p in unique(d$pop)) {
  ii <- which(d$pop == p)
  group <- d$super_pop[ii[1]]
  rate <- if (group %in% names(rates)) rates[[group]] else .5
  ncase <- max(1L, min(length(ii)-1L, round(length(ii)*rate)))
  d$PHENO[sample(ii, ncase, replace=FALSE)] <- 2L
}
tab <- table(d$pop, d$PHENO)
if (ncol(tab) != 2L || any(tab == 0L)) stop("Cada población debe contener casos y controles")
write.table(d, output_file, sep="\t", quote=FALSE, row.names=FALSE)
cat("Semilla:", seed, "\n"); print(tab)
