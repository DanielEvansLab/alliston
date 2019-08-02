
#filename <- "../../data/GWAS/Biobank2-British-FracA-As-C-Gwas-SumStats.txt.gz"
filename <- "../../data/GWAS/Biobank2-British-Bmd-As-C-Gwas-SumStats.txt.gz"

format_magma_snploc <- function(filename){
  #select SNPID, chromosome, and base pair position
  require(data.table)
  if(grepl("Frac", filename)){
    dat <- fread(paste0("gunzip -c ", filename))
    snp_loc <- dat[,.(SNP, CHR, BP)]
    fwrite(snp_loc, file = "../data/snp_loc_frac.txt", sep = "\t", na = "NA", quote = F)
  } else {
    dat <- fread(paste0("gunzip -c ", filename))
    snp_loc <- dat[,.(RSID, CHR, BP)]
    setnames(snp_loc, c("SNP", "CHR", "BP"))
    fwrite(snp_loc, file = "../data/snp_loc_bmd.txt", sep = "\t", na = "NA", quote = F)
  }
}

format_magma_snploc(filename)

format_magma_pvalfile <- function(filename){
  #select SNPID, p values, and N
  require(data.table)
  if(grepl("Frac", filename)){
    dat <- fread(paste0("gunzip -c ", filename))
    pval_file <- dat[,.(SNP, P.I, N)]
    setnames(pval_file, c("SNP", "P", "N"))
    fwrite(pval_file, file = "../data/pvalfile_frac.txt", sep = "\t", na = "NA", quote = F)
  } else {
    dat <- fread(paste0("gunzip -c ", filename))
    pval_file <- dat[,.(RSID, P.NI, N)]
    setnames(pval_file, c("SNP", "P", "N") )
    fwrite(pval_file, file = "../data/pvalfile_bmd.txt", sep = "\t", na = "NA", quote = F)
  }
}

format_magma_pvalfile(filename)

#creating geneset file
format_magma_geneset <- function(file1){
  require(data.table)
  dat <- fread(paste0("../../data/serra_tables/", file1))
  paste(c(sub(".csv", "", file1), dat$Entrez_Id), collapse = " ")
}

all.files <- list.files("../../data/serra_tables/")
result <- lapply(all.files, format_magma_geneset)
result <- do.call(rbind, result)
write.table(result, file = "~/bigdata/alliston/HOBeSNP/magma/data/serra_sets.txt", col.names = F, row.names = F, quote = F)

#format gsa results and plot
library(data.table)
library(ggplot2)

bmd <- fread("../results/bmd.gsa.out", skip = 3)
bmd[, outcome := "bmd"]
str(bmd)

frac <- fread("../results/frac.gsa.out", skip = 3)
frac[, outcome := "frac"]
str(frac)

dat <- rbind(bmd, frac)
dat[, log10p := -log10(P)]
ggplot(data = dat) +
	geom_col(mapping = aes(x = VARIABLE, y = log10p, fill = outcome), position = "dodge") +
	geom_hline(yintercept = -log10(0.05), linetype = 2) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
	labs(x = "gene sets", y = "-log10 P-value")
ggsave("pvals.png")


