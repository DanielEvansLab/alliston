library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(data.table)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#one way to extract from db is to retrieve all genes into an object, then 
#subset the object by gene id
#ghs <- genes(txdb)
#extract 1 gene
#ghs[ghs$gene_id == 1,]
#extract multiple genes
#ghs[ghs$gene_id %in% c(1,10,100),]

#list available meta columns
columns(txdb)
#the other way is to select columns and row filters from the gene function
#genes(txdb, columns = c("GENEID", "TXSTART", "TXEND", "TXID"),
#      filter = list(gene_id = c(1,10,100)))

#return start/stop for all gene ids, then use that to annotate the GWAS results,
#then do split/apply/combine by gene id with the GWAS results

#this returns for all genes
ghs <- as.data.frame(genes(txdb, columns = c("GENEID")), stringsAsFactors = F)
ghs$GENEID <- as.integer(unlist(ghs$GENEID))

gwas <- fread("gunzip -c ../data/GWAS/Biobank2-British-FracA-As-C-Gwas-SumStats.txt.gz")
gwas[,gene := NA]

#annotate GWAS file with gene names
for ( i in seq_along(ghs$GENEID)){
  mychr <- ghs$chr[i]
  mystart <- ghs$start[i]
  myend <- ghs$end[i]
  mygene <- ghs$GENEID[i]
  if( gwas[CHR == mychr][BP >= mystart & BP <= myend][, .N] > 0 ){
    gwas[CHR == mychr & BP >= mystart & BP <= myend, gene := mygene]
  } 
}

gwas <- gwas[!is.na(gene)]

gwas[, min(P.NI, na.rm = T), by = gene] # gene-based score



tab <- read.csv("~/bigdata/alliston/HOBeSNP/data/serra_tables/12 TGFb_FDR 0_1.csv", header= T, stringsAsFactors=F)
tab <- read.csv("~/bigdata/alliston/HOBeSNP/data/serra_tables/2 Common DE Aging TGFb.csv", header= T, stringsAsFactors=F)
ghs <- genes(txdb, columns = c("GENEID"),
      filter = list(gene_id = tab$Entrez_Id))
ghs <- as.data.frame(ghs, stringsAsFactors=F)
ghs$GENEID <- as.integer(unlist(ghs$GENEID))
str(ghs)
ghs[1:5,]

geneids <- as.data.frame(org.Hs.egSYMBOL)
#geneids <- select(org.Hs.eg.db, keys = as.character(ghs$GENEID), c("ENTREZID", "GENENAME"), keytype)
ghs <- merge(ghs, geneids, by.x = "GENEID", by.y = "gene_id")
str(ghs)

#there are 343 genes in Serra's table, but only 334 in annotation database
#how many Serra table genes are not found in annotation database?
tab <- merge(tab, ghs, by.x = "Entrez_Id", by.y = "GENEID", all.x = T)
#gene id 199 is on chr 6, but it is absent from annotation file. I'll try a few other files later.

#start is always smaller than end
sum(ghs$start > ghs$end)
ghs$start <- ghs$start - 150000L
ghs$end <- ghs$end + 150000L
ghs$chr <- as.character(ghs$seqnames)
ghs$chr <- gsub("chr", "", ghs$chr)
ghs$chr[ghs$chr == "X"] <- "23"
ghs$chr <- as.integer(ghs$chr)

gwas <- fread("gunzip -c ../data/GWAS/Biobank2-British-FracA-As-C-Gwas-SumStats.txt.gz")
gwas <- fread("gunzip -c ../data/GWAS/Biobank2-British-Bmd-As-C-Gwas-SumStats.txt.gz")
result <- vector(mode = "list", length = length(ghs$GENEID))

#select most significant SNP from each gene region
for ( i in seq_along(ghs$GENEID)){
  mychr <- ghs$chr[i]
  mystart <- ghs$start[i]
  myend <- ghs$end[i]
  if( gwas[CHR == mychr][BP >= mystart & BP <= myend][, .N] > 0 ){
    gwassub <- gwas[CHR == mychr][BP >= mystart & BP <= myend]
    minP <- min(gwassub$P.NI, na.rm=T)
    gwassub2 <- gwassub[P.NI == minP][1]
  } else {
    gwassub <- gwas[1]
    gwassub2 <- as.data.table(lapply(gwassub, function(x) x <- NA))
  }
  result[[i]] <- gwassub2

}
result2 <- rbindlist(result)
ghs <- as.data.table(ghs)
result2 <- cbind(ghs, result2)
fwrite(result2, file = "../results/table12_GWAS_BMD.csv", na = "NA")






