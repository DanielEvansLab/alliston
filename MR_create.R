#merge eSNPs with MrOS GWAS data to create score
library(data.table)

vitD <- data.frame(SNP = c("rs3755967","rs12785878","rs10741657","rs17216707","rs10745742","rs8018720"),
		   beta = c(0.089, 0.036, 0.031, 0.026, 0.017, 0.017),
		   se = c(0.0023, 0.0022, 0.0022, 0.0027, 0.0022, 0.0029),
		   effect_allele = c("C","T","A","T","T","G"),
		   other_allele = c("T","G","G","C","C","C"),
		   eaf = c(0.72, 0.75, 0.4, 0.79, 0.4, 0.18),
		   gene = c("GC","NADSYN1/DHCR7","CYP2R1","CYP24A1","AMDHD1","SEC23A"),
		   chr_b37 = c(4L, 11L, 11L, 20L, 12L, 14L),
		   pos_b37 = c(72609398L, 71167449L, 14914878L, 52732362L, 96358529L, 39556185L),
		   samplesize = rep(79366L,6)
		   )

MRinst<-vitD
MRinst$rsID <- paste(MRinst$chr_b37,MRinst$pos_b37,sep=":")
MRinst <- as.data.table(MRinst)

MRinst <- fread("~/bigdata/alliston/HOBeSNP/data/LD_result")
#there are 3 rows with expression as blank
#MRinst[, table(expression)]
#MRinst[expression!="down" & expression!="up"]
#MRinst[symbol_human == "ANKH"]
#MRinst[symbol_human == "NDST4"]
mychr <- base::strsplit(MRinst$snp, ":", fixed = TRUE)
mychr <- as.integer(sapply(mychr, function(x) x[1]))
MRinst[ , chr := mychr ] 
MRinst[ , alleleSame := 0 ]
MRinst[a1 == tempA & a2 == tempB , alleleSame := 1]
MRinst <- MRinst[table_num == 5L ][order(chr)]

#merge to MrOS 1KG map files
mychr <- unique(MRinst$chr)
result <- vector("list", length(mychr))
for (i in seq_along(mychr)){
  map <- fread(paste0("~/U24/GWAScohorts/MrOS/1KGinfo/SNP_info_chr",mychr[i],".txt"))
  MRinst2 <- MRinst[chr == mychr[i]]
  result[[i]] <- merge(MRinst2, map, by.x = "snp", by.y = "SNPID", all.x = TRUE)
}

all_map <- rbindlist(result)
setnames(all_map, c("a1","a2"), c("effect_allele","non_effect_allele"))
fwrite(all_map,file="../data/table5_QC.csv", na=NA, quote = FALSE)

all_map <- fread("../data/table5_QC.csv")
###########################

#check that alleles match
alleleCheck <- function(mydat){
  #MR instrument alleles are effect_allele and non_effect_allele
  #Cohort alleles are coded_all and noncoded_all
  mydat[ ,effect_allele:=toupper(effect_allele ) ]
  mydat[ ,non_effect_allele:=toupper(non_effect_allele) ]
  mydat[ ,coded_all:=toupper(coded_all) ]
  mydat[ ,noncoded_all:=toupper(noncoded_all) ]

  #test for same alleles 
  mydat[,noMatch:=0]
  mydat[(effect_allele!=coded_all | non_effect_allele!=noncoded_all) & !is.na(non_effect_allele), noMatch:=1]
  
  #test if switching coding and non-coding alleles solves non-matching.
  mydat[ , allele_switch:=0]
  mydat[noMatch==1 & effect_allele==noncoded_all & non_effect_allele==coded_all, allele_switch:=1  ]

  #recheck SNPs with NA in non_effect_allele 
  mydat[is.na(non_effect_allele) & effect_allele!=coded_all & effect_allele!=noncoded_all , noMatch:=1 ]
  mydat[is.na(non_effect_allele) & effect_allele!=coded_all & effect_allele==noncoded_all , allele_switch:=1 ]

  return(list(noMatchNoFlip = mydat[noMatch==1 & allele_switch==0, .N ], noMatchYesFlip = mydat[noMatch==1 & allele_switch==1,.N] ))
}


#setnames(MRinst,c("other_allele","codedAll","noncodedAll"),c("non_effect_allele","coded_all","noncoded_all"))
MRinst <- all_map
setnames(MRinst, "chr.x", "chr")
alleleCheck(MRinst)
MRinst <- MRinst[order(chr,SNPorder),]
fwrite(MRinst, file="../data/table5_snp_info.csv",na=NA)

########################
#now that map file is made, extract the SNPs

#extract SNPs from dosage file
my_chr<-unique(MRinst[,chr] )
for(i in 1:length(my_chr)){
  snp_annot<-MRinst[chr==my_chr[i],] 
  for(j in 1:snp_annot[,.N]) {
    mySNPorder<-snp_annot[j,SNPorder ]
    #this works because header is deleted from dose files, so SNPorder
    #same as line number
    SNPline<-fread(paste0("~/bigdata/U24/GWAScohorts/MrOS/plink1KG/chr",my_chr[i],".plink.dose"),skip=mySNPorder-1,nrows=1,header=F,colClasses=c("character"))
    SNPline <- unname(unlist(SNPline))
    if(i==1 & j==1  ){
      dose_out<-SNPline
    } else {
      dose_out<-rbind(dose_out,SNPline)
    }    
  }
}

dose_out <- as.data.table(dose_out)
fwrite(dose_out, file = "../data/dose.csv", col.names = F, quote = F)

#####
