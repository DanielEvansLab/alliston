library(data.table)
library(stringr)

#####
#read in dose file, flip dosage data as needed, multiply betas, sum across SNPs
doseHeader<-scan("~/bigdata/U24/GWAScohorts/MrOS/vcf1KG/headerDosage.txt",sep="\t",what="character")
doseHeader <- sapply(str_split(doseHeader,"_" ),function(x) x[1] )
dose<-fread("../data/dose.csv",header=F,col.names=doseHeader)
dose<-as.data.frame(dose)
map<-fread("../data/table5_snp_info.csv")
sum(map$coded_all != dose$A1 | map$noncoded_all != dose$A2)
#alleles in info files and dose files match up
#allele_score collects the allele I want to count.
map[,allele_score := vector("character", length = map[,.N])]
map[expression=="up"   & beta>=0 , allele_score := effect_allele ]
map[expression=="up"   & beta<0  , allele_score := non_effect_allele ]
map[expression=="down" & beta<0  , allele_score := effect_allele ]
map[expression=="down" & beta>=0 , allele_score := non_effect_allele]

#identify which rows contain SNPs that need to be flipped, then loop through those rows and flip dosage
map[, allele_switch := -99]
map[allele_score == coded_all   , allele_switch := 0]
map[allele_score == noncoded_all, allele_switch := 1]
map[,table(allele_switch)]
myrows<-which(map[,allele_switch]==1)
for(i in myrows){
  temp<-unlist(dose[i,4:4618])
  temp<-2-temp
  dose[i,4:4618]<-temp
}
#multiply beta by dosage for each snp
#will do that later

#make new score that counts the number of SNPs that are in the same direction as the profile
#If there are 12 SNPs that represent 2 genes up and 10 genes down, create a count for the SNPs that are concordant for each gene
row_up <- which(map[,expression]=="up")
row_down <- which(map[,expression]=="down")

result <- vector(mode = "list", length = length(4:4618))
for(i in 4:4618){
  dose_id <- dose[,i]
  dose_up <- dose_id[row_up]
  dose_down <- dose_id[row_down]
  score <- sum(dose_up) 

}


#sum allele dosage per person for genes up
row_up <- which(map[,expression]=="up")
IVup <- sapply(dose[row_up, 4:4618], sum ) #4615 length vector
IVup <- unname(IVup)
row_down <- which(map[,expression]=="down")
IVdown <- sapply(dose[row_down, 4:4618], sum ) #4615 length vector
IVdown <- unname(IVdown)
ID <- names(dose)[4:4618]
IV_predictor<-data.table(ID = ID, IVup = IVup, IVdown = IVdown)
summary(IV_predictor$IVup)
cut75 <- unname(quantile(IV_predictor$IVup)["75%"])
cut25 <- unname(quantile(IV_predictor$IVup)["25%"])
#mybreaks <- quantile(IV_predictor$IVup, names=T)
#IV_predictor[,IVup := cut(IVup, breaks = mybreaks)]
IV_predictor[,IVupcat := ifelse(IVup>=cut75,1,0)]
IV_predictor[IVup <= cut25, IVupcat := -1]
cut75 <- unname(quantile(IV_predictor$IVdown)["75%"])
cut25 <- unname(quantile(IV_predictor$IVdown)["25%"])
IV_predictor[,IVdowncat := ifelse(IVdown>=cut75,1,0)]
IV_predictor[IVdown <= cut25, IVdowncat := -1]
IV_predictor[,table(IVupcat, IVdowncat)]
IV_predictor[IVdowncat == -1 & IVupcat == -1, topVbottom := 0 ]
IV_predictor[IVdowncat == 1 & IVupcat == 1  , topVbottom := 1 ]
IV_predictor[IVdowncat == 0 & IVupcat == 0  , topVmid    := 0 ]
IV_predictor[IVdowncat == 1 & IVupcat == 1  , topVmid    := 1 ]
IV_predictor[,table(IVupcat, IVdowncat)]
IV_predictor[,table(topVbottom)]
IV_predictor[,table(topVmid)]

##Merge phenotype data
v1 <- fread("~/bigdata/U24/GWAScohorts/MrOS/pheno/U24/V1FEB14.CSV")
b1 <- fread("~/bigdata/U24/GWAScohorts/MrOS/pheno/U24/B1AUG16.CSV")
fa <- fread("~/bigdata/U24/GWAScohorts/MrOS/pheno/U24/fafeb18.csv")
m1 <- fread("~/bigdata/U24/GWAScohorts/MrOS/pheno/U24/M1AUG16.CSV")
v1_2 <- fread("~/bigdata/U24/GWAScohorts/MrOS/pheno/U24/mros.csv")
setnames(v1_2,"MROSID","ID")

v1 <- v1[,c("ID","GIAGE1"),with=FALSE]
b1 <- b1[,c("ID","B1FND","B1LSD"),with=FALSE]
fa <- fa[,c("ID","FAANYMOS","FAMOSFV1","FAANYHIP","FAHIPFV1"),with=FALSE]
m1 <- m1[,c("ID","M1MEDOST","M1CORTO"),with=FALSE]

pheno <- merge(v1_2,v1,by="ID")
pheno <- merge(pheno,b1,by="ID")
pheno <- merge(pheno,fa,by="ID")
pheno <- merge(pheno,m1,by="ID")

pheno <- pheno[M1MEDOST==0 & M1CORTO==0,]

pheno<-merge(pheno,IV_predictor,by="ID")

#fracture
library(survival)
L1 <- coxph(Surv(FAMOSFV1,FAANYMOS) ~ topVbottom + GIAGE1 + factor(SITE) + EV1 + EV2 + EV3 + EV4, data=pheno)
summary(L1)
L1 <- coxph(Surv(FAMOSFV1,FAANYMOS) ~ topVmid + GIAGE1 + factor(SITE) + EV1 + EV2 + EV3 + EV4, data=pheno)
summary(L1)


###continuous phenotypes
myfit_topVbottom <- function(outcome){
  outcome1 <- pheno[[outcome]]
  lm1<-lm(outcome1 ~ topVbottom + GIAGE1 + factor(SITE) + EV1 + EV2 + EV3 + EV4,data=pheno)
  beta<-summary(lm1)$coefficients[2,1]
  se<-summary(lm1)$coefficients[2,2]
  p<-summary(lm1)$coefficients[2,4]
  return(c(beta,se,p))
}

myfit_topVmid <- function(outcome){
  outcome1 <- pheno[[outcome]]
  lm1<-lm(outcome1 ~ topVmid + GIAGE1 + factor(SITE) + EV1 + EV2 + EV3 + EV4,data=pheno)
  beta<-summary(lm1)$coefficients[2,1]
  se<-summary(lm1)$coefficients[2,2]
  p<-summary(lm1)$coefficients[2,4]
  return(c(beta,se,p))
}


myfit_topVbottom("B1FND")
myfit_topVbottom("B1LSD")
myfit_topVmid("B1FND")
myfit_topVmid("B1LSD")



