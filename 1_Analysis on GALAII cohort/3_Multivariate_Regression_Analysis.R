library(ggplot2)
library(cowplot)
library(ggpubr)
library(rms)
library(openxlsx)

#================#
# load phenotype #
#================#

phen <- read.table("695.GALA.final.phen.with.add.txt", header=T, strings=F, row.names=1, sep='\t')
rownames(phen) <- phen$SubjectID
phen$asthma_status <- factor(plyr::mapvalues(phen$control, from=c(0,1), to=c("asthmatic","HC")), levels=c("HC", "asthmatic"))

# load type2 and interferon status
type2_interferon.df <- read.table("interferon_type2.status.txt", header=T, strings=F, row.names=1, sep='\t')

phen$type2_status <- factor(type2_interferon.df[rownames(phen), "type2_status"], levels=c("type2_low", "type2_high"))
phen$interferon_status <- factor(type2_interferon.df[rownames(phen), "interferon_status"], levels=c("interferon_low", "interferon_high"))
phen$Male <- factor(phen$Male, levels=c("Female", "Male"))

#=================#
# load expression #
#=================#

vcf <- "/Seibold/proj/GALA_720donors/WGS_DATA_EQTL/DATA/EQTL_PREP/GENOTYPE_DATA/GALA_WGS_GENO.vcf.gz"
vcf.samples <- system(paste0("bcftools query -l ", vcf), intern=T)
expr <- read.table("GALA_681_FASTQTL_EXPRM.bed.gz", header=F)
colnames(expr) <- c("chr", "chr_start", "chr_end", "gene", vcf.samples)
rownames(expr) <- expr[,4]
expr <- expr[,-c(1:4)]

X_vcf <- "/Seibold/proj/GALA_720donors/WGS_DATA_EQTL/DATA/XCHR/681_X.vcf.gz"
vcf.samples <- system(paste0("bcftools query -l ", X_vcf), intern=T)
expr_X <- read.table("GALA_681_CHRX_FASTQTL_EXPRM.bed.gz", header=F)
colnames(expr_X) <- c("chr", "chr_start", "chr_end", "gene", vcf.samples)
rownames(expr_X) <- expr_X[,4]
expr_X <- expr_X[,-c(1:4)]

phen$ACE2_eQTL_norm <- NA
phen[colnames(expr_X), "ACE2_eQTL_norm"] <-  t(expr_X["ACE2", ])
phen$TMPRSS2_eQTL_norm <- NA
phen[colnames(expr), "TMPRSS2_eQTL_norm"] <-  t(expr["TMPRSS2", ])

#====================#
# Load genotype data #
#====================#

library(biomaRt)

snp <- "rs181603331"
vcf.file <- "GALA_WGS_GENO.vcf.gz"
vcf <- X_vcf.file <- "/Seibold/proj/GALA_720donors/WGS_DATA_EQTL/DATA/XCHR/681_X.vcf.gz"

snpmart = useMart(biomart = "ENSEMBL_MART_SNP", dataset="hsapiens_snp")

pos <- getBM(attributes = c('refsnp_id','chr_name','chrom_start'),
  filters = c('snp_filter'),
		values = snp,
		mart = snpmart)
if(dim(pos)[1]==0) stop("Could not find SNP ID!")
snp.region <- paste0("chr", pos$chr_name,":", pos$chrom_start)[1]

vcf.snp.df <- system(paste0("bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' -r ", snp.region, " ", vcf), intern=T)
vcf.snp.df <- read.table(text=vcf.snp.df, header=F, colClasses="character", stringsAsFactors=F, sep='\t')

vcf.samples <- system(paste0("bcftools query -l ", vcf), intern=T)

colnames(vcf.snp.df) <- c("CHROM", "POS", "ID", "REF", "ALT", vcf.samples)

A1 <- as.character(vcf.snp.df[1, "ALT"])
A2 <- as.character(vcf.snp.df[1, "REF"])

snp.df <- data.frame(t(vcf.snp.df[,vcf.samples]), stringsAsFactors=F)
colnames(snp.df) <- snp

phen_geno <- cbind(phen[vcf.samples,], snp.df)
phen_geno <- phen_geno[phen_geno[,snp]!="./.",]
phen_geno$rs181603331 <- as.numeric(plyr::mapvalues(phen_geno$rs181603331, from=c("0/0", "0/1", "1/1"), to=c(0,1,2)))

#TMPRSS2
snp <- "rs1475908"
snp2 <- "rs74659079"
snp3 <- "rs2838057"
vcf <- vcf.file <- "/Seibold/proj/GALA_720donors/WGS_DATA_EQTL/DATA/EQTL_PREP/GENOTYPE_DATA/GALA_WGS_GENO.vcf.gz"

snpmart = useMart(biomart = "ENSEMBL_MART_SNP", dataset="hsapiens_snp")

pos <- getBM(attributes = c('refsnp_id','chr_name','chrom_start'),
  filters = c('snp_filter'),
        values = c(snp,snp2,snp3),
        mart = snpmart)
if(dim(pos)[1]==0) stop("Could not find SNP ID!")
snp.region <- paste0(paste0(pos$chr_name,":", pos$chrom_start), collapse=",")

vcf.snp.df <- system(paste0("bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' -r ", snp.region, " ", vcf), intern=T)
vcf.snp.df <- read.table(text=vcf.snp.df, header=F, colClasses="character", stringsAsFactors=F, sep='\t')

vcf.samples <- system(paste0("bcftools query -l ", vcf), intern=T)

colnames(vcf.snp.df) <- c("CHROM", "POS", "ID", "REF", "ALT", vcf.samples)

A1 <- as.character(vcf.snp.df[1, "ALT"])
A2 <- as.character(vcf.snp.df[1, "REF"])

snp.df <- data.frame(t(vcf.snp.df[,vcf.samples]), stringsAsFactors=F)
colnames(snp.df) <- c(snp3, snp, snp2)

phen_geno2 <- cbind(phen[vcf.samples,], snp.df)
#phen_geno2 <- phen_geno2[phen_geno2[,snp]!="./.",]
phen_geno2$rs1475908 <- as.numeric(plyr::mapvalues(phen_geno2$rs1475908, from=c("0/0", "0/1", "1/1"), to=c(0,1,2)))
phen_geno2$rs74659079 <- as.numeric(plyr::mapvalues(phen_geno2$rs74659079, from=c("0/0", "0/1", "1/1"), to=c(0,1,2)))
phen_geno2$rs2838057 <- as.numeric(plyr::mapvalues(phen_geno2$rs2838057, from=c("0/0", "0/1", "1/1"), to=c(2,1,0)))

#=================#
# load covariates #
#=================#

cov <- read.table("COVARIATES_FILE_GALA_681.txt", sep='\t', header=T, row.names=1)
cov <- data.frame(t(cov))

cov$ACE2 <- t(expr_X["ACE2",rownames(cov)])
cov$TMPRSS2 <- t(expr["TMPRSS2",rownames(cov)])
cov$rs181603331 <- NA
cov[rownames(phen_geno), "rs181603331"] <- phen_geno$rs181603331
cov$rs1475908 <- NA
cov[rownames(phen_geno2), "rs1475908"] <- phen_geno2$rs1475908
cov[rownames(phen_geno2), "rs74659079"] <- phen_geno2$rs74659079
cov[rownames(phen_geno2), "rs2838057"] <- phen_geno2$rs2838057
cov$Male <- phen[rownames(cov), "Male"]

#interferon status and type2 status
cov$interferon_status <- phen[rownames(cov), "interferon_status"]
cov$type2_status <- phen[rownames(cov), "type2_status"]

#==================================================#
# Model expression of ACE2 with all the covariates #
#==================================================#


### univariate model
traits_to_test <- c("age", "interferon_status", "type2_status", "sex", "Asthma_Status", "rs181603331")
ACE2_uni_cov.df <- lapply(traits_to_test, function(x) {
    fit_lm <- lm(as.formula(paste0("ACE2 ~ ", x )), data=cov)
    coef <- summary(fit_lm)$coefficients[2,]
    R2 <- summary(fit_lm)$r.squared
    names(coef) <- c("Estimate", "SE", "t", "p-value")
    data.frame(predictor=x, R2=R2*100, t(coef))
})

ACE2_uni_cov.df <- do.call(rbind, ACE2_uni_cov.df)

### multivariate model without genotype
ACE2_form <- as.formula(paste0("ACE2 ~ age + interferon_status + type2_status + sex + Asthma_Status "))
ACE2_multi_cov_fit <- ols(ACE2_form, data=cov)
ACE2_multi_cov_fit2 <- lm(ACE2_form, data=cov)

var_of_interest <- rownames(summary(ACE2_multi_cov_fit2)$coefficients)[-1]
var_of_interest[2:3] <- c("interferon_status", "type2_status")

plt <- plot(anova(ACE2_multi_cov_fit), what='partial R2')
ACE2_multi_cov.df <- data.frame(predictor=var_of_interest, partial_R2=plt[var_of_interest]*100, summary(ACE2_multi_cov_fit2)$coef[-1,])
colnames(ACE2_multi_cov.df) <- c("predictor", "partial_R2", "Estimate", "SE", "t", "p-value")

### multivariate model with genotype
ACE2_geno_form <- as.formula(paste0("ACE2 ~ age + interferon_status + type2_status + sex + Asthma_Status + rs181603331 "))
ACE2_multi_cov_geno_fit <- ols(ACE2_geno_form, data=cov)
ACE2_multi_cov_geno_fit2 <- lm(ACE2_geno_form, data=cov)

var_of_interest <- rownames(summary(ACE2_multi_cov_geno_fit2)$coefficients)[-1]
var_of_interest[2:3] <- c("interferon_status", "type2_status")

plt <- plot(anova(ACE2_multi_cov_geno_fit), what='partial R2')
ACE2_multi_cov_geno.df <- data.frame(predictor=var_of_interest, partial_R2=plt[var_of_interest]*100, summary(ACE2_multi_cov_geno_fit2)$coef[-1,])
colnames(ACE2_multi_cov_geno.df) <- c("predictor", "partial_R2", "Estimate", "SE", "t", "p-value")

wb <- createWorkbook()
addWorksheet(wb, sheetName="Uni. ACE2")
addWorksheet(wb, sheetName="Multi. ACE2")
addWorksheet(wb, sheetName="Multi. ACE2 with geno")
writeData(wb, sheet="Uni. ACE2", x=ACE2_uni_cov.df)
writeData(wb, sheet="Multi. ACE2", x=ACE2_multi_cov.df)
writeData(wb, sheet="Multi. ACE2 with geno", x=ACE2_multi_cov_geno.df)
saveWorkbook(wb, file="ACE2_noPEER_new_prediction.xlsx", overwrite=T)

#==========================================================================#
# Model expression of TMPRSS2 with all the covariates (with all ind. hits) #
#==========================================================================#

### univariate model
traits_to_test <- c("age", "interferon_status", "type2_status", "sex", "Asthma_Status", "rs1475908", "rs74659079", "rs2838057")
TMPRSS2_uni_cov.df <- lapply(traits_to_test, function(x) {
    fit_lm <- lm(as.formula(paste0("TMPRSS2 ~ ", x )), data=cov)
    coef <- summary(fit_lm)$coefficients[2,]
    R2 <- summary(fit_lm)$r.squared
    names(coef) <- c("Estimate", "SE", "t", "p-value")
    data.frame(predictor=x, R2=R2*100, t(coef))
})

TMPRSS2_uni_cov.df <- do.call(rbind, TMPRSS2_uni_cov.df)

### multivariate model without genotype
TMPRSS2_form <- as.formula(paste0("TMPRSS2 ~ age + interferon_status + type2_status + sex + Asthma_Status "))
TMPRSS2_multi_cov_fit <- ols(TMPRSS2_form, data=cov)
TMPRSS2_multi_cov_fit2 <- lm(TMPRSS2_form, data=cov)

var_of_interest <- rownames(summary(TMPRSS2_multi_cov_fit2)$coefficients)[-1]
var_of_interest[2:3] <- c("interferon_status", "type2_status")

plt <- plot(anova(TMPRSS2_multi_cov_fit), what='partial R2')
TMPRSS2_multi_cov.df <- data.frame(predictor=var_of_interest, partial_R2=plt[var_of_interest]*100, summary(TMPRSS2_multi_cov_fit2)$coef[-1,])
colnames(TMPRSS2_multi_cov.df) <- c("predictor", "partial_R2", "Estimate", "SE", "t", "p-value")

### multivariate model with genotype
TMPRSS2_geno_form <- as.formula(paste0("TMPRSS2 ~ age + interferon_status + type2_status + sex + Asthma_Status + rs1475908 + rs74659079 + rs2838057"))
TMPRSS2_multi_cov_geno_fit <- ols(TMPRSS2_geno_form, data=cov)
TMPRSS2_multi_cov_geno_fit2 <- lm(TMPRSS2_geno_form, data=cov)

var_of_interest <- rownames(summary(TMPRSS2_multi_cov_geno_fit2)$coefficients)[-1]
var_of_interest[2:3] <- c("interferon_status", "type2_status")

plt <- plot(anova(TMPRSS2_multi_cov_geno_fit), what='partial R2')
TMPRSS2_multi_cov_geno.df <- data.frame(predictor=var_of_interest, partial_R2=plt[var_of_interest]*100, summary(TMPRSS2_multi_cov_geno_fit2)$coef[-1,])
colnames(TMPRSS2_multi_cov_geno.df) <- c("predictor", "partial_R2", "Estimate", "SE", "t", "p-value")

wb <- createWorkbook()
addWorksheet(wb, sheetName="Uni. TMPRSS2")
addWorksheet(wb, sheetName="Multi. TMPRSS2")
addWorksheet(wb, sheetName="Multi. TMPRSS2 with geno")
writeData(wb, sheet="Uni. TMPRSS2", x=TMPRSS2_uni_cov.df)
writeData(wb, sheet="Multi. TMPRSS2", x=TMPRSS2_multi_cov.df)
writeData(wb, sheet="Multi. TMPRSS2 with geno", x=TMPRSS2_multi_cov_geno.df)
saveWorkbook(wb, file="TMPRSS2_noPEER_new_preds._all_ind._hits.xlsx", overwrite=T)

