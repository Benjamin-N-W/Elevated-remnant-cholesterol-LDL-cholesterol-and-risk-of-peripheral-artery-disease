
rm(list=ls())


# #Setting working directory and packages ----------------------------------------------


setwd("H:/Projekt remnant cholesterol/Artikel 6 - Mendelian randomization PAD in CGPS/Working folder")

library(data.table)
library(TwoSampleMR)
library(ieugwasr)
library(mygene)
library(ggplot2)
library(eoffice)
library(cowplot)
library(forestploter)
library(grid)
library(stringr)
library(officer)
library(data.table)


#Formatting P-values function for forestploter
format_p<-function(plot,pvals,plot_cols){
  
  for (q in 1:length(pvals))
    for(i in 1:length(pvals[[q]])){
      
      p_i<-if (grepl( "e",pvals[[q]][i])){
        bquote(.(paste0(word(pvals[[q]][i], 1, sep="e")," x 10"))^.(paste0("-",sub("^0+", "", word(pvals[[q]][i], 2, sep="-")))))} else {
          str_pad(pvals[[q]][i], width=4, side="right", pad="0")}
      
      plot <- add_text(plot, p_i, part="body",just ="left", row=i,col=plot_cols[q],gp=gpar(fontsize=9))
    }
  return(plot)
}

# Filtering GWAS by P-values to get smaller work files --------------------


# Extracting exposures from local GWAS summary files
ukbb_remnant_c<-fread("GWAS_remnant_cholesterol_MAF0.01")
ukbb_ldl_c<-fread("GWAS_ldl_direct_MAF0.01")
ukbb_hdl_c<-fread("GWAS_hdl_cholesterol_MAF0.01")

ukbb_remnant_c$TRAIT<-"Remnant cholesterol"
ukbb_ldl_c$TRAIT<-"LDL cholesterol"
ukbb_hdl_c$TRAIT<-"HDL cholesterol"

# Getting SNPs with low P-values for smaller files (computational efficiency)
snps_remnant_C<-subset(ukbb_remnant_c,pval.exposure<5e-6)$ID
snps_ldl_C<-subset(ukbb_ldl_c,P<5e-6)$ID
snps_ukbb_hdl_c<-subset(ukbb_hdl_c,P<5e-6)$ID

# Filtering data sets to only retain SNPS with P values <1e-8 for either remnant or LDL cholesterol
ukbb_remnant_c<-ukbb_remnant_c[match(unique(c(snps_remnant_C,snps_ldl_C,snps_ukbb_hdl_c)), ukbb_remnant_c$ID, nomatch=0),]
ukbb_ldl_c<-ukbb_ldl_c[match(unique(c(snps_remnant_C,snps_ldl_C,snps_ukbb_hdl_c)), ukbb_ldl_c$ID, nomatch=0),]
ukbb_hdl_c<-ukbb_hdl_c[match(unique(c(snps_remnant_C,snps_ldl_C,snps_ukbb_hdl_c)), ukbb_hdl_c$ID, nomatch=0),]

fwrite(ukbb_remnant_c,"GWAS_remnant_cholesterol_lowpval",sep = "\t")
fwrite(ukbb_ldl_c,"GWAS_ldl_direct_lowpval",sep = "\t")
fwrite(ukbb_hdl_c,"GWAS_hdl_cholesterol_lowpval",sep = "\t")


# And for no statins GWAS

# Extracting exposures from local GWAS summary files
ukbb_nostatins_remnant_c<-fread("GWAS_nostatins_remnant_cholesterol_MAF0.01")
ukbb_nostatins_ldl_c<-fread("GWAS_nostatins_ldl_direct_MAF0.01")

ukbb_nostatins_remnant_c$TRAIT<-"Remnant cholesterol"
ukbb_nostatins_ldl_c$TRAIT<-"LDL cholesterol"

# Getting SNPs with low P-values for smaller files (computational efficiency)
snps_remnant_C<-subset(ukbb_nostatins_remnant_c,P<5e-6)$ID
snps_ldl_C<-subset(ukbb_nostatins_ldl_c,P<5e-6)$ID

# Filtering data sets to only retain SNPS with P values <1e-8 for either remnant or LDL cholesterol
ukbb_nostatins_remnant_c<-ukbb_nostatins_remnant_c[match(unique(c(snps_remnant_C,snps_ldl_C)), ukbb_nostatins_remnant_c$ID, nomatch=0),]
ukbb_nostatins_ldl_c<-ukbb_nostatins_ldl_c[match(unique(c(snps_remnant_C,snps_ldl_C)), ukbb_nostatins_ldl_c$ID, nomatch=0),]


fwrite(ukbb_nostatins_remnant_c,"GWAS_nostatins_remnant_cholesterol_lowpval",sep = "\t")
fwrite(ukbb_nostatins_ldl_c,"GWAS_nostatins_ldl_direct_lowpval",sep = "\t")




# Remnant cholesterol gene-based score ---------------------------


# Reading remnant cholesterol data
ukbb_remnant_c<-read_exposure_data("GWAS_remnant_cholesterol_lowpval",sep="\t",
                                   snp_col = "ID",
                                   beta_col = "BETA",
                                   se_col = "SE",
                                   effect_allele_col = "ALLELE1",
                                   other_allele_col = "ALLELE0",
                                   eaf_col = "A1FREQ",
                                   pval_col = "P",
                                   units_col = "",
                                   gene_col = "",
                                   samplesize_col = "",
                                   phenotype_col="TRAIT",
                                   chr_col="CHROM",
                                   pos_col="GENPOS")

# Getting positions of genes
a<-mygene::queryMany(c('LPL','LIPC','GPIHBP1','APOC3','APOA5','ANGPTL3','ANGPTL4','ANGPTL8' # Triglyceride metabolism
               #'CETP',"APOA1"				# HDL metabolism
               ), scopes="symbol", fields="genomic_pos_hg19", species="human")[,c("genomic_pos_hg19.chr","genomic_pos_hg19.start","genomic_pos_hg19.end")]

# Adding 100,000 base pairs of each side of gene 
a<-data.frame(chr=a$genomic_pos_hg19.chr,
              start=a$genomic_pos_hg19.start-1e5,
              end=a$genomic_pos_hg19.end+1e5)


# Selecting variants within gene regions by chromosome
ukbb_remnant_c_gene<-data.frame()
for (i in unique(a$chr)){
chr_snps<-setDT(subset(ukbb_remnant_c,chr.exposure==i))[pos.exposure %inrange% a[a$chr==i, c("start","end")]  ]
ukbb_remnant_c_gene<-rbind(ukbb_remnant_c_gene,chr_snps)
}




# Clumping data
ukbb_remnant_c_gene<-clump_data(ukbb_remnant_c_gene,clump_p1=5e-6,clump_r2=0.1,clump_kb = 10000)

# Reading in only the SNPs which are selected for remnant cholesterol
ukbb_remnant_c_ldl<-data.frame(setDT(read_exposure_data("GWAS_ldl_direct_lowpval",sep="\t",
                                   snp_col = "ID",
                                   beta_col = "BETA",
                                   se_col = "SE",
                                   effect_allele_col = "ALLELE1",
                                   other_allele_col = "ALLELE0",
                                   eaf_col = "A1FREQ",
                                   pval_col = "P",
                                   units_col = "",
                                   gene_col = "",
                                   samplesize_col = "",
                                   phenotype_col="TRAIT",
                                   chr_col="CHROM",
                                   pos_col="GENPOS"))[SNP %in% ukbb_remnant_c_gene$SNP  ])

# Combining remnant SNPs effects on remnant and ldl
ukbb_remnant_c_gene<-rbind(ukbb_remnant_c_gene,ukbb_remnant_c_ldl)


ukbb_remnant_c_gene<-ukbb_remnant_c_gene[,c("SNP", "exposure", "id.exposure",  "effect_allele.exposure", "other_allele.exposure",
                  "eaf.exposure", "beta.exposure","se.exposure","pval.exposure"  )]

# Removing redunant objects
rm(ukbb_remnant_c_ldl,chr_snps,a)





# LDL cholesterol gene-based score ----------------------------------------

# Reading remnant cholesterol data
ukbb_ldl_c<-read_exposure_data("GWAS_ldl_direct_lowpval",sep="\t",
                                   snp_col = "ID",
                                   beta_col = "BETA",
                                   se_col = "SE",
                                   effect_allele_col = "ALLELE1",
                                   other_allele_col = "ALLELE0",
                                   eaf_col = "A1FREQ",
                                   pval_col = "P",
                                   units_col = "",
                                   gene_col = "",
                                   samplesize_col = "",
                                   phenotype_col="TRAIT",
                                   chr_col="CHROM",
                                   pos_col="GENPOS")

# Getting positions of genes
b<-mygene::queryMany(c('PCSK9','HMGCR','LDLR','NPC1L1','APOB' 	# LDL metabolism 
                       #'CETP',"APOA1"				# HDL metabolism
), scopes="symbol", fields="genomic_pos_hg19", species="human")[,c("genomic_pos_hg19.chr","genomic_pos_hg19.start","genomic_pos_hg19.end")]


# Adding 100,000 base pairs of each side of gene 
a<-data.frame(chr=b$genomic_pos_hg19.chr,
              start=b$genomic_pos_hg19.start-1e5,
              end=b$genomic_pos_hg19.end+1e5)

# ACLY gne had patch and  had to be formatted differently
a<-rbind(a,mygene::queryMany(c('ACLY' 
), scopes="symbol", fields="genomic_pos_hg19", species="human")[,"genomic_pos_hg19"][[1]][2,c("chr","start","end")])


# Selecting variants within gene regions by chromosome
ukbb_ldl_c_gene<-data.frame()
for (i in unique(a$chr)){
  chr_snps<-setDT(subset(ukbb_ldl_c,chr.exposure==i))[pos.exposure %inrange% a[a$chr==i, c("start","end")]  ]
  ukbb_ldl_c_gene<-rbind(ukbb_ldl_c_gene,chr_snps)
}

# Clumping data
ukbb_ldl_c_gene<-clump_data(ukbb_ldl_c_gene,clump_p1=5e-6,clump_r2=0.1,clump_kb = 10000)

# Reading in only the SNPs which are selected for remnant cholesterol
ukbb_ldl_c_remnant<-data.frame(setDT(read_exposure_data("GWAS_remnant_cholesterol_lowpval",sep="\t",
                                                        snp_col = "ID",
                                                        beta_col = "BETA",
                                                        se_col = "SE",
                                                        effect_allele_col = "ALLELE1",
                                                        other_allele_col = "ALLELE0",
                                                        eaf_col = "A1FREQ",
                                                        pval_col = "P",
                                                        units_col = "",
                                                        gene_col = "",
                                                        samplesize_col = "",
                                                        phenotype_col="TRAIT",
                                                        chr_col="CHROM",
                                                        pos_col="GENPOS"))[SNP %in% ukbb_ldl_c_gene$SNP  ])

# Combining ldl SNPs effects on remnant and ldl
ukbb_ldl_c_gene<-rbind(ukbb_ldl_c_gene,ukbb_ldl_c_remnant)

#removing redundant objects
rm(ukbb_ldl_c_remnant,a,chr_snps)

# Correcting header
ukbb_ldl_c_gene<-ukbb_ldl_c_gene[,c("SNP", "exposure", "id.exposure",  "effect_allele.exposure", "other_allele.exposure",
                                            "eaf.exposure", "beta.exposure","se.exposure","pval.exposure"  )]


# Remnant cholesterol genome-wide score -----------------------------------



ukbb_remnant_c$samplesize.exposure<-77231 


# Clumping data
ukbb_remnant_c<-clump_data(ukbb_remnant_c,clump_p1=5e-8,clump_r2=0.001,clump_kb = 10000)

# Removing SNPs associated with BMI
ukbb_remnant_c<-steiger_filtering(harmonise_data(ukbb_remnant_c,add_metadata(extract_outcome_data(ukbb_remnant_c$SNP, outcomes = 'ukb-b-2303'),cols="sample_size"))) # 'ukb-b-2303 is BMI GWAS


# Reading in only the SNPs which are selected for remnant cholesterol
ukbb_remnant_c_ldl<-data.frame(setDT(read_exposure_data("GWAS_ldl_direct_lowpval",sep="\t",
                                                        snp_col = "ID",
                                                        beta_col = "BETA",
                                                        se_col = "SE",
                                                        effect_allele_col = "ALLELE1",
                                                        other_allele_col = "ALLELE0",
                                                        eaf_col = "A1FREQ",
                                                        pval_col = "P",
                                                        units_col = "",
                                                        gene_col = "",
                                                        samplesize_col = "",
                                                        phenotype_col="TRAIT",
                                                        chr_col="CHROM",
                                                        pos_col="GENPOS"))[SNP %in% ukbb_remnant_c$SNP  ])[,
                              c("SNP", "exposure", "id.exposure",  "effect_allele.exposure", "other_allele.exposure",
                              "eaf.exposure", "beta.exposure","se.exposure","pval.exposure"  )]

ukbb_remnant_c<-ukbb_remnant_c[,c("SNP", "exposure", "id.exposure",  "effect_allele.exposure", "other_allele.exposure",
                                  "eaf.exposure", "beta.exposure","se.exposure","pval.exposure"  )]


# Formatting data frame for multivariable MR
ukbb_remnant_c<-rbind(ukbb_remnant_c,ukbb_remnant_c_ldl)
rm(ukbb_remnant_c_ldl)


# LDL cholesterol genome-wide score ---------------------------------------


ukbb_ldl_c$samplesize.exposure<-95206  

# Clumping data
ukbb_ldl_c<-clump_data(ukbb_ldl_c,clump_p1=5e-8,clump_r2=0.001,clump_kb = 10000)

# Removing SNPs associated with BMI
ukbb_ldl_c<-steiger_filtering(harmonise_data(ukbb_ldl_c,add_metadata(extract_outcome_data(ukbb_ldl_c$SNP, outcomes = 'ukb-b-2303'),cols="sample_size"))) # 'ukb-b-2303 is BMI GWAS


# Reading in only the SNPs which are selected for remnant cholesterol
ukbb_ldl_c_remnant<-data.frame(setDT(read_exposure_data("GWAS_remnant_cholesterol_lowpval",sep="\t",
                                                        snp_col = "ID",
                                                        beta_col = "BETA",
                                                        se_col = "SE",
                                                        effect_allele_col = "ALLELE1",
                                                        other_allele_col = "ALLELE0",
                                                        eaf_col = "A1FREQ",
                                                        pval_col = "P",
                                                        units_col = "",
                                                        gene_col = "",
                                                        samplesize_col = "",
                                                        phenotype_col="TRAIT",
                                                        chr_col="CHROM",
                                                        pos_col="GENPOS"))[SNP %in% ukbb_ldl_c$SNP  ])[,
                                                                                                           c("SNP", "exposure", "id.exposure",  "effect_allele.exposure", "other_allele.exposure",
                                                                                                             "eaf.exposure", "beta.exposure","se.exposure","pval.exposure"  )]

ukbb_ldl_c<-ukbb_ldl_c[,c("SNP", "exposure", "id.exposure",  "effect_allele.exposure", "other_allele.exposure",
                                  "eaf.exposure", "beta.exposure","se.exposure","pval.exposure"  )]


# Formatting data frame for multivariable MR
ukbb_ldl_c<-rbind(ukbb_ldl_c,ukbb_ldl_c_remnant)
rm(ukbb_ldl_c_remnant)


# Reading outcome data ----------------------------------------



snps<-unique(unlist(lapply(exposures, '[[', 'SNP')))

# All SNPs to be extracted
snps<-unique(c(ukbb_remnant_c_gene$SNP,ukbb_ldl_c_gene$SNP,
               ukbb_remnant_c$SNP,ukbb_remnant_c$SNP))


#Reading Million Veteran Program summary data

#PAD in MVP, european american
MVP_PAD_EUR <- read_outcome_data(
  snps = snps,
  filename = "H:/Projekt remnant cholesterol/GWAS summary/dbGaP/GWAS results/CLEANED.MVP.EUR.PAD.results.anno.nodup.txt",
  sep = "\t",
  snp_col = "SNPID",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "EFFECT_ALLELE",
  other_allele_col = "OTHER_ALLELE",
  eaf_col = "EAF",
  pval_col = "PVAL",
  units_col = "",
  samplesize_col = "N_TOTAL",
  ncase_col="N_CASES",
  ncontrol_col="N_CONTROLS",
  chr_col="CHROM",
  pos_col="POS")

#PAD in MVP, african american
MVP_PAD_AFR <- read_outcome_data(
  snps = snps,
  filename = "H:/Projekt remnant cholesterol/GWAS summary/dbGaP/GWAS results/CLEANED.MVP.AFR.PAD.results.anno.nodup.txt",
  sep = "\t",
  snp_col = "SNPID",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "EFFECT_ALLELE",
  other_allele_col = "OTHER_ALLELE",
  eaf_col = "EAF",
  pval_col = "PVAL",
  units_col = "",
  samplesize_col = "N_TOTAL",
  ncase_col="N_CASES",
  ncontrol_col="N_CONTROLS",
  chr_col="CHROM",
  pos_col="POS")

#PAD in MVP, hispanic
MVP_PAD_HIS <- read_outcome_data(
  snps = snps,
  filename = "H:/Projekt remnant cholesterol/GWAS summary/dbGaP/GWAS results/CLEANED.MVP.HIS.PAD.results.anno.nodup.txt",
  sep = "\t",
  snp_col = "SNPID",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "EFFECT_ALLELE",
  other_allele_col = "OTHER_ALLELE",
  eaf_col = "EAF",
  pval_col = "PVAL",
  units_col = "",
  samplesize_col = "N_TOTAL",
  ncase_col="N_CASES",
  ncontrol_col="N_CONTROLS",
  chr_col="CHROM",
  pos_col="POS")

#Coronary artery disease in MVP
MVP_CAD_EUR <- read_outcome_data(
  snps = snps,
  filename = "H:/Projekt remnant cholesterol/GWAS summary/dbGaP/GWAS results/MVP_CAD_Meta_analysis_EUR.withrsID.txt",
  sep = "\t",
  snp_col = "rsmid",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "EA",
  other_allele_col = "NEA",
  eaf_col = "EAF",
  pval_col = "pvalue",
  units_col = "",
  samplesize_col = "TotalSampleSize",
  chr_col="CHR",
  pos_col="POS_b37")


#coronary artery disease data
IHD_C4D <- extract_outcome_data(snps = snps, outcomes = 'ieu-a-7')
IHD_ukbb <- extract_outcome_data(snps = snps, outcomes = 'ebi-a-GCST005194')
IHD_finngen <- extract_outcome_data(snps = snps, outcomes = 'finn-b-I9_ISCHHEART')


#Peripheral artery disease data
PAD_finngen <- extract_outcome_data(snps = snps, outcomes = 'finn-b-I9_PAD')
PAD_bbj <- extract_outcome_data(snps = snps, outcomes = 'bbj-a-144')
PAD_ukbb <- read_outcome_data(snps = snps, filename = "GWAS_PAD_MAF0.01", sep = "\t",
                              snp_col = "ID",
                              beta_col = "BETA",
                              se_col = "SE",
                              effect_allele_col = "ALLELE1",
                              other_allele_col = "ALLELE0",
                              eaf_col = "A1FREQ",
                              pval_col = "P",
                              units_col = "",
                              gene_col = "",
                              samplesize_col = "",
                              phenotype_col="TRAIT",
                              chr_col="CHROM",
                              pos_col="GENPOS")


# Meta-analysis for PAD and IHD and storing results in lists ----------------------------------------

#Data frame with all aNPS to harmonize 
dat_harmonize<-read_exposure_data("GWAS_remnant_cholesterol_lowpval",sep="\t",
                                   snp_col = "ID",
                                   beta_col = "BETA",
                                   se_col = "SE",
                                   effect_allele_col = "ALLELE1",
                                   other_allele_col = "ALLELE0",
                                   eaf_col = "A1FREQ",
                                   pval_col = "P",
                                   units_col = "",
                                   gene_col = "",
                                   samplesize_col = "",
                                   phenotype_col="TRAIT",
                                   chr_col="CHROM",
                                   pos_col="GENPOS")

#Harmonizing to have all alleles on same strand
PAD_ukbb_harm <- harmonise_data(dat_harmonize, PAD_ukbb)
PAD_finngen_harm <- harmonise_data(dat_harmonize, PAD_finngen)
PAD_MVP_EUR_harm <- harmonise_data(dat_harmonize, MVP_PAD_EUR)

library(metafor)

colnames<-colnames(PAD_finngen_harm)[colnames(PAD_finngen_harm) %in%colnames(PAD_ukbb_harm)  ]
# Fixed-effect meta-analysis of each SNP
PAD_meta<-PAD_ukbb_harm
PAD_meta[,colnames(PAD_meta)[colnames(PAD_meta) %in%colnames(PAD_ukbb_harm) ]]<-PAD_ukbb_harm[,colnames(PAD_meta)[colnames(PAD_meta) %in%colnames(PAD_ukbb_harm) ]]
for ( i in PAD_ukbb_harm$SNP){

dat<-rbind("UK Biobank"=subset(PAD_ukbb_harm,SNP==i)[,colnames],
           "Finngen"=subset(PAD_finngen_harm,SNP==i)[,colnames],
           "MVP"=subset(PAD_MVP_EUR_harm,SNP==i)[,colnames])

model <- rma(yi=beta.outcome, sei=se.outcome,data=dat, method="FE")

PAD_meta[PAD_meta$SNP==i,"beta.outcome"]<-model$beta[1]
PAD_meta[PAD_meta$SNP==i,"se.outcome"]<-model$se[1]
PAD_meta[PAD_meta$SNP==i,"pval.outcome"]<-model$pval[1]
}

rm(PAD_finngen_harm,PAD_bbj_harm,PAD_MVP_EUR_harm)


#Harmonizing to have all alleles on same strand
CAD_ukbb_harm <- harmonise_data(dat_harmonize, IHD_ukbb)
CAD_finngen_harm <- harmonise_data(dat_harmonize, IHD_finngen)
CAD_MVP_EUR_harm <- harmonise_data(dat_harmonize, MVP_CAD_EUR)
CAD_C4D_harm <- harmonise_data(dat_harmonize, IHD_C4D)

library(metafor)

colnames<-colnames(CAD_ukbb_harm)[colnames(CAD_ukbb_harm) %in%colnames(CAD_MVP_EUR_harm)  ]

# Fixed-effect meta-analysis of each SNP
IHD_meta<-CAD_ukbb_harm
IHD_meta[,colnames(IHD_meta)[colnames(IHD_meta) %in%colnames(CAD_ukbb_harm) ]]<-CAD_ukbb_harm[,colnames(IHD_meta)[colnames(IHD_meta) %in%colnames(CAD_ukbb_harm) ]]
for ( i in CAD_ukbb_harm$SNP){
  
  dat<-rbind("UK Biobank"=subset(CAD_ukbb_harm,SNP==i)[,colnames],
             "Finngen"=subset(CAD_finngen_harm,SNP==i)[,colnames],
             "MVP"=subset(CAD_MVP_EUR_harm,SNP==i)[colnames],
             "CARDIoGRAMplusC4D"=subset(CAD_C4D_harm,SNP==i)[,colnames])
  
  model <- rma(yi=beta.outcome, sei=se.outcome,data=dat, method="FE")
  
  IHD_meta[IHD_meta$SNP==i,"beta.outcome"]<-model$beta[1]
  IHD_meta[IHD_meta$SNP==i,"se.outcome"]<-model$se[1]
  IHD_meta[IHD_meta$SNP==i,"pval.outcome"]<-model$pval[1]
}

rm(CAD_ukbb_harm,CAD_finngen_harm,CAD_MVP_EUR_harm,CAD_C4D_harm)

PAD_meta<-PAD_meta[,colnames(PAD_meta)[colnames(PAD_meta) %in% colnames(IHD_ukbb)]]
IHD_meta<-IHD_meta[,colnames(IHD_meta)[colnames(PAD_meta) %in% colnames(IHD_ukbb)]]


# Storing all data sets in list
exposures<-list("remnant_c_gene"=ukbb_remnant_c_gene,"ldl_c_gene"=ukbb_ldl_c_gene,
                "remnant_c_wide"=ukbb_remnant_c,"ldl_c_wide"=ukbb_ldl_c)

rm(ukbb_remnant_c_gene,ukbb_ldl_c_gene,ukbb_remnant_c,ukbb_ldl_c)

PAD<-list("ukbb"=PAD_ukbb,"finngen"=PAD_finngen,"MVP_EUR"=MVP_PAD_EUR,"MVP_AFR"=MVP_PAD_AFR,"MVP_HIS"=MVP_PAD_HIS,"bbj"=PAD_bbj,"meta"=PAD_meta)

IHD<-list("ukbb"=IHD_ukbb,"finngen"=IHD_finngen,"MVP_EUR"=MVP_CAD_EUR,"C4D"=IHD_C4D,"meta"=IHD_meta)

IHD_PAD<-list("PAD"=PAD_meta,"IHD"=IHD_meta)

rm(PAD_ukbb,PAD_finngen,PAD_bbj,PAD_meta)
rm(IHD_C4D)


# Saving datasets
saveRDS(exposures,"exposures.RDS")
saveRDS(PAD,"PAD.RDS")
saveRDS(IHD,"IHD.RDS")
saveRDS(IHD_PAD,"IHD_PAD.RDS")


# Preparing data from IEU database -------------------------------------------------------


# Extrracting exposure data from IEU database
exposures_mrbase<-list()

remnant_genes<-c('LPL','LIPC','GPIHBP1','APOC3','APOA5','ANGPTL3','ANGPTL4','ANGPTL8')
ldl_genes<-c('PCSK9','HMGCR','LDLR','NPC1L1','APOB',"ACLY")

# Extracting SNPs from genes for online dataset
a<-fread("genes_key.txt")

# For triglyceride gene-specific SNPs
r<-c()
for(i in remnant_genes){
  r[i]<-paste0(subset(a,gene==i)$chr, ":",
               subset(a,gene==i)$start, "-",
               subset(a,gene==i)$end)
}

#Associations with trig
trig_gene<-ld_clump(associations(r, "ieu-a-302"),clump_r2 = 0.1,clump_p =5e-6 )

#Associations with LDL cholesterol for same genes
trig_gene<-rbind(                                trig_gene[,colnames(trig_gene)[colnames(trig_gene) !=  "pval"]],
                                                 associations(trig_gene$rsid, "ieu-a-300")[,colnames(trig_gene)[colnames(trig_gene)!=  "pval"]])

trig_gene<-format_data(trig_gene,snp_col="rsid",effect_allele_col = "ea",other_allele_col = "nea",
                       phenotype_col ="trait",samplesize_col = "n",pos_col = "position")


exposures_mrbase[["trig_gene"]]<-trig_gene
# For LDL gene-specific SNPs
r<-c()
for(i in ldl_genes){
  r[i]<-paste0(subset(a,gene==i)$chr, ":",
               subset(a,gene==i)$start, "-",
               subset(a,gene==i)$end)
}

#Associations with LDL
ldl_gene<-ld_clump(associations(r, "ieu-a-300"),clump_r2 = 0.1,clump_p =5e-6 )

#Associations with LDL cholesterol for same genes
ldl_gene<-rbind(                                ldl_gene[,colnames(ldl_gene)[colnames(ldl_gene) != "pval"]],
                                                associations(ldl_gene$rsid, "ieu-a-302")[,colnames(ldl_gene)[colnames(ldl_gene) != "pval"]])

ldl_gene<-format_data(ldl_gene,snp_col="rsid",effect_allele_col = "ea",other_allele_col = "nea",
                      phenotype_col ="trait",samplesize_col = "n",pos_col = "position")


exposures_mrbase[["ldl_gene"]]<-ldl_gene

# Reading outcome data

# All SNPs to be extracted
snps<-unique(c(exposures_mrbase[["trig_gene"]]$SNP,exposures_mrbase[["ldl_gene"]]$SNP))

#Data frame with all aNPS to harmonize 
dat_harmonize<-read_exposure_data("GWAS_remnant_cholesterol_lowpval",sep="\t",
                                  snp_col = "ID",
                                  beta_col = "BETA",
                                  se_col = "SE",
                                  effect_allele_col = "ALLELE1",
                                  other_allele_col = "ALLELE0",
                                  eaf_col = "A1FREQ",
                                  pval_col = "P",
                                  units_col = "",
                                  gene_col = "",
                                  samplesize_col = "",
                                  phenotype_col="TRAIT",
                                  chr_col="CHROM",
                                  pos_col="GENPOS")

#Harmonizing to have all alleles on same strand
PAD_ukbb_harm <- harmonise_data(dat_harmonize, PAD_ukbb)
PAD_finngen_harm <- harmonise_data(dat_harmonize, PAD_finngen)
PAD_MVP_EUR_harm <- harmonise_data(dat_harmonize, MVP_PAD_EUR)

library(metafor)

colnames<-colnames(PAD_finngen_harm)[colnames(PAD_finngen_harm) %in%colnames(PAD_ukbb_harm)  ]
# Fixed-effect meta-analysis of each SNP
PAD_meta<-PAD_ukbb_harm
PAD_meta[,colnames(PAD_meta)[colnames(PAD_meta) %in%colnames(PAD_ukbb_harm) ]]<-PAD_ukbb_harm[,colnames(PAD_meta)[colnames(PAD_meta) %in%colnames(PAD_ukbb_harm) ]]
for ( i in PAD_ukbb_harm$SNP){
  
  dat<-rbind("UK Biobank"=subset(PAD_ukbb_harm,SNP==i)[,colnames],
             "Finngen"=subset(PAD_finngen_harm,SNP==i)[,colnames],
             "MVP"=subset(PAD_MVP_EUR_harm,SNP==i)[,colnames])
  
  model <- rma(yi=beta.outcome, sei=se.outcome,data=dat, method="FE")
  
  PAD_meta[PAD_meta$SNP==i,"beta.outcome"]<-model$beta[1]
  PAD_meta[PAD_meta$SNP==i,"se.outcome"]<-model$se[1]
  PAD_meta[PAD_meta$SNP==i,"pval.outcome"]<-model$pval[1]
}

rm(PAD_finngen_harm,PAD_bbj_harm,PAD_MVP_EUR_harm)


#Harmonizing to have all alleles on same strand
CAD_ukbb_harm <- harmonise_data(dat_harmonize, IHD_ukbb)
CAD_finngen_harm <- harmonise_data(dat_harmonize, IHD_finngen)
CAD_MVP_EUR_harm <- harmonise_data(dat_harmonize, MVP_CAD_EUR)
CAD_C4D_harm <- harmonise_data(dat_harmonize, IHD_C4D)

library(metafor)

colnames<-colnames(CAD_ukbb_harm)[colnames(CAD_ukbb_harm) %in%colnames(CAD_MVP_EUR_harm)  ]

# Fixed-effect meta-analysis of each SNP
IHD_meta<-CAD_ukbb_harm
IHD_meta[,colnames(IHD_meta)[colnames(IHD_meta) %in%colnames(CAD_ukbb_harm) ]]<-CAD_ukbb_harm[,colnames(IHD_meta)[colnames(IHD_meta) %in%colnames(CAD_ukbb_harm) ]]
for ( i in CAD_ukbb_harm$SNP){
  
  dat<-rbind("UK Biobank"=subset(CAD_ukbb_harm,SNP==i)[,colnames],
             "Finngen"=subset(CAD_finngen_harm,SNP==i)[,colnames],
             "MVP"=subset(CAD_MVP_EUR_harm,SNP==i)[colnames],
             "CARDIoGRAMplusC4D"=subset(CAD_C4D_harm,SNP==i)[,colnames])
  
  model <- rma(yi=beta.outcome, sei=se.outcome,data=dat, method="FE")
  
  IHD_meta[IHD_meta$SNP==i,"beta.outcome"]<-model$beta[1]
  IHD_meta[IHD_meta$SNP==i,"se.outcome"]<-model$se[1]
  IHD_meta[IHD_meta$SNP==i,"pval.outcome"]<-model$pval[1]
}

rm(CAD_ukbb_harm,CAD_finngen_harm,CAD_MVP_EUR_harm,CAD_C4D_harm)

# Deleting redundant columns
PAD_meta<-PAD_meta[,colnames(PAD_meta)[colnames(PAD_meta) %in% colnames(IHD_ukbb)]]
IHD_meta<-IHD_meta[,colnames(IHD_meta)[colnames(IHD_meta) %in% colnames(IHD_ukbb)]]

outcomes_mrbase<-list("PAD"=PAD_meta,
                      "IHD"=IHD_meta)

#Saving data
saveRDS(exposures_mrbase,"exposures_mrbase.RDS")
saveRDS(outcomes_mrbase,"outcomes_mrbase.RDS")


# Reading prepared data ---------------------------------------------------


# Loading datasets
exposures<-readRDS("exposures.RDS")
PAD<-readRDS("PAD.RDS")
IHD<-readRDS("IHD.RDS")
IHD_PAD<-readRDS("IHD_PAD.RDS")


#Steiger filtering genome-wide scores against outcomes (no SNPs were removed )
PAD_steiger<-PAD[["meta"]]

PAD_steiger$units<-"log odds"
PAD_steiger$prevalence.outcome<-with(PAD_steiger,cases/(controls+cases))
PAD_steiger<-dplyr::rename(PAD_steiger,ncase.outcome=cases,ncontrol.outcome=controls)

IHD_steiger<-IHD[["meta"]]

IHD_steiger$units<-"log odds"
IHD_steiger$prevalence.outcome<-with(IHD_steiger,cases/(controls+cases))
IHD_steiger<-dplyr::rename(IHD_steiger,ncase.outcome=cases,ncontrol.outcome=controls)

# Removing SNPs 
nrow(steiger_filtering(harmonise_data(dplyr::rename(exposures[[3]],samplesize.exposure=sample.size),PAD_steiger))) 
nrow(steiger_filtering(harmonise_data(dplyr::rename(exposures[[3]],samplesize.exposure=sample.size),IHD_steiger)))
nrow(steiger_filtering(harmonise_data(dplyr::rename(exposures[[4]],samplesize.exposure=sample.size),PAD_steiger))) 
nrow(steiger_filtering(harmonise_data(dplyr::rename(exposures[[4]],samplesize.exposure=sample.size),IHD_steiger)))


rem(PAD_steiger,IHD_steiger)


# MR for PAD and IHD  ------------------------------------------

# Plot for SNP effects on remnant and Peripheral artery disease
df<-data.frame()
df_forest_uni<-data.frame()
df_forest_multi<-data.frame()
for (q in c(names(IHD_PAD)))
  for(i in c("remnant_c_gene","ldl_c_gene")){
    
    exp_name<-grep(gsub("\\_.*","",i),unique(exposures[[i]]$exposure),ignore.case = T,value =T)
    
    data<-harmonise_data(action=1,subset(exposures[[i]],exposure==exp_name), IHD_PAD[[q]])
    data<-mr(data, method_list = c("mr_ivw"))
    
    data_mv<-data_mv<-rbind(exposures[["remnant_c_gene"]],exposures[["ldl_c_gene"]])
    data_mv$id.exposure<-data_mv$exposure
    
    data_mv<-mv_harmonise_data(harmonise_strictness =1,data_mv, IHD_PAD[[q]])
    data_mv<-mv_multiple(data_mv,pval_threshold=if(grepl("gene",i)){5e-6} else {5e-8})
    data_mv$result<-subset(data_mv$result,exposure==exp_name)
    
    res<-exp(c(data$b,data$b-(1.959*data$se),data$b+(1.959*data$se)))
    
    res_mv<-exp(c(data_mv$result$b,data_mv$result$b-(1.959*data_mv$result$se),data_mv$result$b+(1.959*data_mv$result$se)))
    
    df[paste(i,q)," "]<-if(grepl("ukbb",q)){"UK Biobank"} else if(
      grepl("fin",q)){"Finngen"}else if(
        grepl("bbj",q)){"Biobank Japan"}else if(
          grepl("meta",q)){"Meta-analysis"}else if(
            grepl("ldl",i)){"artery disease   "}else if(
              grepl("PAD",q) & grepl("remnant",i)){"Peripheral"}else if(
                grepl("IHD",q) & grepl("remnant",i)){"Coronary"}
    
    df[paste(i,q),"\n\nNo. of\ncases"]<-formatC(IHD_PAD[[q]]$cases[1],format="d", big.mark=",")
    df[paste(i,q),"\n\nNo. of\ncontrols"]<-formatC(IHD_PAD[[q]]$controls[1],format="d", big.mark=",")
    
    df[paste(i,q),"\n\nGenetic\nscore"]<-if(grepl("remnant_c",i)){"Remnant"} else {"LDL"}
    
    df[paste(i,q),
       "\n\nUnivariable
odds ratio (95% CI)"]<-sprintf("%.2f (%.2f-%.2f)",  res[1], res[2], res[3])
    
    #Creating column with space for forest plot
    df[paste(i,q),"    "]<-paste(rep(" ", 20), collapse = " ")
    
    df[paste(i,q),"\n\n\n P"]<-"                   "
    df_forest_uni[paste(i,q),"\n\n P"]<-data$pval
    
    
    df[paste(i,q),
       " 
Multivariable
(remnant + LDL)
odds ratio (95% CI)"]<-sprintf("%.2f (%.2f-%.2f)",  res_mv[1], res_mv[2], res_mv[3])
    
    df[paste(i,q),"     "]<-paste(rep(" ", 20), collapse = " ")
    
    df[paste(i,q)," \n\n\n P"]<-"            "
    df_forest_multi[paste(i,q)," \n\n\n P"]<-data_mv$result$pval
    
    df_forest_uni[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-res
    df_forest_multi[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-res_mv
    
    
  }


#Setting theme for plot
tm <- forest_theme(base_size = 9,
                   ci_pch = 20,ci_lwd = 2, ci_Theight = 0.4,
                   xaxis_lwd=2,
                   title_col = "gray40",title_fontface = "plain",
                   core=list(
                     bg_params=list(fill = c("white","white","gray90","gray90"))))

# Fromatting x labels


xlab_1<-"                           Univariable odds ratio (95% CI)"

xlab_2<-"                           Multivariable (remnant + LDL) odds ratio (95% CI)"

p<-forest(df,
          est = list(df_forest_uni$Effect,df_forest_multi$Effect),
          lower =list(df_forest_uni$`Lower 0.95`,df_forest_multi$`Lower 0.95`), 
          upper = list(df_forest_uni$`Upper 0.95`,df_forest_multi$`Upper 0.95`), 
          sizes = 1,
          ref_line = 1,
          ci_column = c(6,9),
          xlab = c(xlab_1,xlab_2),
          xlim = c(0,6),ticks_at=c(1:6),
          title=paste0("Mendelian randomization for peripheral artery disease and coronary artery disease"),
          theme=tm)

#adding line under header
p <- add_border(p, part = "header", row = 1, where = "bottom",gp = grid::gpar(lty=1,lwd=2))

#Bold font for first column
p <- edit_plot(p, col = c(1,4), gp =  grid::gpar(fontface = "bold"))
p <- edit_plot(p, col = c(2,3), which="text" ,hjust = unit(1, "npc"),
               x = unit(0.9, "npc"))

# Changing colors
uneven<-subset(1:200,1:200 %% 2 == 1)
even<-subset(1:200, 1:200 %% 2 == 0)

colors<-c("darkcyan","darkorchid4")

#For text
p <- edit_plot(p, col = c(4),row =uneven , which="text",gp = gpar(col = colors[1]))
p <- edit_plot(p, col = c(4),row =even , which="text",gp = gpar(col = colors[2]))

#For confidence interval
p<-edit_plot(p, col = c(6,9), row = uneven, which = "ci", gp = gpar(col = colors[1], fill = colors[1]))
p<-edit_plot(p, col = c(6,9), row = even, which = "ci", gp = gpar(col = colors[2], fill = colors[2]))

# Deleging repeated text
p <- edit_plot(p, row=c(1:nrow(df))[c(FALSE,TRUE)],col = c(2,3), which="text" ,label="")


#Adding P-values with superscript
df_forest_uni[,1]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest_uni[,1], c(-Inf, 0.10, Inf))], df_forest_uni[,1])
df_forest_multi[,1]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest_multi[,1], c(-Inf, 0.10, Inf))], df_forest_multi[,1])


#Formatting P-values 
p<-format_p(plot=p,pvals=list(df_forest_uni[,1],df_forest_multi[,1]),plot_cols=c(7,10))


p[["widths"]][[6]]<-unit(3,"cm")
p[["widths"]][[9]]<-unit(3,"cm")


doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(p)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")



# MR for PAD by cohorts ------------------------------------------

#Setting sample sizes
exposures[["remnant_c_gene"]]$sample.size<-86872
exposures[["remnant_c_wide"]]$sample.size<-86872
exposures[["ldl_c_gene"]]$sample.size<-95206
exposures[["ldl_c_wide"]]$sample.size<-95206

numbers<-c("cases","controls")
for (i in 1:length(numbers)){
PAD[["ukbb"]][,numbers[i]]<-c(7307,400784)[i]
PAD[["finngen"]][,numbers[i]]<-c(7098,206541)[i]
PAD[["MVP_EUR"]][,numbers[i]]<-c(24009,150983)[i]
PAD[["meta"]][,numbers[i]]<-c(38414,758308)[i]

IHD[["ukbb"]][,numbers[i]]<-c(34541,261984)[i]
IHD[["finngen"]][,numbers[i]]<-c(30952,187840)[i]
IHD[["MVP_EUR"]][,numbers[i]]<-c(95151,197287)[i]
IHD[["C4D"]][,numbers[i]]<-c(60801,123504)[i]
IHD[["meta"]][,numbers[i]]<-c(221445,770615)[i]

IHD_PAD[["IHD"]][,numbers[i]]<-c(221445,770615)[i]
IHD_PAD[["PAD"]][,numbers[i]]<-c(38414,758308)[i]
}



# Plot for SNP effects on remnant and Peripheral artery disease
df<-data.frame()
df_forest_uni<-data.frame()
df_forest_multi<-data.frame()
for (q in c("ukbb","finngen","MVP_EUR","meta"))
  for(i in c("remnant_c_gene","ldl_c_gene")){
    
    exp_name<-grep(gsub("\\_.*","",i),unique(exposures[[i]]$exposure),ignore.case = T,value =T)
    
    data<-harmonise_data(action=1,subset(exposures[[i]],exposure==exp_name), PAD[[q]])
    data<-mr(data, method_list = c("mr_ivw"))
    
    data_mv<-data_mv<-rbind(exposures[["remnant_c_gene"]],exposures[["ldl_c_gene"]])
    data_mv$id.exposure<-data_mv$exposure
    
    data_mv<-mv_harmonise_data(harmonise_strictness =1,data_mv, PAD[[q]])
    data_mv<-mv_multiple(data_mv,pval_threshold=if(grepl("gene",i)){5e-6} else {5e-8})
    data_mv$result<-subset(data_mv$result,exposure==exp_name)
    
    res<-exp(c(data$b,data$b-(1.959*data$se),data$b+(1.959*data$se)))
    
    res_mv<-exp(c(data_mv$result$b,data_mv$result$b-(1.959*data_mv$result$se),data_mv$result$b+(1.959*data_mv$result$se)))
    
    df[paste(i,q)," "]<-if(grepl("ukbb",q)){"UK Biobank"} else if(
                           grepl("fin",q)){"Finngen"}else if(
                           grepl("MVP_EUR",q)){"Million Veteran Program"}else if(
                           grepl("meta",q)){"Meta-analysis"}

    df[paste(i,q),"\nNo. of\ncases"]<-formatC(PAD[[q]]$cases[1],format="d", big.mark=",")
    df[paste(i,q),"\nNo. of\ncontrols"]<-formatC(PAD[[q]]$controls[1],format="d", big.mark=",")
      
    df[paste(i,q),"\nGenetic\nscore"]<-if(grepl("remnant_c",i)){"Remnant"} else {"LDL"}
    
    df[paste(i,q),"\nNo. of\nSNPs"]<-formatC(data$nsnp,format="d", big.mark=",")
    
    df[paste(i,q),
"Univariable odds ratio (95% CI) 
per 1 mmol/L (39 mg/dL)
cholesterol increment"]<-sprintf("%.2f (%.2f-%.2f)",  res[1], res[2], res[3])
    
    #Creating column with space for forest plot
    df[paste(i,q),"    "]<-paste(rep(" ", 20), collapse = " ")
    
    df[paste(i,q),"\n\nP"]<-"                   "
    df_forest_uni[paste(i,q),"\nP"]<-data$pval
    
    
    df[paste(i,q),
"Multivariable (remnant + LDL) odds 
ratio (95% CI) per 1 mmol/L (39 mg/dL)
cholesterol increment"]<-sprintf("%.2f (%.2f-%.2f)",  res_mv[1], res_mv[2], res_mv[3])
    
    df[paste(i,q),"     "]<-paste(rep(" ", 20), collapse = " ")
    
    df[paste(i,q)," \n\nP"]<-"            "
    df_forest_multi[paste(i,q)," \n\nP"]<-data_mv$result$pval
    
    df_forest_uni[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-res
    df_forest_multi[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-res_mv

    
  }


#Setting theme for plot
tm <- forest_theme(base_size = 9,
                   ci_pch = 20,ci_lwd = 2, ci_Theight = 0.4,
                   xaxis_lwd=2,
                   title_col = "gray40",title_fontface = "plain",
                   core=list(
                     bg_params=list(fill = c("white","white","gray90","gray90"))))

# Fromatting x labels


xlab_1<-"\n
                  Univariable odds ratio (95% CI)
                  per 1 mmol/L (39 mg/dL)
                  cholesterol increment"
  
xlab_2<-"\n
                  Multivariable (remnant + LDL) odds ratio (95% CI)
                  per 1 mmol/L (39 mg/dL)
                  cholesterol increment"

p<-forestploter::forest(df,
          est = list(df_forest_uni$Effect,df_forest_multi$Effect),
          lower =list(df_forest_uni$`Lower 0.95`,df_forest_multi$`Lower 0.95`), 
          upper = list(df_forest_uni$`Upper 0.95`,df_forest_multi$`Upper 0.95`), 
          sizes = 1,
          ref_line = 1,
          ci_column = c(7,10),
          xlab = c(xlab_1,xlab_2),
          xlim = c(0,6),ticks_at=c(1:6),
          title=paste0("Mendelian randomization for peripheral artery disease"),
          theme=tm)

#adding line under header
p <- add_border(p, part = "header", row = 1, where = "bottom",gp = grid::gpar(lty=1,lwd=2))

#adding line before last cohort
p <- add_border(p, part = "header", row = 7, where = "bottom",gp = grid::gpar(lty=1,lwd=1))


#Bold font for first column
p <- edit_plot(p, col = c(1,4), gp =  grid::gpar(fontface = "bold"))
p <- edit_plot(p, col = c(2:3), which="text" ,hjust = unit(1, "npc"),
               x = unit(0.9, "npc"))

# Deleging repeated text
p <- edit_plot(p, row=c(1:nrow(df))[c(FALSE,TRUE)],col = c(1:3), which="text" ,label="")


# Changing colors
uneven<-subset(1:200,1:200 %% 2 == 1)
even<-subset(1:200, 1:200 %% 2 == 0)

colors<-c("darkcyan","darkorchid4")

#For text
p <- edit_plot(p, col = c(4),row =uneven , which="text",gp = gpar(col = colors[1]))
p <- edit_plot(p, col = c(4),row =even , which="text",gp = gpar(col = colors[2]))

#For confidence interval
p<-edit_plot(p, col = c(7,10), row = uneven, which = "ci", gp = gpar(col = colors[1], fill = colors[1]))
p<-edit_plot(p, col = c(7,10), row = even, which = "ci", gp = gpar(col = colors[2], fill = colors[2]))


#Adding P-values with superscript
df_forest_uni[,1]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest_uni[,1], c(-Inf, 0.10, Inf))], df_forest_uni[,1])
df_forest_multi[,1]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest_multi[,1], c(-Inf, 0.10, Inf))], df_forest_multi[,1])

p<-format_p(plot=p,pvals=list(df_forest_uni[,1],df_forest_multi[,1]),plot_cols=c(8,11))

p[["widths"]][[2]]<-unit(4.5,"cm")
p[["widths"]][[5]]<-unit(2,"cm")
p[["widths"]][[7]]<-unit(3,"cm")
p[["widths"]][[10]]<-unit(3,"cm")


doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(p)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")



# MR for PAD and IHD by gene   ------------------------------------------


# Reading remnant cholesterol data
snps_location<-read_exposure_data("GWAS_remnant_cholesterol_lowpval",sep="\t",
                                  snp_col = "ID",
                                  beta_col = "BETA",
                                  se_col = "SE",
                                  effect_allele_col = "ALLELE1",
                                  other_allele_col = "ALLELE0",
                                  eaf_col = "A1FREQ",
                                  pval_col = "P",
                                  units_col = "",
                                  gene_col = "",
                                  samplesize_col = "",
                                  phenotype_col="TRAIT",
                                  chr_col="CHROM",
                                  pos_col="GENPOS")[,c("SNP","chr.exposure","pos.exposure")]


# No SNPs in LIPC, GPIHBP1, and ANGPTL4 present in final score, therefore these are excluded
genes_remnant<-c('LPL','LIPC','APOC3','APOA5','ANGPTL3','ANGPTL4','ANGPTL8',"All")
genes_ldl<-c('PCSK9','HMGCR','LDLR','NPC1L1','APOB',"ACLY","All")

# Plot for SNP effects on remnant and Peripheral artery disease and coronary artery disease
df<-data.frame()
df_forest_uni<-data.frame()
df_forest_multi<-data.frame()
for (q in c(names(IHD_PAD)))
  for(i in c("remnant_c_gene","ldl_c_gene"))
    for (j in if (i=="remnant_c_gene") {genes_remnant} else {genes_ldl} ){
      
      if (j=="APOC3"){ # APOC3/APOA5 are too close to differ between
        j<-c("APOC3/APOA5")
      } else if (j=="APOA5"){
        next}
      
      
      gene_select<-if(j=="All" & i=="remnant_c_gene"){genes_remnant[!genes_remnant %in% "All"]} else if(
        j=="All" & i=="ldl_c_gene"){genes_ldl[!genes_ldl %in% "All"]} else if (j=="APOC3/APOA5"){c("APOC3","APOA5")}else {j}
      
      a<-data.frame()
      
      if(gene_select[1]!= "ACLY" ){
        a<-mygene::queryMany(gene_select[!gene_select %in% "ACLY"], scopes="symbol", fields="genomic_pos_hg19", species="human")[,c("genomic_pos_hg19.chr","genomic_pos_hg19.start","genomic_pos_hg19.end")]
      }
      
      if("ACLY" %in% gene_select ){
        ACLY<-mygene::queryMany(c('ACLY' 
        ), scopes="symbol", fields="genomic_pos_hg19", species="human")[,"genomic_pos_hg19"][[1]]
        
        a<-rbind(a,data.frame(genomic_pos_hg19.chr=ACLY[1,"chr"],genomic_pos_hg19.start=as.numeric(ACLY[2,"start"]),genomic_pos_hg19.end=as.numeric(ACLY[2,"end"])))
        
      }
      
      
      # Adding 100,000 base pairs of each side of gene 
      a<-data.frame(chr=a$genomic_pos_hg19.chr,
                    start=a$genomic_pos_hg19.start-1e5,
                    end=a$genomic_pos_hg19.end+1e5)
      
      
      # Selecting variants within gene regions by chromosome
      data_gene<-data.frame()
      for (b in unique(a$chr)){
        dat<-setDT(subset(snps_location,chr.exposure==b))
        chr_snps<-dat[dat$pos.exposure %inrange% a[a$chr==b, c("start","end")]  ]
        data_gene<-rbind(data_gene,chr_snps)
      }
      
      exp_name<-grep(gsub("\\_.*","",i),unique(exposures[[i]]$exposure),ignore.case = T,value =T)
      
      data<-subset(exposures[[i]],exposure==exp_name)
      
      data<-harmonise_data(action=1,data[data$SNP %in% data_gene$SNP,], IHD_PAD[[q]])
      if(length(data$SNP)==0) {
        next} else if(length(data$SNP)>1) {
          data<- mr(data, method_list = c("mr_ivw"))} else {
            data<- mr(data)}
      
      res<-exp(c(data$b[1],data$b[1]-(1.959*data$se[1]),data$b[1]+(1.959*data$se[1])))
      
      summary_row<-(match(gene_select,genes_remnant,nomatch =F)==1 | match(gene_select,genes_ldl,nomatch =F)==1)[1] & j!="All"
      
      df[paste(i,q,j)," "]<-if(summary_row){if(grepl("ukbb",q)){"UK Biobank"} else if(
        grepl("fin",q)){"Finngen"}else if(
          grepl("bbj",q)){"Biobank Japan"}else if(
            grepl("meta",q)){"Meta-analysis"}else if(
              grepl("PAD",q)){"Peripheral artery disease"}else if(
                grepl("IHD",q)){"Coronary artery disease"}} else {""}
      
      df[paste(i,q,j),"\n\nNo. of\ncases"]<-if(summary_row){formatC(IHD_PAD[[q]]$cases[1],format="d", big.mark=",")} else {""}
      
      df[paste(i,q,j),"\n\nGenetic\nscore"]<-if(summary_row){if(grepl("remnant_c",i)){"Remnant"} else {"LDL"}} else {""}
      
      df[paste(i,q,j),"\n\n\nGene"]<-j
      
      df[paste(i,q,j),"\n\nNo. of\nSNPs"]<-data$nsnp
      
      if(nrow(data)==0) next
      
      df[paste(i,q,j),
         "\n\Univariablen\nodds ratio (95% CI)"]<-sprintf("%.2f (%.2f-%.2f)",  res[1], res[2], res[3])
      
      #Creating column with space for forest plot
      df[paste(i,q,j),"    "]<-paste(rep(" ", 40), collapse = " ")
      
      df[paste(i,q,j),"\n\n\nP"]<-"                   "
      
      
      df_forest_uni[paste(i,q,j),"\n\nP"]<-data$pval
      df_forest_uni[paste(i,q,j),c("Effect","Lower 0.95","Upper 0.95")]<-res
      
      
      
    }


#Setting theme for plot
tm <- forest_theme(base_size = 9,
                   ci_pch = 20,ci_lwd = 2, ci_Theight = 0.4,
                   xaxis_lwd=2,
                   title_col = "gray40",title_fontface = "plain",
                   core=list(
                     bg_params=list(fill = c(rep("white",1),rep("gray90",1)))))

# Fromatting x labels


xlab_1<-"
                                    Univariable odds ratio (95% CI)"


# Only using univariable model for now
p<-forest(df,
          est = list(df_forest_uni$Effect),
          lower =list(df_forest_uni$`Lower 0.95`), 
          upper = list(df_forest_uni$`Upper 0.95`), 
          sizes = 1,
          ref_line = 1,
          ci_column = c(7),
          xlab = c(xlab_1),
          xlim = c(0,10),ticks_at=c(1:10),
          title=paste0("Mendelian randomization for peripheral artery disease and coronary artery disease by gene-specific SNPs"),
          theme=tm)



#adding line under header and etween categories
p <- add_border(p, part = "header", row = c(1,14), where = "bottom",gp = grid::gpar(lty=1,lwd=2))

p <- add_border(p, part = "header", col=c(3:8),row = c(8,14,21), where = "bottom",gp = grid::gpar(lty=1,lwd=2))

#adding line between categories
p <- add_border(p, part = "header",col=c(4:8), row = c(7,13,20,26), where = "bottom",gp = grid::gpar(lty=1,lwd=1))

# White background first columns
p <- edit_plot(p, col = c(1:3), which="background" ,gp=gpar(fill = "white"))

#Bold font for first column
p <- edit_plot(p, col = c(1,3,4), gp =  grid::gpar(fontface = "bold"))

# Deleging repeated text
p <- edit_plot(p, row=c(8,21),col = c(1:2), which="text" ,label="")

# Changing colors

colors<-rep(c(rep("Darkcyan",7),rep("Darkorchid4",6)),99)

for (i in 1:nrow(df)){
  #For text
  p <- edit_plot(p, col = c(3),row=i , which="text",gp = gpar(col = colors[i]))
  
  #For confidence interval
  p<-edit_plot(p, col = c(7),row=i, which = "ci", gp = gpar(col = colors[i], fill = colors[i]))
}


#Adding P-values with superscript
df_forest_uni[,1]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest_uni[,1], c(-Inf, 0.10, Inf))], df_forest_uni[,1])

#Formatting P-values 
p<-format_p(plot=p,pvals=list(df_forest_uni[,1]),plot_cols=c(8))

p[["widths"]][[2]]<-unit(4.5,"cm")
p[["widths"]][[7]]<-unit(4,"cm")
p

doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(p)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")

# Tables for SNPs in genetic scores ----------------------------------------------

# Table for gene-based scores

## first getting genes for all SNPS


# Reading remnant cholesterol data
data_genes<-rbind(read_exposure_data("GWAS_remnant_cholesterol_lowpval",sep="\t",
                                     snp_col = "ID",
                                     chr_col="CHROM",
                                     pos_col="GENPOS"),
                  read_exposure_data("GWAS_ldl_direct_lowpval",sep="\t",
                                     snp_col = "ID",
                                     chr_col="CHROM",
                                     pos_col="GENPOS"))[,c("SNP","chr.exposure","pos.exposure")]

# Getting positions of genes
a<-mygene::queryMany(c('LPL','LIPC','GPIHBP1','APOC3','APOA5','ANGPTL3','ANGPTL4','ANGPTL8',
                       'PCSK9','HMGCR','LDLR','NPC1L1','APOB'# Triglyceride metabolism
                       #'CETP',"APOA1"				# HDL metabolism
), scopes="symbol", fields="genomic_pos_hg19", species="human")[,c("query","genomic_pos_hg19.chr","genomic_pos_hg19.start","genomic_pos_hg19.end")]

# Adding 100,000 base pairs of each side of gene 
a<-data.frame(gene=a$query,
              chr=a$genomic_pos_hg19.chr,
              start=a$genomic_pos_hg19.start-1e5,
              end=a$genomic_pos_hg19.end+1e5)

# ACLY gne had patch and  had to be formatted differently
a<-rbind(a,cbind("gene"="ACLY",mygene::queryMany(c('ACLY' 
), scopes="symbol", fields="genomic_pos_hg19", species="human")[,"genomic_pos_hg19"][[1]][2,c("chr","start","end")]))
a[a$gene=="ACLY","chr"]<-17



for( i in a$gene){ # for each gene, find which SNPS are in region
  
  dat<-subset(a,gene==i)
  
  data_genes[data_genes$chr.exposure==dat$chr & data_genes$pos.exposure %between% dat[,c("start","end")] ,i]<-i
  
}

data_genes$gene<-apply(data_genes[,genes_key$gene],1,paste0,collapse=", ") # formatting genes into single column
data_genes$gene<-str_remove_all(data_genes$gene,"NA, ") #removing NA genes
data_genes$gene<-str_remove_all(data_genes$gene,", NA")
data_genes<-data_genes[!duplicated(data_genes$SNP),c("SNP","gene")] # removing duplicated SNPs


#Remnant genetic score
dat<-harmonise_data(action=1,subset(exposures[["remnant_c_gene"]],exposure=="Remnant cholesterol"),IHD_PAD[["PAD"]])

dat_ldl<-harmonise_data(action=1,subset(exposures[["remnant_c_gene"]],exposure=="LDL cholesterol"),IHD_PAD[["PAD"]])

table_remant_gene<-with(dat, #Creating table with relevant information from data frames
                        data.table("Genetic score"="Remnant","SNP"=SNP,"Gene"=data_genes[match(dat$SNP,data_genes$SNP),]$gene,"Effect allele frequency"=eaf.exposure,
                                   "Beta.Remnant cholesterol"=beta.exposure,"SE.Remnant cholesterol"=se.exposure,"Pval.Remnant cholesterol"=pval.exposure,
                                   "Beta.LDL cholesterol"=dat_ldl$beta.exposure,"SE.LDL cholesterol"=dat_ldl$se.exposure,"Pval.LDL cholesterol"=dat_ldl$pval.exposure))


dat_IHD<-harmonise_data(action=1,subset(exposures[["remnant_c_gene"]],exposure=="Remnant cholesterol"),IHD_PAD[["IHD"]])[,c("SNP","beta.outcome","se.outcome","pval.outcome")]

dat_IHD<-merge(dat[,c("SNP","beta.exposure")],dat_IHD,by="SNP",all=T)

table_remant_gene<-cbind(table_remant_gene, #Adding estimates for outcomes
                         "Beta.PAD"=dat$beta.outcome,"SE.PAD"=dat$se.outcome,"Pval.PAD"=dat$pval.outcome,
                         "Beta.CAD"=dat_IHD$beta.outcome,"SE.CAD"=dat_IHD$se.outcome,"Pval.CAD"=dat_IHD$pval.outcome)


# same operation again but for LDL genetic score
dat<-harmonise_data(action=1,subset(exposures[["ldl_c_gene"]],exposure=="Remnant cholesterol"),IHD_PAD[["PAD"]])

dat_ldl<-harmonise_data(action=1,subset(exposures[["ldl_c_gene"]],exposure=="LDL cholesterol"),IHD_PAD[["PAD"]])

table_ldl_gene<-with(dat,
                     data.table("Genetic score"="LDL","SNP"=SNP,"Gene"=data_genes[match(dat$SNP,data_genes$SNP),]$gene,"Effect allele frequency"=eaf.exposure,
                                "Beta.Remnant cholesterol"=beta.exposure,"SE.Remnant cholesterol"=se.exposure,"Pval.Remnant cholesterol"=pval.exposure,
                                "Beta.LDL cholesterol"=dat_ldl$beta.exposure,"SE.LDL cholesterol"=dat_ldl$se.exposure,"Pval.LDL cholesterol"=dat_ldl$pval.exposure))



dat_IHD<-harmonise_data(action=1,subset(exposures[["ldl_c_gene"]],exposure=="Remnant cholesterol"),IHD_PAD[["IHD"]])[,c("SNP","beta.outcome","se.outcome","pval.outcome")]

dat_IHD<-merge(dat[,c("SNP","beta.exposure")],dat_IHD[dat_IHD$SNP %in% dat$SNP ,],by="SNP",all=T)


table_ldl_gene<-cbind(table_ldl_gene,
                      "Beta.PAD"=dat$beta.outcome,"SE.PAD"=dat$se.outcome,"Pval.PAD"=dat$pval.outcome,
                      "Beta.CAD"=dat_IHD$beta.outcome,"SE.CAD"=dat_IHD$se.outcome,"Pval.CAD"=dat_IHD$pval.outcome)

table_gene<-rbind(table_remant_gene,table_ldl_gene)


table_gene$`Genetic score`<-factor(table_gene$`Genetic score`,ordered = T,levels = c("Remnant","LDL"))

table_gene$Gene<-factor(table_gene$Gene,ordered = T,levels = c("LPL","LIPC","APOA5","APOC3","APOC3, APOA5", "ANGPTL3","ANGPTL4","ANGPTL8, LDLR",
                                                               "PCSK9","HMGCR","LDLR","NPC1L1","APOB"))


table_gene<-dplyr::arrange(table_gene,`Genetic score`,Gene)
table_gene

writexl::write_xlsx(table_gene,"Temp.xlsx")


### Same operations again but for genome-wide scores, now without genes
dat<-harmonise_data(action=1,subset(exposures[["remnant_c_wide"]],exposure=="Remnant cholesterol"),IHD_PAD[["PAD"]])

dat_ldl<-harmonise_data(action=1,subset(exposures[["remnant_c_wide"]],exposure=="LDL cholesterol"),IHD_PAD[["PAD"]])

table_remant<-with(dat,
                   data.table("Genetic score"="Remnant","SNP"=SNP,"Effect allele frequency"=eaf.exposure,
                              "Beta.Remnant cholesterol"=beta.exposure,"SE.Remnant cholesterol"=se.exposure,"Pval.Remnant cholesterol"=pval.exposure,
                              "Beta.LDL cholesterol"=dat_ldl$beta.exposure,"SE.LDL cholesterol"=dat_ldl$se.exposure,"Pval.LDL cholesterol"=dat_ldl$pval.exposure))


dat_IHD<-harmonise_data(action=1,subset(exposures[["remnant_c_wide"]],exposure=="Remnant cholesterol"),IHD_PAD[["IHD"]])[,c("SNP","beta.outcome","se.outcome","pval.outcome")]

dat_IHD<-merge(dat[,c("SNP","beta.exposure")],dat_IHD,by="SNP",all=T)

table_remant<-cbind(table_remant,
                    "Beta.PAD"=dat$beta.outcome,"SE.PAD"=dat$se.outcome,"Pval.PAD"=dat$pval.outcome,
                    "Beta.CAD"=dat_IHD$beta.outcome,"SE.CAD"=dat_IHD$se.outcome,"Pval.CAD"=dat_IHD$pval.outcome)


# LDL genetic score
dat<-harmonise_data(action=1,subset(exposures[["ldl_c_wide"]],exposure=="Remnant cholesterol"),IHD_PAD[["PAD"]])

dat_ldl<-harmonise_data(action=1,subset(exposures[["ldl_c_wide"]],exposure=="LDL cholesterol"),IHD_PAD[["PAD"]])

dat<-merge(dat_ldl[,c("SNP","mr_keep")],dat,by="SNP",all=T)

table_ldl<-with(dat,
                data.table("Genetic score"="LDL","SNP"=SNP,"Effect allele frequency"=eaf.exposure,
                           "Beta.Remnant cholesterol"=beta.exposure,"SE.Remnant cholesterol"=se.exposure,"Pval.Remnant cholesterol"=pval.exposure,
                           "Beta.LDL cholesterol"=dat_ldl$beta.exposure,"SE.LDL cholesterol"=dat_ldl$se.exposure,"Pval.LDL cholesterol"=dat_ldl$pval.exposure))



dat_IHD<-harmonise_data(action=1,subset(exposures[["ldl_c_wide"]],exposure=="Remnant cholesterol"),IHD_PAD[["IHD"]])[,c("SNP","beta.outcome","se.outcome","pval.outcome")]

dat_IHD<-merge(dat[,c("SNP","beta.exposure")],dat_IHD[dat_IHD$SNP %in% dat$SNP ,],by="SNP",all=T)


table_ldl<-cbind(table_ldl,
                 "Beta.PAD"=dat$beta.outcome,"SE.PAD"=dat$se.outcome,"Pval.PAD"=dat$pval.outcome,
                 "Beta.CAD"=dat_IHD$beta.outcome,"SE.CAD"=dat_IHD$se.outcome,"Pval.CAD"=dat_IHD$pval.outcome)

table_wide<-rbind(table_remant,table_ldl)


table_wide$`Genetic score`<-factor(table_wide$`Genetic score`,ordered = T,levels = c("Remnant","LDL"))

table_wide<-dplyr::arrange(table_wide,`Genetic score`)

writexl::write_xlsx(table_wide,"Temp.xlsx")




# MR for IHD different cohorts  ------------------------------------------

# Plot for SNP effects on remnant and Peripheral artery disease
df<-data.frame()
df_forest_uni<-data.frame()
df_forest_multi<-data.frame()
for (q in c(names(IHD)))
  for(i in c("remnant_c_gene","ldl_c_gene")){
    
    exp_name<-grep(gsub("\\_.*","",i),unique(exposures[[i]]$exposure),ignore.case = T,value =T)
    
    data<-harmonise_data(action=1,subset(exposures[[i]],exposure==exp_name), IHD[[q]])
    data<-mr(data, method_list = c("mr_ivw"))
    
    data_mv<-data_mv<-rbind(exposures[["remnant_c_gene"]],exposures[["ldl_c_gene"]])
    data_mv$id.exposure<-data_mv$exposure
    
    data_mv<-mv_harmonise_data(harmonise_strictness =1,data_mv, IHD[[q]])
    data_mv<-mv_multiple(data_mv,pval_threshold=if(grepl("gene",i)){5e-6} else {5e-8})
    data_mv$result<-subset(data_mv$result,exposure==exp_name)
    
    res<-exp(c(data$b,data$b-(1.959*data$se),data$b+(1.959*data$se)))
    
    res_mv<-exp(c(data_mv$result$b,data_mv$result$b-(1.959*data_mv$result$se),data_mv$result$b+(1.959*data_mv$result$se)))
    
    df[paste(i,q)," "]<-if(grepl("ukbb",q)){"UK Biobank"} else if(
      grepl("fin",q)){"Finngen"}else if(
        grepl("bbj",q)){"Biobank Japan"}else if(
          grepl("meta",q)){"Meta-analysis"}else if(
            grepl("PAD",q)){"Peripheral artery disease"}else if(
              grepl("MVP",q)){"Million Veteran Program"} else if (
                grepl("C4D",q)) {"CARDIoGRAMplusC4D"}
    
    df[paste(i,q),"\n\nNo. of\ncases"]<-formatC(IHD[[q]]$cases[1],format="d", big.mark=",")
    df[paste(i,q),"\n\nNo. of\ncontrols"]<-formatC(IHD[[q]]$controls[1],format="d", big.mark=",")
    
    df[paste(i,q),"\n\nGenetic\nscore"]<-if(grepl("remnant_c",i)){"Remnant"} else {"LDL"}
    
    df[paste(i,q),"\n\nNo. of\nSNPs"]<-formatC(data$nsnp,format="d", big.mark=",")
    
    df[paste(i,q),
       "\nUnivariable odds ratio (95% CI)
per 1 mmol/L (39 mg/dL)
cholesterol increment"]<-sprintf("%.2f (%.2f-%.2f)",  res[1], res[2], res[3])
    
    #Creating column with space for forest plot
    df[paste(i,q),"    "]<-paste(rep(" ", 20), collapse = " ")
    
    df[paste(i,q),"\n\n\n P"]<-"                   "
    df_forest_uni[paste(i,q),"\n\n P"]<-data$pval
    
    
    df[paste(i,q),
       "\nMultivariable (remnant + LDL) odds ratio (95% CI)
per 1 mmol/L (39 mg/dL)
cholesterol increment"]<-sprintf("%.2f (%.2f-%.2f)",  res_mv[1], res_mv[2], res_mv[3])
    
    df[paste(i,q),"     "]<-paste(rep(" ", 20), collapse = " ")
    
    df[paste(i,q)," \n\n\n P"]<-"            "
    df_forest_multi[paste(i,q)," \n\n\n P"]<-data_mv$result$pval
    
    df_forest_uni[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-res
    df_forest_multi[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-res_mv
    
    
  }


#Setting theme for plot
tm <- forest_theme(base_size = 9,
                   ci_pch = 20,ci_lwd = 2, ci_Theight = 0.4,
                   xaxis_lwd=2,
                   title_col = "gray40",title_fontface = "plain",
                   core=list(
                     bg_params=list(fill = c("white","white","gray90","gray90"))))

# Fromatting x labels


xlab_1<-"
                  Univariable odds ratio (95% CI) per
                  1 mmol/L (39 mg/dL) cholesterol increment"

xlab_2<-"
                  Multivariable (remnant + LDL) odds ratio (95% CI) per
                  1 mmol/L (39 mg/dL) cholesterol increment"

p<-forest(df,
          est = list(df_forest_uni$Effect,df_forest_multi$Effect),
          lower =list(df_forest_uni$`Lower 0.95`,df_forest_multi$`Lower 0.95`), 
          upper = list(df_forest_uni$`Upper 0.95`,df_forest_multi$`Upper 0.95`), 
          sizes = 1,
          ref_line = 1,
          ci_column = c(7,10),
          xlab = c(xlab_1,xlab_2),
          xlim = c(0,6),ticks_at=c(1:6),
          title=paste0("Mendelian randomization for coronary artery disease"),
          theme=tm)

#adding line under header
p <- add_border(p, part = "header", row = 1, where = "bottom",gp = grid::gpar(lty=1,lwd=2))

#adding line before last cohort
p <- add_border(p, part = "header", row = 9, where = "bottom",gp = grid::gpar(lty=1,lwd=1))


#Bold font for first column
p <- edit_plot(p, col = c(1,4), gp =  grid::gpar(fontface = "bold"))
p <- edit_plot(p, col = c(2,3), which="text" ,hjust = unit(1, "npc"),
               x = unit(0.9, "npc"))

# Changing colors
uneven<-subset(1:200,1:200 %% 2 == 1)
even<-subset(1:200, 1:200 %% 2 == 0)

colors<-c("darkcyan","darkorchid4")

#For text
p <- edit_plot(p, col = c(4),row =uneven , which="text",gp = gpar(col = colors[1]))
p <- edit_plot(p, col = c(4),row =even , which="text",gp = gpar(col = colors[2]))

#For confidence interval
p<-edit_plot(p, col = c(7,10), row = uneven, which = "ci", gp = gpar(col = colors[1], fill = colors[1]))
p<-edit_plot(p, col = c(7,10), row = even, which = "ci", gp = gpar(col = colors[2], fill = colors[2]))

# Deleging repeated text
p <- edit_plot(p, row=c(1:nrow(df))[c(F,T)],col = c(1:2), which="text" ,label="")


#Adding P-values with superscript
df_forest_uni[,1]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest_uni[,1], c(-Inf, 0.10, Inf))], df_forest_uni[,1])
df_forest_multi[,1]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest_multi[,1], c(-Inf, 0.10, Inf))], df_forest_multi[,1])

#Formatting P-values 
p<-format_p(plot=p,pvals=list(df_forest_uni[,1],df_forest_multi[,1]),plot_cols=c(8,11))

p[["widths"]][[2]]<-unit(4.5,"cm")
p[["widths"]][[6]]<-unit(3,"cm")
p[["widths"]][[9]]<-unit(3,"cm")

doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(p)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")




# MR for PAD by cohorts Egger ------------------------------------------


# Plot for SNP effects on remnant and Peripheral artery disease
df<-data.frame()
df_forest_uni<-data.frame()
df_forest_multi<-data.frame()
for (q in c("ukbb","finngen","MVP_EUR","meta"))
  for(i in c("remnant_c_gene","ldl_c_gene")){
    
    exp_name<-grep(gsub("\\_.*","",i),unique(exposures[[i]]$exposure),ignore.case = T,value =T)
    
    data<-harmonise_data(action=1,subset(exposures[[i]],exposure==exp_name), PAD[[q]])
    data<-mr(data, method_list = c("mr_egger_regression"))
    
    data_mv<-data_mv<-rbind(exposures[["remnant_c_gene"]],exposures[["ldl_c_gene"]])
    data_mv$id.exposure<-data_mv$exposure
    
    data_mv<-mv_harmonise_data(harmonise_strictness =1,data_mv, PAD[[q]])
    
    num<-grep(data_mv$expname[data_mv$expname$exposure==exp_name,"id.exposure"],colnames(data_mv$exposure_beta)) # getting integer of place of exposure in multivariable data
    
    data_mv<-MendelianRandomization::mr_mvegger(
      with(data_mv,MendelianRandomization::mr_mvinput(
        bx = exposure_beta, bxse = exposure_se, by = outcome_beta, byse = outcome_se,snps = names(exposure_se[,1]))))
    
    res<-exp(c(data$b[1],data$b[1]-(1.959*data$se[1]),data$b[1]+(1.959*data$se[1])))
    
    res_mv<-exp(c(data_mv@Estimate[num],data_mv@Estimate[num]-(1.959*data_mv@StdError.Est[num]),data_mv@Estimate[num]+(1.959*data_mv@StdError.Est[num])))
    
    df[paste(i,q)," "]<-if(grepl("ukbb",q)){"UK Biobank"} else if(
      grepl("fin",q)){"Finngen"}else if(
        grepl("bbj",q)){"Biobank Japan"}else if(
          grepl("MVP",q)){"Million Veteran Program"}else if(
          grepl("meta",q)){"Meta-analysis"}
    
    df[paste(i,q),"\n\nNo. of\ncases"]<-formatC(PAD[[q]]$cases[1],format="d", big.mark=",")
    df[paste(i,q),"\n\nNo. of\ncontrols"]<-formatC(PAD[[q]]$controls[2],format="d", big.mark=",")
    
    df[paste(i,q),"\n\nGenetic\nscore"]<-if(grepl("remnant_c",i)){"Remnant"} else {"LDL"}
    
    df[paste(i,q),"\n\nNo. of\nSNPs"]<-formatC(data$nsnp,format="d", big.mark=",")
    
    df[paste(i,q),
       "\nUnivariable Egger odds ratio (95% CI)  
per 1 mmol/L (39 mg/dL) 
cholesterol increment"]<-sprintf("%.2f (%.2f-%.2f)",  res[1], res[2], res[3])
    
    #Creating column with space for forest plot
    df[paste(i,q),"    "]<-paste(rep(" ", 20), collapse = " ")
    
    df[paste(i,q),"\n\n\nP"]<-"                   "
    df_forest_uni[paste(i,q),"\n\nP"]<-data$pval
    
    
    df[paste(i,q),
       "\nMultivariable (remnant + LDL) Egger odds ratio (95% CI) 
per 1 mmol/L (39 mg/dL) 
       cholesterol increment"]<-sprintf("%.2f (%.2f-%.2f)",  res_mv[1], res_mv[2], res_mv[3])
    
    df[paste(i,q),"     "]<-paste(rep(" ", 20), collapse = " ")
    
    df[paste(i,q)," \n\n\nP"]<-"            "
    df_forest_multi[paste(i,q)," \n\n\nP"]<-data_mv@Pvalue.Est[num]
    
    df_forest_uni[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-res
    df_forest_multi[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-res_mv
    
    
  }


#Setting theme for plot
tm <- forest_theme(base_size = 9,
                   ci_pch = 20,ci_lwd = 2, ci_Theight = 0.4,
                   xaxis_lwd=2,
                   title_col = "gray40",title_fontface = "plain",
                   core=list(
                     bg_params=list(fill = c("white","white","gray90","gray90"))))

# Fromatting x labels


xlab_1<-"
                  Univariable Egger odds ratio (95% CI)
                  per 1 mmol/L (39 mg/dL) cholesterol increment"

xlab_2<-"
                  Multivariable (remnant + LDL) Egger odds ratio (95% CI)
                  per 1 mmol/L (39 mg/dL) cholesterol increment"

p<-forest(df,
          est = list(df_forest_uni$Effect,df_forest_multi$Effect),
          lower =list(df_forest_uni$`Lower 0.95`,df_forest_multi$`Lower 0.95`), 
          upper = list(df_forest_uni$`Upper 0.95`,df_forest_multi$`Upper 0.95`), 
          sizes = 1,
          ref_line = 1,
          ci_column = c(7,10),
          xlab = c(xlab_1,xlab_2),
          xlim = c(0,6),ticks_at=c(1:6),
          title=paste0("Mendelian randomization for peripheral artery disease with Egger regression"),
          theme=tm)

#adding line under header
p <- add_border(p, part = "header", row = 1, where = "bottom",gp = grid::gpar(lty=1,lwd=2))

#adding line before last cohort
p <- add_border(p, part = "header", row = 7, where = "bottom",gp = grid::gpar(lty=1,lwd=1))


#Bold font for first column
p <- edit_plot(p, col = c(1,4), gp =  grid::gpar(fontface = "bold"))
p <- edit_plot(p, col = c(2), which="text" ,hjust = unit(1, "npc"),
               x = unit(0.9, "npc"))

# Changing colors
uneven<-subset(1:200,1:200 %% 2 == 1)
even<-subset(1:200, 1:200 %% 2 == 0)

colors<-c("darkcyan","darkorchid4")

#For text
p <- edit_plot(p, col = c(4),row =uneven , which="text",gp = gpar(col = colors[1]))
p <- edit_plot(p, col = c(4),row =even , which="text",gp = gpar(col = colors[2]))

#For confidence interval
p<-edit_plot(p, col = c(7,10), row = uneven, which = "ci", gp = gpar(col = colors[1], fill = colors[1]))
p<-edit_plot(p, col = c(7,10), row = even, which = "ci", gp = gpar(col = colors[2], fill = colors[2]))

# Deleging repeated text
p <- edit_plot(p, row=c(1:nrow(df))[c(F,T)],col = c(1:3), which="text" ,label="")


#Adding P-values with superscript
df_forest_uni[,1]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest_uni[,1], c(-Inf, 0.10, Inf))], df_forest_uni[,1])
df_forest_multi[,1]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest_multi[,1], c(-Inf, 0.10, Inf))], df_forest_multi[,1])


#Formatting P-values 
p<-format_p(plot=p,pvals=list(df_forest_uni[,1],df_forest_multi[,1]),plot_cols=c(8,11))

p[["widths"]][[2]]<-unit(4.5,"cm")
p[["widths"]][[7]]<-unit(3,"cm")
p[["widths"]][[10]]<-unit(3,"cm")


doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(p)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")


# MR for PAD and IHD without ANGPTL8 / LDLR ------------------------------------------


#Using table from from "Tables for SNPs in genetic scores" with information on genes close to SNPs
snps<-subset(table_gene,Gene=="ANGPTL8, LDLR")$SNP # Selecting SNPs in ANGTPL8 and LDLR regions


# Plot for SNP effects on remnant and Peripheral artery disease
df<-data.frame()
df_forest_uni<-data.frame()
df_forest_multi<-data.frame()
for (q in c(names(IHD_PAD)))
  for(i in c("remnant_c_gene","ldl_c_gene")){
    
    exp_name<-grep(gsub("\\_.*","",i),unique(exposures[[i]]$exposure),ignore.case = T,value =T)
    
    data<-harmonise_data(action=1,subset(exposures[[i]],exposure==exp_name), IHD_PAD[[q]])
    data<-data[data$SNP %in% snps==F, ] #Removing selected SNPs
    
    data<-mr(data, method_list = c("mr_ivw"))
    
    data_mv<-data_mv<-rbind(exposures[["remnant_c_gene"]],exposures[["ldl_c_gene"]])
    data_mv$id.exposure<-data_mv$exposure
    data_mv<-data_mv[data_mv$SNP %in% snps==F, ] #Removing selected SNPs
    
    data_mv<-mv_harmonise_data(harmonise_strictness =1,data_mv, IHD_PAD[[q]])
    data_mv<-mv_multiple(data_mv,pval_threshold=if(grepl("gene",i)){5e-6} else {5e-8})
    data_mv$result<-subset(data_mv$result,exposure==exp_name)
    
    res<-exp(c(data$b,data$b-(1.959*data$se),data$b+(1.959*data$se)))
    
    res_mv<-exp(c(data_mv$result$b,data_mv$result$b-(1.959*data_mv$result$se),data_mv$result$b+(1.959*data_mv$result$se)))
    
    df[paste(i,q)," "]<-if(grepl("ukbb",q)){"UK Biobank"} else if(
      grepl("fin",q)){"Finngen"}else if(
        grepl("bbj",q)){"Biobank Japan"}else if(
          grepl("meta",q)){"Meta-analysis"}else if(
            grepl("PAD",q)){"Peripheral artery disease   "}else if(
              grepl("IHD",q)){"Coronary artery disease   "}
    
    df[paste(i,q),"\n\nNo. of\n cases"]<-formatC(IHD_PAD[[q]]$cases[1],format="d", big.mark=",")
    
    df[paste(i,q),"\n\nGenetic\nscore"]<-if(grepl("remnant_c",i)){"Remnant"} else {"LDL"}
    
    df[paste(i,q),"\n\nNo. of\nSNPs"]<-formatC(data$nsnp,format="d", big.mark=",")
    
    df[paste(i,q),
       "
Univariable odds ratio (95% CI)
per 1 mmol/L (39 mg/dL)
cholesterol increment"]<-sprintf("%.2f (%.2f-%.2f)",  res[1], res[2], res[3])
    
    #Creating column with space for forest plot
    df[paste(i,q),"    "]<-paste(rep(" ", 20), collapse = " ")
    
    df[paste(i,q),"\n\n\n P"]<-"                   "
    df_forest_uni[paste(i,q),"\n\n P"]<-data$pval
    
    
    df[paste(i,q),
       " 
Multivariable (remnant + LDL) odds ratio (95% CI)
per 1 mmol/L (39 mg/dL)
cholesterol increment"]<-sprintf("%.2f (%.2f-%.2f)",  res_mv[1], res_mv[2], res_mv[3])
    
    df[paste(i,q),"     "]<-paste(rep(" ", 20), collapse = " ")
    
    df[paste(i,q)," \n\n\n P"]<-"            "
    df_forest_multi[paste(i,q)," \n\n\n P"]<-data_mv$result$pval
    
    df_forest_uni[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-res
    df_forest_multi[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-res_mv
    
    
  }


#Setting theme for plot
tm <- forest_theme(base_size = 9,
                   ci_pch = 20,ci_lwd = 2, ci_Theight = 0.4,
                   xaxis_lwd=2,
                   title_col = "gray40",title_fontface = "plain",
                   core=list(
                     bg_params=list(fill = c("white","white","gray90","gray90"))))

# Fromatting x labels


xlab_1<-"                           Univariable odds ratio (95% CI) per
                                    1 mmol/L (39 mg/dL) cholesterol increment"

xlab_2<-"                           Multivariable (remnant + LDL) odds ratio (95% CI) per
                                    1 mmol/L (39 mg/dL) cholesterol increment"

p<-forest(df,
          est = list(df_forest_uni$Effect,df_forest_multi$Effect),
          lower =list(df_forest_uni$`Lower 0.95`,df_forest_multi$`Lower 0.95`), 
          upper = list(df_forest_uni$`Upper 0.95`,df_forest_multi$`Upper 0.95`), 
          sizes = 1,
          ref_line = 1,
          ci_column = c(6,9),
          xlab = c(xlab_1,xlab_2),
          xlim = c(0,6),ticks_at=c(1:6),
          title=paste0("Mendelian randomization for peripheral artery disease and coronary artery disease without ANGPTL8/LDLR variants"),
          theme=tm)

#adding line under header
p <- add_border(p, part = "header", row = 1, where = "bottom",gp = grid::gpar(lty=1,lwd=2))

#Bold font for first column
p <- edit_plot(p, col = c(1,3), gp =  grid::gpar(fontface = "bold"))
p <- edit_plot(p, col = c(2), which="text" ,hjust = unit(1, "npc"),
               x = unit(0.9, "npc"))

# Changing colors
uneven<-subset(1:200,1:200 %% 2 == 1)
even<-subset(1:200, 1:200 %% 2 == 0)

colors<-c("darkcyan","darkorchid4")

#For text
p <- edit_plot(p, col = c(3),row =uneven , which="text",gp = gpar(col = colors[1]))
p <- edit_plot(p, col = c(3),row =even , which="text",gp = gpar(col = colors[2]))

#For confidence interval
p<-edit_plot(p, col = c(6,9), row = uneven, which = "ci", gp = gpar(col = colors[1], fill = colors[1]))
p<-edit_plot(p, col = c(6,9), row = even, which = "ci", gp = gpar(col = colors[2], fill = colors[2]))

# Deleging repeated text
p <- edit_plot(p, row=c(1:nrow(df))[c(F,T)],col = c(1:2), which="text" ,label="")


#Adding P-values with superscript
df_forest_uni[,1]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest_uni[,1], c(-Inf, 0.10, Inf))], df_forest_uni[,1])
df_forest_multi[,1]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest_multi[,1], c(-Inf, 0.10, Inf))], df_forest_multi[,1])


#Formatting P-values 
p<-format_p(plot=p,pvals=list(df_forest_uni[,1],df_forest_multi[,1]),plot_cols=c(7,10))


p[["widths"]][[6]]<-unit(3,"cm")
p[["widths"]][[9]]<-unit(3,"cm")


doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(p)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")


# MR for PAD and IHD different MR models   ------------------------------------------
# Plot for SNP effects on remnant and Peripheral artery disease and coronary artery disease

df<-data.frame()
df_forest_uni<-data.frame()
for (q in names(IHD_PAD))
  for(i in c("remnant_c_gene","ldl_c_gene")){
    
    
    exp_name<-grep(gsub("\\_.*","",i),unique(exposures[[i]]$exposure),ignore.case = T,value =T)
    
    data<-harmonise_data(action=1,subset(exposures[[i]],exposure==exp_name), IHD_PAD[[q]])
    data<-mr(data,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
    
    res<-exp(data.frame(data$b,data$b-(1.959*data$se),data$b+(1.959*data$se)))
    
    df[paste(i,q,data$method)," "]<-if(grepl("ukbb",q)){"UK Biobank"} else if(
      grepl("fin",q)){"Finngen"}else if(
        grepl("bbj",q)){"Biobank Japan"}else if(
          grepl("meta",q)){"Meta-analysis"}else if(
            grepl("PAD",q)){"Peripheral artery disease"}else if(
              grepl("IHD",q)){"Coronary artery disease"}
    
    df[paste(i,q,data$method),"\n\nNo. of\ncases"]<-formatC(IHD_PAD[[q]]$cases[1],format="d", big.mark=",")
    df[paste(i,q,data$method),"\n\nNo. of\ncontrols"]<-formatC(IHD_PAD[[q]]$controls[2],format="d", big.mark=",")
    
    df[paste(i,q,data$method),"\n\nGenetic\nscore"]<-if(grepl("remnant_c",i)){"Remnant"} else {"LDL"}
    
    df[paste(i,q,data$method),"\n\nNo. of\nSNPs"]<-data$nsnp
    
    
    df[paste(i,q,data$method),"\n\n\nMethod"]<-c("Original (inverse variance weighted)",data$method[2:4])
    
    if(nrow(data)==0) next
    
    df[paste(i,q,data$method),
       "\nUnivariable odds ratio (95% CI)
per 1 mmol/L (39 mg/dL)
cholesterol increment"]<-sprintf("%.2f (%.2f-%.2f)",  res[,1], res[,2], res[,3])
    
    #Creating column with space for forest plot
    df[paste(i,q,data$method),"    "]<-paste(rep(" ", 40), collapse = " ")
    
    df[paste(i,q,data$method),"\n\n\nP"]<-"                   "
    
    df_forest_uni[paste(i,q,data$method),"\n\nP"]<-data$pval
    df_forest_uni[paste(i,q,data$method),c("Effect","Lower 0.95","Upper 0.95")]<-res
    
    
  }


#Setting theme for plot
tm <- forest_theme(base_size = 9,
                   ci_pch = 20,ci_lwd = 2, ci_Theight = 0.4,
                   xaxis_lwd=2,
                   title_col = "gray40",title_fontface = "plain",
                   core=list(
                     bg_params=list(fill = c(rep("white",1),rep("gray90",1)))))

# Fromatting x labels
xlab_1<-"
                  Univariable odds ratio (95% CI) per
                  1 mmol/L (39 mg/dL) cholesterol increment"

xlab_2<-"
                  Multivariable (remnant + LDL) odds ratio (95% CI) per
                  1 mmol/L (39 mg/dL) cholesterol increment"


# Only using univariable model for now
p<-forest(df,
          est = list(df_forest_uni$Effect),
          lower =list(df_forest_uni$`Lower 0.95`), 
          upper = list(df_forest_uni$`Upper 0.95`), 
          sizes = 1,
          ref_line = 1,
          ci_column = c(8),
          xlab = c(xlab_1),
          xlim = c(0,4),ticks_at=c(1:4),
          title=paste0("Mendelian randomization for peripheral artery disease and coronary artery disease by Mendelian randomization method"),
          theme=tm)

#adding line under header and between categories
p <- add_border(p, part = "header", row = c(1,9), where = "bottom",gp = grid::gpar(lty=1,lwd=2))

p <- add_border(p, part = "header", row = c(1,5,9,13),col=4:ncol(df), where = "bottom",gp = grid::gpar(lty=1,lwd=1))


# White background first columns
p <- edit_plot(p, col = c(1:6), which="background" ,gp=gpar(fill = "white"))

#Bold font for first column
p <- edit_plot(p, col = c(1,4,6), gp =  grid::gpar(fontface = "bold"))
p <- edit_plot(p, col = c(2:4), which="text" ,hjust = unit(1, "npc"),
               x = unit(0.9, "npc"))

#deleting repeated text first columns
p <- edit_plot(p, row=c(1:nrow(df))[-c(1,5,9,13)],col = c(1:5), which="text" ,label="")
p <- edit_plot(p, row=c(1:nrow(df))[-c(1,9)],col = c(1:3), which="text" ,label="")

# Changing colors

colors<-rep(c(rep("Darkcyan",4),rep("Darkorchid4",4)),2)

for (i in 1:nrow(df)){
  #For text
  p <- edit_plot(p, col = c(4),row=i , which="text",gp = gpar(col = colors[i]))
  
  #For confidence interval
  p<-edit_plot(p, col = c(8),row=i, which = "ci", gp = gpar(col = colors[i], fill = colors[i]))
}


#Adding P-values with superscript
df_forest_uni[,1]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest_uni[,1], c(-Inf, 0.10, Inf))], df_forest_uni[,1])

#Formatting P-values 
p<-format_p(plot=p,pvals=list(df_forest_uni[,1]),plot_cols=c(9))

p[["widths"]][[7]]<-unit(7,"cm")
p[["widths"]][[8]]<-unit(3,"cm")


doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(p)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")


# MR for PAD and IHD with SNP correlation   ------------------------------------------

# Plot for SNP effects on remnant and Peripheral artery disease and coronary artery disease

df<-data.frame()
df_forest_uni<-data.frame()
for (q in names(IHD_PAD))
  for(i in c("remnant_c_gene","ldl_c_gene")){
    
    exp_name<-grep(gsub("\\_.*","",i),unique(exposures[[i]]$exposure),ignore.case = T,value =T)
    
    data<-harmonise_data(action=1,subset(exposures[[i]],exposure==exp_name), IHD_PAD[[q]])
    
    
    data_ivw<-mr(data,method_list = c("mr_ivw"))
    data_corr<-MendelianRandomization::mr_ivw(dat_to_MRInput(data,get_correlations=T)[[1]],correl=T)
    
    
    data<-data.frame(method=c("Original model","Correlation-accounted"),b=c(data_ivw$b,data_corr@Estimate),se=c(data_ivw$se,data_corr@StdError),
                     pval=c(data_ivw$pval,data_corr@Pvalue),nsnp=c(data_ivw$nsnp,data_corr@SNPs))
    
    
    res<-exp(data.frame(data$b,data$b-(1.959*data$se),data$b+(1.959*data$se)))
    
    df[paste(i,q,data$method)," "]<-if(grepl("ukbb",q)){"UK Biobank"} else if(
      grepl("fin",q)){"Finngen"}else if(
        grepl("bbj",q)){"Biobank Japan"}else if(
          grepl("meta",q)){"Meta-analysis"}else if(
            grepl("PAD",q)){"Peripheral artery disease"}else if(
              grepl("IHD",q)){"coronary artery disease"}
    
    df[paste(i,q,data$method),"\n\nNo. of\ncases"]<-formatC(IHD_PAD[[q]]$cases[1],format="d", big.mark=",")
    df[paste(i,q,data$method),"\n\nNo. of\ncontrols"]<-formatC(IHD_PAD[[q]]$controls[2],format="d", big.mark=",")
    
    df[paste(i,q,data$method),"\n\nGenetic\nscore"]<-if(grepl("remnant_c",i)){"Remnant"} else {"LDL"}
    
    df[paste(i,q,data$method),"\n\nNo. of\nSNPs"]<-data$nsnp
    
    df[paste(i,q,data$method),"\n\n\nMethod"]<-data$method
    
    if(nrow(data)==0) next
    
    df[paste(i,q,data$method),
       "\nUnivariable odds ratio (95% CI)
per 1 mmol/L (39 mg/dL)
cholesterol increment"]<-sprintf("%.2f (%.2f-%.2f)",  res[,1], res[,2], res[,3])
    
    #Creating column with space for forest plot
    df[paste(i,q,data$method),"    "]<-paste(rep(" ", 40), collapse = " ")
    
    df[paste(i,q,data$method),"\n\n\nP"]<-"                   "
    
    df_forest_uni[paste(i,q,data$method),"\n\nP"]<-data$pval
    df_forest_uni[paste(i,q,data$method),c("Effect","Lower 0.95","Upper 0.95")]<-res
    
    
  }


#Setting theme for plot
tm <- forest_theme(base_size = 9,
                   ci_pch = 20,ci_lwd = 2, ci_Theight = 0.4,
                   xaxis_lwd=2,
                   title_col = "gray40",title_fontface = "plain",
                   core=list(
                     bg_params=list(fill = c(rep("white",1),rep("gray90",1)))))

# Fromatting x labels
xlab_1<-"
                  Univariable odds ratio (95% CI) per
                  1 mmol/L (39 mg/dL) cholesterol increment"


# Only using univariable model for now
p<-forest(df,
          est = list(df_forest_uni$Effect),
          lower =list(df_forest_uni$`Lower 0.95`), 
          upper = list(df_forest_uni$`Upper 0.95`), 
          sizes = 1,
          ref_line = 1,
          ci_column = c(8),
          xlab = c(xlab_1),
          xlim = c(0,4),ticks_at=c(1:4),
          title=paste0("Mendelian randomization for peripheral artery disease and coronary artery disease by Mendelian randomization method"),
          theme=tm)

#adding line under header and between categories
p <- add_border(p, part = "header", row = c(1,5), where = "bottom",gp = grid::gpar(lty=1,lwd=2))

p <- add_border(p, part = "header", row = c(1,3,5,7),col=4:ncol(df), where = "bottom",gp = grid::gpar(lty=1,lwd=1))

# White background first columns
p <- edit_plot(p, col = c(1:5), which="background" ,gp=gpar(fill = "white"))

#Bold font for first column
p <- edit_plot(p, col = c(1,4,6), gp =  grid::gpar(fontface = "bold"))
p <- edit_plot(p, col = c(2:4), which="text" ,hjust = unit(1, "npc"),
               x = unit(0.9, "npc"))

#deleting repeated text first columns
p <- edit_plot(p, row=c(1:nrow(df))[-c(1,3,5,7)],col = c(1:5), which="text" ,label="")
p <- edit_plot(p, row=c(1:nrow(df))[-c(1,5)],col = c(1:3), which="text" ,label="")


# Changing colors
colors<-rep(c(rep("Darkcyan",2),rep("Darkorchid4",2)),2)

for (i in 1:nrow(df)){
  #For text
  p <- edit_plot(p, col = c(4),row=i , which="text",gp = gpar(col = colors[i]))
  
  #For confidence interval
  p<-edit_plot(p, col = c(8),row=i, which = "ci", gp = gpar(col = colors[i], fill = colors[i]))
}


#Adding P-values with superscript
df_forest_uni[,1]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest_uni[,1], c(-Inf, 0.10, Inf))], df_forest_uni[,1])

#Formatting P-values 

p<-format_p(plot=p,pvals=list(df_forest_uni[,1]),plot_cols=c(9))


p[["widths"]][[8]]<-unit(3,"cm")


doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(p)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")


# MR for PAD and IHD without statin users ------------------------------------------



#Reading data files
dat_remnant<-read_exposure_data("GWAS_nostatins_remnant_cholesterol_MAF0.01",sep="\t",
                                snp_col = "ID",
                                beta_col = "BETA",
                                se_col = "SE",
                                effect_allele_col = "ALLELE1",
                                other_allele_col = "ALLELE0",
                                eaf_col = "A1FREQ",
                                pval_col = "P",
                                units_col = "",
                                gene_col = "",
                                samplesize_col = "",
                                phenotype_col="TRAIT",
                                chr_col="CHROM",
                                pos_col="GENPOS")  

dat_ldl<-read_exposure_data("GWAS_nostatins_ldl_direct_MAF0.01",sep="\t",
                            snp_col = "ID",
                            beta_col = "BETA",
                            se_col = "SE",
                            effect_allele_col = "ALLELE1",
                            other_allele_col = "ALLELE0",
                            eaf_col = "A1FREQ",
                            pval_col = "P",
                            units_col = "",
                            gene_col = "",
                            samplesize_col = "",
                            phenotype_col="TRAIT",
                            chr_col="CHROM",
                            pos_col="GENPOS")  

dat_remnant$exposure<-"Remnant cholesterol" #Setting exposure names
dat_ldl$exposure<-"LDL cholesterol"

# reading in and saving data
exposures_nostatins<-list()

# Combining ldl SNPs effects on remnant and ldl
exposures_nostatins[["remnant_c_gene"]]<-rbind(dat_remnant[dat_remnant$SNP %in% exposures[["remnant_c_gene"]]$SNP,],
                                               dat_ldl[dat_ldl$SNP %in% exposures[["remnant_c_gene"]]$SNP,])


exposures_nostatins[["ldl_c_gene"]]<-rbind(dat_ldl[dat_ldl$SNP %in% exposures[["ldl_c_gene"]]$SNP,],
                                           dat_remnant[dat_remnant$SNP %in% exposures[["ldl_c_gene"]]$SNP,])

saveRDS(exposures_nostatins,"exposures_nostatins.RDS") #saving file

#reading data
exposures_nostatins<-readRDS("exposures_nostatins.RDS")


# Plot for SNP effects on remnant and Peripheral artery disease
df<-data.frame()
df_forest_uni<-data.frame()
df_forest_multi<-data.frame()
for (q in c(names(IHD_PAD)))
  for(i in c("remnant_c_gene","ldl_c_gene")){
    
    exp_name<-grep(gsub("\\_.*","",i),unique(exposures_nostatins[[i]]$exposure),ignore.case = T,value =T)
    
    data<-harmonise_data(action=1,subset(exposures_nostatins[[i]],exposure==exp_name), IHD_PAD[[q]])
    data<-mr(data, method_list = c("mr_ivw"))
    
    data_mv<-data_mv<-rbind(exposures_nostatins[["remnant_c_gene"]],exposures_nostatins[["ldl_c_gene"]])
    data_mv$id.exposure<-data_mv$exposure
    
    data_mv<-mv_harmonise_data(harmonise_strictness =1,data_mv, IHD_PAD[[q]])
    data_mv<-mv_multiple(data_mv,pval_threshold=if(grepl("gene",i)){5e-6} else {5e-8})
    data_mv$result<-subset(data_mv$result,exposure==exp_name)
    
    res<-exp(c(data$b,data$b-(1.959*data$se),data$b+(1.959*data$se)))
    
    res_mv<-exp(c(data_mv$result$b,data_mv$result$b-(1.959*data_mv$result$se),data_mv$result$b+(1.959*data_mv$result$se)))
    
    df[paste(i,q)," "]<-if(grepl("ukbb",q)){"UK Biobank"} else if(
      grepl("fin",q)){"Finngen"}else if(
        grepl("bbj",q)){"Biobank Japan"}else if(
          grepl("meta",q)){"Meta-analysis"}else if(
            grepl("PAD",q)){"Peripheral artery disease"}else if(
              grepl("IHD",q)){"coronary artery disease"}
    
    df[paste(i,q),"\n\nNo. of\n cases"]<-formatC(IHD_PAD[[q]]$cases[1],format="d", big.mark=",")
    df[paste(i,q),"\n\nNo. of\n controls"]<-formatC(IHD_PAD[[q]]$controls[2],format="d", big.mark=",")
    
    df[paste(i,q),"\n\nGenetic\nscore"]<-if(grepl("remnant_c",i)){"Remnant"} else {"LDL"}
    
    df[paste(i,q),"\n\nNo. of\n SNPs"]<-formatC(data$nsnp,format="d", big.mark=",")
    
    df[paste(i,q),
       "\nUnivariable odds ratio (95% CI) 
per 1 mmol/L (39 mg/dL)
cholesterol increment"]<-sprintf("%.2f (%.2f-%.2f)",  res[1], res[2], res[3])
    
    #Creating column with space for forest plot
    df[paste(i,q),"    "]<-paste(rep(" ", 20), collapse = " ")
    
    df[paste(i,q),"\n\n\n P"]<-"                   "
    df_forest_uni[paste(i,q),"\n\n P"]<-data$pval
    
    
    df[paste(i,q),
       "\nMultivariable (remnant + LDL) odds ratio (95% CI)
per 1 mmol/L (39 mg/dL)
cholesterol increment"]<-sprintf("%.2f (%.2f-%.2f)",  res_mv[1], res_mv[2], res_mv[3])
    
    df[paste(i,q),"     "]<-paste(rep(" ", 20), collapse = " ")
    
    df[paste(i,q)," \n\n\n P"]<-"            "
    df_forest_multi[paste(i,q)," \n\n\n P"]<-data_mv$result$pval
    
    df_forest_uni[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-res
    df_forest_multi[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-res_mv
    
    
  }


#Setting theme for plot
tm <- forest_theme(base_size = 9,
                   ci_pch = 20,ci_lwd = 2, ci_Theight = 0.4,
                   xaxis_lwd=2,
                   title_col = "gray40",title_fontface = "plain",
                   core=list(
                     bg_params=list(fill = c("white","white","gray90","gray90"))))

# Fromatting x labels


xlab_1<-"
                  Univariable odds ratio (95% CI) per
                  1 mmol/L (39 mg/dL) cholesterol increment"

xlab_2<-"
                  Multivariable (remnant + LDL) odds ratio (95% CI) per
                  1 mmol/L (39 mg/dL) cholesterol increment"

p<-forest(df,
          est = list(df_forest_uni$Effect,df_forest_multi$Effect),
          lower =list(df_forest_uni$`Lower 0.95`,df_forest_multi$`Lower 0.95`), 
          upper = list(df_forest_uni$`Upper 0.95`,df_forest_multi$`Upper 0.95`), 
          sizes = 1,
          ref_line = 1,
          ci_column = c(7,10),
          xlab = c(xlab_1,xlab_2),
          xlim = c(0,6),ticks_at=c(1:6),
          title=paste0("Mendelian randomization for peripheral artery disease and coronary artery disease excluding statin users"),
          theme=tm)

#adding line under header
p <- add_border(p, part = "header", row = 1, where = "bottom",gp = grid::gpar(lty=1,lwd=2))

#Bold font for first column
p <- edit_plot(p, col = c(1,4), gp =  grid::gpar(fontface = "bold"))
p <- edit_plot(p, col = c(2:5), which="text" ,hjust = unit(1, "npc"),
               x = unit(0.9, "npc"))

# Changing colors
uneven<-subset(1:200,1:200 %% 2 == 1)
even<-subset(1:200, 1:200 %% 2 == 0)

colors<-c("darkcyan","darkorchid4")

#For text
p <- edit_plot(p, col = c(4),row =uneven , which="text",gp = gpar(col = colors[1]))
p <- edit_plot(p, col = c(4),row =even , which="text",gp = gpar(col = colors[2]))

#For confidence interval
p<-edit_plot(p, col = c(7,10), row = uneven, which = "ci", gp = gpar(col = colors[1], fill = colors[1]))
p<-edit_plot(p, col = c(7,10), row = even, which = "ci", gp = gpar(col = colors[2], fill = colors[2]))

# Deleging repeated text
p <- edit_plot(p, row=c(1:nrow(df))[c(F,T)],col = c(1:3), which="text" ,label="")


#Adding P-values with superscript
df_forest_uni[,1]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest_uni[,1], c(-Inf, 0.10, Inf))], df_forest_uni[,1])
df_forest_multi[,1]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest_multi[,1], c(-Inf, 0.10, Inf))], df_forest_multi[,1])

#Formatting P-values 
p<-format_p(plot=p,pvals=list(df_forest_uni[,1],df_forest_multi[,1]),plot_cols=c(8,11))


p[["widths"]][[7]]<-unit(3,"cm")
p[["widths"]][[10]]<-unit(3,"cm")


doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(p)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")



# MR for PAD and IHD  with gene and genome-wide SNPs ------------------------------------------



# Plot for SNP effects on remnant and Peripheral artery disease
df<-data.frame()
df_forest_uni<-data.frame()
df_forest_multi<-data.frame()
for (q in names(IHD_PAD))
  for(i in names(exposures)){
    
    exp_name<-grep(gsub("\\_.*","",i),unique(exposures[[i]]$exposure),ignore.case = T,value =T)
    
    data<-harmonise_data(action=1,subset(exposures[[i]],exposure==exp_name), IHD_PAD[[q]])
    data<-mr(data, method_list = c("mr_ivw"))
    
    data_mv<-data_mv<-if(grepl("gene",i)) {rbind(exposures[["remnant_c_gene"]],exposures[["ldl_c_gene"]])} else {
      rbind(exposures[["remnant_c_wide"]],exposures[["ldl_c_wide"]])}
    
    data_mv$id.exposure<-data_mv$exposure
    
    data_mv<-mv_harmonise_data(harmonise_strictness =1,data_mv, IHD_PAD[[q]])
    data_mv<-mv_multiple(data_mv,pval_threshold=if(grepl("gene",i)){5e-6} else {5e-8})
    data_mv<-subset(data_mv$result,exposure==exp_name)
    
    res<-exp(c(data$b[1],data$b[1]-(1.959*data$se[1]),data$b[1]+(1.959*data$se[1])))
    
    res_mv<-exp(c(data_mv$b[1],data_mv$b[1]-(1.959*data_mv$se[1]),data_mv$b[1]+(1.959*data_mv$se[1])))
    
    df[paste(i,q)," "]<-if(grepl("ukbb",q)){"UK Biobank"} else if(
      grepl("fin",q)){"Finngen"}else if(
        grepl("bbj",q)){"Biobank Japan"}else if(
          grepl("meta",q)){"Meta-analysis"}else if(
            grepl("PAD",q)){"Peripheral artery disease"}else if(
              grepl("IHD",q)){"Coronary artery disease"}
    
    df[paste(i,q),"\n\nNo. of\ncases"]<-formatC(IHD_PAD[[q]]$cases[1],format="d", big.mark=",")
    df[paste(i,q),"\n\nNo. of\ncontrols"]<-formatC(IHD_PAD[[q]]$controls[2],format="d", big.mark=",")
    
    df[paste(i,q),"\n\nSelection\nof SNPs"]<-if(grepl("gene",i)){"Original (gene-specific)"} else {"Genome-wide"}
    
    df[paste(i,q),"\n\nGenetic\nscore"]<-if(grepl("remnant_c",i)){"Remnant"} else {"LDL"}
    
    df[paste(i,q),"\n\nNo. of\nSNPs"]<-formatC(data$nsnp,format="d", big.mark=",")
    
    
    df[paste(i,q),
       "
Univariable odds ratio (95% CI)
per 1 mmol/L (39 mg/dL)
cholesterol increment"]<-sprintf("%.2f (%.2f-%.2f)",  res[1], res[2], res[3])
    
    #Creating column with space for forest plot
    df[paste(i,q),"    "]<-paste(rep(" ", 20), collapse = " ")
    
    df[paste(i,q),"\n\n\n P"]<-"                   "
    df_forest_uni[paste(i,q),"\n\n P"]<-data$pval
    
    
    df[paste(i,q),
       " 
Multivariable (remnant + LDL) odds ratio (95% CI)
per 1 mmol/L (39 mg/dL)
cholesterol increment"]<-sprintf("%.2f (%.2f-%.2f)",  res_mv[1], res_mv[2], res_mv[3])
    
    df[paste(i,q),"     "]<-paste(rep(" ", 20), collapse = " ")
    
    df[paste(i,q)," \n\n\nP"]<-"            "
    df_forest_multi[paste(i,q)," \n\n\nP"]<-data_mv$pval[1]
    
    df_forest_uni[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-res
    df_forest_multi[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-res_mv
    
    
  }


#Setting theme for plot
tm <- forest_theme(base_size = 9,
                   ci_pch = 20,ci_lwd = 2, ci_Theight = 0.4,
                   xaxis_lwd=2,
                   title_col = "gray40",title_fontface = "plain",
                   core=list(
                     bg_params=list(fill = c("white","gray90"))))

# Fromatting x labels


xlab_1<-"
                          Univariable odds ratio (95% CI) per
                          1 mmol/L (39 mg/dL) cholesterol increment"

xlab_2<-"
                          Multivariable (remnant + LDL) odds ratio (95% CI) per
                          1 mmol/L (39 mg/dL) cholesterol increment"

p<-forest(df,
          est = list(df_forest_uni$Effect,df_forest_multi$Effect),
          lower =list(df_forest_uni$`Lower 0.95`,df_forest_multi$`Lower 0.95`), 
          upper = list(df_forest_uni$`Upper 0.95`,df_forest_multi$`Upper 0.95`), 
          sizes = 1,
          ref_line = 1,
          ci_column = c(8,11),
          xlab = c(xlab_1,xlab_2),
          xlim = c(0,6),ticks_at=c(1:6),
          title=paste0("Mendelian randomization for peripheral artery disease and coronary artery disease by SNP selection"),
          theme=tm)

#adding line under header
p <- add_border(p, part = "header", row = 1, where = "bottom",gp = grid::gpar(lty=1,lwd=2))

#adding line between categories
p <- add_border(p, part = "header",col=c(4:12), row = c(3,7), where = "bottom",gp = grid::gpar(lty=1,lwd=1))
p <- add_border(p, part = "header", row = 5, where = "bottom",gp = grid::gpar(lty=1,lwd=2))

# White background first columns
p <- edit_plot(p, col = c(1:4), which="background" ,gp=gpar(fill = "white"))

#deleting repeated text first columns
p <- edit_plot(p, row=c(2:4,6:9),col = c(1:3), which="text" ,label="")
p <- edit_plot(p, row=c(2,4,6,8),col = c(4), which="text" ,label="")

#Bold font for first column
p <- edit_plot(p, col = c(1,4,5), gp =  grid::gpar(fontface = "bold"))
p <- edit_plot(p, col = c(2:3), which="text" ,hjust = unit(1, "npc"),
               x = unit(0.9, "npc"))

# Changing colors
uneven<-subset(1:200,1:200 %% 2 == 1)
even<-subset(1:200, 1:200 %% 2 == 0)

colors<-c("darkcyan","darkorchid4")

#For text
p <- edit_plot(p, col = c(5),row =uneven , which="text",gp = gpar(col = colors[1]))
p <- edit_plot(p, col = c(5),row =even , which="text",gp = gpar(col = colors[2]))

#For confidence interval
p<-edit_plot(p, col = c(8,11), row = uneven, which = "ci", gp = gpar(col = colors[1], fill = colors[1]))
p<-edit_plot(p, col = c(8,11), row = even, which = "ci", gp = gpar(col = colors[2], fill = colors[2]))

#Adding P-values with superscript
df_forest_uni[,1]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest_uni[,1], c(-Inf, 0.10, Inf))], df_forest_uni[,1])
df_forest_multi[,1]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest_multi[,1], c(-Inf, 0.10, Inf))], df_forest_multi[,1])


#Formatting P-values 
p<-format_p(plot=p,pvals=list(df_forest_uni[,1],df_forest_multi[,1]),plot_cols=c(9,12))


#Change widths of some columns
p[["widths"]][[5]]<-unit(4,"cm")
p[["widths"]][[8]]<-unit(3,"cm")
p[["widths"]][[11]]<-unit(3,"cm")


doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(p)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")




# MR for LDL and triglycerides different cohorts  ------------------------------------------


snps<-unique(unlist(lapply(exposures, '[[', 'SNP')))

#Reading Million Veteran Program summary data
MVP_LDL_EUR <- read_outcome_data(
  snps = snps,
  filename = "H:/Projekt remnant cholesterol/GWAS summary/dbGaP/GWAS results/MVP.EUR.LDL.gwas.dbGAP.txt",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "Effect",
  se_col = "SE",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "EAF",
  pval_col = "PValue",
  units_col = "",
  samplesize_col = "SampleSize",
  chr_col="Chromosome",
  pos_col="Position")


MVP_TG_EUR<-read.table("H:/Projekt remnant cholesterol/GWAS summary/dbGaP/GWAS results/MVP.EUR.TG.gwas.dbGAP.txt",fill=T,comment.char = "") #Problems with table formatting so had to read in using read.table and then convert
names(MVP_TG_EUR) <- MVP_TG_EUR[1,]
MVP_TG_EUR <- MVP_TG_EUR[-1,]

MVP_TG_EUR<-MVP_TG_EUR[MVP_TG_EUR$SNP_ID %in% snps,]

MVP_TG_EUR<-format_data(MVP_TG_EUR,type="outcome", snp_col = "SNP_ID",
                        beta_col = "Effect",
                        se_col = "SE",
                        effect_allele_col = "Allele1",
                        other_allele_col = "Allele2",
                        eaf_col = "EAF",
                        pval_col = "PValue",
                        units_col = "",
                        samplesize_col = "SampleSize",
                        chr_col="Chromosome",
                        pos_col="Position")



# Reading in dataa from published estimates for LDL cholesterol and Triglycerides from UK Biobank
UKBB_TG<-associations(variants=snps,id="ieu-b-111")
UKBB_TG<-format_data(UKBB_TG,type="outcome",snp_col="rsid",beta_col="beta",se_col="se",eaf="eaf",effect_allele_col = "ea",other_allele_col = "nea",
                     pval_col="p", id_col="id",chr_col = "chr",pos_col="position")

UKBB_LDL<-associations(variants=snps,id="ieu-b-110")
UKBB_LDL<-format_data(UKBB_LDL,type="outcome",snp_col="rsid",beta_col="beta",se_col="se",eaf="eaf",effect_allele_col = "ea",other_allele_col = "nea",
                      pval_col="p", id_col="id",chr_col = "chr",pos_col="position")


lipids<-list("MVP_TG"=MVP_TG_EUR,"ukbb_TG"=UKBB_TG,
             "MVP_LDL"=MVP_LDL_EUR,"ukbb_LDL"=UKBB_LDL)


# Plot for SNP effects on remnant and Peripheral artery disease
df<-data.frame()
df_forest_uni<-data.frame()
df_forest_multi<-data.frame()
for (q in c(names(lipids)))
  for(i in c("remnant_c_gene","ldl_c_gene")){
    
    exp_name<-grep(gsub("\\_.*","",i),unique(exposures[[i]]$exposure),ignore.case = T,value =T)
    
    data<-harmonise_data(action=1,subset(exposures[[i]],exposure==exp_name), lipids[[q]])
    data<-mr(data, method_list = c("mr_ivw"))
    
    res<-(c(data$b,data$b-(1.959*data$se),data$b+(1.959*data$se)))
    

    df[paste(i,q),"\n\nCohort"]<-if(grepl("ukbb",q)){"UK Biobank"} else if(
      grepl("fin",q)){"Finngen"}else if(
        grepl("bbj",q)){"Biobank Japan"}else if(
          grepl("meta",q)){"Meta-analysis"}else if(
            grepl("PAD",q)){"Peripheral artery disease"}else if(
              grepl("MVP",q)){"Million Veteran Program"} else if (
                grepl("IHD",q)) {"CARDIoGRAMplusC4D consortium"}
    
    df[paste(i,q),"\n\nLipid"]<-if(grepl("TG",q)){"Triglycerides"} else if(
      grepl("LDL",q)){"LDL cholesterol"}
    
    df[paste(i,q),"   "]<-" "
    
    df[paste(i,q),"\nGenetic\nscore"]<-if(grepl("remnant_c",i)){"Remnant"} else {"LDL"}
    
    df[paste(i,q),"\nNo. of\nSNPs"]<-formatC(data$nsnp,format="d", big.mark=",")
    
    df[paste(i,q),
       "\nChange in lipid,
standard deviations (95% CI)"]<-sprintf("%.2f (%.2f-%.2f)",  res[1], res[2], res[3])
    
    #Creating column with space for forest plot
    df[paste(i,q),"    "]<-paste(rep(" ", 20), collapse = " ")
    
    df[paste(i,q),"\n\n P"]<-"                   "
    df_forest_uni[paste(i,q),"\n\n P"]<-data$pval
    
    
    df_forest_uni[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-res

    
    
  }


#Setting theme for plot
tm <- forest_theme(base_size = 9,
                   ci_pch = 20,ci_lwd = 2, ci_Theight = 0.4,
                   xaxis_lwd=2,
                   title_col = "gray40",title_fontface = "plain",
                   core=list(
                     bg_params=list(fill = c("white","white","gray90","gray90"))))

# Fromatting x labels


xlab_1<-"
                  
                  Change in lipid,
                  standard deviations (95% CI)"



p<-forest(df,
          est = df_forest_uni$Effect,
          lower =df_forest_uni$`Lower 0.95`, 
          upper = df_forest_uni$`Upper 0.95`, 
          sizes = 1,
          ref_line = 0,
          ci_column = c(7),
          xlab = c(xlab_1),
          xlim = c(-2,5),ticks_at=c(-2:5),
          title=paste0("Mendelian randomization for triglycerides and LDL cholesterol"),
          theme=tm)

#adding line under header
p <- add_border(p, part = "header", row = 1, where = "bottom",gp = grid::gpar(lty=1,lwd=2))

p <- add_border(p, part = "header", row = 5, where = "bottom",gp = grid::gpar(lty=1,lwd=1))


#Bold font for first column
p <- edit_plot(p, col = c(1:4), gp =  grid::gpar(fontface = "bold"))

# Changing colors
uneven<-subset(1:200,1:200 %% 2 == 1)
even<-subset(1:200, 1:200 %% 2 == 0)

colors<-c("darkcyan","darkorchid4")

#For text
p <- edit_plot(p, col = c(4),row =uneven , which="text",gp = gpar(col = colors[1]))
p <- edit_plot(p, col = c(4),row =even , which="text",gp = gpar(col = colors[2]))

#For confidence interval
p<-edit_plot(p, col = c(7), row = uneven, which = "ci", gp = gpar(col = colors[1], fill = colors[1]))
p<-edit_plot(p, col = c(7), row = even, which = "ci", gp = gpar(col = colors[2], fill = colors[2]))

# Deleging repeated text
p <- edit_plot(p, row=c(1:nrow(df))[c(F,T)],col = c(1:3), which="text" ,label="")


#Adding P-values with superscript
df_forest_uni[,1]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest_uni[,1], c(-Inf, 0.10, Inf))], df_forest_uni[,1])


#Formatting P-values 
p<-format_p(plot=p,pvals=list(df_forest_uni[,1]),plot_cols=c(8))


p[["widths"]][[2]]<-unit(4.5,"cm")
p[["widths"]][[7]]<-unit(5,"cm")
p[["widths"]][[8]]<-unit(7,"cm")


doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(p)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")




# MR for PAD non-European ancestry ------------------------------------------

#Setting sample sizes
exposures[["remnant_c_gene"]]$sample.size<-77231
exposures[["remnant_c_wide"]]$sample.size<-77231
exposures[["ldl_c_gene"]]$sample.size<-95206
exposures[["ldl_c_wide"]]$sample.size<-95206

numbers<-c("cases","controls")
for (i in 1:length(numbers)){
  PAD[["meta"]][,numbers[i]]<-c(38414,758308)[i]
  PAD[["MVP_AFR"]][,numbers[i]]<-c(5373,42485)[i]
  PAD[["MVP_HIS"]][,numbers[i]]<-c(1925,18285)[i]
  PAD[["bbj"]][,numbers[i]]<-c(3593,208860)[i]
}


# Plot for SNP effects on remnant and Peripheral artery disease
df<-data.frame()
df_forest_uni<-data.frame()
df_forest_multi<-data.frame()
for (q in c("MVP_AFR","MVP_HIS","bbj"))
  for(i in c("remnant_c_gene","ldl_c_gene")){
    
    exp_name<-grep(gsub("\\_.*","",i),unique(exposures[[i]]$exposure),ignore.case = T,value =T)
    
    data<-harmonise_data(action=1,subset(exposures[[i]],exposure==exp_name), PAD[[q]])
    data<-mr(data, method_list = c("mr_ivw"))
    
    data_mv<-data_mv<-rbind(exposures[["remnant_c_gene"]],exposures[["ldl_c_gene"]])
    data_mv$id.exposure<-data_mv$exposure
    
    data_mv<-mv_harmonise_data(harmonise_strictness =1,data_mv, PAD[[q]])
    data_mv<-mv_multiple(data_mv,pval_threshold=if(grepl("gene",i)){5e-6} else {5e-8})
    data_mv$result<-subset(data_mv$result,exposure==exp_name)
    
    res<-exp(c(data$b,data$b-(1.959*data$se),data$b+(1.959*data$se)))
    
    res_mv<-exp(c(data_mv$result$b,data_mv$result$b-(1.959*data_mv$result$se),data_mv$result$b+(1.959*data_mv$result$se)))
    
    df[paste(i,q)," "]<-if(grepl("MVP_AFR",q)){"Million Veteran Program, African American ancestry"} else if(
      grepl("MVP_HIS",q)){"Million Veteran Program, Hispanic American ancestry"}else if(
        grepl("bbj",q)){"Biobank Japan"}else if(
          grepl("meta",q)){"Meta-analysis of European American ancestry"}
    
    df[paste(i,q),"\n\nNo. of\ncases"]<-formatC(PAD[[q]]$cases[1],format="d", big.mark=",")
    df[paste(i,q),"\n\nNo. of\ncontrols"]<-formatC(PAD[[q]]$controls[2],format="d", big.mark=",")
    
    df[paste(i,q),"\n\nGenetic\nscore"]<-if(grepl("remnant_c",i)){"Remnant"} else {"LDL"}
    
    df[paste(i,q),"\n\nNo. of\nSNPs"]<-formatC(data$nsnp,format="d", big.mark=",")
    
    df[paste(i,q),
       " 
Univariable odds ratio (95% CI) 
per 1 mmol/L (39 mg/dL)
cholesterol increment"]<-sprintf("%.2f (%.2f-%.2f)",  res[1], res[2], res[3])
    
    #Creating column with space for forest plot
    df[paste(i,q),"    "]<-paste(rep(" ", 20), collapse = " ")
    
    df[paste(i,q),"\n\n\nP"]<-"                   "
    df_forest_uni[paste(i,q),"\n\nP"]<-data$pval
    
    
    df[paste(i,q),
       "
Multivariable (remnant + LDL) odds ratio (95% CI)
per 1 mmol/L (39 mg/dL) 
cholesterol increment"]<-sprintf("%.2f (%.2f-%.2f)",  res_mv[1], res_mv[2], res_mv[3])
    
    df[paste(i,q),"     "]<-paste(rep(" ", 20), collapse = " ")
    
    df[paste(i,q)," \n\n\nP"]<-"            "
    df_forest_multi[paste(i,q)," \n\n\nP"]<-data_mv$result$pval
    
    df_forest_uni[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-res
    df_forest_multi[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-res_mv
    
    
  }


#Setting theme for plot
tm <- forest_theme(base_size = 9,
                   ci_pch = 20,ci_lwd = 2, ci_Theight = 0.4,
                   xaxis_lwd=2,
                   title_col = "gray40",title_fontface = "plain",
                   core=list(
                     bg_params=list(fill = c("white","white","gray90","gray90"))))

# Fromatting x labels


xlab_1<-"\n
                  Univariable odds ratio (95% CI) per
                  1 mmol/L (39 mg/dL) cholesterol increment"

xlab_2<-"\n
                  Multivariable (remnant + LDL) odds ratio (95% CI) per
                  1 mmol/L (39 mg/dL) cholesterol increment"

p<-forestploter::forest(df,
                        est = list(df_forest_uni$Effect,df_forest_multi$Effect),
                        lower =list(df_forest_uni$`Lower 0.95`,df_forest_multi$`Lower 0.95`), 
                        upper = list(df_forest_uni$`Upper 0.95`,df_forest_multi$`Upper 0.95`), 
                        sizes = 1,
                        ref_line = 1,
                        ci_column = c(7,10),
                        xlab = c(xlab_1,xlab_2),
                        xlim = c(0,6),ticks_at=c(1:6),
                        title=paste0("Mendelian randomization for peripheral artery disease"),
                        theme=tm)

#adding line under header
p <- add_border(p, part = "header", row = 1, where = "bottom",gp = grid::gpar(lty=1,lwd=2))

#Bold font for first column
p <- edit_plot(p, col = c(1,4), gp =  grid::gpar(fontface = "bold"))
p <- edit_plot(p, col = c(2:3), which="text" ,hjust = unit(1, "npc"),
               x = unit(0.9, "npc"))

# Deleging repeated text
p <- edit_plot(p, row=c(1:nrow(df))[c(F,T)],col = c(1:3), which="text" ,label="")


# Changing colors
uneven<-subset(1:200,1:200 %% 2 == 1)
even<-subset(1:200, 1:200 %% 2 == 0)

colors<-c("darkcyan","darkorchid4")

#For text
p <- edit_plot(p, col = c(4),row =uneven , which="text",gp = gpar(col = colors[1]))
p <- edit_plot(p, col = c(4),row =even , which="text",gp = gpar(col = colors[2]))

#For confidence interval
p<-edit_plot(p, col = c(7,10), row = uneven, which = "ci", gp = gpar(col = colors[1], fill = colors[1]))
p<-edit_plot(p, col = c(7,10), row = even, which = "ci", gp = gpar(col = colors[2], fill = colors[2]))


#Adding P-values with superscript
df_forest_uni[,1]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest_uni[,1], c(-Inf, 0.10, Inf))], df_forest_uni[,1])
df_forest_multi[,1]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest_multi[,1], c(-Inf, 0.10, Inf))], df_forest_multi[,1])

p<-format_p(plot=p,pvals=list(df_forest_uni[,1],df_forest_multi[,1]),plot_cols=c(8,11))

p[["widths"]][[2]]<-unit(9,"cm")
p[["widths"]][[7]]<-unit(3,"cm")
p[["widths"]][[10]]<-unit(3,"cm")


doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(p)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")



# MR for LDL and triglycerides non-European ancestry  ------------------------------------------


snps<-unique(unlist(lapply(exposures, '[[', 'SNP')))

#Reading Million Veteran Program summary data
MVP_LDL_AFR <- read_outcome_data(
  snps = snps,
  filename = "H:/Projekt remnant cholesterol/GWAS summary/dbGaP/GWAS results/MVP.AFR.LDL.gwas.dbGAP.txt",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "Effect",
  se_col = "SE",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "EAF",
  pval_col = "PValue",
  units_col = "",
  samplesize_col = "SampleSize",
  chr_col="Chromosome",
  pos_col="Position")


MVP_TG_AFR <- read_outcome_data(
  snps = snps,
  filename = "H:/Projekt remnant cholesterol/GWAS summary/dbGaP/GWAS results/MVP.AFR.TG.gwas.dbGAP.txt",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "Effect",
  se_col = "SE",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "EAF",
  pval_col = "PValue",
  units_col = "",
  samplesize_col = "SampleSize",
  chr_col="Chromosome",
  pos_col="Position")


MVP_LDL_HIS <- read_outcome_data(
  snps = snps,
  filename = "H:/Projekt remnant cholesterol/GWAS summary/dbGaP/GWAS results/MVP.HIS.LDL.gwas.dbGAP.txt",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "Effect",
  se_col = "SE",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "EAF",
  pval_col = "PValue",
  units_col = "",
  samplesize_col = "SampleSize",
  chr_col="Chromosome",
  pos_col="Position")


MVP_TG_HIS <- read_outcome_data(
  snps = snps,
  filename = "H:/Projekt remnant cholesterol/GWAS summary/dbGaP/GWAS results/MVP.HIS.TG.gwas.dbGAP.txt",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "Effect",
  se_col = "SE",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "EAF",
  pval_col = "PValue",
  units_col = "",
  samplesize_col = "SampleSize",
  chr_col="Chromosome",
  pos_col="Position")



MVP_TG_AFR<-read.table("H:/Projekt remnant cholesterol/GWAS summary/dbGaP/GWAS results/MVP.AFR.TG.gwas.dbGAP.txt",fill=T,comment.char = "") #Problems with table formatting so had to read in using read.table and then convert
names(MVP_TG_EUR) <- MVP_TG_EUR[1,]
MVP_TG_EUR <- MVP_TG_EUR[-1,]

MVP_TG_EUR<-MVP_TG_EUR[MVP_TG_EUR$SNP_ID %in% snps,]

MVP_TG_EUR<-format_data(MVP_TG_EUR,type="outcome", snp_col = "SNP_ID",
                        beta_col = "Effect",
                        se_col = "SE",
                        effect_allele_col = "Allele1",
                        other_allele_col = "Allele2",
                        eaf_col = "EAF",
                        pval_col = "PValue",
                        units_col = "",
                        samplesize_col = "SampleSize",
                        chr_col="Chromosome",
                        pos_col="Position")



BBJ_TG<-associations(variants=snps,id="bbj-a-55")
BBJ_TG<-format_data(BBJ_TG,type="outcome",snp_col="rsid",beta_col="beta",se_col="se",eaf="eaf",effect_allele_col = "ea",other_allele_col = "nea",
                     pval_col="p", id_col="id",chr_col = "chr",pos_col="position")

BBJ_LDL<-associations(variants=snps,id="bbj-a-31")
BBJ_LDL<-format_data(BBJ_LDL,type="outcome",snp_col="rsid",beta_col="beta",se_col="se",eaf="eaf",effect_allele_col = "ea",other_allele_col = "nea",
                      pval_col="p", id_col="id",chr_col = "chr",pos_col="position")


# Reading in dataa from published estimates for LDL cholesterol and Triglycerides from UK Biobank
UKBB_TG<-associations(variants=snps,id="ieu-b-111")
UKBB_TG<-format_data(UKBB_TG,type="outcome",snp_col="rsid",beta_col="beta",se_col="se",eaf="eaf",effect_allele_col = "ea",other_allele_col = "nea",
                     pval_col="p", id_col="id",chr_col = "chr",pos_col="position")

UKBB_LDL<-associations(variants=snps,id="ieu-b-110")
UKBB_LDL<-format_data(UKBB_LDL,type="outcome",snp_col="rsid",beta_col="beta",se_col="se",eaf="eaf",effect_allele_col = "ea",other_allele_col = "nea",
                      pval_col="p", id_col="id",chr_col = "chr",pos_col="position")


lipids<-list("MVP_AFR_TG"=MVP_TG_AFR,"MVP_HIS_TG"=MVP_TG_HIS,"bbj_TG"=BBJ_TG,"ukbb_TG"=UKBB_TG,
             "MVP_AFR_LDL"=MVP_LDL_AFR,"MVP_HIS_LDL"=MVP_LDL_HIS,"bbj_LDL"=BBJ_LDL,"ukbb_LDL"=UKBB_LDL)


# Plot for SNP effects on remnant and Peripheral artery disease
df<-data.frame()
df_forest_uni<-data.frame()
df_forest_multi<-data.frame()
for (q in c(names(lipids)))
  for(i in c("remnant_c_gene","ldl_c_gene")){
    
    exp_name<-grep(gsub("\\_.*","",i),unique(exposures[[i]]$exposure),ignore.case = T,value =T)
    
    data<-harmonise_data(action=1,subset(exposures[[i]],exposure==exp_name), lipids[[q]])
    data<-mr(data, method_list = c("mr_ivw"))
    
    res<-(c(data$b,data$b-(1.959*data$se),data$b+(1.959*data$se)))
    
    
    df[paste(i,q),"\n\nCohort"]<-if(grepl("ukbb",q)){"UK Biobank"} else if(
      grepl("fin",q)){"Finngen"}else if(
        grepl("bbj",q)){"Biobank Japan"}else if(
          grepl("meta",q)){"Meta-analysis"}else if(
            grepl("PAD",q)){"Peripheral artery disease"}else if(
              grepl("MVP_AFR",q)){"Million Veteran Program: African American"} else if (
                grepl("MVP_HIS",q)){"Million Veteran Program: Hispanic American"} else if (
                grepl("IHD",q)) {"CARDIoGRAMplusC4D consortium"}
    
    df[paste(i,q),"\n\nLipid"]<-if(grepl("TG",q)){"Triglycerides"} else if(
      grepl("LDL",q)){"LDL cholesterol"}
    
    df[paste(i,q),"   "]<-" "
    
    df[paste(i,q),"\nGenetic\nscore"]<-if(grepl("remnant_c",i)){"Remnant"} else {"LDL"}
    
    df[paste(i,q),"\nNo. of\nSNPs"]<-formatC(data$nsnp,format="d", big.mark=",")
    
    df[paste(i,q),
       "\nChange in lipid,
standard deviations (95% CI)"]<-sprintf("%.2f (%.2f-%.2f)",  res[1], res[2], res[3])
    
    #Creating column with space for forest plot
    df[paste(i,q),"    "]<-paste(rep(" ", 20), collapse = " ")
    
    df[paste(i,q),"\n\n P"]<-"                   "
    df_forest_uni[paste(i,q),"\n\n P"]<-data$pval
    
    
    df_forest_uni[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-res
    
    
    
  }


#Setting theme for plot
tm <- forest_theme(base_size = 9,
                   ci_pch = 20,ci_lwd = 2, ci_Theight = 0.4,
                   xaxis_lwd=2,
                   title_col = "gray40",title_fontface = "plain",
                   core=list(
                     bg_params=list(fill = c("white","white","gray90","gray90"))))

# Fromatting x labels


xlab_1<-"
                  
                  Change in lipid,
                  standard deviations (95% CI)"



p<-forest(df,
          est = df_forest_uni$Effect,
          lower =df_forest_uni$`Lower 0.95`, 
          upper = df_forest_uni$`Upper 0.95`, 
          sizes = 1,
          ref_line = 0,
          ci_column = c(7),
          xlab = c(xlab_1),
          xlim = c(-2,5),ticks_at=c(-2:5),
          title=paste0("Mendelian randomization for triglycerides and LDL cholesterol"),
          theme=tm)

#adding line under header
p <- add_border(p, part = "header", row = 1, where = "bottom",gp = grid::gpar(lty=1,lwd=2))

p <- add_border(p, part = "header", row = 9, where = "bottom",gp = grid::gpar(lty=1,lwd=1))


#Bold font for first column
p <- edit_plot(p, col = c(1:4), gp =  grid::gpar(fontface = "bold"))

# Changing colors
uneven<-subset(1:200,1:200 %% 2 == 1)
even<-subset(1:200, 1:200 %% 2 == 0)

colors<-c("darkcyan","darkorchid4")

#For text
p <- edit_plot(p, col = c(4),row =uneven , which="text",gp = gpar(col = colors[1]))
p <- edit_plot(p, col = c(4),row =even , which="text",gp = gpar(col = colors[2]))

#For confidence interval
p<-edit_plot(p, col = c(7), row = uneven, which = "ci", gp = gpar(col = colors[1], fill = colors[1]))
p<-edit_plot(p, col = c(7), row = even, which = "ci", gp = gpar(col = colors[2], fill = colors[2]))

# Deleging repeated text
p <- edit_plot(p, row=c(1:nrow(df))[c(F,T)],col = c(1:3), which="text" ,label="")


#Adding P-values with superscript
df_forest_uni[,1]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest_uni[,1], c(-Inf, 0.10, Inf))], df_forest_uni[,1])


#Formatting P-values 
p<-format_p(plot=p,pvals=list(df_forest_uni[,1]),plot_cols=c(8))


p[["widths"]][[2]]<-unit(9,"cm")
p[["widths"]][[7]]<-unit(5,"cm")
p[["widths"]][[8]]<-unit(7,"cm")



doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(p)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")




# MR for PAD cohorts and IHD with Global lipid genetics consortium data------------------------------------------

exposures_mrbase<-readRDS("exposures_mrbase.RDS")
outcomes_mrbase<-readRDS("outcomes_mrbase.RDS")

#Setting  sizes
for( i in 1:length(exposures_mrbase)){
exposures_mrbase[[i]][,"sample.size"]<-min(exposures_mrbase[[i]][,"samplesize.exposure"])
}

numbers<-c("cases","controls")
for (i in 1:length(numbers)){
  outcomes_mrbase[["IHD"]][,numbers[i]]<-c(221445,770615)[i]
  outcomes_mrbase[["PAD"]][,numbers[i]]<-c(38414,758308)[i]
}



# Plot for SNP effects on remnant and Peripheral artery disease
df<-data.frame()
df_forest_uni<-data.frame()
df_forest_multi<-data.frame()
for (q in names(outcomes_mrbase))
  for(i in names(exposures_mrbase)){
    
    exp_name<-grep(gsub("\\_.*","",i),unique(exposures_mrbase[[i]]$exposure),ignore.case = T,value =T)
             
    data<-harmonise_data(action=1,subset(exposures_mrbase[[i]],exposure==exp_name), outcomes_mrbase[[q]])
    data<-mr(data, method_list = c("mr_ivw"))
    
    data_mv<-rbind(exposures_mrbase[["trig_gene"]],exposures_mrbase[["ldl_gene"]])
    data_mv$id.exposure<-data_mv$exposure
    
    data_mv<-mv_harmonise_data(harmonise_strictness =1,data_mv, IHD_PAD[[q]])
    data_mv<-mv_multiple(data_mv,pval_threshold=if(grepl("gene",i)){5e-6} else {5e-8})
    data_mv$result<-subset(data_mv$result,exposure==exp_name)
    
    res<-exp(c(data$b,data$b-(1.959*data$se),data$b+(1.959*data$se)))
    
    res_mv<-exp(c(data_mv$result$b,data_mv$result$b-(1.959*data_mv$result$se),data_mv$result$b+(1.959*data_mv$result$se)))
    
    df[paste(i,q)," "]<-if(grepl("ukbb",q)){"UK Biobank"} else if(
      grepl("fin",q)){"Finngen"}else if(
        grepl("bbj",q)){"Biobank Japan"}else if(
          grepl("PAD",q)){"Peripheral artery disease"}else if(
            grepl("IHD",q)){"Coronary artery disease"}
    
    df[paste(i,q),"\nNo. of\ncases"]<-formatC(outcomes_mrbase[[q]]$cases[1],format="d", big.mark=",")
    df[paste(i,q),"\nNo. of\ncontrols"]<-formatC(outcomes_mrbase[[q]]$controls[2],format="d", big.mark=",")
    
    df[paste(i,q),paste(rep(" ",5),collapse=" ")]<-if(grepl("trig",i)){"Triglycerides"} else {"LDL cholesterol"}
    
    df[paste(i,q),"\nNo. of\nSNPs"]<-formatC(data$nsnp,format="d", big.mark=",")
    
    df[paste(i,q),
       "Univariable 
odds ratio (95% CI) per
standard deviation increment"]<-sprintf("%.2f (%.2f-%.2f)",  res[1], res[2], res[3])
    
    #Creating column with space for forest plot
    df[paste(i,q),"    "]<-paste(rep(" ", 20), collapse = " ")
    
    df[paste(i,q),"\n\nP"]<-"                   "
    df_forest_uni[paste(i,q),"\nP"]<-data$pval
    
    
    df[paste(i,q),
       "Multivariable (remnant + LDL) 
odds ratio (95% CI) per
standard deviation increment"]<-sprintf("%.2f (%.2f-%.2f)",  res_mv[1], res_mv[2], res_mv[3])
    
    df[paste(i,q),"     "]<-paste(rep(" ", 20), collapse = " ")
    
    df[paste(i,q)," \n\nP"]<-"            "
    df_forest_multi[paste(i,q)," \n\nP"]<-data_mv$result$pval[1]
    
    df_forest_uni[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-res
    df_forest_multi[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-res_mv
    
    
  }


# theme for plot
tm <- forest_theme(base_size = 9,
                   ci_pch = 20,ci_lwd = 2, ci_Theight = 0.4,
                   xaxis_lwd=2,
                   title_col = "gray40",title_fontface = "plain",
                   core=list(
                     bg_params=list(fill = c("white","white","gray90","gray90"))))

# Fromatting x labels


xlab_1<-"\n
                  Univariable odds ratio (95% CI)
                  per standard deviation increment"

xlab_2<-"\n
                  Multivariable (remnant + LDL) odds ratio (95% CI)
                  per standard deviation increment"

p<-forest(df,
          est = list(df_forest_uni$Effect,df_forest_multi$Effect),
          lower =list(df_forest_uni$`Lower 0.95`,df_forest_multi$`Lower 0.95`), 
          upper = list(df_forest_uni$`Upper 0.95`,df_forest_multi$`Upper 0.95`), 
          sizes = 1,
          ref_line = 1,
          ci_column = c(7,10),
          xlab = c(xlab_1,xlab_2),
          xlim = c(0.5,2.2),ticks_at=c(0,0.5,1,1.5,2),
          title=paste0("Mendelian randomization for peripheral artery disease with exposure data from The Global Lipids Genetics Consortium"),
          theme=tm)

#adding line under header
p <- add_border(p, part = "header", row = 1, where = "bottom",gp = grid::gpar(lty=1,lwd=2))


#Bold font for first column
p <- edit_plot(p, col = c(1,4), gp =  grid::gpar(fontface = "bold"))


# Deleging repeated text
p <- edit_plot(p, row=c(1:nrow(df))[c(F,T)],col = c(1:3), which="text" ,label="")

# Changing colors
uneven<-subset(1:200,1:200 %% 2 == 1)
even<-subset(1:200, 1:200 %% 2 == 0)

colors<-c("darkcyan","darkorchid4")

#For text
p <- edit_plot(p, col = c(4),row =uneven , which="text",gp = grid::gpar(col = colors[1]))
p <- edit_plot(p, col = c(4),row =even , which="text",gp = grid::gpar(col = colors[2]))

#For confidence interval
p<-edit_plot(p, col = c(7,10), row = uneven, which = "ci", gp = grid::gpar(col = colors[1], fill = colors[1]))
p<-edit_plot(p, col = c(7,10), row = even, which = "ci", gp = grid::gpar(col = colors[2], fill = colors[2]))


#Adding P-values with superscript
df_forest_uni[,1]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest_uni[,1], c(-Inf, 0.10, Inf))], df_forest_uni[,1])
df_forest_multi[,1]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest_multi[,1], c(-Inf, 0.10, Inf))], df_forest_multi[,1])

p<-format_p(plot=p,pvals=list(df_forest_uni[,1],df_forest_multi[,1]),plot_cols=c(8,11))

p[["widths"]][[2]]<-unit(5,"cm")
p[["widths"]][[5]]<-unit(3,"cm")
p[["widths"]][[7]]<-unit(3,"cm")
p[["widths"]][[10]]<-unit(3,"cm")

doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(p)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")



# Genes and Phenoscanner of SNPs  ------------------------------------------------

genes_key<-fread("genes_key.txt") # reading gene locations in hg19

pheno_snps<-list()
for (i in c("remnant_c_gene","ldl_c_gene")){

if (length(exposures[[i]]$SNP) >100) {
  
data<-MendelianRandomization::phenoscanner(snpquery=exposures[[i]]$SNP[1:100],
                                           pvalue=5e-6) 
data2<-MendelianRandomization::phenoscanner(snpquery=exposures[[i]]$SNP[101:length(exposures[[i]]$SNP)],
                                           pvalue=5e-6) # finding GWAS associations

data_genes<-rbind(data$snps,data2$snps)# data for SNP locations
data<-rbind(data$results,data2$results) # data for associations

} else {
  
  data<-MendelianRandomization::phenoscanner(snpquery=exposures[[i]]$SNP,
                                             pvalue=5e-6) # finding GWAS associations
  
  data_genes<-data$snps # data for SNP locations
  data<-data$results # data for associations
}


data<-data[!duplicated(data[,c("snp","trait")]),] # removing duplicated GWAS hits

 #removing repeated SNP row labels

for( i in genes_key$gene){ # for each gene, find which SNPS are in region
  
  dat<-subset(genes_key,gene==i)
  
  data_genes[data_genes$chr==dat$chr & data_genes$pos_hg19 %between% dat[,c("start","end")] ,i]<-i
  
}

data_genes$gene<-apply(data_genes[,genes_key$gene],1,paste0,collapse=", ") # formatting genes into single column
data_genes$gene<-str_remove_all(data_genes$gene,"NA, ") #removing NA genes
data_genes$gene<-str_remove_all(data_genes$gene,", NA")
data_genes<-data_genes[!duplicated(data_genes$snp),] # removing duplicated SNPs

for(i in grep("rs",unique(data$snp),value=T)){
  data[data$snp==i,"gene"]<-data_genes[data_genes$snp==i,"gene"] # for each SNP, find matching genes
}


data<-data[,c("snp","gene","trait","pmid")] # removing redunant columns

data<-rbind(data,data.frame(snp=data_genes$snp[! data_genes$snp %in% data$snp], #adding SNPs without hits in GWAS
                            gene=data_genes$gene[! data_genes$snp %in% data$snp],
                            trait="-",pmid="-"))

data<-data[order(data$snp), ] #ordering by snp
data[duplicated(data$snp),c("gene","snp")]<-" " # removing repeated labels

pheno_snps[[i]]<-flextable::flextable(data)
}

pheno_snps[[2]]


# Exporting SNPs to one-sample MR -----------------------------------------

  
  # Defining which SNP is for which PRS
exposures[["remnant_c_gene"]]$PRS<-"Remnant cholesterol gene-based"
exposures[["ldl_c_gene"]]$PRS<-"LDL cholesterol gene-based"
exposures[["remnant_c_wide"]]$PRS<-"Remnant cholesterol genome-wide"
exposures[["ldl_c_wide"]]$PRS<-"LDL cholesterol genome-wide"

  # Formatting SNPs into one file
  PRS_remnant_LDL<-rbind(exposures[["remnant_c_gene"]],exposures[["ldl_c_gene"]],
                         exposures[["remnant_c_wide"]],exposures[["ldl_c_wide"]],fill=TRUE)
  
  # Exporting SNPs by chromosome for individual-level data extraction via PLINK
  for (i in unique(PRS_remnant_LDL$chr.exposure)){
    fwrite(PRS_remnant_LDL[!duplicated(PRS_remnant_LDL$SNP) & PRS_remnant_LDL$chr.exposure==i ,"SNP"],paste0("snps_c",i,".txt"),sep = "\t",col.names=F)
  }  
  
  # Writing table with coefficients for risk score for calculation of individual PRS
  snps_location<-read_exposure_data("GWAS_remnant_cholesterol_lowpval",sep="\t",
                                    snp_col = "ID",
                                    beta_col = "BETA",
                                    se_col = "SE",
                                    effect_allele_col = "ALLELE1",
                                    other_allele_col = "ALLELE0",
                                    eaf_col = "A1FREQ",
                                    pval_col = "P",
                                    units_col = "",
                                    gene_col = "",
                                    samplesize_col = "",
                                    phenotype_col="TRAIT",
                                    chr_col="CHROM",
                                    pos_col="GENPOS")[,c("SNP","chr.exposure","pos.exposure")]
  

  # Giving locations to data set
  PRS_remnant_LDL[,c("chr.exposure","pos.exposure")]<- snps_location[ match(PRS_remnant_LDL$SNP,snps_location$SNP),c("chr.exposure","pos.exposure")]
  
  PRS_remnant_LDL[PRS_remnant_LDL$SNP %in% unique(exposures[["remnant_c_gene"]]$SNP),"Remnant cholesterol gene-based"]<-1
  PRS_remnant_LDL[PRS_remnant_LDL$SNP %in% unique(exposures[["ldl_c_gene"]]$SNP),"LDL cholesterol gene-based"]<-1
  PRS_remnant_LDL[PRS_remnant_LDL$SNP %in% unique(exposures[["remnant_c_wide"]]$SNP),"Remnant cholesterol genome-wide"]<-1
  PRS_remnant_LDL[PRS_remnant_LDL$SNP %in% unique(exposures[["ldl_c_wide"]]$SNP),"LDL cholesterol genome-wide"]<-1
  
  fwrite(PRS_remnant_LDL,"snps_key.txt",sep = "\t")

  PRS_remnant_LDL[PRS_remnant_LDL$SNP %in% unique(exposures[["remnant_c_gene"]]$SNP),]
  
  #Creating file with locations for genes (mygene package does not work in UK Biobank RAP)
  genes<-c('LPL','LIPC','GPIHBP1','APOC3','APOA5','ANGPTL3','ANGPTL4','ANGPTL8', # Triglyceride metabolism
           'PCSK9','HMGCR','LDLR','NPC1L1','APOB') 			# LDL metabolism
           
  a<-mygene::queryMany(genes, 
      scopes="symbol", fields="genomic_pos_hg19", species="human")[,c("genomic_pos_hg19.chr","genomic_pos_hg19.start","genomic_pos_hg19.end")]
  
  a<-data.frame(gene=genes,
                chr=a$genomic_pos_hg19.chr,
                start=a$genomic_pos_hg19.start,
                end=a$genomic_pos_hg19.end)

  #ACLY gene has corrected position and therfore different output format, adding after
  a<-rbind(a,data.frame(c(gene="ACLY",mygene::queryMany(c('ACLY' ),
      scopes="symbol", fields="genomic_pos_hg19", species="human")[,"genomic_pos_hg19"][[1]][2,c("chr","start","end")])))
  
  
  a<-data.frame(gene=a$gene,
                chr=a$chr,
                start=a$start-1e5,
                end=a$end+1e5)
  
  
  fwrite(a,"genes_key.txt")
  
  
  

# SNP plots ---------------------------------------------------------------
  
  # Plot for SNP effects on remnant and Peripheral artery disease
  fig<-list()
  for (q in c(names(PAD)))
    for(i in c("remnant_c_gene","ldl_c_gene","remnant_c_wide","ldl_c_wide")){
      
      harm<-harmonise_data(action=1,exposures[[i]][1:(nrow(exposures[[i]])/2),], PAD[[q]])
      data<-mr(harm, method_list = c("mr_ivw"))
      
      fig[[paste(i,q)]]<-mr_scatter_plot(mr_results=data,dat=harm)[[1]]+
        theme(legend.position = "none")+
        geom_point(size=2)+
        labs(y=paste(q,"PAD beta"),x=paste(i,"beta"))+
        coord_cartesian(ylim=c(-0.1,0.3))
      
      
    }
  fargs  <- formals(plot_grid)
  fargs$ncol           <- 4
  formals(plot_grid) <- fargs
  
  fig<-do.call(plot_grid,fig)
  
  toffice(figure = fig, format = "pptx", 
          filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,
          append = FALSE, width = 10, height = 10,devsize = FALSE, units = "in")
  
  
  # Plot for SNP effects on remnant and Peripheral artery disease versus IHD
  
  fig<-list()
  for (q in c(names(IHD_PAD)))
    for(i in c("remnant_c_gene","ldl_c_gene","remnant_c_wide","ldl_c_wide")){
      harm<-harmonise_data(action=1,exposures[[i]][1:(nrow(exposures[[i]])/2),], IHD_PAD[[q]])
      data<-mr(harm, method_list = c("mr_ivw"))
      
      fig[[paste(i,q)]]<-mr_scatter_plot(mr_results=data,dat=harm)[[1]]+
        theme(legend.position = "none")+
        geom_point(size=2)+
        labs(y=paste(q,"beta"),x=paste(i,"beta"))+
        coord_cartesian(ylim=c(-0.1,0.3))
      
      
    }
  fargs  <- formals(plot_grid)
  fargs$ncol           <- 4
  formals(plot_grid) <- fargs
  
  
  fig<-do.call(plot_grid,fig)
  
  
  toffice(figure = fig, format = "pptx", 
          filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,
          append = FALSE, width = 10, height = 4,devsize = FALSE, units = "in")
  
  
  
  
  fwrite(subset(exposures[["remnant_c_gene"]]),"remn_temp")
  fwrite(subset(exposures[["ldl_c_gene"]]),"ldl_temp")
  
  dat<-mv_extract_exposures_local(c("remn_temp","ldl_temp"),sep=",",
  snp_col = "SNP",
  phenotype_col="exposure",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure",
  pval_col = "pval.exposure",
  samplesize_col = "sample.size",pval_threshold=5e-6,clump_r2 = 0.1)
  
  
  dat<-extract_instruments(outcomes="ieu-a-302")
  data_mv<-harmonise_data(action=1,dat,extract_outcome_data(outcomes="ieu-a-7",snps=unique(dat$SNP)))
  mr(data_mv)
  
  
  dat<-mv_extract_exposures(id_exposure=c("met-d-VLDL_C","ieu-a-300","ieu-a-302"))
  data_mv<-mv_harmonise_data(harmonise_strictness =1,dat,extract_outcome_data(outcomes="ieu-a-7",snps=unique(dat$SNP)))
  mv_multiple(data_mv,pval_threshold=5e-8)
  
  
  
  # No SNPs in LIPC, GPIHBP1, and ANGPTL4 present in final score, therefore these are excluded
  genes_remnant<-c('LPL','LIPC','APOC3','APOA5','ANGPTL3','ANGPTL4','ANGPTL8',"All")
  genes_ldl<-c('PCSK9','HMGCR','LDLR','NPC1L1','APOB',"ACLY","All")
  
  for(i in c("remnant_c_gene","ldl_c_gene"))
    for (j in if (i=="remnant_c_gene") {genes_remnant} else {genes_ldl} ){
      
      
      gene_select<-if(j=="All" & i=="remnant_c_gene"){genes_remnant[!genes_remnant %in% "All"]} else if(
        j=="All" & i=="ldl_c_gene"){genes_ldl[!genes_ldl %in% "All"]} else {j}
      
      a<-data.frame()
      
      if(gene_select[1]!= "ACLY" ){
      a<-mygene::queryMany(gene_select[!gene_select %in% "ACLY"], scopes="symbol", fields="genomic_pos_hg19", species="human")[,c("genomic_pos_hg19.chr","genomic_pos_hg19.start","genomic_pos_hg19.end")]
}
      
      # Adding 100,000 base pairs of each side of gene 
      a<-data.frame(chr=a$genomic_pos_hg19.chr,
                    start=a$genomic_pos_hg19.start-1e5,
                    end=a$genomic_pos_hg19.end+1e5)
      

      if("ACLY" %in% gene_select ){
        ACLY<-mygene::queryMany(c('ACLY' 
        ), scopes="symbol", fields="genomic_pos_hg19", species="human")[,"genomic_pos_hg19"][[1]]
        
        a<-rbind(a,cbind(chr=ACLY[1,"chr"], ACLY[2,c("start","end")]))
      }
      
      
      
      # Selecting variants within gene regions by chromosome
      data_gene<-data.frame()
      for (b in unique(a$chr)){
        dat<-setDT(subset(snps_location,chr.exposure==b))
        chr_snps<-dat[dat$pos.exposure %inrange% a[a$chr==b, c("start","end")]  ]
        data_gene<-rbind(data_gene,chr_snps)
      }
      
      exp_name<-grep(gsub("\\_.*","",i),unique(exposures[[i]]$exposure),ignore.case = T,value =T)
      
      data<-subset(exposures[[i]],exposure==exp_name)
      
      data<-harmonise_data(action=1,data[data$SNP %in% data_gene$SNP,], IHD_PAD[[q]])
      if(length(data$SNP)==0) {
        next} else if(length(data$SNP)>1) {
        data<- mr(data, method_list = c("mr_ivw"))} else {
        data<- mr(data)}

      print(data)
    }

  
  
  
  
# Analysis in progress ----------------------------------------------------


  
results<-list()
for(q in c("met-d-S_LDL_P","met-d-M_LDL_P","met-d-L_LDL_P","met-d-IDL_P",
           "met-d-XS_VLDL_P","met-d-S_VLDL_P","met-d-M_VLDL_P","met-d-L_VLDL_P","met-d-XL_VLDL_P","met-d-XXL_VLDL_P"))
  for(i in c("8:19659228-19924769","19:11000038-11444492")){

    
#Associations with trig
    try(data<-ld_clump(associations(i, "ieu-b-108"),clump_r2 = 0.3,clump_p =5e-2 ))

    if(!length(data$rsid)>1){next}
data<-data[,colnames(data)!="pval"]

#Associations with LDL cholesterol for same genes


data<-format_data(data,snp_col="rsid",effect_allele_col = "ea",other_allele_col = "nea",
                       phenotype_col ="trait",samplesize_col = "n",pos_col = "position")

  data<-harmonise_data(action=1,data, extract_outcome_data(snps = data$SNP, outcomes = q))
  
  data_ivw<-MendelianRandomization::mr_ivw(dat_to_MRInput(data,get_correlations=T)[[1]],correl=T)
  data_egger<-MendelianRandomization::mr_egger(dat_to_MRInput(data,get_correlations=T)[[1]],correl=T)
  
  results[[paste(i,q,"ivw")]]<-data_ivw
  results[[paste(i,q,"egger")]]<-data_egger
  
  }
  

  
#For 
  data<-exposures_mrbase[["ldl_gene"]]
  
  a<-mygene::queryMany("ANGPTL8", scopes="symbol", fields="genomic_pos_hg19", species="human")[,c("genomic_pos_hg19.chr","genomic_pos_hg19.start","genomic_pos_hg19.end")]
  
  # Adding 100,000 base pairs of each side of gene 
  a<-data.frame(chr=a$genomic_pos_hg19.chr,
                start=a$genomic_pos_hg19.start-1e5,
                end=a$genomic_pos_hg19.end+1e5)
  
  
  # Selecting variants within gene regions by chromosome
  data_gene<-data.frame()
  for (b in unique(a$chr)){
    dat<-setDT(subset(data,chr.exposure==b))
    chr_snps<-dat[dat$pos.exposure %inrange% a[a$chr==b, c("start","end")]  ]
    data_gene<-rbind(data_gene,chr_snps)
  }
  
  
harm_data<-harmonise_data(action=1,subset(data,exposure=="LDL cholesterol"), 
                          extract_outcome_data(snps = data$SNP[data$SNP %in% data_gene$SNP], outcomes = "ieu-b-4980"))

mr(harm_data)



trig_gene<-rbind( trig_gene[,colnames(trig_gene)],
                  associations(trig_gene$rsid, "ieu-b-108")[,colnames(trig_gene)])


data_mv<-mv_harmonise_data(harmonise_strictness =1,trig_gene, extract_outcome_data(snps = trig_gene$SNP, outcomes = 'ebi-a-GCST005195'))

mv_multiple(data_mv,pval_threshold=5e-6)


data<-data.frame(SNP=c("rs964184","rs4420638","rs2075650", "rs660240","rs6589564", "rs2075290","rs2980853","rs1260326","rs7254892",
                       "rs562338","rs508487", "rs4803760","rs1178977","rs6976930"),
                 beta.exposure=c(-5.70008,-5.45295,-7.15912,3.45306,5.53037,5.48955,2.6962,-2.60009,-10.9757,-2.65678,-4.95643,
                                  2.73614,2.55226,-2.54050),
                 se.exposure=c(0.511384,0.490271,0.748784,0.432418,0.700672,0.701319,0.353669,0.366214, 1.78507,0.464554, 0.868883,0.49648,0.466442,0.464764),
                 effect_allele.exposure=c("C","A","A","C", "C","C","A","C","A","A","C","C","A","A"),
                 other_allele.exposure=c("G","G","G","T","G","T","C","T","G","G","T","T","G","G"),
                 eaf.exposure=c(0.8588,0.8264,0.8611,0.7873, 0.0723, 0.0716,0.5315,0.5856,0.0534,0.1811,0.9402,0.7692, 0.8000,0.2001),
                 pval.exposure=c(7.46e-29,9.77e-29, 1.17e-21,  1.40e-15, 2.95e-15,  4.98e-15, 8.42e-14,1.25e-12,  7.82e-10,1.07e-08,1.17e-08,  3.57e-08,4.46e-08,
                                 4.60e-08),
                
                 
                 id.exposure="sd_LDLC",
                 exposure="sd_LDLC")


data<-extract_instruments("ebi-a-GCST90092943")


library(ieugwasr)
data_ldl<-associations(variants=data$SNP,id="ieu-b-110")

data_ldl<-format_data(data_ldl,snp_col ="rsid",effect_allele_col ="ea",other_allele_col = "nea",pval_col = "p",phenotype_col ="trait")

data_ldl<-data_ldl[,colnames(data)]


data_remnant<-associations(variants=data$SNP,id="ebi-a-GCST90092943")
data_remnant<-format_data(data_remnant,snp_col ="rsid",effect_allele_col ="ea",other_allele_col = "nea",pval_col = "p",phenotype_col ="trait")
data_remnant<-data_remnant[,colnames(data)]

data_apob<-associations(variants=data$SNP,id="ebi-a-GCST90025952")
data_apob<-format_data(data_apob,snp_col ="rsid",effect_allele_col ="ea",other_allele_col = "nea",pval_col = "p",phenotype_col ="trait")
data_apob<-data_apob[,colnames(data)]


data_mv<-rbind(data,data_ldl)
data_mv2<-rbind(data,data_remnant)
data_mv3<-rbind(data,data_apob)

harm_data0<-harmonise_data(data, 
                           extract_outcome_data(snps = data$SNP, outcomes = "finn-b-I9_PAD"))

harm_data<-mv_harmonise_data(data_mv, 
                          extract_outcome_data(snps = data$SNP, outcomes = "ebi-a-GCST005195"))

harm_data2<-mv_harmonise_data(data_mv2, 
                             extract_outcome_data(snps = data$SNP, outcomes = "ebi-a-GCST005195"))
harm_data3<-mv_harmonise_data(data_mv3, 
                              extract_outcome_data(snps = data$SNP, outcomes = "ebi-a-GCST005195"))

mr(harm_data0)
mv_multiple(harm_data)
mv_multiple(harm_data2)
mv_multiple(harm_data3)



library(MRInstruments)

data_tissue<-list()
for (s in c("ieu-a-302","ieu-a-300")) 
  for (i in c("Artery Coronary","Liver") ){
  
irak1bp1_exp_dat <-
  format_gtex_eqtl(subset(
    gtex_eqtl,
    tissue == i
  ))


snps<-unique(irak1bp1_exp_dat$SNP)

intervals<-split(1:length(snps), ceiling(seq_along(1:length(snps))/1000))

data_tissue[[paste(s,i)]]<-associations(snps[1:1000], s)

for (k in 2:length(intervals)){
data_tissue[[paste(s,i)]]<-rbind(data_tissue[[paste(s,i)]],associations(snps[intervals[[k]]], s))
}


}


library(TwoSampleMR)
s<-"ieu-a-300"
i<-"Liver"

data<-
  format_gtex_eqtl(subset(
    gtex_eqtl,gene_name=="GIP!"
  ))


data<-ld_clump(data,clump_r2 = 0.01,clump_p =5e-6 )

data<-format_data(data,snp_col="rsid",effect_allele_col = "ea",other_allele_col = "nea",
                  phenotype_col ="trait",samplesize_col = "n",pos_col = "position")

data$id.exposure<-"GLP1R"
try(data<-harmonise_data(action=1,data, extract_outcome_data("finn-b-DEATH",snps = data$SNP)))

print(i)
print(mr(data))








MVP_CAD_EUR <- read_outcome_data(
  snps = snps,
  filename = "H:/Projekt remnant cholesterol/GWAS summary/dbGaP/GWAS results/MVP_CAD_Meta_analysis_EUR.withrsID.txt",
  sep = "\t",
  snp_col = "rsmid",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "EA",
  other_allele_col = "NEA",
  eaf_col = "EAF",
  pval_col = "pvalue",
  units_col = "",
  samplesize_col = "TotalSampleSize",
  chr_col="CHR",
  pos_col="POS_b37")
                
MVP_TG_EUR<-convert_outcome_to_exposure(MVP_TG_EUR)
MVP_LDL_EUR<-convert_outcome_to_exposure(MVP_LDL_EUR)



MVP_TG_EUR$exposure<-"Triglycerides"
MVP_LDL_EUR$exposure<-"LDL cholesterol"

harm_data<-harmonise_data(action=1,exposures[["remnant_c_gene"]], MVP_PAD_EUR)

mr(harm_data)

data_mv<-rbind(exposures[["remnant_c_gene"]],exposures[["ldl_c_gene"]])
data_mv$id.exposure<-data_mv$exposure

data_mv<-mv_harmonise_data(harmonise_strictness =1,data_mv, IHD_PAD[["PAD"]])

mv_multiple(data_mv,pval_threshold=5e-6)

# Conditional F-statistic
MVMR::strhet_mvmr(with(data_mv,MendelianRandomization::mr_mvinput(
  bx = exposure_beta, bxse = exposure_se, by = outcome_beta, byse = outcome_se,snps = names(exposure_se[,1]))))


# F-statistic

F   = with(exposures[["ldl_c_gene"]],beta.exposure^2/se.exposure^2)
mean(F)

