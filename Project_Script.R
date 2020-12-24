#2 Linear Regression 
library("data.table")

# read in data
phenotype_data <- fread("phenotype.txt")

# function for doing linear regression on a single snp
snp_lin <- function(rnum = 1, df = snpdata, phen_dat = phenotype_data$V3){
  #get all the allele columns from data
  snp_one <- as.character(df[rnum, -(1:4)])
  #paste0: pair up the alleles 
  #as.numeric(factor()): encode the alleles, aka 0,1 & 2
  coded_pairs <- as.numeric(factor(paste0(snp_one[c(T, F)], snp_one[c(F, T)])))
  #return p-val
  summary(lm(phen_dat ~ coded_pairs))$coefficients[2,4]
}

#create output vector
lin_outputs <- numeric(828325)

# nrows_read = number of rows read at a time
nrows_read <- 1000
for (i in 0:(828325/nrows_read)){
  # skip and nrows to load each part of the data
  snpdata <- fread("plink.tped", nrows = nrows_read, skip = i * nrows_read)
  
  # Store outputs of linear regressions to a list
  # set nrows_read to the amount of rows read at a time
  lin_outputs[(i * nrows_read + 1):(i * nrows_read + nrow(snpdata))] <- sapply(1:nrow(snpdata), snp_lin)
}

length(lin_outputs)
head(lin_outputs, 50)
tail(lin_outputs, 50)
 
#------------------------------------------------------------------------------------------------------------
#3 Draw Manhattan Plot 

# From qqman.r Script
# manhattan plot using base graphics
manhattan <- function(dataframe, colors=c("gray10", "gray50"), ymax="max", limitchromosomes=1:23, suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8), annotate=NULL, ...) {
  
  d=dataframe
  if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
  
  if (any(limitchromosomes)) d=d[d$CHR %in% limitchromosomes, ]
  d=subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
  
  numRow = dim(dataframe)[1]
  genomewideline = -log10(0.05/numRow)
  
  d$logp = -log10(d$P)
  d$pos=NA
  ticks=NULL
  lastbase=0
  colors <- rep(colors,max(d$CHR))[1:max(d$CHR)]
  if (ymax=="max") ymax<-ceiling(max(d$logp))
  if (ymax<8) ymax<-8
  
  print(max(d$logp))
  
  numchroms=length(unique(d$CHR))
  if (numchroms==1) {
    d$pos=d$BP
    ticks=floor(length(d$pos))/2+1
  } else {
    for (i in unique(d$CHR)) {
      if (i==1) {
        d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
      } else {
        lastbase=lastbase+tail(subset(d,CHR==i-1)$BP, 1)
        d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
      }
      ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
    }
  }
  
  if (numchroms==1) {
    with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",unique(d$CHR),"position"), ...))
  }	else {
    with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab="Chromosome", xaxt="n", type="n", ...))
    axis(1, at=ticks, lab=unique(d$CHR), ...)
    icol=1
    for (i in unique(d$CHR)) {
      with(d[d$CHR==i, ],points(pos, logp, col=colors[icol], ...))
      icol=icol+1
    }
  }
  
  if (!is.null(annotate)) {
    d.annotate=d[which(d$SNP %in% annotate), ]
    with(d.annotate, points(pos, logp, col="green3", ...)) 
  }
  
  #if (suggestiveline) abline(h=suggestiveline, col="blue")
  if (genomewideline) abline(h=genomewideline, col="red")
}


# Make a pretty QQ plot of p-values
qq = function(pvector, ...) {
  if (!is.numeric(pvector)) stop("D'oh! P value vector is not numeric.")
  pvector <- pvector[!is.na(pvector) & pvector<1 & pvector>0]
  o = -log10(sort(pvector,decreasing=F))
  e = -log10( ppoints(length(pvector) ))
  plot(e,o,pch=19,cex=1, xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))), xlim=c(0,max(e)), ylim=c(0,max(o)), ...)
  abline(0,1,col="red")
}


#My Code
library("data.table")
plink_linear <- fread("plink.assoc.linear")
plink_linear <- plink_linear[, !c("SNP", "A1", "TEST", "NMISS", "BETA", "STAT")]
plink_linear <- as.data.frame(plink_linear)
manhattan(plink_linear)

#------------------------------------------------------------------------------------------------------------
#4 Minimum P-Values
library("data.table")
linear <- fread("plink.assoc.linear")

#vector of unique chromosomes
chrs <- unique(linear$CHR)

#split data into chromosomes
eval(str2expression(paste0("chr", chrs, " <- linear[linear$CHR==", chrs ,"]")))
chr_list <- lapply(chrs, function(a) eval(str2expression(paste0("linear[linear$CHR==", a ,"]"))))

#find minimum p-value SNP in each chromosome
eval(str2expression(paste0("chr", chrs, "_minp <- chr", chrs, "[P == min(chr", chrs, "$P)]")))
chr_minp_list <- lapply(chrs, function(a) eval(str2expression(paste0("chr", a, "[P == min(chr", a, "$P)]"))))

#extract p-values from each SNP
eval(str2expression(paste0("min_p <- as.numeric(c(", paste0("chr", chrs, "_minp$P", sep = ",", collapse = " "), '""',"))[-(length(chrs)+1)]")))
min_p <- sapply(chr_minp_list, function(a) a$P)

#store SNP names
SNPs <- sapply(chr_minp_list, function(a) a$SNP)

R.utils::gunzip("30760_irnt.gwas.imputed_v3.both_sexes.tsv.bgz", "UK.txt")

UK <- fread("UK.txt")
#find SNP's p-value in UK Biobank
chr1_UK <- UK[grep("230301574", UK$variant),]$pval
chr2_UK = UK[grep("21206183", UK$variant),]$pval
chr3_UK = UK[grep("65454372", UK$variant),]$pval
chr4_UK = UK[grep("88909033", UK$variant),]$pval
chr5_UK = UK[grep("173583412", UK$variant),]$pval
chr6_UK = UK[grep("95301301", UK$variant),]$pval
chr7_UK = UK[grep("135863198", UK$variant),]$pval
chr8_UK = UK[grep("19870263", UK$variant),]$pval
chr9_UK = UK[grep("107661742", UK$variant),]$pval
chr10_UK = UK[grep("55127010", UK$variant),]$pval
chr11_UK = UK[grep("116664040", UK$variant),]$pval
chr12_UK = UK[grep("125412376", UK$variant),]$pval
chr13_UK = UK[grep("106835388", UK$variant),]$pval
chr14_UK = UK[grep("24194045", UK$variant),]$pval
chr15_UK = UK[grep("58723426", UK$variant),]$pval
chr16_UK = UK[grep("57006378", UK$variant),]$pval[1]
chr17_UK = UK[grep("1087242", UK$variant),]$pval[7]
chr18_UK = UK[grep("47167214", UK$variant),]$pval
chr19_UK = UK[grep("45411941", UK$variant),]$pval
chr20_UK = UK[grep("44551855", UK$variant),]$pval
chr21_UK = UK[grep("42586130", UK$variant),]$pval
chr22_UK = UK[grep("40051166", UK$variant),]$pval

#store in vector
Uk_minp <-c(chr1_UK, chr2_UK, chr3_UK, chr4_UK, chr5_UK, chr6_UK, chr7_UK, chr8_UK, chr9_UK, chr10_UK, chr11_UK, chr12_UK, chr13_UK, chr14_UK, chr15_UK, chr16_UK, chr17_UK, chr18_UK, chr19_UK, chr20_UK, chr21_UK,chr22_UK)

#create table
minp_table <- cbind(SNPs, min_p, Uk_minp)
snpdata <- as.data.frame(minp_table)
snpdata

