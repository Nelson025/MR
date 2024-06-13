#install.packages("devtools")
#devtools::install_github("mrcieu/gwasglue", force = TRUE)

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("VariantAnnotation")

#install.packages("dplyr")
#install.packages("tidyr")
#install.packages("CMplot")


#引用包
library(VariantAnnotation)
library(gwasglue)
library(dplyr)
library(tidyr)
library(CMplot)
library(TwoSampleMR)

setwd("F:\\阿尔兹海默症\\先兆子痫_AD\\01.exposure")      #设置工作目录

#############################################04.exposure----
#读取输入文件, 并对输入文件进行格式转换
#### 1.加载结局数据 ----
load('finndataR9_Preeclampsia_exposure.rdata')
data=finndata

#根据pvalue<5e-08对结果进行过滤
outTab<-subset(data, pval.exposure<5e-08)
write.csv(outTab, file="exposure.pvalue.csv", row.names=F)

#准备绘制曼哈顿图的数据
data=data[,c("SNP", "chr.exposure", "pos.exposure", "pval.exposure")]
colnames(data)=c("SNP","CHR","BP","pvalue")

#绘制线性的曼哈顿图
CMplot(data,  plot.type="m",
       LOG10=TRUE, threshold=5e-08, threshold.lwd=3, threshold.lty=1, signal.cex=0.2,
       chr.den.col=NULL, cex=0.2, bin.size=1e5, ylim=c(0,50),
       file="pdf", file.output=TRUE, width=15, height=9, verbose=TRUE)

#绘制圈图
CMplot(data,  plot.type="c",
       LOG10=TRUE, threshold=5e-08, threshold.lwd=3, threshold.lty=1, signal.cex=0.2,
       chr.den.col=NULL, cex=0.2, bin.size=1e5, ylim=c(0,100),
       file="pdf", file.output=TRUE, width=7, height=7, verbose=TRUE)

#############################################05.DL----
#去除连锁不平衡的SNP
exposureFile="exposure.pvalue.csv"     #输入文件

#读取输入文件
exposure_dat<-read_exposure_data(filename=exposureFile,
                                 sep = ",",
                                 snp_col = "SNP",
                                 beta_col = "beta.exposure",
                                 se_col = "se.exposure",
                                 effect_allele_col = "effect_allele.exposure",
                                 other_allele_col = "other_allele.exposure",
                                 eaf_col = "eaf.exposure",
                                 samplesize_col = "samplesize.exposure",
                                 clump = F)

#去除连锁不平衡的SNP
exposure_dat_clumped <- clump_data(exposure_dat, clump_kb=10000, clump_r2=0.001)
write.csv(exposure_dat_clumped, file="exposure.LD.csv", row.names=F)

#############################################06.F----
Ffilter=10        #F值过滤条件
inputFile="exposure.LD.csv"      #输入文件

#读取输入文件
dat<-read.csv(inputFile, header=T, sep=",", check.names=F)

#计算F检验值
N=dat[1,"samplesize.exposure"]     #获取样品的数目
dat=transform(dat,R2=2*((beta.exposure)^2)*eaf.exposure*(1-eaf.exposure))     #计算R2
dat=transform(dat,F=(N-2)*R2/(1-R2))      #计算F检验值

#根据F值>10进行过滤, 删除弱工具变量
outTab=dat[dat$F>Ffilter,]
write.csv(outTab, "exposure.F.csv", row.names=F)

