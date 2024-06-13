setwd("F:\\糖尿病_先兆子痫\\R10\\T1D\\02.TwoSampleMR")      #设置工作目录

#### 1.加载结局数据 ----
load('finndataR10_preeclampsia_outcome.rdata')

#### 2.加载暴露数据 ----
#引用包
library(VariantAnnotation)
library(gwasglue)
library(TwoSampleMR)

exposureName="Type1 diabetes"                         #图形中展示暴露数据的名称
outcomeName="preeclampsia"       #图形中展示结局数据的名称

exposureFile="exposure.F.csv"     #暴露数据输入文件


#读取暴露数据
exposure_dat<-read_exposure_data(filename=exposureFile,
                                 sep = ",",
                                 snp_col = "SNP",
                                 beta_col = "beta.exposure",
                                 se_col = "se.exposure",
                                 effect_allele_col = "effect_allele.exposure",
                                 other_allele_col = "other_allele.exposure",
                                 eaf_col = "eaf.exposure",
                                 clump = F)

#### 3.从结局数据中提取工具变量 ----
outcomeTab<-merge(exposure_dat, finndata, by.x="SNP", by.y="SNP")
write.csv(outcomeTab[,-(2:ncol(exposure_dat))], file="outcome.csv") 

#### 4.读取整理好的结局数据  ----
outcome_dat<-read_outcome_data(snps=exposure_dat$SNP,
                               filename="outcome.csv", sep = ",",
                               snp_col = "SNP",
                               beta_col = "beta.outcome",
                               se_col = "se.outcome",
                               effect_allele_col = "effect_allele.outcome",
                               other_allele_col = "other_allele.outcome",
                               pval_col = "pval.outcome",
                               eaf_col = "eaf.outcome")

#### 5.将暴露数据和结局数据合并 ----
exposure_dat$exposure=exposureName
outcome_dat$outcome=outcomeName
dat<-harmonise_data(exposure_dat=exposure_dat,
                    outcome_dat=outcome_dat)

#输出用于孟德尔随机化的工具变量
outTab=dat[dat$mr_keep=="TRUE",]
write.csv(outTab, file="table.SNP.csv", row.names=F)

#MR-PRESSO异常值检测(偏倚的SNP)
presso=run_mr_presso(dat)
write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`, file="table.MR-PRESSO.csv")

#孟德尔随机化分析
mrResult=mr(dat)
#选择孟德尔随机化的方法
#mr_method_list()$obj
#mrResult1=mr(dat, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode"))
#对结果进行OR值的计算
mrTab=generate_odds_ratios(mrResult)
write.csv(mrTab, file="table.MRresult.csv", row.names=F)

#异质性分析
heterTab=mr_heterogeneity(dat)
write.csv(heterTab, file="table.heterogeneity.csv", row.names=F)

#多效性检验
pleioTab=mr_pleiotropy_test(dat)
write.csv(pleioTab, file="table.pleiotropy.csv", row.names=F)

#绘制散点图
pdf(file="pic.scatter_plot.pdf", width=7.5, height=7)
mr_scatter_plot(mrResult, dat)
dev.off()

#森林图
res_single=mr_singlesnp(dat)      #得到每个工具变量对结局的影响
pdf(file="pic.forest.pdf", width=7, height=6.5)
mr_forest_plot(res_single)
dev.off()

#漏斗图
pdf(file="pic.funnel_plot.pdf", width=7, height=6.5)
mr_funnel_plot(singlesnp_results = res_single)
dev.off()

#留一法敏感性分析
pdf(file="pic.leaveoneout.pdf", width=7, height=6.5)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
dev.off()

