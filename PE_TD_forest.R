######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("grid")
#install.packages("readr")
#install.packages("forestploter")


#引用包
library(grid)
library(readr)
library(forestploter)

setwd("F:\\先兆子痫_糖尿病\\TD_forest")       #设置工作目录
files=dir()                           #获取目录下所有文件
files=grep("csv$", files, value=T)    #提取csv结尾的文件

#读取孟德尔随机化分析的结果
data=data.frame()
for(i in files){
  rt=read.csv(i, header=T, sep=",", check.names=F)
  data=rbind(data, rt)
}
lineVec=cumsum(c(1,table(data[,c('exposure','outcome')])))


#对数据进行整理
data$' ' <- paste(rep(" ", 10), collapse = " ")
data$'OR(95% CI)'=ifelse(is.na(data$or), "", sprintf("%.3f (%.3f to %.3f)", data$or, data$or_lci95, data$or_uci95))
data$pval = ifelse(data$pval<0.001, "<0.001", sprintf("%.3f", data$pval))
data$exposure = ifelse(is.na(data$exposure), "", data$exposure)
data$nsnp = ifelse(is.na(data$nsnp), "", data$nsnp)
data2=data[,c('exposure','outcome')]
data[duplicated(data2),]$exposure=""
data[duplicated(data2),]$outcome=""

#准备图形参数
tm <- forest_theme(base_size = 15,   #图形整体的大小
                   #可信区间的形状、线条类型、宽度、颜色、两端竖线高度
                   ci_pch = 16, ci_lty = 1, ci_lwd = 1.5, ci_col = "black", ci_Theight = 0.2, 
                   #参考线条的形状、宽度、颜色
                   refline_lty="dashed", refline_lwd=1, refline_col="grey20",
                   #x轴刻度字体的大小
                   xaxis_cex=0.8,
                   #脚注大小、颜色
                   footnote_cex = 0.6, footnote_col = "blue")

#绘制图形
plot <- forestploter::forest(data[, c("exposure","outcome","nsnp","method","pval"," ","OR(95% CI)")],
                             est = data$or,
                             lower = data$or_lci95,
                             upper = data$or_uci95,
                             ci_column = 6,     #可信区间所在的列
                             ref_line = 1,      #参考线条的位置
                             xlim = c(0, 2),    #X轴的范围
                             ticks_at = c(0,0.5,1,1.5,2),  #ticks_at参数设置坐标轴的刻度位置,大小对应xaxis_cex
                             theme = tm,        #图形的参数
)

#修改图形中可信区间的颜色
boxcolor = c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000","#7E6148")
boxcolor = boxcolor[as.numeric(as.factor(data$method))]
for(i in 1:nrow(data)){
  plot <- edit_plot(plot, col=6,row = i, which = "ci", gp = gpar(fill = boxcolor[i],fontsize=25)) # 改col，box的列
}

#设置pvalue的字体
pos_bold_pval = which(as.numeric(gsub('<',"",data$pval))<0.05)
if(length(pos_bold_pval)>0){
  for(i in pos_bold_pval){
    plot <- edit_plot(plot, col=5,row = i, which = "text", gp = gpar(fontface="bold"))  # 改col pvalue的列
  }
}

#在图形中增加线段
plot <- add_border(plot, part = "header", row =1,where = "top",gp = gpar(lwd =2))
plot <- add_border(plot, part = "header", row = lineVec, gp = gpar(lwd =1)) #上下边框线
#设置字体大小, 并且将文字居中
plot <- edit_plot(plot, col=1:ncol(data),row = 1:nrow(data), which = "text", gp = gpar(fontsize=12))
plot <- edit_plot(plot, col = 1:ncol(data), which = "text",hjust = unit(0.5, "npc"),part="header",
                  x = unit(0.5, "npc"))
plot <- edit_plot(plot, col = 1:ncol(data), which = "text",hjust = unit(0.5, "npc"),
                  x = unit(0.5, "npc"))

#输出图形
pdf("forest5.pdf", width=18, heigh=8)
print(plot)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

