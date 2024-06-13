pFilter=1       #pvalue过滤条件
setwd("F:\\糖尿病_先兆子痫\\R10\\T2D\\02.TwoSampleMR")     #设置工作目录

############定义森林图函数############
bioForest=function(inputFile=null, forestFile=null, forestCol=null){
  #读取输入文件
  rt=read.csv(inputFile, header=T, sep=",", check.names=F)
  row.names(rt)=rt$method
  rt=rt[rt$pval<pFilter,]
  method <- rownames(rt)
  or <- sprintf("%.3f",rt$"or")
  orLow  <- sprintf("%.3f",rt$"or_lci95")
  orHigh <- sprintf("%.3f",rt$"or_uci95")
  OR <- paste0(or,"(",orLow,"-",orHigh,")")
  pVal <- ifelse(rt$pval<0.001, "<0.001", sprintf("%.3f", rt$pval))
  
  #输出图形
  pdf(file=forestFile, width=7, height=4.6)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3.5,2))
  
  #绘制森林图左边的信息
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,method,adj=0,cex=text.cex)
  text(1.9,n:1,pVal,adj=1,cex=text.cex);text(1.9,n+1,'pvalue',cex=1,font=2,adj=1)
  text(3.1,n:1,OR,adj=1,cex=text.cex);text(2.7,n+1,'OR',cex=1,font=2,adj=1)
  
  #绘制右边的森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(min(as.numeric(orLow)*0.975,as.numeric(orHigh)*0.975,0.9),max(as.numeric(orLow),as.numeric(orHigh))*1.025)
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="OR")
  arrows(as.numeric(orLow),n:1,as.numeric(orHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=3)
  abline(v=1, col="black", lty=2, lwd=2)
  boxcolor = ifelse(as.numeric(or)>1, forestCol, forestCol)
  points(as.numeric(or), n:1, pch = 15, col = boxcolor, cex=2)
  axis(1)
  dev.off()
}

#调用函数,绘制森林图
bioForest(inputFile="table.MRresult.csv", forestFile="forest.pdf", forestCol="red")



