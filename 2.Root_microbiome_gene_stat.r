
##heatmap

SAshort <- read.table(file.choose(), header = TRUE)
SAnew = SAshort[,2:7]
SAnew.log=log2(SAnew+1)
SAnew.log.hr=hclust(as.dist(1-cor(t(SAnew.log))))
SAnew.log.range<-quantile(as.matrix(SAnew.log) ,prob = seq(0,1,0.01))
SAnew.log.range.breaks <- seq(SAnew.log.range["5%"],SAnew.log.range["95%"],0.1)
SAnew.log.range.palette <- colorRampPalette(c("blue","black","yellow"))(length(SAnew.log.range.breaks)-1)
heatmap(as.matrix(SAnew.log),Colv=NA,Rowv=NA,col=SAnew.log.range.palette,labRow="")


##plot
SA <- read.table(file.choose(), header = TRUE)
SA.log=log2(SA+1)
library("car")
SA.pca <- prcomp(SA.log[,1:7], scale.=TRUE)
round(SA.pca$rotation[,1:2] ,2)
library("ggplot2")
scores <- as.data.frame(SA.pca$x, col = SA.log[,1:7])
wiltest <- function(x){
res1 = wilcox.test(as.numeric(x[1:4]),as.numeric(x[5:7]), alternative="greater")
res2 = wilcox.test(as.numeric(x[1:4]),as.numeric(x[5:7]), alternative="less")
if(res1$p.value <= 0.05){
return ("blue")
}
else if(res2$p.value <= 0.05){
return ("green")
}
else{
return ("black")
}
}
wiltestres = apply(SA.log,1,wiltest)
plot(scores[,c(1,2)],col=wiltestres,pch=20,cex=0.5)



##volcano

SA <- read.table(file.choose(), header = TRUE)
SA.log = log2(SA+1)
mean1 = apply(SA.log[1:3],1,mean)
mean2 = apply(SA.log[4:6],1,mean)
M = mean2-mean1
wiltest <- function(x){
    p = t.test(as.numeric(x[1:3]),as.numeric(x[4:6]))
  p$p.value
 }
p = apply(SA.log,1,wiltest)
fdr = p.adjust(p, method = "BH")
color = rep("black",length(fdr))
color[(fdr <=0.05 & M >=1)] = "purple" 
color[(fdr <=0.05 & M <= -1)] = "green"
plot(M,-log(p), pch=19,cex=0.5,col=color, ylab = "-log10(p)")

write.table(fdr,"~/fdr",sep="\t")
