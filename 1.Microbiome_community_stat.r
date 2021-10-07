

## Hpyergeometric test
summary(genes)
Length Class Mode   
#[1,]  110   table numeric
#[2,] 1321   table numeric
#[3,] 1907   table numeric

sum(apply(data,1,function(x){x[1]>0 & x[2]>0} ))
# 81 Overlap genes


phyper(80, 110, 1907-110, 1321, lower.tail = FALSE, log.p = FALSE)

# 0.1803999

# phyper(80, 1321, 1907-1321, 110, lower.tail = FALSE, log.p = FALSE) 

#Same answer  0.1803999

# p-value: 0.18039

## Modified hypergeometric test

*Input of dataAll looks like (.csv):
gene
Rhizophagus
Suillus
Rhizopogen
Tuber
Hebeloma
Terfezia boudier
Cenococcum
Laccaria
Sebacina
Leptosphaeria
Cadophora
Fusarium
Pleurotus
Glonium
Mortierella
Ilyonectria
Atractiellales
Umbelopsis
*Input of dataV1 looks like (.csv):
gene
Atractiellales
Cadophora
Cadophora
Cadophora
Cadophora
Cadophora
Cadophora
Cadophora
Cenococcum
Cenococcum
*Input of dataV2 looks like (.csv):
gene
Atractiellales
Atractiellales
Atractiellales
Atractiellales
Atractiellales
Cadophora
Cadophora
Cadophora
Cadophora
Cadophora
Cadophora
Cadophora
Cadophora
Cadophor

dataV1 <- read.table(file.choose(), header = TRUE)
dataV2 <- read.table(file.choose(), header = TRUE)
dataAll<- read.table(file.choose(), header = TRUE)

genes<- list()
genes[[1]]<- table(dataV1$gene)
genes[[2]]<- table(dataV2$gene)
genes[[3]]<- table(dataAll$gene)

## V1 vs V2
# allgenes<- sort(as.character(unique(dataAll$gene)))
allgenes<- sort(unique( c(names(genes[[1]]), names(genes[[2]]) ) ))
data<-data.frame(matrix(0, nrow=length(allgenes), ncol=2))
rownames(data)<- allgenes

for(i in 1:2){
    data[ rownames(data) %in% names(genes[[i]]), i ] <- genes[[i]]
}
# sum(apply(data,1, function(x){all(is.na(x))} )) # sanity check == 0
# sum(apply(data,1, function(x){sum(x)==0} ))  # sanity check == 0

head(data)
# X1   X2
#Atractiellales  1    5
#Cadophora       7 1060
#Cenococcum     20   79
#Fusarium        1   76
#Glonium         2   63
#Hebeloma        1   15

# write.csv(data, file="DomainsV1V2.csv", row.names=F)


## Hpyergeometric test
summary(domains)
# [1,]  727   table numeric
# [2,]  459   table numeric
# [3,] 2465   table numeric
sum(apply(data,1,function(x){x[1]>0 & x[2]>0} ))
# 324 Overlap domains
phyper(323, 727, 2465-727, 459, lower.tail = FALSE, log.p = FALSE)
# phyper(323, 459, 2465-459, 727, lower.tail = FALSE, log.p = FALSE) #Same answer
# p-value: 1.831952e-93

##
chisq.test(x=data[,1],y=data[,2], simulate=T, B=1000000)
# p-value = 1e-05

cor.test(x=data[,1],y=data[,2])
# p-value < 2.2e-16










##hyper and fisher##

library(readxl)

data_raw <- data.frame(read_xlsx("genes assign to taxa.xlsx", skip=2))

index <- grep("taxa", names(data_raw))
data2 <- data_raw[,(index-1):(index+8)]
data2 <- data_raw[,(index+1):(index+8)]

nrow <- NROW(data2)
result <- matrix(nrow=nrow-1, ncol=8, dimnames=list(data_raw[-nrow, index], colnames(data2)))
fisher_list <- list()
groups <- list(1:2, 3:4, 5:6, 7:8)
for(g in groups){
  data_sub <- data2[,g]

  m <- data_sub[nrow, 1]
  n <- data_sub[nrow, 2]
  data_sub <- data_sub[-nrow, ]
  result[,g] <- t(apply(data_sub, 1, function(x){
      total <- sum(x)
      p1 <- phyper(x[1], m, n, total, lower.tail=F)
      p2 <- phyper(x[2], n, m, total, lower.tail=F)
      return(c(p1, p2))
  }))

  fisher_list[[g[1]]] <- apply(data_sub, 1, function(x){
    ctable <- matrix(c(x, m-x[1], n-x[2]), nrow=2)
    fisher.test(ctable)$p.value
  })

}

cname <- c('A. FunGene_P.tri (9988)_P.tri (P.tri x P.tri)_gene counts_vs_P. tri (Ptri x P. taeda)_gene_counts',
'B. FunGene_P.tri (4771)_P. taeda (P.tri x P.taeda)_gene counts_vs_P. taeda (P.taedax P.taeda)_gene counts',
'C. FunGene_P.taeda (1527)_P.tri (P.tri x P.tri)_gene counts_vs_P. tri (Ptri x P. taeda)_gene_counts',
'D. FunGene_P.taeda (2312)_P. taeda (P.tri x P.taeda)_gene counts_vs_P. taeda (P.taedax P.taeda)_gene counts')

result_fisher <- do.call("cbind", fisher_list)
dimnames(result_fisher) <- list(rownames(result), cname)

write.csv(result_fisher, "genes_to_taxa_fisher.csv")
write.csv(result, "genes_to_taxa_hyper.csv")








##Permanova

# taxa_name$group <- gsub("Endophye", "endophyte", taxa_name$group)
library(readxl)
library(vegan)
library(RColorBrewer)
library(gplots)


# librdata_header <- read.csv("Neighbor_Dataset_2C.csv", nrow=1)
# data <- read.csv("Neighbor_Dataset_2C.csv", skip=1)
data_header <- data.frame(read_excel("Neighbor_Dataset_2A&B.xlsx", n_max=2, col_names=F, skip=1, sheet="2C"))
data <- data.frame(read_excel("Neighbor_Dataset_2A&B.xlsx", col_names=F, skip=3, sheet="2C"))
head(data)

if(any(dim(data)!=c(51,18))){
  stop("Check data import!!!!!!!!!!!!!!!!")
}

taxa_name <- setNames(data[,2:3], data_header[1, 2:3])
taxa_name$group <- taxa_name[[data_header[1,2]]]
index <- names(which(table(taxa_name$group) <= 2))
taxa_name$group[taxa_name$group %in% index] <- "Other"
taxa_name$group <- factor(taxa_name$group)
# taxa_name$group <- factor(gsub("[_ ].+", "", taxa_name$group))
# taxa_name$group5 <- factor(gsub("/.+", "", taxa_name$group))
# taxa_name2 <- apply(taxa_name, 1, function(x){paste(x, collapse="_")})

h2 <- data_header[2,]
# h2 <- gsub("P. taeda", "Taead", h2)
h2 <- gsub("Ptri", "P\\. tri", h2)
h2 <- na.omit(h2)

data2 <- data[,-(1:3)]
colnames(data2) <- h2
# rownames(data) <- group_name2
data[is.na(data)]<- 0

dataExp <- t(apply(data2, 2, function(x){x/sum(x, na.rm=T)}))
apply(dataExp,1,sum, na.rm=T)
site <- data.frame(ID=gsub("Sum of(.+)_\\%read", "\\1", data_header[1,-(1:3)]), label=factor(rownames(dataExp)))


cPal8 <- brewer.pal(8, "Dark2")


# site <- rep(c("A","B"), each=4)

# dataTrt <- data.frame(site=site, rep=paste0(site, 1:4))
# cPal8 <- brewer.pal(8, "Dark2")

# for(sheet in sheets){
#   pdf(sprintf("overallResult_%s.pdf", sheet), width=8, height=8)
#   data <- data.frame(read_xlsx("V10_Fig3_statistics.xlsx", skip=1, sheet=sheet))
#   names(data) <- gsub("ecological.function", "group", names(data))
#   summary(data)
#   data[["group"]] <- gsub("Sap .+", "Sap", data[["group"]])
#   data[["group"]] <- gsub("pathogen", "Pathogen", data[["group"]])
#   data[["group"]] <- gsub("Ecto-,.+", "Ecto", data[["group"]])
#   data[["group"]] <- factor(data[["group"]])
#
#   dataExp <- t(data[, grep("[A|B][1-4]", names(data), value=T)])
#   colnames(dataExp) <- data$ref


all_labels <- levels(site$label)
resultAll <- list()
for(cc in combn(4,2, simplify=FALSE)){
  index <- site$label %in% all_labels[cc]
  dataSub <- dataExp[index, ]
  site_cc <- site$label[index]
  dd <- rep(NA, length=NCOL(dataSub))
  for(i in 1:NCOL(dataSub)){
     try({
       abDist <- vegdist(dataSub[,i], method="bray")
       # zero-adjusted Bray–Curtis
       abDist[is.na(abDist)] <- 0
       result <- adonis(abDist ~ site_cc ) #, data=site)
       # result <- adonis(dataSub[,i] ~ site_cc, method="bray" ) #, data=site)
       dd[i] <-   result$aov.tab[["Pr(>F)"]][1]
     }, silent=T)
  }
  resultAll[[paste(all_labels[cc], collapse=" : ")]] <- dd

}
result_output <- cbind(taxa_name, do.call("cbind", resultAll))
write.csv(result_output, file="2C_permanova_each_ref.csv")


pdf("2C_permanova_each_sites.pdf")
all_labels <- levels(site$label)
for(cc in combn(4,2, simplify=FALSE)){
  index <- site$label %in% all_labels[cc]
  dataSub <- dataExp[index, ]
  abDist <- vegdist(dataSub, method="bray")
  result <- adonis(abDist ~ site$label[index] )#, data=site)
  textplot(c(all_labels[cc], capture.output(result)), cex=0.7)
}
dev.off()


abDist <- vegdist(dataExp, method="bray")
result <- adonis(abDist ~ site$label )#, data=site)
pdf("2C_permanova_all_sites.pdf")
textplot(capture.output(result))
dev.off()

mm <- metaMDS((dataExp), k=2, distance="euclidean", trymax=100, autotransform=F, noshare=F, wascores=T)

stress <- vector()
for(k in 1:10){
  mm <- metaMDS((dataExp), k=k, distance="bray", trymax=100, autotransform=F, noshare=F, wascores=T)
  stress[k] <- mm$stress
}
plot(stress)

k <- 2

png("2C_Site.png", res=300, width=2000, height=2000)
tiff("2C_Site.tiff", res=300, width=2000, height=2000)

cPalSet1 <- brewer.pal(8, "Set1")
customColour <- as.numeric(site$label)
customColour[customColour==1] <- cPalSet1[4]
customColour[customColour==2] <- cPalSet1[3]
customColour[customColour==3] <- cPalSet1[2]
customColour[customColour==4] <- cPalSet1[1]

mm <- metaMDS((dataExp), k=k, distance="bray", trymax=100, autotransform=F, noshare=F, wascores=T)
choice <- c(1, 2)
plot(mm, type="n", display = "sites", cex=2, choice=choice, xlim=c(-1.2, 0.9))
  # main=sprintf("%s   p-value:%.3f", "", result$aov.tab[["Pr(>F)"]][1]), )
# text(mm, display="sites", labels = site$ID, col=cPal8[as.numeric(site$label)], choice=choice )
# points(mm, col=cPal8[as.numeric(site$label)], cex=1.2)
points(mm, pch=as.numeric(site$label)-1, cex=1.2, col=customColour)
with(site, ordiellipse(mm, label, scaling = "symmetric", kind = "ehull", col = cPalSet1[4:1], lwd=2, label=T))
legend("bottomleft", legend=levels(site$label), pch=1:nlevels(site$label)-1, bty="n", col=cPalSet1[4:1])

dev.off()

# png("2C_Site_BW.png", res=300, width=2000, height=2000)
tiff("2C_Site_BW.tiff", res=300, width=2000, height=2000)

plot(mm, type="n", display = "sites", cex=2, choice=choice, xlim=c(-1.2, 0.9))
points(mm, pch=as.numeric(site$label)-1, cex=1.2)
with(site, ordiellipse(mm, label, scaling = "symmetric", kind = "ehull", col = cPalSet1[4:1], lwd=2, label=T))
legend("bottomleft", legend=levels(site$label), pch=1:nlevels(site$label)-1, bty="n")
dev.off()

  # plot(mm, type="p", display = "sites", cex=2, choice=c(2,3),
  #   main=sprintf("%s   p-value:%.3f", "", result$aov.tab[["Pr(>F)"]][1]))
  # text(mm, display="sites", labels = site$ID, col=cPal8[as.numeric(site$label)], choice=c(2,3) )
  # with(site, ordiellipse(mm, label, scaling = "symmetric", kind = "ehull", col = cPal8[1:nlevels(label)], lwd=3, label=T, choice=c(2,3)))




png("2C_Species.png", res=300, width=2000, height=2000)
tiff("2C_Species.tiff", res=300, width=2000, height=2000)
# dev.off()
# png("2C_Species.png", width=600, height=600)
new_legend <- levels(taxa_name$group)
# new_legend[3] <- "Endophyte or/and EM"
# new_legend[4] <- "Pathogen or/and endophyte"
plot(mm, main="", type="n")
# text(mm, display="sites", labels = site$ID, col=cPal8[8])
points(mm, pch=as.numeric(site$label)-1, col=1, cex=1.2)
text(mm, display="species", labels = taxa_name[["Row Labels"]], col=cPal8[taxa_name$group], cex=0.7)
legend("bottomleft", legend=new_legend, fill=cPal8[1:nlevels(taxa_name$group)], bty="n")
dev.off()

write.csv(taxa_name, file="2C_taxa.csv")

  # , col=cPal8[as.numeric(group)], cex=0.7))
  # with(dataExp, text(mm, display="species", labels = as.character(ref), col=cPal8[as.numeric(group)], cex=0.7))
  # #
  # # dev.off()
  #
}

# TODO: Analysis 2C
# TODO: 2A. 44 by 7 (all + pairwise) PERANOVA p-value table.
# TODO: 2A. by groups (EM, AM..etc) by 7

