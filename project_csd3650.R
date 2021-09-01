gds4879.raw <- read.table("GDS4879.clean", h = T)
gds4879 <- log2(gds4879.raw[,-c(1,2)])
##plotting data
pdf("boxplot.pdf")
boxplot(gds4879)
dev.off()

##get genes names
geneNames4879 <- as.character(gds4879.raw[,2])
##test to samples , save greater
myttest <- function(v, a, b){
  t.test(v[a], v[b], alternative = "greater")$p.value
}
##comparing alcoholic's and control's greater genes via myttest func
gds4879alcohol.pvalues <- apply(gds4879, 1, myttest, c(1:6,13:26),c(7:12,27:39))

##comparing male's and female's greater genes via myttest func
gds4879gender.pvalues <- apply(gds4879, 1, myttest, c(1:12),c(13:39))

##apply fdr on alcoholic and print their sum
gds4879alcohol.cor <- p.adjust(gds4879alcohol.pvalues, "fdr")
sum(gds4879alcohol.cor < 0.1)
##apply fdr on genger and print their sum
gds4879gender.cor <- p.adjust(gds4879gender.pvalues, "fdr")
sum(gds4879alcohol.cor < 0.1)

##apply bonferroni on alcoholic and print their sum
gds4879alcohol.cor.bon <- p.adjust(gds4879alcohol.pvalues, "bonferroni")
sum(gds4879alcohol.cor.bon < 0.1)

##apply bonferroni on gender and print their sum
gds4879gender.cor.bon <- p.adjust(gds4879gender.pvalues, "bonferroni")
sum(gds4879alcohol.cor.bon < 0.1)

geneNames4879new <- geneNames4879[gds4879alcohol.pvalues < 0.1]
geneNames4879new <- geneNames4879[gds4879gender.cor < 0.1]

##find names/id based on p-values

geneNames4879gprof.alc <- geneNames4879[gds4879alcohol.pvalues< 0.1]
geneNames4879gprof.gen <- geneNames4879[gds4879gender.pvalues< 0.1]

##extracting csv -> input to gprofiler and get plot
write.csv(c(geneNames4879gprof.alc,geneNames4879gprof.gen), file = "ex06.csv", row.names=FALSE, quote=FALSE)

mat <-gds4879


mypca<-prcomp(t(mat), center = TRUE, scale. = TRUE)
##create 12 red and 27 blue colored objects
colvec = c(rep("red",12),rep("blue",27))

##create cross= id 13 circle = id 19 for each different category
symbol = c(rep(3,6),rep(19,6),rep(3,14),rep(19,13))

pdf("pca.pdf")
##plotting 
plot(mypca$x[,1:39], col=colvec, pch=symbol)
dev.off()

pdf("boxplot.pdf")
boxplot(mat)
dev.off()
##author csd3650