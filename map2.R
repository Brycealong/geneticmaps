install.packages("qtl")
library(qtl)
https://storage.cloud.google.com/bucket-quickstart_mindful-linker-430503-t8/hq.vcf.gz
df <- read.csv("hq.vcf.gz.csv", header = FALSE)
df[2,1] <- ""
write.table(df, file = "hq.vcf.gz.rpl.csv", sep = ",",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
mapthis <- read.cross("csv", "", "hq.vcf.gz.rpl.csv",
                      estimate.map = FALSE, crosstype = "riself")
summary(mapthis)
plotMissing(mapthis)

par(mfrow=c(1,2), las=1)
plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual") 
plot(ntyped(mapthis, "mar"), ylab="No. typed individuals", 
     main="No. genotypes by marker")
mapthis <- subset(mapthis, ind=(ntyped(mapthis)>60000))

nt.bymar <- ntyped(mapthis, "mar")
todrop <- names(nt.bymar[nt.bymar < 250])
mapthis <- drop.markers(mapthis, todrop)

# Identify duplicate individuals
cg <- comparegeno(mapthis)
par(mfrow=c(1,1), las = 1)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes") 
rug(cg[lower.tri(cg)])

wh <- which(cg > 0.9, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
wh

# g <- pull.geno(mapthis)
# for(i in 1:nrow(wh)) {
#   tozero <- !is.na(g[wh[i,1],]) & !is.na(g[wh[i,2],]) & g[wh[i,1],] != g[wh[i,2],]
#   mapthis$geno[[1]]$data[wh[i,1],tozero] <- NA
#   }
mapthis <- subset(mapthis, ind=-wh[,2])

print(dup <- findDupMarkers(mapthis, exact.only=FALSE))
mapthis <- drop.markers(mapthis, unlist(dup))

# Look for markers with distorted segregation patterns
gt <- geno.table(mapthis)
gt[gt$P.value < 0.05/totmar(mapthis),]

todrop <- rownames(gt[gt$P.value < 1e-10,])
mapthis <- drop.markers(mapthis, todrop)

# Study individuals’ genotype frequencies
g <- pull.geno(mapthis)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:2)))
gfreq <- t(t(gfreq) / colSums(gfreq))
par(mfrow=c(1,2), las=1)
for(i in 1:2)
  {plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "BB")[i], ylim=c(0,1))}

# Study pairwise marker linkages; look for switched alleles
mapthis <- est.rf(mapthis)
checkAlleles(mapthis, threshold=5)

rf <- pull.rf(mapthis)
lod <- pull.rf(mapthis, what="lod")
pdf()
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")
