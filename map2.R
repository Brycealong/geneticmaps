#install.packages("qtl")
library(qtl)
options(timeout = 240)
# https://storage.cloud.google.com/bucket-quickstart_mindful-linker-430503-t8/hq.vcf.gz
# df <- read.csv("hq.vcf.gz.csv", header = FALSE)
# df[2,1] <- ""
# write.table(df, file = "hq.vcf.gz.rpl.csv", sep = ",",
#             quote = FALSE, row.names = FALSE, col.names = FALSE)
mapthis <- read.cross("csv", "https://storage.googleapis.com/bucket-quickstart_mindful-linker-430503-t8/datasets", "hq.filt.vcf.gz.rpl.csv",
                      estimate.map = FALSE, crosstype = "riself")

#prepearing a vcf for setting the alleles 
tx  <- readLines("hq.filt.vcf.gz")
tx2  <- gsub(pattern = "#CHROM", replace = "CHROM", x = tx)
writeLines(tx2, con="cross.vcf.wohash.vcf")
rm(tx,tx2)

vcf <- read.table("cross.vcf.wohash.vcf", header = T, stringsAsFactors = F)
vcf[,c(10:ncol(vcf))] <- sapply(vcf[,c(10:ncol(vcf))], function(x)(
  gsub("|", "/", substr(x, 1, 3), fixed = TRUE)
))

#setting the maternal genotype which will be set to As
parent.A = "NM9"
parent.B = "Y158"

#replacing vcf encoded alleles to As, Bs, and Hs
vcf[,10:ncol(vcf)] <-  t(apply(vcf[,10:ncol(vcf)], 1, function(x)(
  ifelse(x == x[names(x) == parent.A], "A",
         ifelse(x == "0/1" | x  == "1/0", "H", 
                ifelse(x == "./.", "-", "B")))
)))

#creating qtl formatted file 
Rqtl <- vcf[,c(3,1,2,10:ncol(vcf))]
Rqtl$ID <- ifelse(Rqtl$ID == ".", paste0("S", Rqtl$CHROM, "_", Rqtl$POS), Rqtl$ID)
# Rqtl <- subset(Rqtl, select = -POS)
Rqtl <- t(Rqtl)
colnames(Rqtl) <- Rqtl[1,]
Rqtl <- Rqtl[-1,]

# #approximation of physical distance to genetic distance (1cM ~ 1Mbp)
# Rqtl[1,] <- as.numeric(sub(".+_", "", colnames(Rqtl))) / 1000000
# 
# #adding chromosome numbers
# Rqtl <- rbind(sub("_.+", "", sub("S", "", colnames(Rqtl))), Rqtl)

# this is optional step needed if the map is going to be reconstructed de novo without prioir information on markers
#Rqtl[1,] <- 1

#formatting to qtl object
Rqtl <-  as.data.frame(Rqtl)
Rqtl <- data.frame(id = rownames(Rqtl), Rqtl, row.names = NULL)

Rqtl[c(1,2),1] <- ""

write.table(Rqtl, "cross.qtl.csv", quote = F, sep = ",", row.names = F, na = "-")

### ----
summary(mapthis)
if (!dir.exists("images")){
  dir.create("images")
}
png("images/plotMissing.png", width = 1200, height = 1200, pointsize = 20)
plotMissing(mapthis)
dev.off()

png("images/geno-markers-individuals.png", width = 1200, height = 1200, pointsize = 20)
par(mfrow=c(1,2), las=1)
plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual") 
plot(ntyped(mapthis, "mar"), ylab="No. typed individuals", 
     main="No. genotypes by marker")
dev.off()

mapthis <- subset(mapthis, ind=(ntyped(mapthis)>60000))
summary(mapthis)
threshold <- 0.1
nt.bymar <- ntyped(mapthis, "mar")
todrop <- names(nt.bymar[nt.bymar < (1-threshold)*totmar(mapthis)])
mapthis <- drop.markers(mapthis, todrop)

# Identify duplicate individuals
cg <- comparegeno(mapthis)
png("images/dup-individuals.png", width = 1200, height = 1200, pointsize = 20)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes", main = "") 
rug(cg[lower.tri(cg)])
dev.off()

wh <- which(cg > 0.9, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
# wh

# g <- pull.geno(mapthis)
# for(i in 1:nrow(wh)) {
#   tozero <- !is.na(g[wh[i,1],]) & !is.na(g[wh[i,2],]) & g[wh[i,1],] != g[wh[i,2],]
#   mapthis$geno[[1]]$data[wh[i,1],tozero] <- NA
#   }
mapthis <- subset(mapthis, ind=-wh[,2])

# Identify duplicate markers
print(dup <- findDupMarkers(mapthis, exact.only=FALSE))
mapthis <- drop.markers(mapthis, unlist(dup))

# Look for markers with distorted segregation patterns
gt <- geno.table(mapthis)
# gt[gt$P.value < 0.05/totmar(mapthis),]
# use bonferroni p.adjust
todrop <- rownames(gt[p.adjust(gt$P.value, method = "bonferroni") < 0.05,])
mapthis <- drop.markers(mapthis, todrop)

# Study individuals’ genotype frequencies
g <- pull.geno(mapthis)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:2)))
gfreq <- t(t(gfreq) / colSums(gfreq))
png("images/geno-freq.png", width = 1200, height = 1200, pointsize = 20)
par(mfrow=c(1,2), las=1)
for(i in 1:2){
  plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "BB")[i], ylim=c(0,1))
  }
dev.off()

# Study pairwise marker linkages; look for switched alleles
mapthis <- est.rf(mapthis)
checkAlleles(mapthis, threshold=5)

rf <- pull.rf(mapthis)
lod <- pull.rf(mapthis, what="lod")
png("images/rf-vs-LOD.png", width = 1200, height = 1200, pointsize = 20)
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")
dev.off()


lg <- formLinkageGroups(mapthis, max.rf=0.25, min.lod=6)
table(lg[,2])

mapthis <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=8, reorgMarkers=TRUE)
png("images/plotRF1.png", width = 1200, height = 1200, pointsize = 20)
plotRF(mapthis, alternate.chrid=TRUE)
dev.off()

toswitch <- markernames(mapthis, chr=c(22:27)) 
mapthis <- switchAlleles(mapthis, toswitch)

mapthis <- est.rf(mapthis)
png("images/plotRF1.png", width = 1200, height = 1200, pointsize = 20)
plotRF(mapthis, alternate.chrid=TRUE)
dev.off()


