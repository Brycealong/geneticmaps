## import.R

#install.packages("qtl")
library(qtl)
library(argparse)
# library(LinkageMapView)
options(timeout = 240)

# Create a parser object
parser <- ArgumentParser()

# Define the arguments
parser$add_argument("-i", "--input", type = "character", required = TRUE,
                    help = "Path to the input VCF file. gzipped file is supported.")
parser$add_argument("-c", "--chromInclude", type="character", nargs = "+",
                    help="A list of chromosomes to be INcluded in the analysis, separated by space",
                    metavar="<CHR>")
parser$add_argument("-ce", "--chromExclude", type="character", nargs = "+",
                    help="A list of chromosomes to be EXcluded in the analysis, separated by space",
                    metavar="<CHR>")
parser$add_argument("-pA", "--parentA", type = "character", required = TRUE,
                    help = "Name of parent A in the vcf file.")
parser$add_argument("-pB", "--parentB", type = "character", required = TRUE,
                    help = "Name of parent B in the vcf file.")
# cross type
parser$add_argument("--crosstype", type = "character", required = T,
                    choices = c("bc", "f2", "riself", "risib"), help = "Cross type for the analysis, choose from [bc, f2, riself, risib]")

# parser$add_argument("-o", "--output", type = "character", required = TRUE,
#                     help = "Output direcotory.")

# Parse the arguments
args <- parser$parse_args()

# Create output directory if it doesn't exist
if (!dir.exists("output/import")) {
  dir.create("output/import", recursive = T)
}

# echo copy
cat("Copying original file into output/import...\n")
if (file.copy(args$input, to = "output/import/org.vcf.gz", overwrite = T)){
  cat("Copy complete.\n")
}

# vcf2rqtl csv
tx  <- readLines("output/import/org.vcf.gz")
tx2  <- gsub(pattern = "#CHROM", replace = "CHROM", x = tx)
writeLines(tx2, con="output/import/org.temp.vcf")
rm(tx,tx2)

vcf <- read.table("output/import/org.temp.vcf", header = T, stringsAsFactors = F)
if(!is.null(args$chromInclude)){
  vcf <- vcf[vcf$CHROM %in% args$chromInclude,]
  if(!is.null(args$chromExclude)){
    warnings("proving chromInclude and chromExclude at the same time, using only chromInclude...")
  }
} else if(!is.null(args$chromExclude)){
  vcf <- vcf[!vcf$CHROM %in% args$chromExclude,]
} else{
  vcf <- vcf
}

vcf[,c(10:ncol(vcf))] <- sapply(vcf[,c(10:ncol(vcf))], function(x)(
  gsub("|", "/", substr(x, 1, 3), fixed = TRUE)
))

#setting the maternal genotype which will be set to As
parent.A = args$parentA
parent.B = args$parentB

#replacing vcf encoded alleles to As, Bs, and Hs
vcf[,10:ncol(vcf)] <-  t(apply(vcf[,10:ncol(vcf)], 1, function(x)(
  ifelse(x == x[names(x) == parent.A], "A",
         ifelse(x == "0/1" | x  == "1/0", "H", 
                ifelse(x == "./.", "-", "B")))
)))

#creating qtl formatted file 
Rqtl <- vcf[,c(3,1,2,10:ncol(vcf))]
Rqtl <- Rqtl[, !(colnames(Rqtl) %in% c(parent.A, parent.B))]
Rqtl$ID <- ifelse(Rqtl$ID == ".", paste0("S", Rqtl$CHROM, "_", Rqtl$POS), Rqtl$ID)
# Rqtl <- subset(Rqtl, select = -POS)
#approximation of physical distance to genetic distance (1cM ~ 1Mbp)
Rqtl$POS <- as.numeric(Rqtl$POS) / 1e06
Rqtl <- t(Rqtl)
colnames(Rqtl) <- Rqtl[1,]
Rqtl <- Rqtl[-1,]

# this is optional step needed if the map is going to be reconstructed de novo without prioir information on markers
#Rqtl[1,] <- 1

#formatting to qtl object
Rqtl <- as.data.frame(Rqtl)
Rqtl <- data.frame(id = rownames(Rqtl), Rqtl, row.names = NULL)

Rqtl[c(1,2),1] <- ""

write.table(Rqtl, "output/import/rqtl.csv", quote = F, sep = ",", row.names = F, na = "-")

mapthis <- read.cross("csv", "output/import", "rqtl.csv",
                      estimate.map = FALSE, crosstype = args$crosstype)

print(summaryMap(mapthis))
saveRDS(mapthis, file = file.path("output", "import", "mapthis.RDS"))

png(file.path("output", "import", "org_map.png"), width = 1200, height = 1200, pointsize = 20)
plotMap(mapthis)
dev.off()
