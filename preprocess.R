# preprocess.R
#install.packages("qtl")
library(qtl)
library(argparse)
options(timeout = 240)
set.seed(61777369)

# crosstype:
# "bc" for backcross
# "f2" for intercross
# "riself" for 2-way RIL by selfing
# "risib" for 2-way RIL by sib-mating

# Create a parser object
parser <- ArgumentParser()

# cross type
parser$add_argument("--crosstype", type = "character", required = T,
                    choices = c("bc", "f2", "riself", "risib"), help = "Cross type for the analysis, choose from [bc, f2, riself, risib]")

# Filtering options
parser$add_argument("--filterMissingMarkers", action = "store_true",
                    help = "Filter out markers with missing data")
parser$add_argument("--filterMissingMarkersThres", type = "double", default = 0.05, metavar = "<FLOAT>",
                    help = "Filter out markers with missing data. For example, if <Thres>=0.05, then any marker with more than 5%% missing observations will be removed. (default: %(default)s)" )
parser$add_argument("--filterDupMarkers", action = "store_true",
                    help = "Filter out markers that have the same genotypes.")
parser$add_argument("--filterCloseMarkers", action = "store_true",
                    help = "Filter out markers that are closer than a specified threshold")
parser$add_argument("--filterCloseMarkersThres", type = "double", default = 0.5, metavar = "<FLOAT>",
                    help = "Identify the largest subset of markers for which no two adjacent markers are separated by less than the specified distance (in cM). (default: %(default)s)")

parser$add_argument("--filterSegregDistMarkers", action = "store_true",
                    help = "Remove markers that do not pass the segregation test")

parser$add_argument("--filterMatchingIndividuals", action = "store_true",
                    help = "Omit individuals with a high proportion of matching genotypes")
parser$add_argument("--filterMatchingIndividualsThres", type = "double", default = 0.9, metavar = "<FLOAT>",
                    help = "Threshold for removing one of a pair of individuals with more than <Thres> proportion of matching genotypes across all markers. (default: %(default)s)")

# Parse the arguments
args <- parser$parse_args()

# Create output directory if it doesn't exist
if (!dir.exists("output/preprocess")) {
  dir.create("output/preprocess", recursive = T)
}

mapthis <- read.cross("csv", "output/import", "rqtl.csv",
                      estimate.map = FALSE, crosstype = args$crosstype)

### Filtering ----
# 1. Filter markers with a lot of missing data
if (args$filterMissingMarkers) {
  cat(paste0(
    "Filtering markers that have more than ",
    sprintf("%1.0f%%", args$filterMissingMarkersThres * 100),
    " missing observations...",
    "\n"
  ))
  org_count <- totmar(mapthis)
  nt.bymar <- ntyped(mapthis, "mar")
  todrop <- names(nt.bymar[nt.bymar < (1 - args$filterMissingMarkersThres) * nind(mapthis)])
  mapthis <- drop.markers(mapthis, todrop)
  count <- totmar(mapthis)
  cat(paste0(
    "Original marker number: ",
    org_count,
    ", Filtered: ",
    org_count - count,
    ", Remaining: ",
    count,
    "\n"
  ))
}

# 2. Identify and drop duplicate markers
if (args$filterDupMarkers){
  cat(paste0(
    "Filtering markers that are duplicates...",
    "\n"
  ))
  org_count <- totmar(mapthis)
  dup <- findDupMarkers(mapthis, exact.only = FALSE)
  mapthis <- drop.markers(mapthis, unlist(dup))
  count <- totmar(mapthis)
  cat(paste0(
    "Original marker number: ",
    org_count,
    ", Filtered: ",
    org_count - count,
    ", Remaining: ",
    count,
    "\n"
  ))
}


# 3. Drop markers that are too close to each other
if (args$filterCloseMarkers) {
  cat(paste0(
    "Identifying the largest subset of markers where no two adjacent markers are separated by less than ",
    args$filterCloseMarkersThres,
    " cM...",
    "\n"
  ))
  org_count <- totmar(mapthis)
  for (chr in chrnames(mapthis)) {
    mark.cur <- markernames(mapthis, chr)
    tokeep <- pickMarkerSubset(pull.map(mapthis)[[chr]], args$filterCloseMarkersThres)
    todrop <- setdiff(mark.cur, tokeep)
    mapthis <- drop.markers(mapthis, todrop)
  }
  count <- totmar(mapthis)
  cat(paste0(
    "Original marker number: ",
    org_count,
    ", Filtered: ",
    org_count - count,
    ", Remaining: ",
    count,
    "\n"
  ))
}

# 4. Remove markers that do not pass the segregation test
if (args$filterSegregDistMarkers) {
  cat(paste0(
    "Filtering markers that have segregation distortion...",
    "\n"
  ))
  org_count <- totmar(mapthis)
  gt <- geno.table(mapthis)
  todrop <- rownames(gt[p.adjust(gt$P.value, method = "bonferroni") < 0.05,])
  mapthis <- drop.markers(mapthis, todrop)
  count <- totmar(mapthis)
  cat(paste0(
    "Original marker number: ",
    org_count,
    ", Filtered: ",
    org_count - count,
    ", Remaining: ",
    count,
    "\n"
  ))
}

## Filter individuals
#Identify and omit individuals with high matching genotypes
if (args$filterMatchingIndividuals) {
  cat(paste0(
    "Filtering individuals that have more than ",
    sprintf("%1.0f%%", args$filterMatchingIndividualsThres * 100),
    " matching genotypes across all markers...",
    "\n"
  ))
  org_count <- nind(mapthis)
  cg <- comparegeno(mapthis)
  wh <- which(cg > args$filterMatchingIndividualsThres, arr = TRUE)
  wh <- wh[wh[,1] < wh[,2],]
  
  for (i in 1:nrow(wh)) {
    for (chr in chrnames(mapthis)) {
      g <- pull.geno(mapthis, chr)
      tozero <- !is.na(g[wh[i,1],]) & !is.na(g[wh[i,2],]) & g[wh[i,1],] != g[wh[i,2],]
      mapthis$geno[[chr]]$data[wh[i,1], tozero] <- NA
    }
  }
  mapthis <- subset(mapthis, ind = -wh[,2])
  count <- nind(mapthis)
  cat(paste0(
    "Original individual number: ",
    org_count,
    ", Filtered: ",
    org_count - count,
    ", Remaining: ",
    count,
    "\n"
  ))
}

print(summaryMap(mapthis))
saveRDS(mapthis, file = file.path("output", "preprocess", "mapthis.RDS"))

