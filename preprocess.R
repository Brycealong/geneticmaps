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

mapthis <- readRDS("output/import/mapthis.RDS")

# Initialize an empty data frame to store the differences in marker counts
marker_diff <- data.frame(Chromosome = chrnames(mapthis))
rownames(marker_diff) <- marker_diff$Chromosome
marker_diff$Chromosome <- NULL  # Remove the Chromosome column as it's redundant

# print(summaryMap(mapthis))
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
  marker_before <- nmar(mapthis)
  ###
  nt.bymar <- ntyped(mapthis, "mar")
  todrop <- names(nt.bymar[nt.bymar < (1 - args$filterMissingMarkersThres) * nind(mapthis)])
  mapthis <- drop.markers(mapthis, todrop)
  ###
  marker_after <- nmar(mapthis)
  marker_diff$MissingMarkers <- marker_before - marker_after
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
  marker_before <- nmar(mapthis)
  ###
  dup <- findDupMarkers(mapthis, exact.only = FALSE)
  mapthis <- drop.markers(mapthis, unlist(dup))
  ###
  marker_after <- nmar(mapthis)
  marker_diff$DupMarkers <- marker_before - marker_after
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
  marker_before <- nmar(mapthis)
  ###
  for (chr in chrnames(mapthis)) {
    mark.cur <- markernames(mapthis, chr)
    tokeep <- pickMarkerSubset(pull.map(mapthis)[[chr]], args$filterCloseMarkersThres)
    todrop <- setdiff(mark.cur, tokeep)
    mapthis <- drop.markers(mapthis, todrop)
  }
  ###
  marker_after <- nmar(mapthis)
  marker_diff$CloseMarkers <- marker_before - marker_after
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
  marker_before <- nmar(mapthis)
  ###
  gt <- geno.table(mapthis)
  p <- gt$P.value
  todrop <- rownames(gt[p.adjust(p, method = "bonferroni") < 0.05 & !is.na(p),])
  mapthis <- drop.markers(mapthis, todrop)
  ###
  marker_after <- nmar(mapthis)
  marker_diff$SegregDistMarkers <- marker_before - marker_after
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

# left markers
marker_diff$MarkersLeft <- nmar(mapthis)

# Transpose the data frame to get a suitable format for plotting
marker_diff_matrix <- t(marker_diff)

# Define colors for each filtering step
cs <- c("MissingMarkers" = "#993404",
        "DupMarkers" = "#FB6A4A", 
        "CloseMarkers" = "#FED976", 
        "SegregDistMarkers" = "#FFFFCC",
        "MarkersLeft" = "white")

png(file.path("output", "preprocess", "marker_diff.png"), width = 1200, height = 1200, pointsize = 20)
# par(mar = c(5, 4, 4, 6) + 0.1)
barplot(marker_diff_matrix,
        col = cs[rownames(marker_diff_matrix)],
        legend.text = rownames(marker_diff_matrix), 
        args.legend = list(x = "topleft", cex = 0.8),
        main = "Stacked Barplot for Filtered Markers",
        xlab = "Chromosome", 
        ylab = "Number of Markers"
)
dev.off()

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

png(file.path("output", "preprocess", "map.png"), width = 1200, height = 1200, pointsize = 20)
plotMap(mapthis, show.marker.names = F)
dev.off()

cat("preprocess complete.\n")