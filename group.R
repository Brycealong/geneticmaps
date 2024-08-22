# group.R
#install.packages("qtl")
library(qtl)
library(argparse)
options(timeout = 240)

# Create a parser object
parser <- ArgumentParser()

# Grouping options
parser$add_argument("--by", type = "character", choices = c("obs", "infer"), default = "obs",
                    help = "Specify 'obs' to use predefined groups (observation) or 'infer' to use max_rf and min_lod for forming linkage groups (inference)")
parser$add_argument("--max_rf", type = "double", default = 0.25,
                    help = "Maximum recombination fraction for forming linkage groups (ignored if --group is 'byRef'). (default: %(default)s)")
parser$add_argument("--min_lod", type = "double", default = 3,
                    help = "Minimum LOD score for forming linkage groups (ignored if --group is 'byRef'). (default: %(default)s)")

# Parse the arguments
args <- parser$parse_args()

# Create output directory if it doesn't exist
if (!dir.exists("output/group")) {
  dir.create("output/group", recursive = T)
}

mapthis <- readRDS(file = file.path("output", "preprocess", "mapthis.RDS"))

### Grouping ----
# Study pairwise marker linkages and estimate recombination fraction
mapthis <- est.rf(mapthis)

# rf <- pull.rf(mapthis)
# lod <- pull.rf(mapthis, what = "lod")
# 
# # Save recombination fraction vs. LOD plot
# png("output/group/rf-vs-LOD.png", width = 1200, height = 1200, pointsize = 20)
# plot(as.numeric(rf), as.numeric(lod), xlab = "Recombination fraction", ylab = "LOD score")
# dev.off()

# Form linkage groups
if (args$by == "obs") {
  # show map
  cat("Using the input groups...\n")
  print(summaryMap(mapthis))
} else if (args$by == "infer"){
  cat(paste0(
    "Two markers are placed in the same linkage group if the estimated recombination fraction between them is <= ",
    args$max_rf,
    " and the LOD score (for the test of the rec. frac. = 1/2) is >= ",
    args$min_lod,
    ". The transitive property (if A is linked to B and B is linked to C then A is linked to C) is used to close the groups.",
    "\n"
  ))
  # show map
  lg <- formLinkageGroups(mapthis, max.rf = args$max_rf, min.lod = args$min_lod)
  write.table(lg, file = "output/group/compare_sum.csv", quote = F, 
              sep = ",")
  # reorganize
  mapthis <- formLinkageGroups(mapthis, max.rf = args$max_rf, min.lod = args$min_lod, reorgMarkers = TRUE)
  print(summaryMap(mapthis))
} else {
  stop("Please specify one method to group the markers.")
}

map <- pull.map(mapthis)
maptbl <- map2table(map)
write.table(maptbl, file = "output/group/sum.csv", quote = F, 
            sep = ",")
saveRDS(mapthis, file = file.path("output", "group", "mapthis.RDS"))

png(file.path("output", "group", "map.png"), width = 1200, height = 1200, pointsize = 20)
plotMap(mapthis, show.marker.names = F)
dev.off()
