# output.R
library(qtl)
library(argparse)
# library(LinkageMapView)
# options(timeout = 240)

parser <- ArgumentParser()
parser$add_argument("--dropone", action = "store_true",
                    help = "If specified, will output the results after drop one marker.")

# Parse the arguments
args <- parser$parse_args()

# Create output directory if it doesn't exist
if (!dir.exists("output/output")) {
  dir.create("output/output", recursive = T)
}

### output----
if (args$dropone){
  mapthis <- readRDS("output/order/mapthis_dropone.RDS")
} else if (file.exists("output/ripple/mapthis.RDS")){
  mapthis <- readRDS("output/ripple/mapthis.RDS")
} else {
  mapthis <- readRDS("output/order/mapthis.RDS")
}

map <- pull.map(mapthis)
maptbl <- map2table(map)
write.table(maptbl, file = "output/output/sum.csv", quote = F, 
            sep = ",")
# write.table(pull.rf(mapthis), file = "output/output/rf.csv", quote = F,
#             sep = ",")
# write.table(pull.rf(mapthis, what = "lod"), file = "output/output/lod.csv", quote = F,
#             sep = ",")
# maxpos <- floor(max(map2table(pull.map(mapthis))$pos))
# at.axis <- seq(0, maxpos)
# ## put labels on ruler at every 50 cM
# axlab <- vector()
# for (lab in 0:maxpos) {
#   if (!lab %% 50) {
#     axlab <- c(axlab, lab)
#   }
#   else {
#     axlab <- c(axlab, NA)
#   }
# }
# 
# outfile <- file.path("output", "output", "lmv-denmap.pdf")
# lmv.linkage.plot(mapthis, outfile = outfile, denmap = TRUE,
#                  cex.axis = 1, at.axis = at.axis, labels.axis = axlab,
#                  pdf.height = 20, pdf.width = 20, pdf.pointsize = 0.5)

png(file.path("output", "output", "map.png"), width = 1200, height = 1200, pointsize = 20)
plotMap(mapthis, show.marker.names = F)
dev.off()

for (i in chrnames(mapthis)){
  png(file.path("output", "output", paste0(i, ".png")), width = 1200, height = 1200, pointsize = 20)
  plotMap(mapthis, chr = i, show.marker.names = T, shift = F)
  dev.off()
}

cat("output complete.\n")