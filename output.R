# output.R
library(qtl)
library(argparse)
# library(LinkageMapView)
# options(timeout = 240)

parser <- ArgumentParser()
parser$add_argument("-m", "--show_marker_names", action = "store_false",
                    help = "If specified, marker names are included.")

# Parse the arguments
args <- parser$parse_args()

# Create output directory if it doesn't exist
if (!dir.exists("output/output")) {
  dir.create("output/output", recursive = T)
}

### output----
if (file.exists("output/ripple/mapthis.RDS")){
  mapthis <- readRDS("output/ripple/mapthis.RDS")
} else {
  mapthis <- readRDS("output/order/mapthis.RDS")
}

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
plotMap(mapthis, show.marker.names = args$show_marker_names)
dev.off()

