# dropone.R
#install.packages("qtl")
suppressPackageStartupMessages(
  {
    library(qtl)
    library(argparse)
    library(parallel)
    library(snow)
  }
)
options(timeout = 240)
set.seed(61777369)

# Create a parser object
parser <- ArgumentParser()

# Define the arguments
parser$add_argument("-e", "--error_prob", type = "double", default = 0.0001, 
                    help = "Assumed genotyping error rate used in the final estimated map. (default: %(default)s)")
parser$add_argument("-f", "--map_function", type = "character", choices = c("haldane","kosambi","c-f","morgan"), default = "haldane", 
                    help = "map function to use. default is 'haldane'.")

# Parse the arguments
args <- parser$parse_args()

### ----
# Create output directory if it doesn't exist
if (!dir.exists("output/order")) {
  stop("You should run order.R first...")
}

mapthis <- readRDS("output/order/mapthis.RDS")

dropone <- droponemarker(mapthis, 
                         error.prob = args$error_prob, 
                         map.function = args$map_function,
                         tol = 1e-4)

png(file.path("output", "order", "dropone_result.png"), width = 1200, height = 1200, pointsize = 20)
par(mfrow=c(2,1))
plot(dropone, lod=1)
plot(dropone, lod=2, ylab="Change in chromosome length")
dev.off()

sum_df <- summary(dropone, lodcolumn=2)
badmar <- rownames(sum_df)[sum_df$LOD > 0]

mapthis <- drop.markers(mapthis, badmar)

# re-estimate map
nm <- est.map(mapthis, error.prob=args$error_prob, map.function=args$map_function, 
              tol = 1e-4, verbose=T, n.cluster = detectCores())
mapthis <- replace.map(mapthis, nm)

print(summaryMap(mapthis))

saveRDS(mapthis, file = file.path("output", "order", "mapthis_dropone.RDS"))

png(file.path("output", "order", "map_dropone.png"), width = 1200, height = 1200, pointsize = 20)
plotMap(mapthis, show.marker.names = F)
dev.off()
