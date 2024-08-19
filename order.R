# order.R
#install.packages("qtl")
library(qtl)
library(argparse)
options(timeout = 240)

# Create a parser object
parser <- ArgumentParser()

# Define the arguments
parser$add_argument("--by", type = "character", choices = c("obs", "infer"), default = "obs",
                    help = "Specify 'obs' to use predefined orders (observation) but re-estimate the genetic distances or 'infer' to reorder markers (inference)")
parser$add_argument("-e", "--error_prob", type = "double", default = 0.0001, 
                    help = "Assumed genotyping error rate used in the final estimated map. (default: %(default)s)")
parser$add_argument("-f", "--map_function", type = "character", choices = c("haldane","kosambi","c-f","morgan"), default = "haldane", 
                    help = "map function to use. default is 'haldane'.")
parser$add_argument("-w", "--window", type = "integer", default = 3,
                    help = "window size used to ripple. (default: %(default)s)")


# Parse the arguments
args <- parser$parse_args()

### ----
# Create output directory if it doesn't exist
if (!dir.exists("output/order")) {
  dir.create("output/order", recursive = T)
}

### Ordering ----
mapthis <- readRDS("output/group/mapthis.RDS")
if (args$by == "obs") {
  cat("Using the input orders...\n")
  newmap <- est.map(mapthis, error.prob = args$error_prob, map.function = args$map_function)
  mapthis <- replace.map(mapthis, newmap)
  # show map
  print(summaryMap(mapthis))
} else if (args$by == "infer"){
  for (chr in chrnames(mapthis)) {
    mapthis <- orderMarkers(mapthis, chr = chr, window = args$window, error.prob = args$error_prob, map.function = args$map_function)
  }
  print(summaryMap(mapthis))
} else {
  stop("Please specify one method to order the markers.")
}

map <- pull.map(mapthis)
maptbl <- map2table(map)
write.table(maptbl, file = "output/order/sum.csv", quote = F, 
            sep = ",")
saveRDS(mapthis, file = file.path("output", "order", "mapthis.RDS"))

