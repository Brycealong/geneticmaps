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
  
  ## cleaning large gap
  # for (chr in chrnames(mapthis)){
  #   maptbl <- map2table(pull.map(mapthis, chr))
  #   iqr <- IQR(maptbl$pos)
  #   lower_bound <- quantile(maptbl$pos, 0.25) - 1.5 * iqr
  #   upper_bound <- quantile(maptbl$pos, 0.75) + 1.5 * iqr
  #   badmar <- rownames(maptbl)[which(maptbl$pos > upper_bound | maptbl$pos < lower_bound)]
  #   mapthis <- drop.markers(mapthis, badmar)
  # }
  # # re-estimate map and rf
  # newmap <- est.map(mapthis, error.prob = args$error_prob, map.function = args$map_function)
  # mapthis <- replace.map(mapthis, newmap)
  # mapthis <- est.rf(mapthis)

  # show map
  print(summaryMap(mapthis))
} else if (args$by == "infer"){
  for (chr in chrnames(mapthis)) {
    mapthis <- orderMarkers(mapthis, chr = chr, window = args$window, error.prob = args$error_prob, map.function = args$map_function)
  }
  ## cleaning large gap
  # for (chr in chrnames(mapthis)){
  #   maptbl <- map2table(pull.map(mapthis, chr))
  #   badmar <- rownames(maptbl)[which(maptbl$pos > quantile(maptbl$pos, 0.9))]
  #   mapthis <- drop.markers(mapthis, badmar)
  # }
  # # re-estimate map and rf
  # newmap <- est.map(mapthis, error.prob = args$error_prob, map.function = args$map_function)
  # mapthis <- replace.map(mapthis, newmap)
  # mapthis <- est.rf(mapthis)

  print(summaryMap(mapthis))
} else {
  stop("Please specify one method to order the markers.")
}

iteration_count <- 0  # Initialize the counter

cat("Start removing bad markers...")
repeat {
  iteration_count <- iteration_count + 1  # Increment the counter
  cat("Iteration: ", iteration_count, "\n")  # Print the current iteration
  
  dropone <- droponemarker(mapthis, error.prob = args$error_prob, map.function = args$map_function, verbose = F)
  sum_df <- summary(dropone, lodcolumn=2)
  badmar <- rownames(sum_df)[sum_df$LOD > 0]
  
  if (length(badmar) == 0) {
    cat("No more bad markers found. Exiting loop.\n")
    break
  }
  
  cat("Dropping bad markers:", paste(badmar, collapse=", "), "\n")
  
  mapthis <- drop.markers(mapthis, badmar)
  
  # re-estimate map and rf
  newmap <- est.map(mapthis, error.prob = args$error_prob, map.function = args$map_function)
  mapthis <- replace.map(mapthis, newmap)
}

mapthis <- est.rf(mapthis)

map <- pull.map(mapthis)
maptbl <- map2table(map)
write.table(maptbl, file = "output/order/sum.csv", quote = F, 
            sep = ",")
saveRDS(mapthis, file = file.path("output", "order", "mapthis.RDS"))

