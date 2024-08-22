# order.R
#install.packages("qtl")
library(qtl)
library(argparse)
library(parallel)
library(snow)
options(timeout = 240)
set.seed(61777369)

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

  # show map
  print(summaryMap(mapthis))
} else if (args$by == "infer"){
  # use rippling
  rip <- vector("list", nchr(mapthis))
  names(rip) <- chrnames(mapthis)
  for(i in chrnames(mapthis)){
    rip[[i]] <- ripple(fill.geno(mapthis, error.prob=args$error_prob), i, window=args$window,
                       error.prob=args$error_prob, map.function=args$map_function, 
                       tol = 1e-4, verbose=T, n.cluster = detectCores())
  }
  dif.nxo <- sapply(rip, function(a) a[1,ncol(a)]-a[2,ncol(a)])
  for(i in chrnames(mapthis)) {
    if(dif.nxo[i] > 0){
      mapthis <- switch.order(mapthis, i, rip[[i]][2,], 
                              error.prob=args$error_prob, map.function=args$map_function,
                              tol = 1e-4)
    }
  }

  print(summaryMap(mapthis))
} else {
  stop("Please specify one method to order the markers.")
}

nm <- est.map(mapthis, error.prob=args$error_prob, map.function=args$map_function, 
              tol = 1e-4, verbose=T, n.cluster = detectCores())
mapthis <- replace.map(mapthis, nm)

# iteration_count <- 0  # Initialize the counter
# 
# cat("Start removing bad markers...\n")
# repeat {
#   iteration_count <- iteration_count + 1  # Increment the counter
#   cat("Iteration: ", iteration_count, "\n")  # Print the current iteration
#   
#   dropone <- droponemarker(mapthis, error.prob = args$error_prob, map.function = args$map_function)
#   sum_df <- summary(dropone, lodcolumn=2)
#   badmar <- rownames(sum_df)[sum_df$LOD > 0]
#   
#   if (length(badmar) == 0) {
#     cat("No more bad markers found. Exiting loop.\n")
#     break
#   }
#   
#   cat("Dropping bad markers:", paste(badmar, collapse=", "), "\n")
#   
#   mapthis <- drop.markers(mapthis, badmar)
#   
#   # re-estimate map and rf
#   newmap <- est.map(mapthis, error.prob = args$error_prob, map.function = args$map_function)
#   mapthis <- replace.map(mapthis, newmap)
# }
# 
# mapthis <- est.rf(mapthis)

map <- pull.map(mapthis)
maptbl <- map2table(map)
write.table(maptbl, file = "output/order/sum.csv", quote = F, 
            sep = ",")
saveRDS(mapthis, file = file.path("output", "order", "mapthis.RDS"))

png(file.path("output", "order", "map.png"), width = 1200, height = 1200, pointsize = 20)
plotMap(mapthis, show.marker.names = F)
dev.off()