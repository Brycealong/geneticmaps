# ripple.R
#install.packages("qtl")
library(qtl)
library(argparse)
options(timeout = 240)

# Create a parser object
parser <- ArgumentParser()

# Define the arguments
parser$add_argument("-e", "--error_prob", type = "double", default = 0.0001, 
                    help = "Assumed genotyping error rate used in the final estimated map. (default %(default)s)")
parser$add_argument("-f", "--map_function", type = "character", choices = c("haldane","kosambi","c-f","morgan"), default = "haldane", 
                    help = "map function to use. default is 'haldane'.")
parser$add_argument("-w", "--window", type = "integer", default = 3,
                    help = "window size used to ripple. (default %(default)s)")

# Parse the arguments
args <- parser$parse_args()


### ----
# Create output directory if it doesn't exist
if (!dir.exists("output/ripple")) {
  dir.create("output/ripple", recursive = T)
}


### Rippling ----
mapthis <- readRDS("output/order/mapthis.RDS")
cat(paste0(
  "Ripple the order using window size ",
  args$window,
  "\n"
))
for (chr in chrnames(mapthis)) {
  result <- tryCatch({
    rip <- ripple(mapthis, chr = chr, window = args$window, method = "likelihood", 
                  error.prob = args$error_prob, map.function = args$map_function,
                  verbose = F)
    mapthis <- switch.order(mapthis, chr, rip[2,])
    TRUE  # Indicating success
  }, error = function(e) {
    message(paste("Error with chromosome", chr, ":", e$message))
    FALSE  # Indicating failure
  })
  
  if (!result) {
    message(paste("Skipping chromosome", chr, "due to an error."))
    next  # Skip to the next iteration if there was an error
  }
}
print(summaryMap(mapthis))
map <- pull.map(mapthis)
maptbl <- map2table(map)
write.table(maptbl, file = "output/ripple/sum.csv", quote = F, 
            sep = ",")
saveRDS(mapthis, file = file.path("output", "ripple", "mapthis.RDS"))

