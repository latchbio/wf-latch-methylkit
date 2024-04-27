library(methylKit)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)

# Assuming the list is the first argument
sample_names <- args[1]
file_paths <- args[2]
output_dir <- args[3]
base_cov_val <- as.numeric(args[4])
tiling_val <- as.numeric(args[5])
tile_coverage <- as.numeric(args[6])
difference_val <- as.numeric(args[7])
q_val <- as.numeric(args[8])
treatments <- args[9]
file_format <- args[10]
genome <- args[11]


treatments <- as.numeric(unlist(strsplit(treatments, ",")))
delimiter <- ","
sample_names <- strsplit(sample_names, delimiter, fixed = TRUE)[[1]]
sample_names <- as.list(sample_names)
file_paths <- strsplit(file_paths, delimiter, fixed = TRUE)[[1]]
file_paths <- as.list(file_paths)


if (file_format == "bedmethyl") {
  print("Processing bedmethyl files")
  myobj=methRead(file_paths,
               sample.id=sample_names,
               assembly=genome,
               treatment=treatments,
               context="CpG",
               pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=3,
                             coverage.col=10,strand.col=6,freqC.col=12)
)
} else if (file_format == "bismark_cov") {
  print("Processing Bismark Coverage files")
  myobj=methRead(file_paths,
                sample.id=sample_names,
                assembly=genome,
                treatment=treatments,
                pipeline="bismarkCoverage",
                context="CpG")
} else {
  print("Processing Bismark Cytosine Report files")
  myobj=methRead(file_paths,
                sample.id=sample_names,
                assembly=genome,
                treatment=treatments,
                pipeline="bismarkCytosineReport",
                context="CpG")

}




for (i in seq_along(myobj)) {
  # Define the file name for the plot
  plot_file <- sprintf("%s/methy_stats_plot_%d.png", output_dir, i)
  png(filename = plot_file, width = 800, height = 600)
  # Generate and plot methylation stats
  getMethylationStats(myobj[[i]], plot = TRUE, both.strands = FALSE)
  dev.off()
  print(getMethylationStats(myobj[[i]], plot = FALSE, both.strands = FALSE))

  plot_file <- sprintf("%s/coverage_stats_plot_%d.png", output_dir, i)
  png(filename = plot_file, width = 800, height = 600)
  # Generate and plot methylation stats
  getCoverageStats(myobj[[i]], plot = TRUE, both.strands = FALSE)
  dev.off()
  print(getCoverageStats(myobj[[i]], plot = FALSE, both.strands = FALSE))

}

print("----------------------------------")
print("MethylKit Object")
print("----------------------------------")
print(myobj)
filtered.myobj=filterByCoverage(myobj,lo.count=base_cov_val,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)
tiles=tileMethylCounts(myobj,win.size=tiling_val,step.size=tiling_val,cov.bases=tile_coverage)
meth=unite(tiles, destrand=TRUE)
print("----------------------------------")
print("Tiles Methylation Data")
print("----------------------------------")
print(meth)

# Print correlation plot
plot_file <- sprintf("%s/correlation_plot.png", output_dir)
png(filename = plot_file, width = 800, height = 600)
getCorrelation(meth,plot=TRUE)
dev.off()

print("----------------------------------")
print("Correlation Data")
print("----------------------------------")
getCorrelation(meth,plot=TRUE)

myDiff=calculateDiffMeth(meth)

print("----------------------------------")
print("Correlation Data")
print("----------------------------------")
myDiff25p=getMethylDiff(myDiff,difference=difference_val,qvalue=q_val)

print("----------------------------------")
print("DMR Regions")
print("----------------------------------")
print(myDiff25p)

# write.csv(myDiff25p@.Data,
#           paste0(output_dir, "/DMR_regions.csv"), row.names=FALSE)
write.csv(myDiff25p,
          paste0(output_dir, "/DMR_regions.csv"), row.names=FALSE)
write.csv(meth,
          paste0(output_dir, "/merged_methylation_data.csv"), row.names=FALSE)

my_data <- read.csv(paste0(output_dir, "/DMR_regions.csv"))
new_column_names <- c("chr", "start", "end", "strand", "p_value", "q_value", "percent_difference")
colnames(my_data) <- new_column_names
write.csv(my_data, paste0(output_dir, "/DMR_regions.csv"), row.names = FALSE)

