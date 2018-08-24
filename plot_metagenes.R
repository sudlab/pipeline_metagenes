### This script plots images from the metagenes pipeline using ggplot and the specifications
# in the pipeline.yaml

library(yaml)
library(ggplot2)
library(dplyr)
library(RSQLite)
library(rlang)

thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}

if (file.exists("pipeline.yml")) {
  options = read_yaml("pipeline.yml")
} else {
  options = read_yaml(file.path(dirname(thisFile()), "pipeline_metagene/pipeline.yml"))
}
print(file.path(getSrcDirectory(function(dummy) {dummy}), "pipeline_metagene/pipeline.yml"))
database <- dbConnect(SQLite(), sub("sqlite:///", "", options$database$url))

metagene_data = dbGetQuery(database, "SELECT * FROM bam2geneprofile")

metagene_data <- mutate(metagene_data, profile = paste(source, condition, replicate, geneset, method))
itr_vals = apply(metagene_data[options$plotting$`iterate-by`], 1, paste, collapse = "_")
metagene_data = split(metagene_data, itr_vals)

cbbPalette =  c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

for (n in names(metagene_data)) {
  this_data = metagene_data[[n]]
  region_info <- this_data %>% 
    group_by(region) %>%
    summarize(end=max(bin),
              midpoint=(max(bin) + min(bin))/2)
  print(region_info)
  tick_points = region_info$midpoint
  labels = region_info$region
  dividing_lines = sort(region_info$end)[1:(length(region_info$end) -1)]
  color_by = syms(options$plotting$`colour-by`)
  p <- ggplot(this_data) +
    aes(y=!!sym(options$bam2geneprofile$profile_normalization), 
        x=bin, 
        color=paste(!!!color_by),
        group=profile) +
    geom_line() +
    geom_vline(xintercept = dividing_lines, lty=2, col="grey75") +
    scale_x_continuous(breaks=tick_points, labels=labels, name = NULL) +
    scale_y_continuous(breaks=c(0), name=NULL, limits=c(0,NA)) +
    theme_bw(base_size = 12) + 
    scale_color_manual(values=cbbPalette,
                       name = paste(options$plotting$`colour-by`, collapse=" ")) +
    theme(legend.position = "bottom")
  
  if (!is.null(options$plotting$`facet-by`)) {
    p <- p + facet_grid(as.formula(options$plotting$`facet-by`), scale="free")
  }
    
  if (!is.null(options$plotting$theme)) {
    p <- p + eval(options$plotting$theme)
  }
    
  output_name = file.path(options$export, "images", paste0(n, ".", options$plotting$format))
  ggsave(output_name, p, width = options$plotting$width, height=options$plotting$height, units="in")
  
}
