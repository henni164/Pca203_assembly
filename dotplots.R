## Eva Henningsen 02/14/2021
# load packages
library(ggplot2)
library(dplyr)
library(purrr)
library(DataCombine)
library(stringr)

# if I wanted to plot numbers, to get rid of scientific notation
options(scipen = 999)

## set custom functions
changenames <- function(element) {
  if (nchar(element) == 5) {
    element <- paste(substr(element, 1,4), "0", substr(element,5,5), sep = "")
    element <- as.character(element)
  } else
    element <- as.character(element)
  
}
changenames_vectorized <- Vectorize(changenames, vectorize.args = "element")

fixnames1 <- function(string) {
  beg <- substr(string, 1, 3)
  end <- substr(string, 9, 11)
  str2 <- paste(beg, end, sep = "")
  str2
}

# load scaffold and contig lengths files
lengths_210 <- read.delim("/Users/evahenningsen/Documents/Research/MS_projects/203_assembly/pgt_scaffold_lengths.txt", sep = "\t", header = FALSE)
lengths_203 <- read.delim("/Users/evahenningsen/Documents/Research/MS_projects/203_assembly/203_assembly_tiglengths.txt", sep = "\t", header = FALSE)
colnames(lengths_210) <- c("ref", "reflen")
colnames(lengths_203) <- c("query", "qlen")
lengths_210["ref"] <- lapply(lengths_210["ref"], changenames_vectorized)

# loop to import data and make plots
coords_files <- list.files(path="/Users/evahenningsen/Documents/Research/MS_projects/203_assembly/203_210_coords", pattern="*.coords", full.names = TRUE, recursive = FALSE)
coords_files <- as.vector(coords_files)
coords_dat = list()
coords_dat <- lapply(coords_files, FUN = read.delim, sep = "\t", header = FALSE)

for (i in 1:length(coords_dat)) {
  colnames(coords_dat[[i]]) <- c("S1", "E1", "S2", "E2", "L1", "L2", "perc_ident", "perc_sim", "perc_stop", "frame1", "frame2", "ref", "query")
  coords_dat[[i]]$query <- as.factor(coords_dat[[i]]$query)
  
}


keep_dat = list()

for (i in 1:length(coords_dat)) {
  keep_dat[[i]] <- as.data.frame(grepl.sub(coords_dat[[i]], "chr", "ref", keep.found = TRUE))
  keep_dat[[i]]["ref"] <- lapply(keep_dat[[i]]["ref"], changenames_vectorized)
  keep_dat[[i]]$ref <- as.factor(keep_dat[[i]]$ref)
}

prdat1 <- keep_dat %>% map(~ left_join(.x, lengths_210, by = "ref"))

prdat2 <- prdat1 %>% map(~ left_join(.x, lengths_203, by = "query"))

ff_dat = list()

for (i in 1:length(prdat2)) {
  ff_dat[[i]] <- prdat2[[i]][,12:15]
  colnames(ff_dat[[i]]) <- c("ref", "query", "S1", "S2")
}

plot_list = list()
plot_list2 = list()
plot_list3 = list()
plot_list4 = list()
final_plotlist = list()

for (i in 1:length(ff_dat)) {
  plot_list[[i]] <- ggplot(data = keep_dat[[i]]) +
    geom_segment(aes(x = S1, xend = E1, y = S2, yend = E2, color = perc_ident)) +
    scale_color_gradient2(low = "#C0060E", mid = "#016E14", high = "#0754B1") +
    facet_grid(query ~ ref, scales = "free", switch = "both", space = "free", labeller = labeller(query = fixnames1)) +
    theme(axis.text.x.bottom = element_blank(),
          axis.text.y.left = element_blank(),
          axis.ticks.x.bottom = element_blank(),
          axis.ticks.y.left = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.spacing.x = unit(0, "lines"),
          panel.spacing.y = unit(0, "lines"),
          strip.text.y.left = element_text(angle = 0, size = 8),
          strip.text.x = element_text(angle = 90, size = 8),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.key.height = unit(0.3, "cm"),
          legend.key.width = unit(0.3, "cm")) +
    labs(x = NULL,
         y = NULL) 

}

for (i in 1:length(plot_list)) {
  plot_list2[[i]] <- plot_list[[i]] + geom_point(data = ff_dat[[i]], x= ff_dat[[i]]$S1, y = ff_dat[[i]]$S2, alpha = 0)
  plot_list3[[i]] <- plot_list2[[i]] + geom_point(dat = ff_dat[[i]], x = 0, y = 0, alpha = 0)
  plot_list4[[i]] <- plot_list3[[i]] + geom_vline(xintercept = 0, color = "#e6e6e6")
  final_plotlist[[i]] <- plot_list4[[i]] + geom_hline(yintercept = 0, color = "#e6e6e6")
}

names_vec = list()
names_vec <- lapply(coords_files, str_match, "203_210_coords/(.*?).coords")
for (w in 1:length(coords_files)) {
  
  names_vec[[w]] <- names_vec[[w]][1,2]
  
}

for (i in 1:length(final_plotlist)) {
  
  ggsave(paste("/Users/evahenningsen/Documents/Research/MS_projects/203_assembly/", names_vec[[i]], ".tiff", sep = ""), final_plotlist[[i]], device = "tiff", width = 5, height = 5, units = "in")
  
}
