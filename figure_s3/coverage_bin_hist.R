library(ggplot2)

# again on cleaned-up assembly
df2 <- read.delim('/Users/evahenningsen/Documents/Research/203_assembly/Old_203_genome_versions/Collapsed_regions_analysis/203_bincoverage_cleaned.txt', sep="\t", header=TRUE, skip=2)
head(df2)

data2 <- subset(df2, select = c("Cov"))

range_to_show = 400

cov_hist <- ggplot(data2, aes(x=Cov)) + geom_histogram(binwidth=1) + 
  xlim(0, range_to_show) +
  labs(x = "Coverage", y = "Count") +
  geom_vline(xintercept = 117, linetype = "longdash") +
  geom_vline(xintercept = 234, linetype = "dotted") +
  geom_vline(xintercept = 175.5, linetype = "solid")

max(cov_hist$data)

ggsave("/Users/evahenningsen/Documents/Research/203_assembly/Old_203_genome_versions/Collapsed_regions_analysis/cov_hist.tiff", cov_hist, device = "tiff", width = 4, height = 4, units = "in")
