library(vcfR)
library(ggplot2)
library(tidyr)

vcf_203 <- read.vcfR("/Users/evahenningsen/Documents/GitHub/Pca203_assembly/figure_s1/203to12sd80.vcf", verbose = FALSE)

#Substitute non-het positions with NA values
vcf_203_gt <- extract.gt(vcf_203, element = 'GT')
vcf_203_hets <- is_het(vcf_203_gt)
is.na(vcf_203@gt[,-1][ !vcf_203_hets]) = TRUE

#Filter het on allele depth
vcf_203_ad = extract.gt(vcf_203, element = 'AD')
vcf_203_allele1 = masplit(vcf_203_ad, record = 1)
vcf_203_allele2 = masplit(vcf_203_ad, record = 2)

# Produce quantiles at 0.15 and 0.95 probabilites
sums_203_ERR = apply(vcf_203_allele1, MARGIN=2, quantile, probs=c(0.15, 0.95), na.rm=TRUE)

# Allele 1
vcf_203_dp2 = sweep(vcf_203_allele1, MARGIN=2, FUN = "-", sums_203_ERR[1,])

# allele1[dp2 < 0] = NA
vcf_203@gt[,-1][ vcf_203_dp2 < 0 & !is.na(vcf_203@gt[,-1]) ] = NA
vcf_203_dp2 = sweep(vcf_203_allele1, MARGIN=2, FUN = "-", sums_203_ERR[2,])

# allele1[dp2 > 0] = NA
vcf_203@gt[,-1][vcf_203_dp2 > 0] = NA

#Allele 2
vcf_203_dp2 = sweep(vcf_203_allele2, MARGIN=2, FUN = "-", sums_203_ERR[1,])
vcf_203@gt[,-1][ vcf_203_dp2 < 0 & !is.na(vcf_203@gt[,-1]) ] = NA
vcf_203_dp2 = sweep(vcf_203_allele2, MARGIN=2, FUN = "-", sums_203_ERR[2,])
vcf_203@gt[,-1][vcf_203_dp2 > 0] = NA

# Calculate allele balance with the new filtered vcf
vcf_203_ad = extract.gt(vcf_203, element = 'AD')
vcf_203_allele1 = masplit(vcf_203_ad, record = 1)
vcf_203_allele2 = masplit(vcf_203_ad, record = 2)
vcf_203_ad1 = vcf_203_allele1 / (vcf_203_allele1 + vcf_203_allele2)
vcf_203_ad2 = vcf_203_allele2 / (vcf_203_allele1 + vcf_203_allele2)

# Make a tidy dataframe to make plot
# NOTE: don't need to remove the reference genome (as there is none here)
vcf_203_ad1t = tidyr::gather(tibble::as_tibble(vcf_203_ad1), "Sample", "Ab")
vcf_203_ad1t$Allele = "ab1"

vcf_203_ad2t = tidyr::gather(tibble::as_tibble(vcf_203_ad2), "Sample", "Ab")
vcf_203_ad2t$Allele = "ab2"
vcf_203_t = rbind(vcf_203_ad1t, vcf_203_ad2t)
vcf_203_t = dplyr::filter(vcf_203_t, !is.na(Ab))

#plot
alfq_203_plot <- ggplot(vcf_203_t, aes(x = Ab)) +
  scale_x_continuous(breaks = c(0,1/2,1), labels = c('0', '1/2', '1')) +
  geom_histogram(data = subset(vcf_203_t, Allele == 'ab1'), fill = "#1f78b4", binwidth = 0.02) +
  geom_histogram(data = subset(vcf_203_t, Allele == 'ab2'), fill = "#a6cee3", binwidth = 0.02) +
  theme_bw() +
  xlab("Allele balance") +
  ylab("Heterozygous positions") +
  theme(panel.grid.major.x=element_line(color = "#A9A9A9", size=0.4),
        panel.grid.major.y=element_line(color = "#A9A9A9", size=0.1, linetype = 'dashed'),
        panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        strip.text.x = element_text(size = 8),
        axis.text.x = element_text(angle = 60, hjust = 1, size=8),
        axis.text.y = element_text(angle = 30, hjust = 1, size=10),
        axis.title.x = element_text(size = 14, angle = 0),
        axis.title.y = element_text(size = 14, angle = 90))

ggsave(filename = "/Users/evahenningsen/Documents/GitHub/Pca203_assembly/figure_s1/203_alfq_plot.tiff", plot = alfq_203_plot, device = "tiff", width = 2, height = 3, units = "in")
