library(ggplot2)
library(ggpubr)

options(scipen = 999)

pca_and_pgt <- read.delim("/Users/evahenningsen/Documents/map_Pgt_201_A1_AssemblyScaffolds_to_Oat_crown_rust_A_cleaned.paf", sep = "\t", header = FALSE)
pca_and_pgt <- pca_and_pgt[,-(13:17)]
colnames(pca_and_pgt) <- c("qname", "qlength", "qstart", "qend", "strand", "tname", "tlength", "tstart", "tend", "nmatches", "alnlength", "mquality")
pca_and_pgt <- pca_and_pgt[grep("chr", pca_and_pgt$qname),]

pca_and_pt <- read.delim("/Users/evahenningsen/Documents/map_chrs_LR1_A_HiFi_to_Oat_crown_rust_A_cleaned.paf", sep = "\t", header = FALSE)
pca_and_pt <- pca_and_pt[,-(13:17)]
colnames(pca_and_pt) <- c("qname", "qlength", "qstart", "qend", "strand", "tname", "tlength", "tstart", "tend", "nmatches", "alnlength", "mquality")

chr_lengths <- read.delim("/Users/evahenningsen/Documents/Research/203_assembly/chr_lengths.txt", header = TRUE, sep = "\t")
CHRNAMES <- as.vector(unique(chr_lengths$chr))
CHRNAMES_203 <- CHRNAMES[37:54]

pca_and_pgt$tname <- factor(pca_and_pgt$tname, levels = CHRNAMES_203[c(18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1)])
pca_and_pgt$qname <- factor(pca_and_pgt$qname, levels = c("chr_1","chr_2","chr_3","chr_4","chr_5","chr_6","chr_7","chr_8","chr_9","chr_10","chr_11","chr_12","chr_13","chr_14","chr_15","chr_16","chr_17","chr_18"))
pca_and_pgt$mquality <- as.numeric(pca_and_pgt$mquality)

pca_and_pt$tname <- factor(pca_and_pt$tname, levels = CHRNAMES_203[c(18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1)])
pca_and_pt$qname <- factor(pca_and_pt$qname, levels = c("chr1_A", "chr2_A", "chr3_A", "chr4_A", "chr5_A", "chr6_A", "chr7_A", "chr8_A", "chr9_A", "chr10_A", "chr11_A", "chr12_A", "chr13_A", "chr14_A", "chr15_A", "chr16_A", "chr17_A", "chr18_A"))
pca_and_pt$mquality <- as.numeric(pca_and_pt$mquality)

levels(pca_and_pgt$tname) <- list("18" = "Pca203_chr18_A", "17" = "Pca203_chr17_A", "16" = "Pca203_chr16_A", "15" = "Pca203_chr15_A", "14" = "Pca203_chr14_A", "13" = "Pca203_chr13_A", "12" = "Pca203_chr12_A", "11" = "Pca203_chr11_A", "10" = "Pca203_chr10_A", "9" = "Pca203_chr9_A", "8" = "Pca203_chr8_A", "7" = "Pca203_chr7_A", "6" = "Pca203_chr6_A", "5" = "Pca203_chr5_A", "4" = "Pca203_chr4_A", "3" = "Pca203_chr3_A", "2" = "Pca203_chr2_A", "1" = "Pca203_chr1_A")
levels(pca_and_pt$tname) <- list("18" = "Pca203_chr18_A", "17" = "Pca203_chr17_A", "16" = "Pca203_chr16_A", "15" = "Pca203_chr15_A", "14" = "Pca203_chr14_A", "13" = "Pca203_chr13_A", "12" = "Pca203_chr12_A", "11" = "Pca203_chr11_A", "10" = "Pca203_chr10_A", "9" = "Pca203_chr9_A", "8" = "Pca203_chr8_A", "7" = "Pca203_chr7_A", "6" = "Pca203_chr6_A", "5" = "Pca203_chr5_A", "4" = "Pca203_chr4_A", "3" = "Pca203_chr3_A", "2" = "Pca203_chr2_A", "1" = "Pca203_chr1_A")
levels(pca_and_pgt$qname) <- list("1" = "chr_1", "2" = "chr_2", "3" = "chr_3", "4" = "chr_4", "5" = "chr_5", "6" = "chr_6", "7" = "chr_7", "8" = "chr_8", "9" = "chr_9", "10" = "chr_10", "11" = "chr_11", "12" = "chr_12", "13" = "chr_13", "14" = "chr_14", "15" = "chr_15", "16" = "chr_16", "17" = "chr_17", "18" = "chr_18")
levels(pca_and_pt$qname) <- list("1" = "chr1_A", "2" = "chr2_A", "3" = "chr3_A", "4" = "chr4_A", "5" = "chr5_A", "6" = "chr6_A", "7" = "chr7_A", "8" = "chr8_A", "9" = "chr9_A", "10" = "chr10_A", "11" = "chr11_A", "12" = "chr12_A", "13" = "chr13_A", "14" = "chr14_A", "15" = "chr15_A", "16" = "chr16_A", "17" = "chr17_A", "18" = "chr18_A")

pca_and_pgt2 <- subset(pca_and_pgt, mquality > 0)
pca_and_pt2 <- subset(pca_and_pt, mquality > 0)


pca_pgt_plot <- ggplot(data = pca_and_pgt2, aes(x = qstart, y = tstart)) +
  geom_point(aes(color = mquality), size = 0.05) +
  facet_grid(rows = vars(tname), cols = vars(qname), scales = "free") +
  scale_color_gradient(low = "#000000", high = "#F36B6B", aesthetics = "color") +
  labs(x = "Pgt21-0", y = "Pca203") +
  theme(axis.text = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.y = element_text(angle = 0, color = "black", size = 10),
        strip.text.x = element_text(color = "black", size = 10),
        axis.title.x = element_text(color = "black", size = 10),
        axis.title.y = element_text(color = "black", size = 10),
        axis.ticks = element_blank(),
        legend.position = "none")

colnames(pca_and_pt2)[12] <- "MAPQ"

pca_pt_plot <- ggplot(data = pca_and_pt2, aes(x = qstart, y = tstart)) +
  geom_point(aes(color = MAPQ), size = 0.05) +
  facet_grid(rows = vars(tname), cols = vars(qname), scales = "free") +
  scale_color_gradient(low = "#000000", high = "#F36B6B", aesthetics = "color") +
  labs(x = "Pt76", y = NULL) +
  theme(axis.text = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.y = element_text(angle = 0, color = "black", size = 10),
        strip.text.x = element_text(color = "black", size = 10),
        axis.title.x = element_text(color = "black", size = 10),
        axis.ticks = element_blank(),
        legend.title = element_text(color = "black", size = 10),
        legend.text = element_text(color = "black", size = 10))

combo_plot <- ggarrange(pca_pgt_plot, pca_pt_plot, ncol = 2, nrow = 1, widths = c(1,1.2), font.label = list(size = 10))

ggsave("/Users/evahenningsen/Documents/Research/203_assembly/full_genome_dotplots_no_multi.tiff", combo_plot, device = "tiff", width = 7.5, height = 3, units = "in")

pca_pgt_chr1 <- subset(pca_and_pgt2, tname == "1")
pca_pgt_chr1 <- subset(pca_pgt_chr1, qname == "1")
pca_pt_chr1 <- subset(pca_and_pt2, tname == "1")
pca_pt_chr1 <- subset(pca_pt_chr1, qname == "1")

pca_pgt_chr1$qstart <- (pca_pgt_chr1$qstart/1000000)

pca_pgt_chr1$tstart <- (pca_pgt_chr1$tstart/1000000)
pca_pt_chr1$qstart <- (pca_pt_chr1$qstart/1000000)
pca_pt_chr1$tstart <- (pca_pt_chr1$tstart/1000000)


pca_pgt_chr1_plot <- ggplot(data = pca_pgt_chr1, aes(x = qstart, y = tstart)) +
  geom_point(aes(color = mquality), size = 0.05) +
  scale_color_gradient(low = "#000000", high = "#F36B6B", aesthetics = "color") +
  labs(x = "Pgt21-0 chromosome 1 (Mbp)", y = "Pca203 chromosome 1 (Mbp)") +
  theme(#axis.text = element_blank(),
    panel.spacing = unit(0, "lines"),
    #axis.ticks = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10))
pca_pt_chr1_plot <- ggplot(data = pca_pt_chr1, aes(x = qstart, y = tstart)) +
  geom_point(aes(color = MAPQ), size = 0.05) +
  scale_color_gradient(low = "#000000", high = "#F36B6B", aesthetics = "color") +
  labs(x = "Pt76 chromosome 1 (Mbp)", y = NULL) +
  theme(#axis.text = element_blank(),
    panel.spacing = unit(0, "lines"),
    strip.text.y = element_text(angle = 0),
    #axis.ticks = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10))

combo_plot_chr1 <- ggarrange(pca_pgt_chr1_plot, pca_pt_chr1_plot, ncol = 2, nrow = 1, widths = c(1,1.3))

ggsave("/Users/evahenningsen/Documents/Research/203_assembly/chr1_dotplots_no_multi.tiff", combo_plot_chr1, device = "tiff", width = 7.5, height = 3, units = "in")

