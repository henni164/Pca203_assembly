library(ggplot2)
library(dplyr)
library(ggchicklet)

coord_dat_203 <- read.delim(file = "/Users/evahenningsen/Documents/Research/203_assembly/full_table_buscov3_Oat_crown_rust_203_chrsA.tsv", header = FALSE, sep = "\t")
coord_dat_pgt <- read.delim(file = "/Users/evahenningsen/Documents/Research/203_assembly/full_table_buscov3_Stem_rust_chrsA.tsv", header = FALSE, sep = "\t")
coord_dat_lr <- read.delim(file = "/Users/evahenningsen/Documents/Research/203_assembly/full_table_buscov3_Leaf_rust_chrsA.tsv", header = FALSE, sep = "\t")
colnames(coord_dat_203) <- c("BUSCOID", "Status", "chr", "start", "end", "s1", "s2")
colnames(coord_dat_pgt) <- c("BUSCOID", "Status", "chr", "start", "end", "s1", "s2")
colnames(coord_dat_lr) <- c("BUSCOID", "Status", "chr", "start", "end", "s1", "s2")

merged_pgt_203 <- merge(coord_dat_pgt, coord_dat_203, by = "BUSCOID", all = FALSE)

merged_pgt_lr <- merge(coord_dat_pgt, coord_dat_lr, by = "BUSCOID", all = FALSE)


full_dat <- rbind(merged_pgt_203, merged_pgt_lr)

chr_lengths <- read.delim("/Users/evahenningsen/Documents/Research/203_assembly/chr_lengths.txt", header = TRUE, sep = "\t")
CHRNAMES_203 <- as.vector(unique(chr_lengths$chr))
CHRNAMES_203 <- CHRNAMES_203[37:54]
colnames(chr_lengths)[1] <- "chr.x"
chr_lengths$chr.y <- chr_lengths$chr.x

with_chrlen_203 <- merge(full_dat, chr_lengths, by = "chr.x", all = FALSE)
colnames(with_chrlen_203)[14] <- "chr_len.x"
with_chrlen_203 <- with_chrlen_203[,-16]
colnames(with_chrlen_203)[9] <- "chr.y"

with_chrlen_203 <- merge(with_chrlen_203, chr_lengths, by = "chr.y", all = FALSE)
colnames(with_chrlen_203)[2] <- "chr.x"
with_chrlen_203 <- with_chrlen_203[,-16]
colnames(with_chrlen_203)[16] <- "chr_len.y"

CHRNAMES_210 <- as.vector(chr_lengths$chr.x)
CHRNAMES_210 <- CHRNAMES_210[1:18]
CHRNAMES_lr <- as.vector(chr_lengths$chr.x)
CHRNAMES_lr <- CHRNAMES_lr[19:36]

with_chrlen_203$chr.x <- factor(with_chrlen_203$chr.x, levels = CHRNAMES_210[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)])


fixed_203 <- with_chrlen_203 %>% 
  group_by(chr.x) %>% 
  summarise(chr_length.x=max(chr_len.x)) %>% 
  mutate(tot.x=cumsum(chr_length.x)-chr_length.x) %>%
  select(-chr_length.x) %>%
  left_join(with_chrlen_203, ., by=c("chr.x"="chr.x")) %>%
  arrange(factor(chr.x, levels = CHRNAMES_210[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)]), start.x) %>%
  mutate(bpcum.x=start.x+tot.x)

fixed_210_203 <- fixed_203[grep("Pca", fixed_203$chr.y),]
fixed_210_lr <- fixed_203[grep("Pt", fixed_203$chr.y),]

fixed_210_203$chr.y <- factor(fixed_210_203$chr.y, levels = CHRNAMES_203[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)])

fixed_210 <- fixed_210_203 %>% 
  group_by(chr.y) %>% 
  summarise(chr_length.y=max(chr_len.y)) %>% 
  mutate(tot.y=cumsum(chr_length.y)-chr_length.y) %>%
  select(-chr_length.y) %>%
  left_join(fixed_210_203, ., by=c("chr.y"="chr.y")) %>%
  arrange(factor(chr.y, levels = CHRNAMES_203[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)]), start.y) %>%
  mutate(bpcum.y=start.y+tot.y)

fixed_210_lr$chr.y <- factor(fixed_210_lr$chr.y, levels = CHRNAMES_lr[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)])

fixed_LR1 <- fixed_210_lr %>% 
  group_by(chr.y) %>% 
  summarise(chr_length.y=max(chr_len.y)) %>% 
  mutate(tot.y=cumsum(chr_length.y)-chr_length.y) %>%
  select(-chr_length.y) %>%
  left_join(fixed_210_lr, ., by=c("chr.y"="chr.y")) %>%
  arrange(factor(chr.y, levels = CHRNAMES_lr[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)]), start.y) %>%
  mutate( bpcum.y=start.y+tot.y)

fixed_210$num.y <- rep("1", length(fixed_210$chr.x))
fixed_LR1$num.y <- rep("3", length(fixed_LR1$chr.x))

really_finished_dat <- rbind(fixed_210, fixed_LR1)
really_finished_dat$num.y <- as.numeric(really_finished_dat$num.y)

plotting_dat <- read.delim(file = "/Users/evahenningsen/Documents/Research/203_assembly/RNAseq/start_end.txt", header = TRUE, sep = "\t")
plotting_dat$color <- as.factor(plotting_dat$color)

plotting_dat2 <- plotting_dat[plotting_dat$color != "2",]

plotting_dat3 <- plotting_dat[plotting_dat$color != "3",]

plotting_dat4 <- plotting_dat[plotting_dat$color != "1",]

# blue is "#619CFF", green is "#00BA38", and red is "#F8766D"

out_FINAL <- really_finished_dat %>%
  group_by(association_number.x) %>%
  mutate(new_bpcum.x = bpcum.x + ((association_number.x-1)*1000000)) %>%
  group_by(association_number.y) %>%
  mutate(new_bpcum.y = bpcum.y + ((association_number.y-1)*1000000))


all_plot <- ggplot() +
  geom_segment(data = out_FINAL, aes(x = 2, xend = num.y, y = new_bpcum.x, yend = new_bpcum.y), size = 0.05) +
  geom_rect(data = plotting_dat, aes(xmin = xstart, xmax = xend, ymin = ystart, ymax = yend, fill = color)) +
  geom_segment(data = plotting_dat, aes(x = xstart, xend = xend, y = centro, yend = centro), color = "red", size = 0.4) + 
  scale_fill_manual(values = c("#00BA38","#F8766D", "#619CFF")) +
  scale_x_continuous(breaks = c(1,2,3), labels = c(expression(paste(italic("Pca"), 203)), expression(paste(italic("Pgt"), 21,'-',0)), expression(paste(italic("Pt"), 76)))) +
  scale_y_continuous(name = NULL, breaks = NULL, labels = NULL, expand = c(0.01,0.01)) +
  labs(x = NULL, y = NULL) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(margin = margin(0,0,0,0, "cm"), size = 10))

ggsave("/Users/evahenningsen/Documents/Research/203_assembly/whole_synteny_plot_buscos.tiff", plot = all_plot, device = "tiff", width = 3, height = 8, unit = "in")
