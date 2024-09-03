library(tidyverse)
library(gridExtra)

pheno = read.table(snakemake@input[["donor_meta_data"]], header = TRUE, sep = "\t")
batches_df = read.table(snakemake@input[["batch_table"]], header = TRUE, sep = "\t")
row.names(batches_df) = batches_df$IID

### Combine all the qc info
missing_qc = read.table(snakemake@input[["missing_qc"]], header = FALSE, sep = "")
names(missing_qc) <- c("FID", "IID", "NMISS", "N", "FRQ_MISS")
sex_qc = read.table(snakemake@input[["sex_qc"]], header = TRUE, sep = "")
ancestry_qc = read.table(snakemake@input[["ancestry_qc"]], header = TRUE, sep = "\t")

sample_qc = merge(batches_df[, c("IID", "Donor", "Batch", "chip_version")], merge(ancestry_qc, merge(missing_qc, sex_qc, by = "IID"), by = "IID"), by = "IID")
sample_qc <- merge(sample_qc, pheno, by = "Donor", all.x = TRUE)

sample_qc <- sample_qc %>%
  mutate(
    Sex = case_match(
      Sex,
      c("Male", "M", "m") ~ "M",
      c("Female", "F", "f") ~ "F",
      c("NA", "nr", "N", "0", " ", "") ~ "0"
    )
  )

### Save sample qc table
write.table(sample_qc, file = "results/geno_qc_summary.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


### Check that duplicates are MZ
duplicated_donors = unique(sample_qc$Donor[duplicated(sample_qc$Donor)])
related_qc = read.table("results/kingrel.kin", header = TRUE, sep = "\t")
MZ = related_qc[related_qc$InfType == "Dup/MZ", ]
for (i in 1:length(duplicated_donors)) {
  if (!duplicated_donors[i] %in% c(MZ$ID1, MZ$ID2)) {
    stop("ERROR: Duplicate donors IDs are not reflected in genotype-based relationship inference.\n")
  }
}
cat("The following donors are duplicated across genotyping batches:\n")
cat(duplicated_donors, sep = "\n")
cat("\n")
cat("Please confirm that relationship inference is consistent:\n")
print(MZ)
cat("\n")


### Choose sample among duplicates
duplicated_df = sample_qc[sample_qc$Donor %in% duplicated_donors, c("FID", "IID", "Donor", "FRQ_MISS")]
duplicated_df = duplicated_df[order(duplicated_df$Donor, duplicated_df$FRQ_MISS, decreasing = TRUE), ]
duplicates_to_drop = duplicated_df[duplicated(duplicated_df$Donor), c("FID", "IID")]
duplicates_to_keep = duplicated_df[!duplicated(duplicated_df$Donor), c("FID", "IID", "FID", "Donor")]
write.table(duplicates_to_drop, file = "results/duplicates_to_drop.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(duplicates_to_keep, file = "results/duplicates_to_update.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


### Sex update file -- can be used to update ped file to reflect inferred sex
sample_qc <- sample_qc[!sample_qc$IID %in% duplicates_to_drop$IID, ]
write.table(sample_qc[, c("FID", "IID", "SNPSEX")], file = "results/sex_update.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


### QC Plots
ancestry_levels = names(c(table(sample_qc$Ancestry)))
ancestry_labels = NULL
for (i in 1:length(ancestry_levels)) {
  ancestry_labels[i] = paste0(ancestry_levels[i], " (N = ", c(table(sample_qc$Ancestry))[i], ")")
}

sex_levels = names(c(table(sample_qc$Sex)))
sex_labels = NULL
for (i in 1:length(sex_levels)) {
  sex_labels[i] = paste0(sex_levels[i], " (N = ", c(table(sample_qc$Sex))[i], ")")
}

chip_levels = names(c(table(sample_qc$chip_version)))
chip_labels = NULL
for (i in 1:length(chip_levels)) {
  chip_labels[i] = paste0(chip_levels[i], " (N = ", c(table(sample_qc$chip_version))[i], ")")
}

pdf("results/geno_qc_plots.pdf", width = 8, height = 3)
p1 = ggplot(sample_qc) +
  geom_point(aes(x = PC1, y = PC2, colour = Ancestry_info)) +
  theme_bw() + ggtitle("Provided Ancestry") + theme(legend.title = element_blank())
p2 = ggplot(sample_qc) +
  geom_point(aes(x = PC1, y = PC2, colour = factor(Ancestry, levels = ancestry_levels, labels = ancestry_labels))) +
  theme_bw() + ggtitle("Inferred Ancestry") + theme(legend.title = element_blank())
grid.arrange(p1, p2, nrow = 1)

p3 = ggplot(sample_qc) +
  geom_point(aes(x = PC1, y = PC2, colour = factor(chip_version, levels = chip_levels, labels = chip_labels))) +
  theme_bw() + ggtitle("Illumina Omni2.5Exome-8 Version") + theme(legend.title = element_blank())
p4 = ggplot(sample_qc) +
  geom_boxplot(aes(x = Sex, y = F, colour = factor(Sex, levels = sex_levels, labels = sex_labels))) +
  theme_bw() + ggtitle("Sex QC") + theme(legend.title = element_blank())
grid.arrange(p3, p4, nrow = 1)

dev.off()