library(tidyverse)
library(gridExtra)
library(optparse)

option_list <- list(
  make_option(
    c("--donor_meta_data"), type = "character", help = " meta "
  ),
  make_option(
    c("--batch_table"), type = "character", help = " batches "
  ),
  make_option(
    c("--smiss"), type = "character", help = " sample missingness "
  ),
  make_option(
    c("--ancestry_qc"), type = "character", help = " ancestry inference by svm"
  ),
  make_option(
    c("--rel_qc"), type = "character", help = " relationship inference by king"
  ),
  make_option(
    c("--outdir"), type = "character", help = " output directory "
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)
print(opts)

### Testing
# opts = list()
# opts$donor_meta_data = "/scratch/scjp_root/scjp0/ccrober/proj-t1d-u01-genotypes/genotype_qc/data/donor_meta_data_2025-04-10.txt"
# opts$batch_table = "/scratch/scjp_root/scjp0/ccrober/proj-t1d-u01-genotypes/genotype_qc/data/genotype_batch_table.txt"
# opts$smiss = "results/filtered_geno_0.05/tmp.smiss"
# opts$ancestry_qc = "results/ancestry_svm/kingsvm_InferredAncestry.txt"
# opts$rel_qc = "results/relationship_inference/kingrel.kin"
# opts$outdir = "results/summary_nosexcheck"

cat("Read meta data...\n")
pheno = read.table(opts$donor_meta_data, header = TRUE, sep = "\t")

cat("Read batches data...\n")
batches_df = read.table(opts$batch_table, header = TRUE)
row.names(batches_df) = batches_df$IID

cat("Read sample missingness data...\n")
missing_qc = read.table(opts$smiss, header = FALSE, sep = "")
names(missing_qc) <- c("FID", "IID", "NMISS", "N", "FRQ_MISS")

cat("Read ancestry inferences data...\n")
ancestry_qc = read.table(opts$ancestry_qc, header = TRUE)

cat("Merge all data...\n")
sample_qc = merge(batches_df[, c("IID", "Donor", "Batch", "chip_version")], merge(ancestry_qc, missing_qc, by = c("FID", "IID")), by = "IID")
sample_qc <- merge(sample_qc, pheno, by = "Donor", all.x = TRUE)

cat("Add dummy variables for SNPSEX and PEDSEX...\n")
sample_qc$SNPSEX = 1 #adding dummy values for downstream analyses
sample_qc$PEDSEX = 1

cat("Harmonize meta data sex variable...\n")
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
cat("Save summary table...\n")
write.table(sample_qc, file = file.path(opts$outdir, "geno_qc_summary.txt"), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


### Check that duplicates are MZ
cat("Check sample duplicates...\n")
duplicated_donors = unique(sample_qc$Donor[duplicated(sample_qc$Donor)])
related_qc = read.table(opts$rel_qc, header = TRUE)
MZ = related_qc[related_qc$InfType == "Dup/MZ", ]

if (length(duplicated_donors)>0) {
for (i in 1:length(duplicated_donors)) {
  if (!duplicated_donors[i] %in% c(MZ$ID1, MZ$ID2)) {
    cat("WARNING: Duplicate donors IDs are not reflected in genotype-based relationship inference.\n")
  }
}
cat("The following donors are duplicated across genotyping batches:\n")
cat(duplicated_donors, sep = "\n")
cat("\n")
cat("Please confirm that relationship inference is consistent:\n")
print(MZ)
cat("\n")
}


### Choose sample among duplicates
cat("Choose sample from each set of duplicate donors to keep...\n")
duplicated_df = sample_qc[sample_qc$Donor %in% duplicated_donors, c("FID", "IID", "Donor", "FRQ_MISS")]
duplicated_df = duplicated_df[order(duplicated_df$Donor, duplicated_df$FRQ_MISS, decreasing = TRUE), ]
duplicates_to_drop = duplicated_df[duplicated(duplicated_df$Donor), c("FID", "IID")]
duplicates_to_keep = duplicated_df[!duplicated(duplicated_df$Donor), c("FID", "IID", "FID", "Donor")]
write.table(duplicates_to_drop, file = file.path(opts$outdir, "duplicates_to_drop.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(duplicates_to_keep, file = file.path(opts$outdir, "duplicates_to_update.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


### Sex update file -- can be used to update ped file to reflect inferred sex
#sample_qc <- sample_qc[!sample_qc$IID %in% duplicates_to_drop$IID, ]
#write.table(sample_qc[, c("FID", "IID", "SNPSEX")], file = file.path(opts$outdir, "sex_update.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


### QC Plots
cat("Generate QC plots...\n")
ancestry_levels = names(c(table(sample_qc$Ancestry)))
ancestry_labels = NULL
for (i in 1:length(ancestry_levels)) {
  ancestry_labels[i] = paste0(ancestry_levels[i], " (N = ", c(table(sample_qc$Ancestry))[i], ")")
}

# sex_levels = names(c(table(sample_qc$Sex)))
# sex_labels = NULL
# for (i in 1:length(sex_levels)) {
#   sex_labels[i] = paste0(sex_levels[i], " (N = ", c(table(sample_qc$Sex))[i], ")")
# }

chip_levels = names(c(table(sample_qc$chip_version)))
chip_labels = NULL
for (i in 1:length(chip_levels)) {
  chip_labels[i] = paste0(chip_levels[i], " (N = ", c(table(sample_qc$chip_version))[i], ")")
}

pdf(file.path(opts$outdir, "geno_qc_plots.pdf"), width = 8, height = 3)

p1 = ggplot(sample_qc) +
  geom_point(aes(x = PC1, y = PC2, colour = Ancestry_info)) +
  theme_bw() + ggtitle("Provided Ancestry") + theme(legend.title = element_blank())
p2 = ggplot(sample_qc) +
  geom_point(aes(x = PC1, y = PC2, colour = factor(Ancestry, levels = ancestry_levels, labels = ancestry_labels))) +
  theme_bw() + ggtitle("Inferred Ancestry") + theme(legend.title = element_blank())
grid.arrange(p1, p2, nrow = 1)

# p3 = ggplot(sample_qc) +
#   geom_point(aes(x = PC1, y = PC2, colour = factor(chip_version, levels = chip_levels, labels = chip_labels))) +
#   theme_bw() + ggtitle("Illumina Omni2.5Exome-8 Version") + theme(legend.title = element_blank())
# p4 = ggplot(sample_qc) +
#   geom_boxplot(aes(x = Sex, y = F, colour = factor(Sex, levels = sex_levels, labels = sex_labels))) +
#   theme_bw() + ggtitle("Sex QC") + theme(legend.title = element_blank())
# grid.arrange(p3, p4, nrow = 1)

p3 = ggplot(sample_qc) +
  geom_point(aes(x = PC1, y = PC2, colour = factor(chip_version, levels = chip_levels, labels = chip_labels))) +
  theme_bw() + ggtitle("Illumina Omni2.5Exome-8 Version") + theme(legend.title = element_blank())
print(p3)

dev.off()
