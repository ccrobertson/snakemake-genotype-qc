library(yaml)
library(tidyverse)
library(gridExtra)

pheno = read.table("/cluster/ifs/projects/collins/users/robertsoncc/islet_genotype_merge/data/donor_meta_data.txt", header = TRUE, sep = "\t")
vcf_yaml = read_yaml("/cluster/ifs/projects/collins/users/robertsoncc/islet_genotype_merge/config.yaml")

### Build sample to batch map
batches_list =lapply(vcf_yaml$batches, function(x) as.data.frame(x))
for (i in 1:length(batches_list)) {
    batches_list[[i]]$Batch = names(batches_list)[i]
}
batches_df = do.call("rbind", batches_list)
batches_df$Donor = batches_df$samples
batches_df$Batch_Donor = paste0(gsub("batch_", "", batches_df$Batch), ":", batches_df$Donor)
batches_df$IID = NA
batches_df$IID[!duplicated(batches_df$Donor)] <- batches_df$samples[!duplicated(batches_df$Donor)]
batches_df$IID[duplicated(batches_df$Donor)] <- batches_df$Batch_Donor[duplicated(batches_df$Donor)]

samples_df = batches_df[,c("IID", "Donor","Batch", "chip_version", "vcf_file")]

### Combine all the qc info
missing_qc = read.table("results/tmp.smiss", header = FALSE, sep = "")
names(missing_qc) <- c("FID", "IID", "NMISS", "N", "FRQ_MISS")
sex_qc = read.table("results/tmp_XY.sexcheck", header = TRUE, sep = "")
ancestry_qc = read.table("results/tmp_king_InferredAncestry.txt", header = TRUE, sep = "\t")

sample_qc = merge(samples_df,
    merge(ancestry_qc,
        merge(missing_qc, sex_qc, 
            by = "IID"), 
        by = "IID"), 
    by = "IID")

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
write.table(sample_qc, file = "results/sample_qc_summary.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


### Check that duplicates are MZ
duplicated_donors = unique(sample_qc$Donor[duplicated(sample_qc$Donor)])
related_qc = read.table("results/tmp_king.kin", header = TRUE, sep = "\t")
MZ = related_qc[related_qc$InfType == "Dup/MZ",]
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
duplicated_df = sample_qc[sample_qc$Donor %in% duplicated_donors,c("FID","IID", "Donor","FRQ_MISS")]
duplicated_df = duplicated_df[order(duplicated_df$Donor, duplicated_df$FRQ_MISS, decreasing = TRUE), ]
duplicates_to_drop = duplicated_df[duplicated(duplicated_df$Donor), c("FID", "IID")]
duplicates_to_keep = duplicated_df[!duplicated(duplicated_df$Donor), c("FID", "IID", "FID", "Donor")]
write.table(duplicates_to_drop, file = "results/duplicates_to_drop.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(duplicates_to_keep, file = "results/duplicates_to_update.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# ### Sex update <-- Currently we are not excluding for sex errors, do we want to update ped file to reflect inferred sex?
# sample_qc <- sample_qc[!sample_qc$IID %in% duplicates_to_drop$IID,]
# #recode sex variable
# sample_qc <- sample_qc %>%
#  mutate(
#         Sex = case_match(
#             Sex,
#             c("Male", "M", "m") ~ "M",
#             c("Female", "F", "f") ~ "F",
#             c("NA", "nr", "N", "0", " ", "") ~ "0"
#         )
#     )
# write.table(sample_qc[,c("FID", "IID", "SNPSEX")], file = "results/sex_update.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


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

pdf("results/pca_plots.pdf", width = 8, height = 3)
p1 = ggplot(sample_qc) +
    geom_point(aes(x = PC1, y = PC2, colour = Ancestry_info)) +
    theme_bw() + ggtitle("Provided Ancestry") + theme(legend.title= element_blank())
p2 = ggplot(sample_qc) +
    geom_point(aes(x = PC1, y = PC2, colour = factor(Ancestry, levels = ancestry_levels, labels = ancestry_labels))) +
    theme_bw() + ggtitle("Inferred Ancestry") + theme(legend.title= element_blank())
grid.arrange(p1, p2, nrow = 1)

p3 = ggplot(sample_qc) +
    geom_point(aes(x = PC1, y = PC2, colour = factor(chip_version, levels = chip_levels, labels = chip_labels))) +
    theme_bw() + ggtitle("Illumina Omni2.5Exome-8 Version") + theme(legend.title= element_blank())
p4 = ggplot(sample_qc) +
    geom_boxplot(aes(x = Sex, y = F, colour = factor(Sex, levels = sex_levels, labels = sex_labels))) +
    theme_bw() + ggtitle("Sex QC") + theme(legend.title= element_blank())
grid.arrange(p3, p4, nrow = 1)

dev.off()