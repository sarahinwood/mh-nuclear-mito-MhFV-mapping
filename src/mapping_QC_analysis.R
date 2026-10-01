#!/usr/bin/env Rscript

#######
# LOG #
#######

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")

#############
# LIBRARIES #
#############

library(data.table)
library(tidyverse)

###########
# GLOBALS #
###########

sample_table_file <- snakemake@input[["sample_table_file"]]

kraken_dir <- snakemake@params[["kraken_dir"]]
bwa_mh_dir <- snakemake@params[["bwa_mh_dir"]]
bwa_mair_dir <- snakemake@params[["bwa_mair_dir"]]
mh_coverage_dir <- snakemake@params[["mh_coverage_dir"]]
bad_samples <- snakemake@input[["bad_samples"]]

########
# MAIN #
########

# read in kraken res and format
kraken_results <- list.files(kraken_dir, "report.txt", full.names = T)
kraken_res <- read_tsv(kraken_results, id="path", col_names=F)
bacteria_fungi_unclassified_taxids <- c(2,4751,0,9605)
kraken_res_main <- data.table(subset(kraken_res, X5 %in% bacteria_fungi_unclassified_taxids))
kraken_res_main$path <- tstrsplit(kraken_res_main$path, "output/04_kraken/reports/kraken_", keep=2)
kraken_res_main$path <- tstrsplit(kraken_res_main$path, "_report", keep=1)
kraken_res_main$classification <- ifelse(grepl("unclassified", kraken_res_main$X6), "kraken_unclassified",
                                         ifelse(grepl("Bacteria", kraken_res_main$X6), "kraken_bacteria",
                                                ifelse(grepl("Fungi", kraken_res_main$X6), "kraken_fungi",
                                                       ifelse(grepl("Homo", kraken_res_main$X6), "kraken_human", "error"))))
kraken_res_main$kraken_pct <- tstrsplit(kraken_res_main$X1, "\t", keep=1)
kraken_res_main <- kraken_res_main[,-c(2)]
kraken_dcast <- dcast(kraken_res_main, path ~ classification)
setnames(kraken_dcast, old="path", new="sample_name")

# merge with sample info
sample_data <- fread(sample_table_file)
sample_info <- sample_data[,c(1,2,3,4,6,7,8)]
sample_data_kraken <- merge(kraken_dcast, sample_info, by="sample_name")

################
## Mh mapping ##
################

## bwa Mh
bwa_mh_results <- list.files(bwa_mh_dir, ".out", full.names = T)
bwa_mh_res <- read_tsv(bwa_mh_results, id="path", col_names=F)
# total reads
bwa_mh_res_reads <- data.table(subset(bwa_mh_res, grepl(" in total", X1), fixed=T))
bwa_mh_res_reads$X1 <- tstrsplit(bwa_mh_res_reads$X1, " +", keep=1)
bwa_mh_res_reads$X1 <- as.numeric(bwa_mh_res_reads$X1)
bwa_mh_res_reads$path <- tstrsplit(bwa_mh_res_reads$path, "/", keep=11)
bwa_mh_res_reads$path <- tstrsplit(bwa_mh_res_reads$path, ".out", fixed=T, keep=1)
setnames(bwa_mh_res_reads, old=c( "path", "X1"), new=c("sample_name", "total_number_reads"))
# mapping
bwa_mh_res_mapping <- data.table(subset(bwa_mh_res, grepl("mapped \\(", X1), fixed=T))
bwa_mh_res_mapping$mapped_reads <-  tstrsplit(bwa_mh_res_mapping$X1, "\\ +", keep=1)
bwa_mh_res_mapping$`Mh_mapping_%` <- tstrsplit(bwa_mh_res_mapping$X1, "\\(", keep=2)
bwa_mh_res_mapping$`Mh_mapping_%` <- tstrsplit(bwa_mh_res_mapping$`Mh_mapping_%`, "\\%", keep=1)
bwa_mh_res_mapping$path <- tstrsplit(bwa_mh_res_mapping$path, "/", keep=11)
bwa_mh_res_mapping$path <- tstrsplit(bwa_mh_res_mapping$path, ".out", fixed=T, keep=1)
setnames(bwa_mh_res_mapping, old=c("path"), new=c("sample_name"))
# reads and mapping
bwa_mh_res_main <- merge(bwa_mh_res_reads, bwa_mh_res_mapping, by="sample_name")
bwa_mh_res_main <- bwa_mh_res_main[,-c("X1")]

## Mh coverage
mh_coverage_results <- list.files(mh_coverage_dir, "coverage.out", full.names = T)
mh_coverage_res <- read_tsv(mh_coverage_results, id="path", col_names=T)
scaffolds <- c("scaffold_1", "scaffold_2", "scaffold_3", "scaffold_4", "scaffold_5", "scaffold_6",
               "scaffold_7", "scaffold_8", "scaffold_9", "scaffold_10", "scaffold_11", "scaffold_12")
mh_coverage_main <- subset(mh_coverage_res, `#rname` %in% scaffolds)

MhFV_coverage <- subset(mh_coverage_res, `#rname`=="MhFV_scaffold")
MhFV_coverage <- data.table(MhFV_coverage)
MhFV_coverage$sample_name <- tstrsplit(MhFV_coverage$path, "/", keep=11, fixed=T)
MhFV_coverage$sample_name <- tstrsplit(MhFV_coverage$sample_name, "_coverage", keep=1, fixed=T)
MhFV_scaff1_cov_depth <- MhFV_coverage[,c(11, 7, 8)]
setnames(MhFV_scaff1_cov_depth, old=c("coverage", "meandepth"), new=c("MhFV_coverage", "MhFV_meandepth"))

mh_scaff1_cov <- subset(mh_coverage_main, `#rname`=="scaffold_1")
mh_scaff1_cov <- data.table(mh_scaff1_cov)
mh_scaff1_cov$sample_name <- tstrsplit(mh_scaff1_cov$path, "/", keep=11, fixed=T)
mh_scaff1_cov$sample_name <- tstrsplit(mh_scaff1_cov$sample_name, "_coverage", keep=1, fixed=T)
mh_scaff1_cov_depth <- mh_scaff1_cov[,c(11, 2, 7, 8)]
setnames(mh_scaff1_cov_depth, c("coverage", "meandepth"), new=c("Mh_scf1_coverage", "Mh_scf1_meandepth"))

## Mh mitochondrial genome mapping & coverage
mh_mito_coverage <- subset(mh_coverage_res, `#rname`=="mh_mh_colemane.1")
mh_mito_coverage <- data.table(mh_mito_coverage)
mh_mito_coverage$sample_name <- tstrsplit(mh_mito_coverage$path, "/", keep=11, fixed=T)
mh_mito_coverage$sample_name <- tstrsplit(mh_mito_coverage$sample_name, "_coverage", keep=1, fixed=T)
Mh_mito_cov_depth <- mh_mito_coverage[,c(11, 7, 8)]
setnames(Mh_mito_cov_depth, old=c("coverage", "meandepth"), new=c("Mh_mito_coverage", "Mh_mito_meandepth"))
mh_mhfv_cov_depth <- merge(mh_scaff1_cov_depth, MhFV_scaff1_cov_depth, by="sample_name")
mh_mito_mhfv_cov_depth <- merge(mh_mhfv_cov_depth, Mh_mito_cov_depth, by="sample_name")
mh_mito_mhfv_cov_depth <- mh_mito_mhfv_cov_depth[,-c("#rname")]

##################
## MaIR mapping ##
##################

## bwa MaIR
bwa_mair_results <- list.files(bwa_mair_dir, ".out", full.names = T)
bwa_mair_res <- read_tsv(bwa_mair_results, id="path", col_names=F)
bwa_mair_res_main <- data.table(subset(bwa_mair_res, grepl("mapped \\(", X1), fixed=T))
bwa_mair_res_main$X1 <- tstrsplit(bwa_mair_res_main$X1, "\\(", keep=2)
bwa_mair_res_main$X1 <- tstrsplit(bwa_mair_res_main$X1, "\\%", keep=1)
bwa_mair_res_main$path <- tstrsplit(bwa_mair_res_main$path, "/", keep=11)
bwa_mair_res_main$path <- tstrsplit(bwa_mair_res_main$path, ".out", fixed=T, keep=1)
setnames(bwa_mair_res_main, old=c("path", "X1"), new=c("sample_name", "MaIR_mapping_%"))

#################
### merge all ###
#################

all_mh_bwa <- merge(bwa_mh_res_main, mh_mito_mhfv_cov_depth, by="sample_name")
all_mapping <- merge(all_mh_bwa, bwa_mair_res_main, by="sample_name")
all_res <- merge(sample_data_kraken, all_mapping, by="sample_name")
all_res$`kraken_bacteria` <- as.numeric(all_res$`kraken_bacteria`)
all_res$`kraken_unclassified` <- as.numeric(all_res$`kraken_unclassified`)
all_res$`Mh_mapping_%` <- as.numeric(all_res$`Mh_mapping_%`)
all_res$`MaIR_mapping_%` <- as.numeric(all_res$`MaIR_mapping_%`)
all_res$mapped_reads <- as.numeric(all_res$mapped_reads)
all_res$MaIR_Mhyp_diff <- all_res$`MaIR_mapping_%`-all_res$`Mh_mapping_%`


# automated identification of bad samples
all_res$concern_label <- ifelse(all_res$kraken_bacteria>10, "a) bacterial contamination above 10%",
                                ifelse(all_res$kraken_unclassified<70, "a) other contamination (unclassified below 70%)",
                                      ifelse(all_res$Species_labelled=="M. aethiopoides Moroccan", "c) Maeth",
                                              ifelse(all_res$Mh_scf1_coverage < 70, "b) Mh scaffold 1 coverage below 70%", # less than 2 mil
                                                     ifelse(all_res$Mh_scf1_meandepth < 9, "b) Mh scaffold 1 depth below 9x",
                                                            ifelse(all_res$mapped_reads <2000000, "b) less than 2mil reads total",
                                                                  ifelse(all_res$`MaIR_mapping_%`>all_res$`Mh_mapping_%`, "d) Maeth mapping better than Mhyp", "e) none")))))))

setorder(all_res, concern_label)

fwrite(all_res, snakemake@output[["mapping_qc_res"]])

# write list of samples with concerns to data folder
exclude_samples_table <- subset(all_res, grepl("a)|b)|d)", all_res$concern_label))
fwrite(list(unique(exclude_samples_table$sample_name)), snakemake@output[["concern_sample_ids"]])

# write log
sessionInfo()