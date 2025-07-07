#!/usr/bin/env Rscript
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
#######################################              first_curation.R in <scellmate first_qc>            #######################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

#---------------Library---------------#
packages <- c("dplyr", "reshape2", "ggplot2", "purrr", "stringr", "data.table", "segmented", "optparse")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

#---------------Command-line options---------------#
option_list <- list(
  make_option(c("--i-featurecount"), type = "character",
              help = "Path to featureCounts filtered table",
              metavar = "file"),

  make_option(c("--i-map"), type = "character",
              help = "Path to mapping stats file",
              metavar = "file"),

  make_option(c("--i-kraken"), type = "character",
              help = "Path to kraken2 stat file",
              metavar = "file"),

  make_option(c("--db"), type = "character",
              help = paste0(
                "Directory with GTDB reference files ",
                "(gtdb_rep_GenomeID.tsv and taxonomy_r220.processed_with_headers.tsv)"
              ),
              metavar = "DIR"),

  make_option(c("-w", "--workdir"), type = "character", default = getwd(),
              help = "Working directory for output files [default: current dir]",
              metavar = "DIR")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

db_dir  <- normalizePath(opt$db)
workdir <- normalizePath(opt$workdir)
#---------------Loading tables---------------#
# feature <- fread("11_featureCounts_v220/featureCounts-chr.filtered.txt")
# mapped_profile <- setNames(
#   read.delim("mapped_stat.txt", sep="|", header=TRUE)[, -2],
#   c("Cell_ID", "total_mapped", "total_mapped_percentage"))
# kraken2_profile <- read.delim('kraken2_stat.txt',header=TRUE,sep='\t')
feature <- fread(opt$`i-featurecount`)

mapped_profile <- setNames(
  read.delim(opt$`i-map`, sep = "|", header = TRUE)[, -2],
  c("Cell_ID", "total_mapped", "total_mapped_percentage")
)

kraken2_profile <- read.delim(opt$`i-kraken`, header = TRUE, sep = '\t')

#---------------Loading database---------------#
mapping_GenomeID <- read.delim(
  file.path(db_dir, "gtdb_rep_GenomeID.tsv"),
  sep = "\t"
)

Genome_database  <- read.delim(
  file.path(db_dir, "taxonomy_r220.processed_with_headers.tsv"),
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)

# mapping_GenomeID <- read.delim('/mnt/md0/wangyanren/database/release220/gtdb_rep_GenomeID.tsv',sep = '\t')
# Genome_database <- read.delim('/mnt/md0/wangyanren/database/release220/taxonomy_r220.processed_with_headers.tsv', sep = '\t', header = TRUE, check.names = FALSE)

################################################################################################################################################
######################################               1. Data clean of marker-gene composition             ######################################
################################################################################################################################################

#---------------pivote in contig---------------#
feature <- feature[,-c(3,4,5,6)]
colnames(feature) <- sub(
  ".*(scDNA_[[:alnum:]_]*_[0-9]{5,}).*",
  "\\1",
  colnames(feature)
)
feature <- data.frame(feature)
system.time({
  library(data.table)
  setDT(feature)
  cols_SAG <- names(feature)[grepl("^scDNA", names(feature))]
  feature_pivoted_contig <- feature[, lapply(.SD, sum, na.rm = TRUE), by = .(Chr), .SDcols = cols_SAG]
  feature_pivoted_contig <- as.data.frame(feature_pivoted_contig)
})


#---------------pivote in GenomeID---------------#
# mapping_GenomeID <- read.delim('/autofs/tong1/wangyanren/database/release220/gtdb_rep_GenomeID.tsv',sep = '\t')
feature_pivoted_contig_with_GenomeID <- merge(feature_pivoted_contig, mapping_GenomeID, by = "Chr", all.x = TRUE)
system.time({
setDT(feature_pivoted_contig_with_GenomeID)
cols_SAG <- names(feature_pivoted_contig_with_GenomeID)[grepl("^scDNA", names(feature_pivoted_contig_with_GenomeID))]
feature_pivoted_GenomeID <- feature_pivoted_contig_with_GenomeID[, lapply(.SD, sum, na.rm = TRUE), by = .(Genome_ID), .SDcols = cols_SAG]
feature_pivoted_GenomeID <- as.data.frame(feature_pivoted_GenomeID)
        })


#---------------remove 0 marker-gene count SAG---------------#
sums <- feature_pivoted_GenomeID %>% summarise(across(starts_with("scDNA"), \(x) sum(x, na.rm = TRUE)))
cols_to_remove <- names(sums)[sums == 0]
if (length(cols_to_remove) > 0) {
  feature_pivoted_GenomeID <- feature_pivoted_GenomeID %>%
    select(-all_of(cols_to_remove))
}


#---------------pivote in Genus---------------#
# Genome_database <- read.delim('/autofs/tong1/wangyanren/database/release220/taxonomy_r220.processed_with_headers.tsv', sep = '\t', header = TRUE, check.names = FALSE)
feature_pivoted_GenomeID_temp <- as.data.frame(feature_pivoted_GenomeID)
subset_Genome_database <- Genome_database[Genome_database$Genome_ID %in% feature_pivoted_GenomeID_temp$Genome_ID, ]
feature3_pivoted_merged_data <- merge(subset_Genome_database, feature_pivoted_GenomeID_temp, by = "Genome_ID")
system.time({
setDT(feature3_pivoted_merged_data)
cols_with_SAG <- names(feature3_pivoted_merged_data)[grepl("^scDNA", names(feature3_pivoted_merged_data))]
feature_pivoted_Genus <- feature3_pivoted_merged_data[, lapply(.SD, sum, na.rm = TRUE), by = .(Genus), .SDcols = cols_with_SAG]
feature_pivoted_Genus <- as.data.frame(feature_pivoted_Genus)
    })
colnames(feature_pivoted_Genus)[1] <- 'taxa'
# write.table(feature_pivoted_Genus, "feature_pivoted_Genus.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


#---------------transfer count to ratio---------------#
ratios <- sweep(feature_pivoted_Genus[, -1], 2, colSums(feature_pivoted_Genus[, -1]), FUN="/")
ratios <- cbind(feature_pivoted_Genus$taxa, ratios)
# write.table(ratios, "ratios-test.tsv", sep = "\t", row.names = FALSE, quote = FALSE)



#---------------Start to rank the Genus in each SAG---------------#
#---------------determine the rank.1 rank.2 and rank.3 taxa---------------#
feature_temp <- as.data.frame(feature_pivoted_Genus)
rownames(feature_temp) <- feature_temp[,1]
feature_temp <- feature_temp[,-1]
feature_temp <- feature_temp[, colSums(feature_temp) != 0]
column_names <- colnames(feature_temp)
char_result <- data.frame(matrix(NA, ncol = length(column_names), nrow = 3), row.names = c("rank.1_taxa", "rank.2_taxa", "rank.3_taxa"))
colnames(char_result) <- column_names
for (col in column_names) {
    current_col <- feature_temp[, col, drop = FALSE]
    sorted_col <- current_col[order(-current_col[, 1]),, drop = FALSE]
    char_result["rank.1_taxa", col] <- if(nrow(sorted_col) > 0) rownames(sorted_col)[1] else NA
    char_result["rank.2_taxa", col] <- if(nrow(sorted_col) > 1) rownames(sorted_col)[2] else NA
    char_result["rank.3_taxa", col] <- if(nrow(sorted_col) > 2) rownames(sorted_col)[3] else NA
}

#---------------give the rank.x taxa with their abundance---------------#
column_names <- colnames(feature_temp)
numeric_result <- data.frame(matrix(0, ncol = length(column_names), nrow = 4))
row_names <- c("rank.1_taxa", "rank.2_taxa", "rank.3_taxa", "Others_taxa")
rownames(numeric_result) <- row_names
colnames(numeric_result) <- column_names
for (col in column_names) {
    current_col <- feature_temp[, col, drop = FALSE]    
    sorted_col <- current_col[order(-current_col[, 1]), , drop = FALSE]
    numeric_result["rank.1_taxa", col] <- if (nrow(sorted_col) >= 1) sorted_col[1, 1] else 0
    numeric_result["rank.2_taxa", col] <- if (nrow(sorted_col) >= 2) sorted_col[2, 1] else 0
    numeric_result["rank.3_taxa", col] <- if (nrow(sorted_col) >= 3) sorted_col[3, 1] else 0
    total_sum   <- sum(feature_temp[, col], na.rm = TRUE)
    top3_sum    <- sum(numeric_result[c("rank.1_taxa", "rank.2_taxa", "rank.3_taxa"), col], na.rm = TRUE)
    numeric_result["Others_taxa", col] <- total_sum - top3_sum
}


# prepare to melt this file
# don't consider rank.3_taxa yet
numeric_result["Others_taxa", ] <- numeric_result["Others_taxa", ] + numeric_result["rank.3_taxa", ]
# remove the rank.3_taxa
numeric_result <- numeric_result[-which(rownames(numeric_result) %in% c("rank.3_taxa")), ]
# clean the data
df <- cbind(rownames(numeric_result) , numeric_result)
df <- df %>% mutate(across(everything(), ~replace(., is.na(.), 0)))
df <- df %>% mutate(across(starts_with("scDNA"), ~ .x / sum(.x)))
# 这里之前不知道为什么要把最后一列删掉？？
colnames(df)[1] <- 'Type'
setDT(df)
melt_df <- melt(df,id.vars = 'Type')
# write.table(melt_df, "melt_df-test.tsv", sep = "\t", row.names = FALSE, quote = FALSE)



################################################################################################################################################
##########################            2. Unknown SAGs check —— kraken2, mapping ratio, marker-gene counts             ##########################
################################################################################################################################################

#----------------------------------------------Segemented regression of mapping ratio---------------------------------------------#
mapped_profile$sample <- "placeholder"
# Initialize empty lists to collect results for each sample
breakpoints_list <- list()
corresponding_values_list <- list()
num_rows_above_corresponding_list <- list()
percentage_num_rows_below_corresponding_list <- list()
plots_list_map <- list()
# Get all unique sample values
unique_samples <- unique(mapped_profile$sample)
setDT(mapped_profile)
# Analyze each sample
for (sample_name in unique_samples) {
  filtered_dt <- mapped_profile[sample == sample_name]
  # Sort by total_mapped_percentage and add a rank column
  filtered_dt <- filtered_dt[order(-total_mapped_percentage)]
  filtered_dt[, rank := .I]  # .I is the data.table row number
  # Select SAG in the 75th to 100th percentile for segmented regression
  start_rank <- ceiling(0.75 * nrow(filtered_dt))
  filtered_dt_limited <- filtered_dt[rank >= start_rank]
  # Fit a simple linear regression as initial model
  fit <- lm(total_mapped_percentage ~ rank, data = filtered_dt_limited)
  # Perform segmented linear regression
  seg_fit <- segmented(fit, seg.Z = ~rank)
  # Extract breakpoint estimates
  breakpoints <- seg_fit$psi[, "Est."]
  breakpoints_int <- round(breakpoints)
  corresponding_values <- filtered_dt_limited[rank %in% breakpoints_int, total_mapped_percentage]
  # Count rows above the maximum breakpoint value and compute percentage
  num_rows_above_corresponding <- filtered_dt[total_mapped_percentage > max(corresponding_values), .N]
  percentage_num_rows_below_corresponding <- num_rows_above_corresponding / nrow(filtered_dt)
  # Store results in the lists
  breakpoints_list[[sample_name]] <- breakpoints
  corresponding_values_list[[sample_name]] <- corresponding_values
  num_rows_above_corresponding_list[[sample_name]] <- num_rows_above_corresponding
  percentage_num_rows_below_corresponding_list[[sample_name]] <- percentage_num_rows_below_corresponding
  # Plot segmented regression results and breakpoints
  p <- ggplot(filtered_dt, aes(x = rank, y = total_mapped_percentage)) +
    geom_point(color = 'blue') +
    geom_line(data = data.frame(rank = seg_fit$rank, total_mapped_percentage = predict(seg_fit)), aes(x = rank, y = total_mapped_percentage), color = 'grey') +
    geom_vline(xintercept = breakpoints, color = "red", linetype = "dashed") +
    geom_point(data = data.frame(rank = breakpoints, total_mapped_percentage = corresponding_values), 
             aes(x = rank, y = total_mapped_percentage), 
             color = "red", shape = 8, size = 3) +
    geom_text(aes(x = breakpoints, y = corresponding_values, label = paste0("(", round(breakpoints, 2), ", ", round(corresponding_values, 2), ")")), 
              color = "red", hjust = 1.2) +
    labs(title = paste("Segmented Regression for", sample_name), x = "Rank", y = "Total Mapped Percentage")
  # Store the plot
  plots_list_map[[sample_name]] <- p
}
# Combine all results into a data.table
results_dt <- data.table(
  sample = unique_samples,
  breakpoints = I(breakpoints_list),
  corresponding_values = I(corresponding_values_list),
  num_rows_above_corresponding = unlist(num_rows_above_corresponding_list),
  percentage_num_rows_below_corresponding = unlist(percentage_num_rows_below_corresponding_list)
)
# Get unique sample values from mapped_profile
sample_source_mapping <- unique(mapped_profile[, .(sample)])
# Merge sample information into results_dt
results_dt <- merge(results_dt, sample_source_mapping, by = "sample", all.x = TRUE)
# Initialize an empty list to store filtered subsets for each sample
filtered_subsets <- list()
# Loop over each row in results_dt
for (i in 1:nrow(results_dt)) {
  sample_name <- results_dt$sample[i]
  cutoff_value <- results_dt$corresponding_values[i]
  # Filter mapped_profile by sample name and cutoff value
  subset_data <- mapped_profile[sample == sample_name & total_mapped_percentage < cutoff_value]
  # Store the filtered subset in the list
  filtered_subsets[[sample_name]] <- subset_data
}
# Combine all filtered subsets into a single data.table
filtered_result <- rbindlist(filtered_subsets)
# Print or inspect the filtered results
variables_mapping_unknown <- filtered_result[sample == "placeholder", unique(Cell_ID)]
# Print
length(variables_mapping_unknown)
print('Unknown-check finish for mapping ratio')

#---------------------------------------------Segemented regression of kraken2 report---------------------------------------------#
setDT(kraken2_profile)
kraken2_profile[, classified_read_count := 100 - unclassified_read_count]
colnames(kraken2_profile)[1] <- 'variable'
kraken2_profile$sample <- 'placeholder'
breakpoints_list <- list()
corresponding_values_list <- list()
num_rows_above_corresponding_list <- list()
percentage_num_rows_below_corresponding_list <- list()
plots_list_kraken <- list()
unique_samples <- unique(kraken2_profile$sample)
setDT(kraken2_profile)
for (sample_name in unique_samples) {
  filtered_dt <- kraken2_profile[sample == sample_name]
  # Sort by classified_read_count and add a rank column
  filtered_dt <- filtered_dt[order(-classified_read_count)]
  filtered_dt[, rank := .I]
  start_rank <- ceiling(0.75 * nrow(filtered_dt))
  filtered_dt_limited <- filtered_dt[rank >= start_rank]
  fit <- lm(classified_read_count ~ rank, data = filtered_dt_limited)
  seg_fit <- segmented(fit, seg.Z = ~rank)
  breakpoints <- seg_fit$psi[, "Est."]
  breakpoints_int <- round(breakpoints)
  corresponding_values <- filtered_dt_limited[rank %in% breakpoints_int, classified_read_count]
  num_rows_above_corresponding <- filtered_dt[classified_read_count > max(corresponding_values), .N]
  percentage_num_rows_below_corresponding <- num_rows_above_corresponding / nrow(filtered_dt)
  breakpoints_list[[sample_name]] <- breakpoints
  corresponding_values_list[[sample_name]] <- corresponding_values
  num_rows_above_corresponding_list[[sample_name]] <- num_rows_above_corresponding
  percentage_num_rows_below_corresponding_list[[sample_name]] <- percentage_num_rows_below_corresponding
  p <- ggplot(filtered_dt, aes(x = rank, y = classified_read_count)) +
    geom_point(color = 'blue') +
    geom_line(data = data.frame(rank = seg_fit$rank, classified_read_count = predict(seg_fit)), aes(x = rank, y = classified_read_count), color = 'grey') +
    geom_vline(xintercept = breakpoints, color = "red", linetype = "dashed") +
geom_point(data = data.frame(rank = breakpoints, total_mapped_percentage = corresponding_values), 
             aes(x = rank, y = total_mapped_percentage), 
             color = "red", shape = 8, size = 3) +    geom_text(aes(x = breakpoints, y = corresponding_values, label = paste0("(", round(breakpoints, 2), ", ", round(corresponding_values, 2), ")")), 
              color = "red", hjust = 1.2) +
    labs(title = paste("Segmented Regression for", sample_name), x = "Rank", y = "Unclassified Read")
    plots_list_kraken[[sample_name]] <- p
}
results_dt <- data.table(
  sample = unique_samples,
  breakpoints = I(breakpoints_list),
  corresponding_values = I(corresponding_values_list),
  num_rows_above_corresponding = unlist(num_rows_above_corresponding_list),
  percentage_num_rows_below_corresponding = unlist(percentage_num_rows_below_corresponding_list)
)
sample_source_kraken2 <- unique(kraken2_profile[, .(sample)])
results_dt <- merge(results_dt, sample_source_kraken2, by = "sample", all.x = TRUE)
filtered_subsets <- list()
for (i in 1:nrow(results_dt)) {
  sample_name <- results_dt$sample[i]
  cutoff_value <- results_dt$corresponding_values[i]
  subset_data <- kraken2_profile[sample == sample_name & classified_read_count < cutoff_value]
    filtered_subsets[[sample_name]] <- subset_data
}
filtered_result <- rbindlist(filtered_subsets)
variables_kraken2_unknown <- filtered_result[sample == "placeholder", unique(variable)]
length(variables_kraken2_unknown)
print('Unknown-check finish for kraken2')


#---------------------------------------------Segemented regression of marker-gene count---------------------------------------------#
union_result <- union(variables_mapping_unknown,variables_kraken2_unknown)
# make marker-gene filtering strict
sums <- feature_pivoted_GenomeID %>% summarise(across(starts_with("scDNA"), \(x) sum(x, na.rm = TRUE)))
t_sums <- as.data.frame(t(sums))
t_sums$ID <- rownames(t_sums)
colnames(t_sums) <- c('SCG_gene_count','ID')
t_sums$SCG_gene_count <- as.numeric(t_sums$SCG_gene_count)
setDT(t_sums)
colnames(t_sums)[2] <- 'variable'
t_sums <- t_sums[!(variable %in% union_result), ]
t_sums$sample <- 'placeholder'
breakpoints_list <- list()
corresponding_values_list <- list()
num_rows_above_corresponding_list <- list()
percentage_num_rows_below_corresponding_list <- list()
plots_list <- list()
unique_samples <- unique(t_sums$sample)
setDT(t_sums)
for (sample_name in unique_samples) {
  filtered_dt <- t_sums[sample == sample_name]
  # Sort by SCG_gene_count and add a rank column
  filtered_dt <- filtered_dt[order(-SCG_gene_count)]
  filtered_dt[, rank := .I]  # .I 是 data.table 的行号
  start_rank <- ceiling(0.75 * nrow(filtered_dt))
  filtered_dt_limited <- filtered_dt[rank >= start_rank]
  fit <- lm(SCG_gene_count ~ rank, data = filtered_dt_limited)
  seg_fit <- segmented(fit, seg.Z = ~rank)
  breakpoints <- seg_fit$psi[, "Est."]
  breakpoints_int <- round(breakpoints)
  corresponding_values <- filtered_dt_limited[rank %in% breakpoints_int, SCG_gene_count]
  num_rows_above_corresponding <- filtered_dt[SCG_gene_count > max(corresponding_values), .N]
  percentage_num_rows_below_corresponding <- num_rows_above_corresponding / nrow(filtered_dt)
  breakpoints_list[[sample_name]] <- breakpoints
  corresponding_values_list[[sample_name]] <- corresponding_values
  num_rows_above_corresponding_list[[sample_name]] <- num_rows_above_corresponding
  percentage_num_rows_below_corresponding_list[[sample_name]] <- percentage_num_rows_below_corresponding
  p <- ggplot(filtered_dt, aes(x = rank, y = SCG_gene_count)) +
    geom_point(color = 'blue') +
    geom_line(data = data.frame(rank = seg_fit$rank, SCG_gene_count = predict(seg_fit)), aes(x = rank, y = SCG_gene_count), color = 'grey') +
    geom_vline(xintercept = breakpoints, color = "red", linetype = "dashed") +
geom_point(data = data.frame(rank = breakpoints, SCG_gene_count = corresponding_values), 
             aes(x = rank, y = SCG_gene_count), 
             color = "red", shape = 8, size = 3) +
    geom_text(aes(x = breakpoints, y = corresponding_values, label = paste0("(", round(breakpoints, 2), ", ", round(corresponding_values, 2), ")")), 
              color = "red", hjust = -0.5) +
    labs(title = paste("Segmented Regression for", sample_name), x = "Rank", y = "SCG gene count")
  plots_list[[sample_name]] <- p
}
results_dt <- data.table(
  sample = unique_samples,
  breakpoints = I(breakpoints_list),
  corresponding_values = I(corresponding_values_list),
  num_rows_above_corresponding = unlist(num_rows_above_corresponding_list),
  percentage_num_rows_below_corresponding = unlist(percentage_num_rows_below_corresponding_list)
)
sample_source_SCG <- unique(t_sums[, .(sample)])
results_dt <- merge(results_dt, sample_source_SCG, by = "sample", all.x = TRUE)
filtered_subsets <- list()
for (i in 1:nrow(results_dt)) {
  sample_name <- results_dt$sample[i]
  cutoff_value <- results_dt$corresponding_values[i]
  subset_data <- t_sums[sample == sample_name & SCG_gene_count < cutoff_value]
  filtered_subsets[[sample_name]] <- subset_data
}
filtered_result <- rbindlist(filtered_subsets)
variables_SCG_unknown <- filtered_result[sample == "placeholder", unique(variable)]
length(variables_SCG_unknown)
print('QC finish by SCG')


#------------------------------------Obtain the union of SAGs filtered by three independent criteria------------------------------------#
union_variables <- union(union(variables_mapping_unknown,variables_kraken2_unknown),variables_SCG_unknown)
length(union_variables)


################################################################################################################################################
#################################            3. first_qc check —— rank.1 taxa based on marker gene             #################################
################################################################################################################################################

#---------------Remove those Unknown SAGs for 1st QC---------------#
setDT(melt_df)
melt_df_filtered <- melt_df[!(variable %in% union_variables), ]
melt_df_filtered$sample <- 'placeholder'
breakpoints_list <- list()
corresponding_values_list <- list()
num_rows_above_corresponding_list <- list()
percentage_num_rows_above_corresponding_list <- list()
plots_list_rank1 <- list()
unique_samples <- unique(melt_df_filtered$sample)
for (sample_name in unique_samples) {
  # Sort by rank.1_taxa and add a rank column
  filtered_dt <- melt_df_filtered[sample == sample_name & Type == "rank.1_taxa"]
  filtered_dt <- filtered_dt[order(-value)]
  filtered_dt[, rank := .I]
  fit <- lm(value ~ rank, data = filtered_dt)
  seg_fit <- segmented(fit, seg.Z = ~rank)
  breakpoints <- seg_fit$psi[, "Est."]
  breakpoints_int <- round(breakpoints)
  corresponding_values <- filtered_dt[rank %in% breakpoints_int, value]
  num_rows_above_corresponding <- filtered_dt[value > max(corresponding_values), .N]
  percentage_num_rows_above_corresponding <- num_rows_above_corresponding / nrow(filtered_dt)
  breakpoints_list[[sample_name]] <- breakpoints
  corresponding_values_list[[sample_name]] <- corresponding_values
  num_rows_above_corresponding_list[[sample_name]] <- num_rows_above_corresponding
  percentage_num_rows_above_corresponding_list[[sample_name]] <- percentage_num_rows_above_corresponding
  p <- ggplot(filtered_dt, aes(x = rank, y = value)) +
    geom_point(color = 'blue') +
    geom_line(data = data.frame(rank = seg_fit$rank, value = predict(seg_fit)), aes(x = rank, y = value), color = 'grey') +
    geom_vline(xintercept = breakpoints, color = "red", linetype = "dashed") +
    geom_point(aes(x = breakpoints, y = corresponding_values), color = "red", shape = 8, size = 3) +
    geom_text(aes(x = breakpoints, y = corresponding_values, label = paste0("(", round(breakpoints, 2), ", ", round(corresponding_values, 2), ")")), 
              color = "red", hjust = -0.5) +
    labs(title = paste("Segmented Regression for", sample_name), x = "Rank", y = "Value")
  plots_list_rank1[[sample_name]] <- p
}
results_dt <- data.table(
  sample = unique_samples,
  breakpoints = I(breakpoints_list),
  corresponding_values = I(corresponding_values_list),
  num_rows_above_corresponding = unlist(num_rows_above_corresponding_list),
  percentage_num_rows_above_corresponding = unlist(percentage_num_rows_above_corresponding_list)
)
sample_source_melt_df <- unique(melt_df_filtered[, .(sample)])
results_dt <- merge(results_dt, sample_source_melt_df, by = "sample", all.x = TRUE)
filtered_subsets <- list()
for (i in 1:nrow(results_dt)) {
  sample_name <- results_dt$sample[i]
  cutoff_value <- results_dt$corresponding_values[i]
  subset_data <- melt_df_filtered[sample == sample_name & Type == "rank.1_taxa" & value > cutoff_value]
  filtered_subsets[[sample_name]] <- subset_data
}
filtered_result <- rbindlist(filtered_subsets)
variables_marker_unknown <- filtered_result[sample == "placeholder", unique(variable)]
length(variables_marker_unknown)
print('QC finish by rank.1_taxa')


# genus-level contaminated SAGs in first-round QC
contaminated_values <- setdiff(unique(melt_df_filtered$variable), filtered_result$variable)
filtered_contaminated_rows <- melt_df_filtered[melt_df_filtered$variable %in% contaminated_values, ]
contaminated_dt <- data.table(variable = contaminated_values)
contaminated_dt$sample <- 'placeholder'
sum(contaminated_dt$sample == "placeholder")
sample_dt <- contaminated_dt[sample == "placeholder", .(variable)]
# fwrite(sample_dt, file = "sample_contaminated.txt", sep = "\t", col.names = FALSE)

# Ready-to-Co-assembly SAG
kraken2_sample_variables <- kraken2_profile[sample == "placeholder", variable]
contaminated_sample_variables <- contaminated_dt[sample == "placeholder", variable]
clean_sample_variables <- setdiff(kraken2_sample_variables, contaminated_sample_variables)
# write.table(clean_sample_variables, file = "clean_sample_variables-test.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


# Unknown SAGs
# writeLines(union_variables, "unknown_id-test.txt")

# genus-level annotation for pass-1st-QC SAGs
temp <- as.data.frame(t(char_result[1,]))
temp$variable <- rownames(temp)
QC_non_mock <- temp[temp$variable %in% filtered_result$variable, ]
length(rownames(QC_non_mock))
write.table(QC_non_mock,
            file = file.path(workdir, "05_first_QC", "QC_non_mock-genus-test.tsv"),
            sep  = "\t",
            row.names = FALSE,
            quote = FALSE)

# output the statistics
rank1_dt <- melt_df[Type == "rank.1_taxa",
                    .(SAG = variable,
                      `rank.1 taxa ratio` = value)]

class_dt <- kraken2_profile[, .(SAG = variable,
                                classified_read_count)]

map_dt   <- mapped_profile[, .(SAG = Cell_ID,
                               total_mapped_percentage)]

marker_dt <- data.table(
  SAG               = colnames(sums),
  marker_gene_count = as.integer(sums[1, ])
)

qc_dt <- Reduce(function(x, y) merge(x, y, by = "SAG", all = TRUE),
                list(rank1_dt, class_dt, map_dt, marker_dt))

genus_dt <- as.data.table(QC_non_mock)[ ,
              .(SAG = variable,
                `Genus annotation` = rank.1_taxa)]

# union of all SAG IDs mentioned anywhere
all_SAGs <- unique(c(qc_dt$SAG,
                     genus_dt$SAG,
                     contaminated_dt$variable,
                     union_variables))

# ensure every SAG appears once, then merge details
qc_full <- merge(data.table(SAG = all_SAGs), qc_dt, by = "SAG", all.x = TRUE)
qc_full <- merge(qc_full, genus_dt, by = "SAG", all.x = TRUE)

# assign QC Status
qc_full[ , Status := fifelse(SAG %in% genus_dt$SAG,             "pass-1st-QC",
                      fifelse(SAG %in% contaminated_dt$variable, "fail-1st-QC",
                              "Unknown"))]

# fill missing Genus annotation
qc_full[is.na(`Genus annotation`), `Genus annotation` := "N/A"]

setcolorder(qc_full,
            c("SAG",
              "Genus annotation",
              "Status",
              "rank.1 taxa ratio",
              "classified_read_count",
              "total_mapped_percentage",
              "marker_gene_count"))

write.table(qc_full,
            file = file.path(workdir,"05_first_QC", "QC_1st.tsv"),
            sep  = "\t",
            row.names = FALSE,
            quote = FALSE)

thres_val <- as.numeric(results_dt$corresponding_values[[i]])  # extract as numeric
writeLines(
  text = format(signif(thres_val, 17), scientific = FALSE, trim = TRUE),
  con  = file.path(workdir, "05_first_QC", "marker_gene_threshold")
)