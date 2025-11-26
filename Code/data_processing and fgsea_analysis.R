# Loading all required datasets using getGEO function
# This is publicaly available dataset
geo_data <- getGEO(filename = "GSE64810_series_matrix.txt")
# Excated metadta from this 
metadata <- pData(geo_data)
# Load normalized count data
fpkm_counts <- read.delim("GSE64810_norm_counts_FPKM_GRCh38.p13_NCBI.tsv", stringsAsFactors = FALSE)
write.csv(fpkm_counts, "norm_counts.csv")

# Cleaning metadata for analysis
# Select specific columns by name
metadata_1 <- metadata %>% select(characteristics_ch1.1:characteristics_ch1.11)
metadata_1[metadata_1 == ""] <- NA
metadata_cleaned <- metadata_1 %>%
  # Rename columns using only the first row of the metadata
  rename_with(~ sapply(seq_along(.), function(i) gsub(":.*", "", metadata_1[50, i]))) %>%
  # Clean the column values by removing everything before ":"
  mutate(across(everything(), ~ gsub(".*: ", "", .))) %>%
  # Remove the first row (used as column names)
  slice(-1)
colnames(metadata_cleaned) <- c("Diagnosis", "pmi", "Age_of_Death", "rin", "mrna-seq-reads", "Age_of_Onset", "Duration", "cag", "vonsattel_grade", "striatal_score", "cortical_score")
write.csv(metadata_cleaned, "metadata.csv")


# Cleaning norm daata for analysis
# Drop the unwanted column if it exists (e.g., row index or unnamed column)
norm_counts_1 <- norm_counts %>% select(-starts_with("X"), -starts_with("..."))
# Reshape the data into long format
formatted_counts <- norm_counts_1 %>%
  pivot_longer(
    cols = -GeneID, # Exclude GeneID column from reshaping
    names_to = "samplenames", # Column name for sample names
    values_to = "counts" # Column name for counts
  )

write.csv(formatted_counts, "formatted_norm_counts.csv", row.names = FALSE)


# Performing GSEA with FGSEA
deseq_results <- read.csv("DESeq2_results.csv", stringsAsFactors = FALSE)

# Rank genes based on log2FoldChange, signifying up/down-regulation
ranks <- deseq_results %>%
  filter(!is.na(log2FoldChange)) %>%
  arrange(desc(log2FoldChange)) %>%
  select(symbol, log2FoldChange) %>%
  deframe()  # Converts data frame to named vector

# Download gene sets for the organism (e.g., Homo sapiens)
gene_sets <- msigdbr(species = "Homo sapiens", category = "H")  # Hallmark gene sets
# Convert gene sets into a list format suitable for fgsea
pathways <- split(gene_sets$gene_symbol, gene_sets$gs_name)

# Perform fgsea
fgsea_results <- fgsea(
  pathways = pathways,  # Gene sets
  stats = ranks,  # Ranked gene list
  minSize = 15,  # Minimum size of a gene set to be included in the analysis
  maxSize = 500,  # Maximum size of a gene set
  nperm = 1000  # Number of permutations
)

# Filter significant pathways
fgsea_significant <- fgsea_results %>%
  filter(padj < 0.05) %>%
  arrange(desc(NES)) 

# Save the results to a CSV file
write.csv(fgsea_significant, "fgsea_results.csv", row.names = FALSE)

# Check which columns are lists
is_list_col <- sapply(fgsea_significant, is.list)

# Flatten the list-columns if any
if (any(is_list_col)) {
  fgsea_significant[is_list_col] <- lapply(fgsea_significant[is_list_col], function(col) {
    sapply(col, function(x) paste(x, collapse = ";")) # Collapse lists into strings
  })
}

fgsea_results$leadingEdge <- sapply(fgsea_results$leadingEdge, function(x) paste(x, collapse = ", "))

# Write the data frame to a .csv file
write.csv(fgsea_results, "fgsea_results_cleaned.csv", row.names = FALSE)
