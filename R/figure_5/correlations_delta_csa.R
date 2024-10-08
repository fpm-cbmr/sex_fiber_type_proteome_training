source(here::here("R/Library.R"))

#load metadata
metadata <- vroom::vroom(here::here("data-raw/Metadata_REX_pooled_fibers.csv")) %>%
    tidyr::unite(col = "group",
                 c(sex, fibertype), sep ="_", remove = FALSE, na.rm = FALSE) %>%
    filter(!subject == "rex08") %>% #rex08 is NA in csa_delta
    filter(!subject == "rex09")

# Filter metadata for pre and post samples
metadata_pre <- metadata %>%
    filter(trial == "pre") %>%
    select(!fibers) #want csa as only numeric variabel

metadata_post <- metadata %>%
    filter(trial == "post") %>%
    select(!fibers)

#use delta csa in post metadata
metadata_delta_csa <- metadata_post %>%
    mutate(csa_delta = metadata_post$CSA - metadata_pre$CSA) %>%
    select(!CSA) %>% #not interested in post csa
column_to_rownames("sample_id")

##CORRELATION OF EXPRESSION PRE AND DELTA CSA##

#load normalized data and filter for pre samples
data_pre <- vroom::vroom(here::here("data/data_log2_normalized.csv")) %>%
    select(!contains("s08")) %>%  #rex08 is NA in csa_delta
    rename_with(~ "gene", colnames(.)[1]) %>%
    column_to_rownames("gene") %>%
    select(contains("pre"))

#Combine pre data and delta csa
data_pre_csa <- data_pre %>%
    t() %>%
    as.data.frame() %>%
    mutate(csa_delta = metadata_delta_csa$csa_delta)

# Extract the protein columns (excluding 'csa_delta')
protein_columns <- data_pre_csa[, !colnames(data_pre_csa) %in% c("csa_delta")]

# Initialize an empty list to store correlation results
cor_results_list <- list()

# Iterate over each protein column
for (protein_column in colnames(protein_columns)) {
    # Extract the protein expression data and csa_delta
    expression_data <- data_pre_csa[[protein_column]]
    csa_delta <- data_pre_csa$csa_delta

    # Perform correlation analysis
    cor_result <- cor.test(expression_data, csa_delta)

    # Store results in a data frame
    correlation_df_pre <- data.frame(
        protein = protein_column,
        cor = cor_result$estimate,
        p_value = cor_result$p.value
    )

    # Append results to the list
    cor_results_list[[protein_column]] <- correlation_df_pre
}

# Combine results into a data frame
correlations_df_pre <- do.call(rbind, cor_results_list)

# Calculate P-values
p_values <- correlations_df_pre$p_value

write.csv(correlations_df_pre, here::here("results/correlations_pre_csa.csv"))
