
##TRANSLATION PRE##

#load data long
long_df <- vroom::vroom(here::here("data/data_log2_long_normalized.csv"))

#create mean dataframe
mean_df <- long_df %>%
    group_by(group, trial, fibertype, sex, gene) %>%
    summarise(mean_expression = mean(expression, na.rm = TRUE)) %>%
    rename(protein = gene)

#Add keywords and GO terms
annotations <- vroom::vroom(here::here("data/keywords.csv")) %>%
    dplyr::rename_with(snakecase::to_snake_case) %>%
    dplyr::select(c("gene_names", "keywords", "gene_ontology_biological_process", "gene_ontology_cellular_component", "gene_ontology_molecular_function")) %>%
    dplyr::rename(gobp = gene_ontology_biological_process,
                  gocc = gene_ontology_cellular_component,
                  gomf = gene_ontology_molecular_function,
                  protein = gene_names) %>%
    dplyr::mutate(protein = gsub("\\ .*","", protein)) %>%
    dplyr::mutate(protein = make.names(protein, unique=TRUE), protein)

#Add annotations to mean dataframe
mean_df <- mean_df %>%
    merge(annotations, by="protein", all.x = TRUE)

#write.csv(mean_df,
#here::here("data/data_mean_keywords.csv"))

#filter mean dataframe to only include pre samples and proteins invovled in translation
translation_pre <- mean_df %>%
    filter(trial == "pre") %>%
    filter(protein %in% c("MRPS9", "GSPT1", "MRPS5", "RPL28", "MRPL49", "RCC1L", "MRPL39", "RPL23A", "NARS1", "RPL15", "RPL18", "FASTKD2", "EIF4G2",
                       "MRPS14", "APEH", "RPS27", "RPS6KA3", "RPS14", "RPS8", "DDX1", "RPS19", "MRPL12", "RPS10", "RPL7A", "DARS1", "TRAP1", "EEF1A1", "EIF3L", "LARS1",
                       "RPS26", "RPL24", "C1QBP", "DAP3", "FARSB", "EEF1G", "MRPL13", "TRMT10C", "RPS23", "EEF2", "MIF4GD", "RPL3L", "RPS6KB2", "RPL6", "MCTS1",
                       "SH3BGRL", "RPL3", "PRMT1", "LSM14B", "PAIP1", "ABCE1", "MRPL4", "RPS16", "EIF3A", "RPL13", "RPL13A", "RPL27A", "LARP4", "NIBAN1", "RPS3", "CALR",
                       "EEF2K", "RPS20", "MRPS17", "RPL10", "RPLP1", "RPL18A", "RPS24", "RPS4X", "EEFSEC", "RPL10A", "EEF1A2", "MARS1", "MRPL15", "RPL5", "RPS13",
                       "MTFMT", "RPL8", "RACK1", "RPS29", "MKNK1", "NGRN", "KRT17", "DHPS", "COPS5", "MRPL37", "EIF3H", "RPS15A", "EIF3C", "PRKDC", "RPL14", "PKM",
                       "EIF4B", "RPSA", "RPL12", "RPS3A", "RPL27", "PTCD3", "RPS2", "EIF4A2", "RPL22", "EEF1B2", "EIF6", "IMPACT", "FARSA", "RPL30", "RPL17", "HSPB1",
                       "FAU", "EARS2", "GAPDH", "GARS1", "IARS2", "SARS2", "RPL38", "FXR1", "EIF3F", "RPS11", "MTOR", "MTPN", "HNRNPR", "SARS1", "MRPL47", "MRPL46",
                       "RPS5", "SHMT2", "ACO1", "RPL23", "NPM1", "TNKS1BP1")) %>%
    dplyr::select(c("protein", "mean_expression", "group"))

#create dataframe for heatmap
df_translation <- translation_pre %>%
    tidyr::pivot_wider(
        id_cols = "protein",
        names_from = "group",
        values_from = "mean_expression"
    ) %>%
    column_to_rownames("protein") %>%
    filter_all(all_vars(!is.na(.)))

#Translation heatmap
heatmap_translation <- pheatmap(df_translation,
                    scale = "row",
                    cluster_rows = T,
                    cluster_cols = T,
                    #annotation_col = dplyr::select(metadata, "group"),
                    #annotation_colors=list(sex=c(male="#ce64c5", female="#080808"),
                    #                       fibertype=c(mhc1="#3c1f51", mhc2="#78b85c")),
                    color=colorRampPalette(c("#3173b2", "white", "#d3352f"))(100),
                    border_color = "white",
                    clustering_distance_cols="euclidean",
                    clustering_distance_rows="euclidean",
                    clustering_method = "ward.D2",
                    show_rownames = F,
                    show_colnames = T,
                    fontsize = 5,
                    annotation_legend = F,
                    legend = T,
                    treeheight_row =15,
                    treeheight_col =15,
                    cutree_rows = 3
)


# Filter for ribosomal proteins and eukaryotic initiation and elon --------

#filter mean dataframe to only include pre samples and proteins involved in translation
translation_filtered <- mean_df %>%
    filter(trial == "pre") %>%
    filter(protein %in% c("RPL28", "RPL23A", "RPL15", "RPL18", "EIF4G2", "RPS27", "RPS6KA3", "RPS14", "RPS8", "RPS19", "RPS10", "RPL7A", "EEF1A1", "EIF3L", "RPS26", "RPL24", "EEF1G", "RPS23", "EEF2",
                          "RPL3L", "RPS6KB2", "RPL6", "RPL3", "RPS16", "EIF3A", "RPL13", "RPL13A", "RPL27A", "RPS3", "EEF2K", "RPS20", "RPL10", "RPLP1", "RPL18A", "RPS24", "RPS4X", "EEFSEC", "RPL10A",
                          "EEF1A2", "RPL5", "RPS13", "RPL8", "RPS29", "EIF3H", "RPS15A", "EIF3C", "RPL14", "EIF4B", "RPSA", "RPL12", "RPS3A", "RPL27", "RPS2", "EIF4A2", "RPL22", "EEF1B2", "EIF6",
                          "RPL30", "RPL17", "RPL38", "EIF3F", "RPS11", "MTOR", "RPS5", "RPL23")) %>%
    dplyr::select(c("protein", "mean_expression", "group"))

#create dataframe for filtered heatmap
df_translation_v2 <- translation_filtered %>%
    tidyr::pivot_wider(
        id_cols = "protein",
        names_from = "group",
        values_from = "mean_expression"
    ) %>%
    column_to_rownames("protein") %>%
    filter_all(all_vars(!is.na(.))) %>%
    t() %>%
    as.data.frame()

#Heatmap of selected proteins in translation
heatmap_translation_v2 <- pheatmap(df_translation_v2,
                                scale = "column",
                                cluster_rows = T,
                                cluster_cols = T,
                                #annotation_col = dplyr::select(metadata, "group"),
                                #annotation_colors=list(sex=c(male="#ce64c5", female="#080808"),
                                #                       fibertype=c(mhc1="#3c1f51", mhc2="#78b85c")),
                                color=colorRampPalette(c("#3173b2", "white", "#d3352f"))(100),
                                border_color = "white",
                                clustering_distance_cols="euclidean",
                                clustering_distance_rows="euclidean",
                                clustering_method = "ward.D2",
                                show_rownames = T,
                                show_colnames = T,
                                fontsize = 6,
                                annotation_legend = F,
                                legend = T,
                                treeheight_row =15,
                                treeheight_col =15,
                                cutree_cols = 2
)

#save heatmap
ggsave(plot = heatmap_translation_v2, here::here('figures/figure_5/heatmap_translation.pdf'), height = 60, width = 150, units = "mm")
