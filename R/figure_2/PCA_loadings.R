source(here::here("R/Library.R"))

#load log2 transformed normalized data and filter for pre samples
data_pre <- vroom::vroom(here::here("data/data_log2_normalized.csv")) %>%
    rename_with(~ "gene", colnames(.)[1]) %>%
    column_to_rownames("gene") %>%
    dplyr::select(contains("pre"))

#load metadata, filter for pre samples, create groups and make sample_id rownames
metadata_pre <- vroom::vroom(here::here("data-raw/Metadata_REX_pooled_fibers.csv")) %>%
    filter(trial == "pre") %>%
    tidyr::unite(col = "group",
                 c(sex, fibertype), sep ="_", remove = FALSE, na.rm = FALSE) %>%
    tibble::column_to_rownames("sample_id")

#set seed for reproducible imputation
set.seed(666)

#transpose dataframe for pca
df_pca_pre <- data_pre %>%
    t() %>%
    as.data.frame()

#Run PCA analysis
pca_pre <- pca(df_pca_pre, method="ppca", nPcs=2)

#merge pca data and metadata
df_pca_pre <- df_pca_pre %>%
    merge(metadata_pre, by=0, all.x=TRUE) %>%
    dplyr::select(-c("trial", "position", "fibers", "CSA" )) %>%
    column_to_rownames("Row.names")

#merge pca with pca dataframe
df_pca_pre <- merge(df_pca_pre, scores(pca_pre), by=0)


#make PCA plot
df_pca_pre %>%
    mutate(group = factor(group, levels = c("female_mhc1", "female_mhc2", "male_mhc1", "male_mhc2"),
                          labels = c("Female type I", "Female type II", "Male type I", "Male type II"))) %>%
    ggplot(aes(x = PC1, y = PC2, col = group)) +
    geom_point(size = 3) +
    geom_text(aes(label = subject), hjust = -0.2, vjust = 0.5, size = 3) +
    theme_classic() +
    xlab(paste("Principal component 1 (", round(pca_pre@R2[1] * 100, digits = 1), "% of variance)", sep="")) +
    ylab(paste("Principal component 2 (", round(pca_pre@R2[2] * 100, digits = 1), "% of variance)", sep=""))

#very clear clustering of fibertype along PC1 and sex along PC2

#calculate contributions of loadings
loadings_pre <- loadings(pca_pre) %>%
    as.data.frame() %>%
    mutate(pc1_contribution = (PC1^2)*100) %>%
    mutate(pc2_contribution = (PC2^2)*100) %>%
    rownames_to_column() %>%
    dplyr::rename(protein=rowname)

#visualize contribution to PC1
loadings_pre %>%
    dplyr::select("protein", "pc1_contribution") %>%
    arrange(desc(pc1_contribution)) %>%
    slice_head(n = 20) %>%
    ggplot(aes(x=reorder(protein, -pc1_contribution), y=pc1_contribution)) +
    geom_col(color="black", fill="gray")+
    theme_classic() +
    xlab("Protein")+
    ylab("Principal component 1 contribution (%)")

#A lot of myosin and troponin isoforms explains PC1

#visualize contribution to PC2
loadings_pre %>%
    dplyr::select("protein", "pc2_contribution") %>%
    arrange(desc(pc2_contribution)) %>%
    slice_head(n = 20) %>%
    ggplot(aes(x=reorder(protein, -pc2_contribution), y=pc2_contribution)) +
    geom_col(color="black", fill="gray")+
    theme_classic() +
    xlab("Protein")+
    ylab("Principal component 2 contribution (%)")

#A lot of proteins in intermediate filament explains PC2. May play a role in mechanotransduction!

# Principal component loadings plot ---------------------------------------

#visualize contribution to PC1
loadings_pre %>%
    dplyr::select(protein,
                  PC1,
                  pc1_contribution) %>%
    arrange(desc(pc1_contribution)) %>%
    slice_head(n = 20) %>%
    dplyr::mutate(
        protein = factor(protein, levels = rev(c(
            "TNNI1",
            "TNNT1",
            "TNNC1",
            "MYL3",
            "MYOZ2",
            "ATP2A2",
            "TPM3",
            "MYL6B",
            "MYL2",
            "MYH7",
            "MYH6",
            "MPO",
            "DNAH7",
            "MYL11",
            "TNNT3",
            "MYBPC2",
            "MYH1",
            "TNNC2",
            "TNNI2",
            "MYH2"
        )))
    ) |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = PC1,
            y = protein,
            fill = pc1_contribution,
            color = pc1_contribution
        )
    ) +
    ggplot2::geom_col(alpha = 0.75,) +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_viridis_c(
        "Contribution\nto PC1 (%)",
        option = "plasma"
    ) +
    ggplot2::ggtitle("PC1 drivers") +
    ggplot2::geom_vline(xintercept = 0, colour="black") +
    ggplot2::scale_color_viridis_c(
        "Contribution\nto PC1 (%)",
        option = "plasma"
    ) +
    ggplot2::theme(
        legend.key.size = ggplot2::unit(2, "mm"),
        text = ggplot2::element_text(face = "bold",size = 7, colour = "black"),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 8, face="bold"),
        axis.title.y = ggplot2::element_blank()
    )

ggplot2::ggsave(here::here("figures/figure_2/PC1_loadings.png"),
                units = "mm",
                height = 60,
                width = 60)

#visualize contribution to PC2
loadings_pre %>%
    dplyr::select(protein,
                  PC2,
                  pc2_contribution) %>%
    arrange(desc(pc2_contribution)) %>%
    slice_head(n = 20) %>%
    dplyr::mutate(
        protein = factor(protein, levels = rev(c(
            "NDUFC1",
            "ARG1",
            "FLG",
            "AZGP1",
            "KRT10",
            "KRT17",
            "KRT1",
            "PIP",
            "S100A7",
            "JUP",
            "KRT9",
            "KRT5",
            "DSP",
            "DSC1",
            "KRT6B",
            "DCD",
            "KRT14",
            "KRT6A",
            "DSG1",
            "KRT16"

        )))
    ) |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = PC2,
            y = protein,
            fill = pc2_contribution,
            color = pc2_contribution
        )
    ) +
    ggplot2::geom_col(alpha = 0.75,) +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_viridis_c(
        "Contribution\nto PC2 (%)",
        option = "plasma"
    ) +
    ggplot2::ggtitle("PC2 drivers") +
    ggplot2::geom_vline(xintercept = 0, colour="black") +
    ggplot2::scale_color_viridis_c(
        "Contribution\nto PC2 (%)",
        option = "plasma"
    ) +
    ggplot2::theme(
        legend.key.size = ggplot2::unit(2, "mm"),
        text = ggplot2::element_text(face = "bold",size = 7, colour = "black"),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 8, face="bold"),
        axis.title.y = ggplot2::element_blank(),
        legend.position = "bottom",
        legend.box.margin = margin(-10, 0, -10, 0)
    )

ggplot2::ggsave(here::here("figures/figure_3/PC2_loadings.png"),
                units = "mm",
                height = 60,
                width = 40)
