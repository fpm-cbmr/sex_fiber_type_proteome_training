source(here::here("R/Library.R"))

#load log2 transformed normalized data and filter for pre samples
data_pre <- vroom::vroom(here::here("data/data_log2_normalized.csv")) %>%
    rename_with(~ "gene", colnames(.)[1]) %>%
    column_to_rownames("gene") %>%
    dplyr::select(contains("pre"))

#load metadata, filter for pre samples, create groups and make sample_id rownames
metadata_pre <- vroom::vroom(here::here("data-raw/Metadata_REX_pooled_fibers.csv")) %>%
    filter(!subject == "rex09") %>%
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

pca <- df_pca_pre %>%
    mutate(group = factor(
        group,
        levels = c("female_mhc1", "female_mhc2", "male_mhc1", "male_mhc2"),
        labels = c("Female type I", "Female type II", "Male type I", "Male type II")
    )) %>%
    dplyr::mutate(sex = dplyr::case_when(
        sex == "male" ~ "Male",
        TRUE ~ "Female"
    ),
    fibertype = dplyr::case_when(
        fibertype == "mhc1" ~ "Type I",
        TRUE ~ "Type II"
    )) |>
    ggplot() +
    ggnewscale::new_scale_colour() +
    ggnewscale::new_scale_fill() +
    geom_point(
        mapping = ggplot2::aes(
            x = PC1,
            y = PC2,
            color = sex,
            fill = fibertype
        ),
        size = 2,
        alpha = 0.75,
        shape = 21,
        stroke = 0.75
    ) +
    # geom_text(aes(label = subject), hjust = -0.2, vjust = 0.5, size = 2) +
    theme_minimal() +
    xlab(paste(
        "PC1 (",
        round(pca_pre@R2[1] * 100, digits = 1),
        "%)",
        sep = ""
    )) +
    ylab(paste(
        "PC2 (",
        round(pca_pre@R2[2] * 100, digits = 1),
        "%)",
        sep = ""
    )) +
    ggplot2::theme(
        text = ggplot2::element_text(size = 8),
        legend.key.size = ggplot2::unit(2, units = "mm"),
        legend.spacing.y = ggplot2::unit(2, units = "mm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, -5, 1, -7)
    ) +
    guides(
        fill = guide_legend(byrow = TRUE),
        color = guide_legend(byrow = TRUE),
        shape = guide_legend(byrow = TRUE)
    ) +
    scale_fill_manual("Fiber type",
                      values = c("#440154FF",
                                 "#67CC5CFF")) +
    ggplot2::scale_color_manual("Sex",
                                values = c("black",
                                           "#FF7518"))

ggsave(plot = pca, here::here('figures/figure_1/pca.pdf'), height = 60, width = 90, units = "mm")


# pca post ----------------------------------------------------------------

#load log2 transformed normalized data and filter for pre samples
data_post <- vroom::vroom(here::here("data/data_log2_normalized.csv")) %>%
    rename_with(~ "gene", colnames(.)[1]) %>%
    column_to_rownames("gene") %>%
    dplyr::select(contains("post"))

#load metadata, filter for post samples, create groups and make sample_id rownames
metadata_post <- vroom::vroom(here::here("data-raw/Metadata_REX_pooled_fibers.csv")) %>%
    filter(!subject == "rex09") %>%
    filter(trial == "post") %>%
    tidyr::unite(col = "group",
                 c(sex, fibertype), sep ="_", remove = FALSE, na.rm = FALSE) %>%
    tibble::column_to_rownames("sample_id")

#set seed for reproducible imputation
set.seed(666)

#transpose dataframe for pca
df_pca_post <- data_post %>%
    t() %>%
    as.data.frame()

#Run PCA analysis
pca_post <- pca(df_pca_post, method="ppca", nPcs=2)

#merge pca data and metadata
df_pca_post <- df_pca_post %>%
    merge(metadata_post, by=0, all.x=TRUE) %>%
    dplyr::select(-c("trial", "position", "fibers", "CSA" )) %>%
    column_to_rownames("Row.names")

#merge pca with pca dataframe
df_pca_post <- merge(df_pca_post, scores(pca_post), by=0)

#make PCA plot

pca <- df_pca_post %>%
    mutate(group = factor(
        group,
        levels = c("female_mhc1", "female_mhc2", "male_mhc1", "male_mhc2"),
        labels = c("Female type I", "Female type II", "Male type I", "Male type II")
    )) %>%
    dplyr::mutate(sex = dplyr::case_when(
        sex == "male" ~ "Male",
        TRUE ~ "Female"
    ),
    fibertype = dplyr::case_when(
        fibertype == "mhc1" ~ "Type I",
        TRUE ~ "Type II"
    )) |>
    ggplot() +
    ggnewscale::new_scale_colour() +
    ggnewscale::new_scale_fill() +
    geom_point(
        mapping = ggplot2::aes(
            x = PC1,
            y = PC2,
            color = sex,
            fill = fibertype
        ),
        size = 2,
        alpha = 0.75,
        shape = 21,
        stroke = 0.75
    ) +
    # geom_text(aes(label = subject), hjust = -0.2, vjust = 0.5, size = 2) +
    theme_minimal() +
    xlab(paste(
        "PC1 (",
        round(pca_post@R2[1] * 100, digits = 1),
        "%)",
        sep = ""
    )) +
    ylab(paste(
        "PC2 (",
        round(pca_post@R2[2] * 100, digits = 1),
        "%)",
        sep = ""
    )) +
    ggplot2::theme(
        text = ggplot2::element_text(size = 8),
        legend.key.size = ggplot2::unit(2, units = "mm"),
        legend.spacing.y = ggplot2::unit(2, units = "mm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, -5, 1, -7)
    ) +
    guides(
        fill = guide_legend(byrow = TRUE),
        color = guide_legend(byrow = TRUE),
        shape = guide_legend(byrow = TRUE)
    ) +
    scale_fill_manual("Fiber type",
                      values = c("#440154FF",
                                 "#67CC5CFF")) +
    ggplot2::scale_color_manual("Sex",
                                values = c("black",
                                           "#FF7518"))

ggsave(plot = pca, here::here('figures/figure_1/pca.pdf'), height = 60, width = 60, units = "mm")
