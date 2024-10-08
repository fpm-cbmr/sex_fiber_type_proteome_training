#FEMALE TYPE I FIBERS PRE TO POST

#load results of female type I fiber differential expression from pre to post
female_mhc1 <- vroom::vroom(here::here("results/differential_expression_female_mhc1.csv")) %>%
    dplyr::rename("genes" = "...1") |>
    mutate(
        regulated = case_when(
            logFC > 0 & q < 0.05 ~ "up",
            logFC < 0 & q < 0.05 ~ "down",
            q > 0.05 ~ "no")) |>
    dplyr::mutate(label = dplyr::case_when(
        genes %in% c(
            "COL3A1",
            "MYH1",
            "SERPINH1",
            "THY1",
            "S100A13",
            "VIM",
            "PTCD3",
            "HSPB6",
            "LIPE"
        )~ genes,
        TRUE ~ ""
    ))

#FEMALE TYPE II FIBERS PRE TO POST

#load results of female type II fiber differential expression from pre to post
female_mhc2 <- vroom::vroom(here::here("results/differential_expression_female_mhc2.csv")) %>%
    dplyr::rename("genes" = "...1") |>
    mutate(
        regulated = case_when(
            logFC > 0 & q < 0.05 ~ "up",
            logFC < 0 & q < 0.05 ~ "down",
            q > 0.05 ~ "no")) |>
    dplyr::mutate(label = dplyr::case_when(
        genes %in% c(
            "HSPB6",
            "S100A13",
            "GBE1",
            "THY1",
            "SERPINH1",
            "CSRP3",
            "FLG2",
            "KRT16",
            "KRT71",
            "ART3",
            "TINAGL1",
            "PHIP",
            "ACBD3",
            "RPL3",
            "RPS6"
        )~ genes,
        TRUE ~ ""
    ))

#MALE TYPE I FIBERS PRE TO POST

#load results of male type I fiber differential expression from pre to post
male_mhc1 <- vroom::vroom(here::here("results/differential_expression_male_mhc1.csv")) %>%
    dplyr::rename("genes" = "...1") |>
    mutate(
        regulated = case_when(
            logFC > 0 & q < 0.05 ~ "up",
            logFC < 0 & q < 0.05 ~ "down",
            q > 0.05 ~ "no")) |>
    dplyr::mutate(label = dplyr::case_when(
        genes %in% c(
            "S100A13",
            "THY1",
            "SERPINH1",
            "BGN",
            "LGALSL",
            "ART3",
            "LAMA4"
        )~ genes,
        TRUE ~ ""
    ))

#MALE TYPE II FIBERS FROM PRE TO POST

#load results of male type II fiber differential expression from pre to post
male_mhc2 <- vroom::vroom(here::here("results/differential_expression_male_mhc2.csv")) %>%
    dplyr::rename("genes" = "...1") |>
    mutate(
        regulated = case_when(
            logFC > 0 & q < 0.05 ~ "up",
            logFC < 0 & q < 0.05 ~ "down",
            q > 0.05 ~ "no")) |>
    dplyr::mutate(label = dplyr::case_when(
        genes %in% c(
            "S100A13",
            "THY1",
            "HSPB6",
            "SERPINH1",
            "CARNS1",
            "CSRP3",
            "DDRGK1",
            "ART3",
            "MYLK2"
        ) ~ genes,
        TRUE ~ ""
    ))

#in each data frame add sex and fibertype column
female_mhc1$fibertype <- "type1"
female_mhc1$sex <- "female"
female_mhc2$fibertype <- "type2"
female_mhc2$sex <- "female"
male_mhc1$fibertype <- "type1"
male_mhc1$sex <- "male"
male_mhc2$fibertype <- "type2"
male_mhc2$sex <- "male"

#combine data frames
merged_df <- rbind(female_mhc1, female_mhc2, male_mhc1, male_mhc2)

#make volcano plot of combined data frame with classifiers
merged <- merged_df %>%
    ggplot(aes(x = logFC, y = `-log10p`, color = regulated)) +
    geom_point(size = 2, alpha = 0.6, stroke = 0) +
    scale_color_manual(
        name = NULL,  # Remove legend title
        values = c("up" = "#b4656f", "down" = "#007ea7", "no" = "gray"),
        labels = c("up" = "Upregulated", "down" = "Downregulated", "no" = "Not regulated")
    ) +
    theme_bw() +
    theme(legend.position = "none",
          text = ggplot2::element_text(size = 6),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          strip.text = element_text(size = 8)) +
    geom_label_repel(
        data = merged_df |>
            dplyr::filter(!label == ""),
        aes(label = label,
            fill = regulated),
        color = "black",
        size = 1.75,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 20,
        max.overlaps = Inf
    ) +
    facet_grid(fibertype ~ sex, labeller = labeller(fibertype = c("type1" = "Type I", "type2" = "Type II"),
                                                    sex = c("female" = "Female", "male" = "Male" ))) +
    ggplot2::scale_fill_manual(values = c(
        "#b3d8e5",
        "#e9d1d4",
        "lightgray"
    )) +
    xlab("Log2fold change (post - pre)") +
    ylab("-log10 p-value")

ggsave(plot = merged, here::here('figures/figure_4/volcano_sex_fibertype.png'), height = 120, width = 120, units = "mm")

#Main effect of training
main_time <- vroom::vroom(here::here("results/trial_main_effect.csv"))

main_time <- main_time %>%
    mutate(
        regulated = case_when(
            logFC > 0 & q < 0.05 ~ "up",
            logFC < 0 & q < 0.05 ~ "down",
            q > 0.05 ~ "no")) |>
    dplyr::mutate(label = dplyr::case_when(
        protein %in% c(
            "THY1",
            "S100A13",
            "KRT71",
            "KRT16",
            "MYBPH",
            "SERPINH1",
            "CAMK2D",
            "GSN",
            "CASP3",
            "CSRP3",
            "ACBD3",
            "ART3",
            "TINAGL1",
            "XIRP1"
        ) ~ protein,
        TRUE ~ ""
    ))

#Volcano of main effect of training
main_volcano <- main_time %>%
    ggplot(aes(x = logFC, y = `-log10p`, color = regulated)) +
    geom_point(size = 2, alpha = 0.6, stroke = 0) +
    scale_color_manual(
        name = NULL,  # Remove legend title
        values = c("up" = "#b4656f", "down" = "#007ea7", "no" = "gray"),
        labels = c("up" = "Upregulated", "down" = "Downregulated"),
        breaks = c("down", "up")
    ) +
    theme_bw() +
    theme(legend.position = "top",
          legend.box.margin = margin(0, 0, -10, 0),
          text = ggplot2::element_text(size = 6),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          legend.key.size = ggplot2::unit(4, units = "mm"),
          legend.text = element_text(size = 6)
          ) +
    geom_label_repel(
        aes(label = label,
            fill = regulated),
        color = "black",
        size = 1.75,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 20,
        max.overlaps = Inf
    ) +
    ggplot2::scale_fill_manual(values = c(
        "#b3d8e5",
        "lightgray",
        "#e9d1d4"
    )) +
    xlab("Log2fold change (post - pre)") +
    ylab("-log10 p-value") +
    guides(fill = FALSE)

ggsave(plot = main_volcano, here::here('figures/figure_4/volcano_main_effect_trial.png'), height = 70, width = 70, units = "mm")
