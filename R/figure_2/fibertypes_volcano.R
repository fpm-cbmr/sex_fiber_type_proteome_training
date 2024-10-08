source(here::here("R/Library.R"))
#FEMALE TYPE I VS TYPE II FIBERS AT BASELINE

#load results of female type I and type II fiber differential expression
female_fibers <- vroom::vroom(here::here("results/differential_expression_female_fibertypes.csv")) %>%
    dplyr::rename("genes" = "...1") |>
    mutate(
        regulated = case_when(
            logFC > 0 & q < 0.05 ~ "type2",
            logFC < 0 & q < 0.05 ~ "type1",
            q > 0.05 ~ "no")) |>
    dplyr::mutate(label = dplyr::case_when(
        genes %in% c(
            "MYBPC2",
            "MYH7",
            "ATP2A2",
            "MYL11",
            "MYH2",
            "PDLIM1",
            "TNNC1",
            "ATP2A1",
            "TPM3",
            "TNNI2",
            "MYH1"
        ) ~ genes,
        TRUE ~ ""
    ))

#too many differentially expressed proteins. Label only top 100 to get a graphical overview.
female_100_de <- female_fibers %>%
    dplyr::filter(!label == "") |>
    dplyr::pull(genes)

#volcano plot of female fiber types
female_fibers %>%
    ggplot(aes(x = logFC, y = `-log10p`, color = regulated)) +
    geom_point(size = 1, alpha = 0.6) +
    scale_color_manual(
        name = NULL,
        values = c("type2" = "#5DC863FF", "type1" = "#440154FF", "no" = "gray"),
        labels = c("type2" = "Type II", "type1" = "Type I", "no" = "Not differentially expressed")
    ) +
    theme_bw() +
    theme(legend.position = "none",
          text = ggplot2::element_text(size = 6)) +
    geom_label_repel(
        data = filter(female_fibers, genes %in% female_100_de),
        aes(label = protein,
            fill = regulated),
        color = "black",
        size = 2,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 20
    ) +
    ggplot2::scale_fill_manual(values = c(
        "#efedf5",
        "#e5f5e0"
    )) +
    xlab("Log2FC (type II - type I)") +
    ylab("-log10(p-value)")


#MALE TYPE I VS TYPE II FIBERS AT BASELINE

#load results of male type I and type II fiber differential expression
male_fibers <- vroom::vroom(here::here("results/differential_expression_male_fibertypes.csv")) %>%
    dplyr::rename("genes" = "...1") |>
    mutate(
        regulated = case_when(
            logFC > 0 & q < 0.05 ~ "type2",
            logFC < 0 & q < 0.05 ~ "type1",
            q > 0.05 ~ "no")) |>
    dplyr::mutate(label = dplyr::case_when(
        genes %in% c(
            "TNNT3",
            "MYH2",
            "TNNI2",
            "ATP2A2",
            "MYH7",
            "TNNC1",
            "MYBPC2",
            "TNNI1",
            "PDLIM1",
            "ATP2A1",
            "TPM3",
            "MYL11",
            "DPYSL3",
            "MYH1"
        ) ~ genes,
        TRUE ~ ""
    ))

male_100_de <- male_fibers %>%
    dplyr::filter(!label == "") |>
    dplyr::pull(genes)

#volcano plot of female fiber types
male_fibers %>%
    ggplot(aes(x = logFC, y = `-log10p`, color = regulated)) +
    geom_point(size = 1, alpha = 0.6) +
    scale_color_manual(
        name = NULL,
        values = c("type2" = "#5DC863FF", "type1" = "#440154FF", "no" = "gray"),
        labels = c("type2" = "Type II", "type1" = "Type I", "no" = "Not differentially expressed")
    ) +
    theme_bw() +
    theme(legend.position = "none",
          text = ggplot2::element_text(size = 6)) +
    geom_label_repel(
        data = filter(male_fibers, genes %in% male_100_de),
        aes(label = protein,
            fill = regulated),
        color = "black",
        size = 2,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 20
    ) +
    ggplot2::scale_fill_manual(values = c(
        "#efedf5",
        "#e5f5e0"
    )) +
    xlab("Log2FC (type II - type I)") +
    ylab("-log10(p-value)")

#add sex as a column in each dataframe
female_fibers$sex <- "female"
male_fibers$sex <- "male"

#combine data frames
merged_df_fibers <- rbind(female_fibers, male_fibers)

#volcano plot of combined data frame of fibertypes at baseline in males and females
merged_volcano <- merged_df_fibers %>%
    ggplot(aes(x = logFC, y = `-log10p`, color = regulated)) +
    geom_point(size = 2, alpha = 0.65, stroke = 0) +
    scale_color_manual(
        name = NULL,
        values = c("type2" = "#67CC5CFF", "type1" = "#440154FF", "no" = "gray"),
        labels = c("type2" = "Greater expression in type II", "type1" = "Greater expression in type I"),
        breaks = c("type2", "type1")
    ) +
    theme_bw() +
    theme(legend.position = "top",
          plot.margin = margin(1, 1, 1, 1),
          legend.box.margin = margin(0, 0, -10, 0),
          text = ggplot2::element_text(size = 6),
          strip.text = ggplot2::element_text(size = 8),
          axis.text.x= element_text(color="black", size = 8),
          axis.text.y= element_text(color="black", size = 8),
          axis.title.x= element_text(color="black", size = 8),
          axis.title.y= element_text(color="black", size = 8),
          legend.text = element_text(size = 6),
          legend.key.size = ggplot2::unit(4, units = "mm")
          ) +
    # ggplot2::scale_y_continuous(expand = c(0,0)) +
    geom_label_repel(
        data = merged_df_fibers |>
            dplyr::filter(!label == ""),
        aes(label = protein,
            fill = regulated),
        color = "black",
        size = 1.75,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 10,
        max.overlaps = Inf
    ) +
    facet_wrap(~ sex, labeller = labeller(sex = c("female" = "Females", "male" = "Males"))) +
    ggplot2::scale_fill_manual(values = c(
        "#efedf5",
        "#e5f5e0"
    )) +
    xlab("Log2fold difference (type II - type I)") +
    ylab("-log10 p-value") +
    guides(fill = FALSE)


ggsave(plot = merged_volcano, here::here('figures/figure_2/figure1_volcano.png'), height = 70, width = 120, units = "mm")
