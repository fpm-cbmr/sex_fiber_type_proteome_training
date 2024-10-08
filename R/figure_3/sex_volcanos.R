source(here::here("R/Library.R"))

#MALE VS FEMALE TYPE I FIBERS AT BASELINE##

#load results of differential expression between male and female type I fibers
type1_fibers <- vroom::vroom(here::here("results/differential_expression_mhc1_sex.csv")) %>%
    dplyr::rename("genes" = 1) |>
    mutate(
        regulated = case_when(
            logFC > 0 & q < 0.05 ~ "male",
            logFC < 0 & q < 0.05 ~ "female",
            q > 0.05 ~ "no")) |>
    dplyr::mutate(label = dplyr::case_when(
        genes %in% c(
            "TNNI2",
            "TNNT3",
            "GBE1",
            "LAP3",
            "ALDH1A1",
            "NAMPT",
            "FTH1",
            "FTL"
        ) ~ genes,
        TRUE ~ ""
    ))

type1_fibers_de <- type1_fibers %>%
    dplyr::filter(!label == "") |>
    dplyr::pull(genes)

#volcano plot of type 1 fibers
type1_fibers %>%
    ggplot(aes(x = logFC, y = `-log10p`, color = regulated)) +
    geom_point(size = 1, alpha = 0.6) +
    scale_color_manual(
        name = NULL,
        values = c("female" = "#000000", "male" = "#d662c4", "no" = "gray"),
        labels = c("female" = "Female", "male" = "Male", "no" = "Not differentially expressed")
    ) +
    theme_bw() +
    theme(legend.position = "none",
          text = ggplot2::element_text(size = 6)) +
    geom_label_repel(
        data = filter(type1_fibers, genes %in% type1_fibers_de),
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
        "white",
        "#f3d0ed"
    )) +
    xlab("Log2FC (Males - Females)") +
    ylab("-log10(p-value)")

#MALE VS FEMALE TYPE II FIBERS AT BASELINE

#load results of differential expression between male and female type II fibers

type2_fibers <- vroom::vroom(here::here("results/differential_expression_mhc2_sex.csv")) %>%
    dplyr::rename("genes" = 1) |>
    mutate(
        regulated = case_when(
            logFC > 0 & q < 0.05 ~ "male",
            logFC < 0 & q < 0.05 ~ "female",
            q > 0.05 ~ "no")) |>
    dplyr::mutate(label = dplyr::case_when(
        genes %in% c(
            "LAP3",
            "ALDH1A1",
            "GPD1",
            # "ETHE1",
            # "VCL",
            # "CYB5R1",
            # "PRR33",
            # "OLA1",
            "FKBP3",
            # "NAMPT",
            "FTL",
            "HACD1",
            # "CP",
            # "SERPINA1",
            # "LAMB1",
            "DCD",
            "IGHM",
            # "DEFA1B",
            # "UFL1",
            #"GPX3"
            "TNNT3",
            # "GAPDH",
            # "GBE1",
            "PGK1"
            # "RPS17",
            # "MRPL4"
        ) ~ genes,
        TRUE ~ ""
    ))

type2_fibers_de <- type2_fibers %>%
    dplyr::filter(!label == "") |>
    dplyr::pull(genes)

#volcano plot of type 2 fibers
type2_fibers %>%
    ggplot(aes(x = logFC, y = `-log10p`, color = regulated)) +
    geom_point(size = 1, alpha = 0.6) +
    scale_color_manual(
        name = NULL,
        values = c("female" = "#000000", "male" = "#FF7518", "no" = "gray"),
        labels = c("female" = "Female", "male" = "Male", "no" = "Not differentially expressed")
    ) +
    theme_bw() +
    theme(legend.position = "none",
          text = ggplot2::element_text(size = 6)) +
    geom_label_repel(
        data = filter(type2_fibers, genes %in% type2_fibers_de),
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
        "white",
        "#f3d0ed"
    )) +
    xlab("Log2FC (Males - Females)") +
    ylab("-log10(p-value)")

#add fibertype as a column in each dataframe
type1_fibers$fibertype <- "type1"
type2_fibers$fibertype <- "type2"

#combine data frames
merged_df_sex <- rbind(type1_fibers, type2_fibers)

#volcano plot of fiber type differences between sexes
merged_sex <- merged_df_sex %>%
    ggplot(aes(x = logFC, y = `-log10p`, color = regulated)) +
    geom_point(size = 2, alpha = 0.6, stroke = 0) +
    scale_color_manual(
        name = NULL,
        values = c("male" = "#FF7518", "female" = "#000000", "no" = "gray"),
        labels = c("male" = "Greater expression in males", "female" = "Greater expression in females", "no" = "Not differentially expressed")
    ) +
    theme_bw() +
    theme(legend.position = "none",
          text = ggplot2::element_text(size = 6),
          strip.text = ggplot2::element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8)) +
    geom_label_repel(
        aes(label = label,
            fill = regulated),
        color = "black",
        size = 1.75,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        max.overlaps = Inf,
        force = 10
    ) +
    ggplot2::scale_fill_manual(values = c(
        "white",
        "#ffe3d1",
        "gray"
    )) +
    facet_wrap(~ fibertype, labeller = labeller(fibertype = c("type1" = "Type I", "type2" = "Type II"))) +
    xlab("Log2fold difference (male - female)") +
    ylab("-log10 p-value")

#ggsave(plot = merged_sex, here::here('R_figures/manuscript_figures/fig2_sex_volcano.pdf'), height = 60, width = 120, units = "mm")


#Main effect of sex
results_male_female <- vroom::vroom(here::here("results/sex_main_effect.csv")) |>
    mutate(
        regulated = case_when(
            logFC > 0 & q < 0.05 ~ "male",
            logFC < 0 & q < 0.05 ~ "female",
            q > 0.05 ~ "no")) |>
    dplyr::mutate(label = dplyr::case_when(
        protein %in% c(
            "FTL",
            "ALDH1A1",
            "FTH1",
            "LAP3",
            "OVCA2",
            # "VCL",
            # "EEF1G",
            # "GBE1",
            # "ATP1B1",
            "SAA1",
            "PIP",
            "DCD",
            "IGHM",
            "AZGP1",
            "CP",
            # "SERPINA1",
            # "XPO1",
            # "PUDP",
            # "DEFA1B",
            # "LAMB1",
            # "EIF3G",
            # "RPS17",
            "MRPL43",
            "COX7C"
            # "LAMA2",
            # "SERPINA6"
            #"SERPINA7"
            #"SOD3"
            #"TF"
        ) ~ protein,
        TRUE ~ ""
    ))

#main effect of sex volcano
main_sex <- results_male_female |>
    ggplot(aes(x = logFC,
               y = -log10(P.Value),
               color = regulated)) +
    geom_point(size = 2, alpha = 0.6, stroke = 0) +
    scale_color_manual(
        name = NULL,
        values = c("male" = "#FF7518", "female" = "#000000", "no" = "gray"),
        labels = c("male" = "Greater expression in males", "female" = "Greater expression in females", "no" = "")
    ) +
    theme_bw() +
    theme(legend.position = "top",
          legend.box.margin = margin(0, 0, -10, 0),
          text = ggplot2::element_text(size = 6),
          strip.text = ggplot2::element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size =8),
          legend.text = element_text(size = 6),
          legend.key.size = ggplot2::unit(4, units = "mm")
          )  +
    geom_label_repel(
        aes(label = label,
            fill = regulated),
        color = "black",
        size = 1.5,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        max.overlaps = Inf,
        force = 20
    ) +
    ggplot2::scale_fill_manual(values = c(
        "white",
        "#ffe3d1",
        "gray"
    )) +
    xlab("Log2fold difference (male - female)") +
    ylab("-log10 p-value") +
    guides(fill = FALSE)

#ggsave(plot = main_sex, here::here('R_figures/manuscript_figures/fig2_main_volcano.pdf'), height = 70, width = 70, units = "mm")


# combined facet ----------------------------------------------------------

tmp <- merged_df_sex |>
    dplyr::select(!c(`regulated_xiao`,
                     `regulated_qiao`,
                     `regulated_q`))

facet_df <- results_male_female |>
    dplyr::mutate(fibertype = "Main effect") |>
    dplyr::rename(genes = 1) |>
    # dplyr::select(colnames(tmp)) |>
    rbind(tmp)

facet_df %>%
    ggplot(aes(x = logFC, y = `-log10p`, color = regulated)) +
    geom_point(size = 2, alpha = 0.6, stroke = 0) +
    scale_color_manual(
        name = NULL,
        values = c("male" = "#FF7518", "female" = "#000000", "no" = "gray"),
        labels = c("male" = "Greater expression in males", "female" = "Greater expression in females"),
        breaks = c("female", "male")
    ) +
    theme_bw() +
    theme(legend.position = "top",
          text = ggplot2::element_text(size = 6),
          plot.margin = ggplot2::margin(0, 0, 0, 0),
          legend.box.margin = ggplot2::margin(-10, 0, -10, 0),
          strip.text = ggplot2::element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8)) +
    geom_label_repel(
        aes(label = label,
            fill = regulated),
        color = "black",
        size = 1.75,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        max.overlaps = Inf,
        force = 10
    ) +
    ggplot2::scale_fill_manual(values = c(
        "white",
        "#ffe3d1",
        "gray"
    ),
    guide = "none") +
    facet_wrap(~ fibertype, labeller = labeller(fibertype = c("Main effect" = "Main effect", "type1" = "Type I", "type2" = "Type II"))) +
    xlab("Log2fold difference (male - female)") +
    ylab("-log10 p-value")

ggsave(here::here('figures/figure_3/sex_volcano.png'), height = 60, width = 150, units = "mm")
