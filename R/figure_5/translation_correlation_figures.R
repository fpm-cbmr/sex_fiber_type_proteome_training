source(here::here("R/Library.R"))
#load data
long_df <- vroom::vroom(here::here("data/data_log2_long_normalized.csv"))

#load metadata
metadata <- vroom::vroom(here::here("data-raw/Metadata_REX_pooled_fibers.csv")) %>%
    tidyr::unite(col = "group",
                 c(sex, fibertype), sep ="_", remove = FALSE, na.rm = FALSE) %>%
    filter(!subject == "rex08") %>% #rex08 is NA in csa_delta
    filter(!subject == "rex09")

#add csa from metadata to data long form
long_csa_df <- long_df %>%
    left_join(metadata %>% select(sample_id, CSA), by = "sample_id") %>%
    rename(csa = CSA)

#add delta csa to pre data
data_pre <- long_csa_df %>%
    filter(trial == "pre")

data_post <- long_csa_df %>%
    filter(trial == "post")

data_pre_delta <- data_pre %>%
    mutate(csa_delta = data_post$csa - data_pre$csa)

#load translation pre data and find summed abundance
translation_df <- data_pre_delta %>%
    filter(gene %in% c("RPL28", "RPL23A", "RPL15", "RPL18", "EIF4G2", "RPS27", "RPS6KA3", "RPS14", "RPS8", "RPS19", "RPS10", "RPL7A", "EEF1A1", "EIF3L", "RPS26", "RPL24", "EEF1G", "RPS23", "EEF2",
                          "RPL3L", "RPS6KB2", "RPL6", "RPL3", "RPS16", "EIF3A", "RPL13", "RPL13A", "RPL27A", "RPS3", "EEF2K", "RPS20", "RPL10", "RPLP1", "RPL18A", "RPS24", "RPS4X", "EEFSEC", "RPL10A",
                          "EEF1A2", "RPL5", "RPS13", "RPL8", "RPS29", "EIF3H", "RPS15A", "EIF3C", "RPL14", "EIF4B", "RPSA", "RPL12", "RPS3A", "RPL27", "RPS2", "EIF4A2", "RPL22", "EEF1B2", "EIF6",
                          "RPL30", "RPL17", "RPL38", "EIF3F", "RPS11", "MTOR", "RPS5", "RPL23")) %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarize(sum = sum(expression, na.rm = T)) %>%
    dplyr::ungroup() %>%
    mutate(csa_delta = data_pre_delta$csa_delta[match(sample_id, data_pre_delta$sample_id)]) %>%
    filter(!grepl("s08", sample_id))

#plot of correlation between summed abundance and delta csa
translation_df %>%
    ggplot(aes(x = sum, y = csa_delta)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    ggpubr::stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
    theme_bw() +
    xlab("summed abundance") +
    ylab("change in fiber CSA") +
    ggtitle("Proteins involved in translation")


#correlation plot ---------------------------------------------------

correlation_results <- cor.test(
x = translation_df$sum,
y = translation_df$csa_delta,
method = "pearson" )


translation_df |>
    dplyr::inner_join(metadata |>
                          dplyr::select(c(fibertype, sex, sample_id))) |>
    ggplot2::ggplot() +
    ggnewscale::new_scale_colour() +
    ggnewscale::new_scale_fill() +
    geom_point(
        mapping = ggplot2::aes(
            x = sum,
            y = csa_delta,
            color = sex,
            fill = fibertype
        ),
        size = 1.5,
        alpha = 0.75,
        shape = 21,
        stroke = 1
    ) +
    theme_minimal() +
    ggplot2::ylab("Change in fiber CSA") +
    ggplot2::xlab("Summed abundance") +
    ggplot2::theme(
        text = ggplot2::element_text(size = 8),
        legend.key.size = ggplot2::unit(2, units = "mm"),
        legend.spacing.y = ggplot2::unit(2, units = "mm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, -5, 1, -7)
    ) +
    ggplot2::geom_smooth(mapping = ggplot2::aes(
        x = sum,
        y = csa_delta),
    method = "lm",
    color = "black",
    linewidth = 0.5)+
    guides(
        fill = guide_legend(byrow = TRUE),
        color = guide_legend(byrow = TRUE),
        shape = guide_legend(byrow = TRUE)
    ) +
    scale_fill_manual("Fiber type",
                      values = c("#440154FF", "#67CC5CFF"),
                      labels = c("mhc1" = "Type I", "mhc2" = "Type II")) +
    ggplot2::scale_color_manual("Sex",
                                values = c("black", "#FF7518"),
                                labels = c("female" = "Female", "male" = "Male")) +
    ggplot2::annotate(
        geom = "text",
        x = 760,
        y = 4000,
        size = 2.5,
        label = paste("r = ", round(correlation_results$estimate, 3), "\np-value = ", round(correlation_results$p.value, 3)
    ))

ggplot2::ggsave(here::here("figures/figure_5/correlation_translation_CSA.png"),
                units = "mm",
                height = 60,
                width = 60)


