source(here::here("R/Library.R"))

#load differential expression between sex with annotations
sex_df <- vroom::vroom(here::here("data/results_sex_keywords.csv"))

#load data long
data_long_pre <- vroom::vroom(here::here("data/data_log2_long_normalized.csv")) %>%
    filter(trial == "pre")


#create dataframe of mean expression
data_mean_pre <- data_long_pre %>%
    group_by(gene, group, fibertype, sex) %>%
    summarise(mean_expression = mean(expression, na.rm = TRUE)) %>%
    ungroup()

#filter for troponin isoforms in data_long
tnn_long <- data_long_pre %>%
    select(!(1)) %>%
    filter(gene %in% c("TNNT3", "TNNC2", "TNNI2"))

#boxplot of troponin isoforms
tnn_plot <- tnn_long %>%
    ggplot(aes(x = gene, y = expression, fill = sex)) +
    geom_boxplot(width = 1, alpha = 0.5) +
    geom_jitter(position = position_dodge(width = 1),
                aes(color = sex), size = 2.5, alpha = 0.5, stroke = 0) +
    scale_fill_manual(values = c("#000000", "#FF7518")) +
    scale_color_manual(values = c("#000000", "#FF7518")) +
    theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        text = element_text(size = 6),
        axis.text.x = element_text(color = "black", size = 6),
        axis.text.y = element_text(color = "black", size = 6),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    ) +
    scale_y_continuous(limits = c(11, 23)) +
    xlab("Protein") +
    ylab("abundance (a.u.)") +
    ggtitle("Fast troponin isoforms") +
    facet_wrap(
        ~fibertype,
        labeller = as_labeller(c(
            mhc1 = "Type I",
            mhc2 = "Type II"
        )))

ggsave(plot = tnn_plot, here::here('figures/figure_3/tnn_figure.pdf'), height = 70, width = 70, units = "mm")

