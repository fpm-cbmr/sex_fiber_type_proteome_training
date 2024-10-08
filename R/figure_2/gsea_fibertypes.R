source(here::here("R/Library.R"))

#open gsea_bp files of fibertype differences
f_fibertypes_bp <- vroom::vroom(here::here("results/gsea_bp_female_fibertypes.csv")) %>%
    mutate(sex = "female") %>%
    mutate(regulated = case_when(
        enrichmentScore > 0 ~ "mhc2",
        enrichmentScore < 0 ~ "mhc1"))

m_fibertypes_bp <- vroom::vroom(here::here("results/gsea_bp_male_fibertypes.csv")) %>%
    mutate(sex = "male") %>%
    mutate(regulated = case_when(
        enrichmentScore > 0 ~ "mhc2",
        enrichmentScore < 0 ~ "mhc1"))



#merge BP data frames
fibertypes_bp <- rbind(f_fibertypes_bp, m_fibertypes_bp)
fibertypes_bp$regulated <- factor(fibertypes_bp$regulated, levels = c("mhc1", "mhc2"))


#Enrichment plot of BP for fibertype differences
gsea <- fibertypes_bp %>%
    dplyr::filter(Description %in% c(
        "striated muscle cell development",
        "proton transmembrane transport",
        "NADH regeneration",
        "mitochondrial gene expression",
        "cytoplasmic translation",
        "ATP metabolic process",
        "lipid catabolic process",
        "fatty acid metabolic process"
    )) |>
    ggplot(aes(x = sex, y = Description, color = qvalue, size = abs(NES))) +
    geom_point() +
    facet_wrap(~regulated, labeller = labeller(regulated = c("mhc1" = "Type I", "mhc2" = "Type II"))) +
    scale_x_discrete(breaks = c("female", "male"),
                     labels = c("Female", "Male" )) +
    theme_bw() +
    theme(text = ggplot2::element_text(size = 7),
          strip.text = ggplot2::element_text(size = 8),
          legend.key.size = ggplot2::unit(4, units = "mm")) +
    xlab("") +
    ylab("") +
    labs(color = "q-value", size = "NES")

ggsave(plot = gsea, here::here('figures/figure_2/gsea_fibertypes.pdf'), height = 60, width = 120, units = "mm")
