source(here::here("R/Library.R"))

##PRE TO POST BP##

#open gsea_bp files of pre to post
f_t1_bp <- vroom::vroom(here::here("results/gsea_bp_female_mhc1.csv")) %>%
    mutate(group = "female_mhc1") %>%
    mutate(regulated = case_when(
        enrichmentScore > 0 ~ "Upregulated",
        enrichmentScore < 0 ~ "Downregulated"))

f_t2_bp <- vroom::vroom(here::here("results/gsea_bp_female_mhc2.csv")) %>%
    mutate(group = "female_mhc2") %>%
    mutate(regulated = case_when(
        enrichmentScore > 0 ~ "Upregulated",
        enrichmentScore < 0 ~ "Downregulated"))

m_t1_bp <- vroom::vroom(here::here("results/gsea_bp_male_mhc1.csv")) %>%
    mutate(group = "male_mhc1") %>%
    mutate(regulated = case_when(
        enrichmentScore > 0 ~ "Upregulated",
        enrichmentScore < 0 ~ "Downregulated"))

m_t2_bp <- vroom::vroom(here::here("results/gsea_bp_male_mhc2.csv")) %>%
    mutate(group = "male_mhc2") %>%
    mutate(regulated = case_when(
        enrichmentScore > 0 ~ "Upregulated",
        enrichmentScore < 0 ~ "Downregulated"))

#merge BP data frames
merged_bp <- rbind(f_t1_bp, f_t2_bp, m_t1_bp, m_t2_bp)
merged_bp$regulated <- factor(merged_bp$regulated, levels = c("Upregulated", "Downregulated"))


#Enrichment plot of BP
gsea <- merged_bp %>%
    dplyr::filter(!Description %in% c(
        "skin development",
        "epidermis development",
        "epidermal cell differentiation"
    )) |>
    ggplot(aes(x = group, y = Description, color = qvalue, size = abs(NES))) +
    geom_point() +
    facet_wrap(~regulated) +
    scale_x_discrete(breaks = c("female_mhc1", "female_mhc2", "male_mhc1", "male_mhc2"),
                     labels = c("Female\nType I", "Female\nType II", "Male\nType I", "Male\nType II" )) +
    theme_bw() +
    xlab("") +
    ylab("") +
    labs(color = "FDR", size = "NES") +
    theme(
        text = ggplot2::element_text(size = 7),
        strip.text = ggplot2::element_text(size = 8),
        legend.key.size = ggplot2::unit(4, units = "mm"),
        legend.spacing.y = ggplot2::unit(0.4, units = "mm"))

ggsave(plot = gsea, here::here('figures/figure_4/training_gsea.pdf'), height = 70, width = 140, units = "mm")
