source(here::here("R/Library.R"))

##SEX DIFFERENCES PRE##

#open gsea_bp files of sex differences

mhc1_sex_bp <- vroom::vroom(here::here("results/gsea_bp_mhc1_sex.csv")) %>%
    mutate(fibertype = "mhc1") %>%
    mutate(regulated = case_when(
        enrichmentScore > 0 ~ "male",
        enrichmentScore < 0 ~ "female"))

mhc2_sex_bp <- vroom::vroom(here::here("results/gsea_bp_mhc2_sex.csv")) %>%
    mutate(fibertype = "mhc2") %>%
    mutate(regulated = case_when(
        enrichmentScore > 0 ~ "male",
        enrichmentScore < 0 ~ "female"))



#merge BP data frames
sex_bp <- rbind(mhc1_sex_bp, mhc2_sex_bp)
sex_bp$regulated <- factor(sex_bp$regulated, levels = c("female", "male"))


#Enrichment plot of BP for sex differences
gsea <- sex_bp %>%
    dplyr::filter(Description %in% c(
        "intermediate filament organization",
        "mitochondrial gene expression",
        "energy derivation by oxidation of organic compounds",
        "long-chain fatty acid transport",
        "skeletal muscle contraction",
        "mRNA destabilization",
        "extracellular matrix organization",
        "tissue remodeling",
        "ribose phosphate metabolic process",
        "cell adhesion",
        "extracellular structure organization",
        "proton transmembrane transport",
        "ATP metabolic process",
        "regulation of fatty acid metabolic process",
        "response to external stimulus"
    )) |>
    ggplot(aes(x = fibertype, y = Description, color = qvalue, size = abs(NES))) +
    geom_point() +
    facet_wrap(~regulated, labeller = labeller(regulated = c("female" = "Enriched in females", "male" = "Enriched in males"))) +
    scale_x_discrete(breaks = c("mhc1", "mhc2"),
                     labels = c("Type I", "Type II" )) +
    theme_bw() +
    theme_bw() +
    theme(text = ggplot2::element_text(size = 7),
          strip.text = ggplot2::element_text(size = 8),
          legend.key.size = ggplot2::unit(2, units = "mm"),
          legend.box.margin = margin(0, -10, 0, -10),
          plot.margin = margin(1,5,1,1)) +
    xlab("") +
    ylab("") +
    labs(color = "FDR", size = "NES")

ggsave(plot = gsea, here::here('figures/figure_3/sex_gsea.pdf'), height = 70, width = 120, units = "mm")
