source(here::here("R/Library.R"))

#load metadata
metadata <- vroom::vroom(here::here("data-raw/Metadata_REX_pooled_fibers.csv")) %>%
    tidyr::unite(col = "group",
                 c(sex, fibertype), sep ="_", remove = FALSE, na.rm = FALSE) %>%
    filter(!subject == "rex08") %>% #rex08 is NA in csa_delta
    filter(!subject == "rex09")

# Filter metadata for pre and post samples
metadata_pre <- metadata %>%
    filter(trial == "pre") %>%
    select(!fibers) #want csa as only numeric variabel

metadata_post <- metadata %>%
    filter(trial == "post") %>%
    select(!fibers)

#use delta csa in post metadata
metadata_delta_csa <- metadata_post %>%
    mutate(csa_delta = metadata_post$CSA - metadata_pre$CSA) %>%
    select(!CSA)

##PRE ABUNDANCE AND DELTA CSA CORRELATIONS##
#load csa and log2fc correlations
correlations_csa_pre <- vroom::vroom(here::here("results/correlations_pre_csa.csv")) %>%
    mutate("log10p" = -log10(.$p_value)) %>%
    mutate(
        correlated = case_when(
            cor > 0 & p_value < 0.05 ~ "positive",
            cor < 0 & p_value < 0.05 ~ "negative",
            p_value > 0.05 ~ "no")) %>%
    dplyr::mutate(label = dplyr::case_when(
        protein %in% c(
            "RPL15",
            "RPL23A",
            "RPL28",
            "SERPINA1",
            "PRKAA2",
            "IGKV2-30",
            "IGHM"
        ) ~ protein,
        TRUE ~ ""
    ))

#plot of all significant correlations
correlations_pre_plot <- correlations_csa_pre %>%
    ggplot(aes(x = cor, y = log10p, color = correlated)) +
    geom_point(size = 2, alpha = 0.6, stroke = 0) +
    geom_hline(yintercept=1.30, linetype="dashed")+
    scale_color_manual(
        name = NULL,  # Remove legend title
        values = c("positive" = "#b4656f", "negative" = "#007ea7", "no" = "gray"),
        labels = c("positive" = "Positive association", "negative" = "Negative association", "no" = "No association")
    ) +
    theme_bw() +
    theme(legend.position = "none",
          text = ggplot2::element_text(size = 6),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          strip.text = element_text(size = 8),
          plot.title = element_text(size = 6, face = "bold", hjust = 0.5)
    ) +
    geom_label_repel(
        data = correlations_csa_pre |>
            dplyr::filter(!label == ""),
        aes(label = label,
            fill = correlated),
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
        "#e9d1d4",
        "lightgray"
    )) +
    xlab("r") +
    ylab("-log10 p-value") +
    ggtitle("Protein abundance at baseline \n and fiber hypertrophy correlations")

ggsave(plot = correlations_pre_plot, here::here('figures/figure_5/correlations_pre_csa.pdf'), height = 60, width = 60, units = "mm")


#enrichment analysis using ranked protein list, ranked based on correlation
#create ranked protein list
gsea_list_r_csa_pre <- as.numeric(correlations_csa_pre$cor)
names(gsea_list_r_csa_pre) = as.character(correlations_csa_pre$protein)
gsea_list_r_csa_pre <- gsea_list_r_csa_pre[!is.na(gsea_list_r_csa_pre)]
gsea_list_r_csa_pre <- sort(gsea_list_r_csa_pre, decreasing = TRUE)

#gene set enrichment analysis biological process of correlations between protein log2fc and delta csa
gsea_r_csa_pre_bp <- gseGO(
    geneList = gsea_list_r_csa_pre,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pvalueCutoff = 0.05,
    keyType = "SYMBOL",
    eps = 0,
)

gsea_r_csa_pre_bp <- clusterProfiler::simplify(gsea_r_csa_pre_bp, cutoff=0.5, by="qvalue", select_fun=min)


#save as data frame
gsea_r_csa_pre_bp <- as.data.frame(gsea_r_csa_pre_bp)

#add -10log p-value to dataframe
gsea_r_csa_pre_bp <- gsea_r_csa_pre_bp %>%
    mutate("log10p" = -log10(.$pvalue))

write.csv(gsea_r_csa_pre_bp, here::here("results/pre_csa_correlations_enrichment.csv"))

#Figure of most enriched GO-terms (ranked protein list)
enriched_r_pre_plot <- gsea_r_csa_pre_bp %>%
    filter(Description %in% c("translation",
                            "negative regulation of peptidase activity",
                            "ribosome biogenesis",
                            "glycolytic process",
                            "pyruvate metabolic process",
                            "cell surface receptor signaling pathway",
                            "intermediate filament organization",
                            "positive regulation of cell adhesion",
                            "post-transcriptional regulation of gene expression",
                            "polyol metabolic process",
                            "regulation of epithelial cell apoptotic process",
                            "cell adhesion",
                            "regulation of macroautophagy",
                            "lipid storage")) %>%
    ggplot(aes(x = log10p, y = reorder(Description, NES), fill = NES)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "#007ea7", high = "#b4656f") +
    theme(
        panel.background = element_rect(color = "black", fill=NA, linewidth = 0.1),
        panel.grid.minor=element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(color = "white"),
        #legend.position = "none",
        text = element_text(size = 6),
        axis.text.x= element_text(color="black", size = 6),
        axis.text.y= element_text(color="black", size = 6),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1)) +
    labs(x = "-log10 p-value", y="", title = "") +
    xlim(0, 14) +
    theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5)
    ) +
    ggtitle("GO:BP terms of baseline protein abundance \n correlated with fiber hypertrophy")

ggsave(plot = enriched_r_pre_plot, here::here('figures/figure_5/enrichment_pre_csa.pdf'), height = 60, width = 120, units = "mm")
