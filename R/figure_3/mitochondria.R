source(here::here("R/Library.R"))

#Load log2fold difference between sexes
log2fd_sex <- vroom::vroom(here::here("data/results_sex_keywords.csv"))

#Add mitocarta
mitocarta <- read_excel(here::here('data-raw/mitocarta.xls'))%>%
    dplyr::select('symbol', 'pathways') %>%
    dplyr::rename(protein=symbol)

log2fd_sex <- log2fd_sex %>%
    merge(mitocarta, by="protein", all.x = T) %>%
    dplyr::rename(mito = pathways)

ci <- log2fd_sex %>%
    dplyr::filter(grepl('CI subunits', mito)) %>%
    dplyr::mutate(complex = "CI")

cii <- log2fd_sex %>%
    dplyr::filter(grepl('CII subunits', mito)) %>%
    dplyr::mutate(complex = "CII")

ciii <- log2fd_sex %>%
    dplyr::filter(grepl('CIII subunits', mito)) %>%
    dplyr::mutate(complex = "CIII")

civ <- log2fd_sex %>%
    dplyr::filter(grepl('CIV subunits', mito)) %>%
    dplyr::mutate(complex = "CIV")

cv <- log2fd_sex %>%
    dplyr::filter(grepl('CV subunits', mito)) %>%
    dplyr::mutate(complex = "CV")

mito <- rbind(ci, cii, ciii, civ, cv)

#mito complexes figure
mito_complexes <- mito %>%
    ggplot(aes(x = complex, y = logFC, fill = fibertype)) +
    geom_violin(trim = TRUE, width = 1, linewidth = 0.5, alpha = 0.5) +
    geom_boxplot(width = 0.25, color = "black", fill = "white", alpha = 0.5) +
    geom_jitter(size=2.5, width=0, alpha = 0.5,  stroke = 0)+
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values=c("#440154FF", "#67CC5CFF"))+
    scale_x_discrete(labels = c("CI" = "Complex \n I", "CII" = "Complex \n II", "CIII" = "Complex \n III", "CIV" = "Complex \n IV", "CV" = "Complex \n V")) +
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
    facet_wrap(~fibertype, labeller = labeller(fibertype = c("mhc1" = "Type I", "mhc2" = "Type II"))) +
    xlab("") +
    ylab("Log2fold difference (male - female)") +
    ggtitle("Mitochondrial complexes")

ggsave(plot = mito_complexes, here::here('figures/figure_3/mito_complexes.pdf'), height = 70, width = 120, units = "mm")

#filter for mitoribosome subunits
mito_ribo <- log2fd_sex %>%
    dplyr::filter(grepl("MRPS|MRPL", protein)) %>%
    dplyr::mutate(subunit = case_when(str_detect(protein, "^MRPS") ~ "MRPS",
                                      TRUE ~ "MRPL"))

#mito complexes figure
mitoribo_fig <- mito_ribo %>%
    ggplot(aes(x = subunit, y = logFC, fill = fibertype)) +
    geom_violin(trim = TRUE, width = 1, linewidth = 0.5, alpha = 0.5) +
    geom_boxplot(width = 0.25, color = "black", fill = "white", alpha = 0.5) +
    geom_jitter(size=2.5, width=0, alpha = 0.5,  stroke = 0)+
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values=c("#440154FF", "#67CC5CFF"))+
    scale_x_discrete(labels = c("MRPL" = "39S \n subunit", "MRPS" = "28S \n subunit")) +
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
    facet_wrap(~fibertype, labeller = labeller(fibertype = c("mhc1" = "Type I", "mhc2" = "Type II"))) +
    xlab("") +
    ylab("Log2fold difference (male - female)") +
    ggtitle("Mitoribosome subunits")

ggsave(plot = mitoribo_fig, here::here('figures/figure_3/mito_ribosome.pdf'), height = 70, width = 70, units = "mm")
