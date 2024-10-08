##INTERMEDIATE FILAMENT PROTEINS##

#load dataframe of protein log2fc with annotations
results_log2fc <- vroom::vroom(here::here("data/results_log2fc_keywords.csv"))

#Intermediate filament organization

#filter for intermediate filament proteins using enriched proteins and gobp keywords
i_filament <-  results_log2fc %>%
    filter(protein %in% c("KRT16", "KRT6A", "KRT6C", "KRT71", "KRT6B", "KRT14", "KRT17", "KRT5", "DSP", "KRT9", "KRT1", "KRT2", "KRT10", "VIM",
                          "NES", "SYNC","AGFG1", "DES")) %>%
    mutate(
        regulated = case_when(
            logFC > 0 & q < 0.05 ~ "up",
            logFC < 0 & q < 0.05 ~ "down",
            q > 0.05 ~ "no"))

#define symbols for figures
symbols_i_filament <- tibble(
    sex = c("female", "female", "female", "female", "female", "female", "female", "female", "female", "female", "female", "female", "female", "female",
            "male", "male", "female"),
    fibertype = c("mhc2", "mhc2", "mhc2", "mhc2", "mhc2", "mhc2", "mhc2", "mhc2", "mhc2", "mhc2", "mhc2", "mhc2", "mhc2", "mhc2",
                  "mhc2", "mhc2", "mhc1"),
    x = c(2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 17, 18, 18),
    y = c(2.25, 2.25, 2.25, 2.25, 2.25, 2.25, 2.25, 2.25, 2.25, 2.25, 2.25, 2.25, 2.25, 2.25, 2.25, 2.25, 2.25),
    label = c("*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*")
)

#plot of log2fc of intermediate filament proteins
i_filament_plot <- i_filament %>%
    ggplot(aes(x = protein, y = logFC, fill = q)) +
    geom_col() +
    facet_grid(fibertype ~ sex, scales = "free_x", space = "free_x",
               labeller = labeller(fibertype = c("mhc1" = "Type I", "mhc2" = "Type II"),
                                   sex = c("female" = "Female", "male" = "Male"))) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
          panel.background = element_rect(color = "black", fill=NA, linewidth = 0.5),
          panel.grid.minor=element_blank(),
          panel.grid.major = element_blank(),
          plot.background = element_blank(),
          strip.text = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          legend.position = "right",
          legend.key.size = ggplot2::unit(4, units = "mm"),
          text = ggplot2::element_text(size = 7),
          plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) +
    geom_hline(yintercept=0, linetype="dashed") +
    xlab("") +
    ylab("Log2fold change (post - pre)") +
    ggtitle("Intermediate filament proteins") +
    labs(fill = "FDR") +
    scale_fill_gradient(breaks = c(0.05, 0.25, 0.50, 0.75)) +
    scale_y_continuous(limits = c(-1, 2.5)) +
    geom_text(data = symbols_i_filament, aes(x = x, y = y, label = label), inherit.aes = FALSE)
    #scale_fill_manual(name = NULL,
                      #values = c("up" = "#b4656f", "down" = "#007ea7", "no" = "gray"),
                      #labels = c("up" = "Upregulated", "down" = "Downregulated", "no" = "Not regulated"))


ggsave(plot = i_filament_plot, here::here('figures/figure_4/i_filament.pdf'), height = 70, width = 190, units = "mm")

##CSRP3##

#Load data of individual log2fc with keywords
data_log2fc <- vroom::vroom(here::here("data/data_log2fc_long_keywords.csv"))

#filter for CSRP3, which  is a positive regulator of myogenesis and plays a role in sensing mechanical stretch
csrp3 <- data_log2fc %>%
    filter(protein == "CSRP3")

# plot of CSRP3
csrp3_plot <- csrp3 %>%
    ggplot(aes(x = fibertype, y = log2fc, fill = fibertype)) +
    geom_violin(trim = TRUE, width=1, linewidth = 0.5, alpha=0.5)+
    geom_boxplot(width=0.25, color="black", fill="white", alpha=0.5)+
    geom_jitter(size=3, width=0, alpha = 0.5,  stroke = 0)+
    geom_hline(yintercept=0, linetype="dashed")+
    scale_fill_manual(values=c("#440154FF", "#67CC5CFF"))+
    scale_x_discrete(labels=c("Type I", "Type II"))+
    ggplot2::theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill=NA, linewidth = 0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title.x = ggplot2::element_blank(),
        text = element_text(size = 8),
        axis.text.x= element_text(color="black", size = 8),
        axis.text.y= element_text(color="black", size = 8),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8)
    )+
    xlab("none") +
    ylab("CSRP3 log2fold change (post - pre)") +
    facet_wrap(~sex,
               labeller=as_labeller(c(male="Males",
                                      female="Females"
               )))

ggsave(plot = csrp3_plot, here::here('figures/figure_4/csrp3.pdf'), height = 60, width = 60, units = "mm")

