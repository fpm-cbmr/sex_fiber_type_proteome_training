source(here::here("R/Library.R"))

##REGULATED PROTEINS##

#Load differential expression results for each group and filter for regulated proteins
#female type I fibers
f_mhc1_de <- vroom::vroom(here::here("results/differential_expression_female_mhc1.csv")) %>%
    mutate(
        regulated = case_when(
            logFC > 0 & q < 0.05 ~ "up",
            logFC < 0 & q < 0.05 ~ "down",
            q > 0.05 ~ "no")) %>%
    filter(regulated %in% c("up", "down")) %>%
    column_to_rownames("protein")

#female type II fibers
f_mhc2_de <- vroom::vroom(here::here("results/differential_expression_female_mhc2.csv")) %>%
    mutate(
        regulated = case_when(
            logFC > 0 & q < 0.05 ~ "up",
            logFC < 0 & q < 0.05 ~ "down",
            q > 0.05 ~ "no")) %>%
    filter(regulated %in% c("up", "down")) %>%
    column_to_rownames("protein")

#male type I fibers
m_mhc1_de <- vroom::vroom(here::here("results/differential_expression_male_mhc1.csv")) %>%
    mutate(
        regulated = case_when(
            logFC > 0 & q < 0.05 ~ "up",
            logFC < 0 & q < 0.05 ~ "down",
            q > 0.05 ~ "no")) %>%
    filter(regulated %in% c("up", "down")) %>%
    column_to_rownames("protein")

#male type II fibers
m_mhc2_de <- vroom::vroom(here::here("results/differential_expression_male_mhc2.csv")) %>%
    mutate(
        regulated = case_when(
            logFC > 0 & q < 0.05 ~ "up",
            logFC < 0 & q < 0.05 ~ "down",
            q > 0.05 ~ "no")) %>%
    filter(regulated %in% c("up", "down")) %>%
    column_to_rownames("protein")

# create list of differential expression data
venn_list_de <- list(f_mhc1_de = as.character(rownames(f_mhc1_de)),
                     f_mhc2_de = as.character(rownames(f_mhc2_de)),
                     m_mhc1_de = as.character(rownames(m_mhc1_de)),
                     m_mhc2_de = as.character(rownames(m_mhc2_de)))

# create mapping vector to rename groups
mapping_vector_de <- c("m_mhc1_de" = "Male type I",
                       "m_mhc2_de" = "Male type II",
                       "f_mhc1_de" = "Female type I",
                       "f_mhc2_de" = "Female type II")


# create UpSet plot of differential expression
upset_plot <- upset(fromList(venn_list_de),
      order.by = "freq",
      mainbar.y.label = "Regulated proteins",
      sets.x.label = "Regulated proteins",
      text.scale = 2,
      point.size = 3,
      line.size = 0.75
)

# Convert UpSet plot to ggplot object
upset_ggplot <- ggplotify::as.ggplot(upset_plot)


ggsave(plot = upset_ggplot, here::here('figures/figure_4/upset_training.pdf'), height = 140, width = 200, units = "mm")

