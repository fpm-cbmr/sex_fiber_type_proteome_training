source(here::here("R/Library.R"))
#PRE DATA##

#load data long format, filter for pre samples and make wide format
data_pre <- vroom::vroom(here::here("data/data_log2_long_normalized.csv")) %>%
    filter(trial == "pre") %>%
    dplyr::select(!1) %>% #column 1 was just numbers so removed
    dplyr::select(!"position") %>%
    dplyr::select(!"subject") %>%
    dplyr::select(!"trial") %>%
    dplyr::select(!"fibertype") %>%
    dplyr::select(!"sex") %>%
    pivot_wider(names_from = "gene", values_from = "expression") %>%
    column_to_rownames("sample_id")

#create metadata for each group
m_mhc1_meta <- vroom::vroom(here::here("data-raw/Metadata_REX_pooled_fibers.csv")) %>%
    filter(!subject == "rex09") %>%
    tidyr::unite(col = "group",
                 c(sex, fibertype), sep ="_", remove = FALSE, na.rm = FALSE) %>%
    filter(group == "male_mhc1") %>%
    filter(trial == "pre")

m_mhc2_meta <- vroom::vroom(here::here("data-raw/Metadata_REX_pooled_fibers.csv")) %>%
    filter(!subject == "rex09") %>%
    tidyr::unite(col = "group",
                 c(sex, fibertype), sep ="_", remove = FALSE, na.rm = FALSE) %>%
    filter(group == "male_mhc2") %>%
    filter(trial == "pre")

f_mhc1_meta <- vroom::vroom(here::here("data-raw/Metadata_REX_pooled_fibers.csv")) %>%
    filter(!subject == "rex09") %>%
    tidyr::unite(col = "group",
                 c(sex, fibertype), sep ="_", remove = FALSE, na.rm = FALSE) %>%
    filter(group == "female_mhc1") %>%
    filter(trial == "pre")

f_mhc2_meta <- vroom::vroom(here::here("data-raw/Metadata_REX_pooled_fibers.csv")) %>%
    filter(!subject == "rex09") %>%
    tidyr::unite(col = "group",
                 c(sex, fibertype), sep ="_", remove = FALSE, na.rm = FALSE) %>%
    filter(group == "female_mhc2") %>%
    filter(trial == "pre")

#create dataframe for each group
m_mhc1 <- data_pre %>%
    filter(group == "male_mhc1") %>%
    dplyr::select(!group) %>%  #remove group variable so this only consists of genes
    t() %>%
    as.data.frame() %>%
    PhosR::selectGrps(grps=m_mhc1_meta$group, percent = 0.7)

m_mhc2 <- data_pre %>%
    filter(group == "male_mhc2") %>%
    dplyr::select(!group) %>%  #remove group variable so this only consists of genes
    t() %>%
    as.data.frame() %>%
    PhosR::selectGrps(grps=m_mhc2_meta$group, percent = 0.7)

f_mhc1 <- data_pre %>%
    filter(group == "female_mhc1") %>%
    dplyr::select(!group) %>%  #remove group variable so this only consists of genes
    t() %>%
    as.data.frame() %>%
    PhosR::selectGrps(grps=f_mhc1_meta$group, percent = 0.7)

f_mhc2 <- data_pre %>%
    filter(group == "female_mhc2") %>%
    dplyr::select(!group) %>%  #remove group variable so this only consists of genes
    t() %>%
    as.data.frame() %>%
    PhosR::selectGrps(grps=f_mhc2_meta$group, percent = 0.7)

# create list of data
venn_list <- list(m_mhc1 = as.character(rownames(m_mhc1)),
                  m_mhc2 = as.character(rownames(m_mhc2)),
                  f_mhc1 = as.character(rownames(f_mhc1)),
                  f_mhc2 = as.character(rownames(f_mhc2)))

# create mapping vector to rename groups
mapping_vector <- c("m_mhc1" = "Male type I",
                    "m_mhc2" = "Male type II",
                    "f_mhc1" = "Female type I",
                    "f_mhc2" = "Female type II")

# Replace the names in the list with labels
names(venn_list) <- mapping_vector[names(venn_list)]

#stroke_colors <- c("#FF7518", "#FF7518", "black", "black")


# create venn diagram of pre
venn <- ggvenn(
    venn_list,
    fill_color = c(ggplot2::alpha("#440154FF", alpha = 0.5),
                   ggplot2::alpha("#67CC5CFF", alpha = 0.5),
                   ggplot2::alpha("#440154FF", alpha = 0.5),
                   ggplot2::alpha("#67CC5CFF", alpha = 0.5)),
    #stroke_color = stroke_colors,
    show_percentage = FALSE,
    stroke_size = 0.5, set_name_size = 2
)


# euler_diagram <-

ggplotify::as.ggplot(
    plot(
        eulerr::venn(venn_list),
        labels = list(cex = 0.5),
        quantities = list(
            font = 1,
            cex = 0.5
        ),
        fill = c(
            ggplot2::alpha("#440154FF", alpha = 0.5),
            ggplot2::alpha("#67CC5CFF", alpha = 0.5),
            ggplot2::alpha("#440154FF", alpha = 0.5),
            ggplot2::alpha("#67CC5CFF", alpha = 0.5)
        ),
        alpha = 0.65,
        edges = list(
            col = c(
                c("#FF7518", "#FF7518", "black", "black")
            ),
            lex = 3
        ),
    )
)



# create UpSet plot of pre
upset_plot <- upset(fromList(venn_list),
      order.by = "freq",
      mainbar.y.label = "Identified proteins",
      sets.x.label = "Identified proteins",
      text.scale = 0.85,
      point.size = 0.5,
      line.size = 0.25,
)

# Convert UpSet plot to ggplot object
upset_ggplot <- ggplotify::as.ggplot(upset_plot) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_x_discrete(
        expand = c(0, 0)
    )


ggsave(plot = upset_ggplot, here::here('figures/figure_1/fig1_upset.pdf'), height = 80, width = 120, units = "mm")


