source(here::here("R/Library.R"))
#load data
data_raw <- vroom::vroom(
    here::here("data-raw/library_based_data_raw.csv")
)

# load Metadata:

metadata <- vroom::vroom(
    here::here("data-raw/Metadata_REX_pooled_fibers.csv")
)


# Setup of dataframe:

data_wrangled <- data_raw |>
    dplyr::select(!Protein.Group) |> #! means "don't"
    dplyr::select(!First.Protein.Description) |>
    dplyr::select(!Protein.Names) |>
    dplyr::select(!Protein.Ids) |>
    dplyr::filter(!duplicated(Genes)) |>
    dplyr::filter(!is.na(Genes) == TRUE) |>
    tibble::column_to_rownames("Genes") |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("column_names") |>
    dplyr::mutate(new_column_names = stringr::str_extract(column_names, "(?<=-)(.*?)(?=_)"))

metadata <- metadata |>
    dplyr::mutate(position = toupper(position)) # toupper changes to upper case

data_wrangled <- data_wrangled |>
    dplyr::rename("position" = new_column_names) |>
    dplyr::group_by(
        metadata |>
            dplyr::select(c(position,
                            sample_id))
    ) |>
    tibble::column_to_rownames("sample_id") |>
    dplyr::select(!position) |>
    dplyr::select(!column_names) |>
    t() |>
    as.data.frame()

#Make sample_id row names in metadata
metadata <- metadata %>%
    tibble::column_to_rownames("sample_id")

#create groups in metadata
metadata <- metadata %>%
    tidyr::unite(col = "group",
                 c(sex, fibertype), sep ="_", remove = FALSE, na.rm = FALSE)

#need to filter for proteins identified in 70% of samples in at least one group
data_filtered <- data_wrangled %>%
    PhosR::selectGrps(grps=metadata$group, percent = 0.7, n = 1)

#setup data in long format
data_long <- data_filtered %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    group_by(metadata %>%
                 rownames_to_column("sample_id") %>%
                 dplyr::select(!fibers) %>%
                 dplyr::select(!CSA)) %>%  #merge with metadata and make sample_id column for common identifyer
    pivot_longer(where(is.numeric), names_to = "gene", values_to = "expression")

#Make pre and post to factors in data_long so they appear in correct order
data_long$trial <- factor(data_long$trial, levels = c("pre", "post"))

#Check how pure fiber pools are with non-log transformed data
myh_df <- data_long %>%
    filter(grepl('MYH', gene)) %>%
    group_by(sample_id) %>%
    mutate(expression = expression / sum(expression, na.rm = TRUE)*100)

#visualize MYH distribution using barplot

myh_df$gene <- factor(myh_df$gene, levels = c("MYH10", "MYH11", "MYH13", "MYH14", "MYH15", "MYH3", "MYH4", "MYH6", "MYH7B", "MYH8", "MYH9", "MYH1", "MYH2", "MYH7"))

myh <- myh_df %>%
    ggplot(aes(x = sample_id,
               y = expression,
               fill = gene)) +
    geom_bar(stat = "identity",
             width = 0.8) +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_manual("Gene", values = c(
        "#000004FF",
        "#150E37FF",
        "#3B0F70FF",
        "#641A80FF",
        "#8C2981FF",
        "#B63679FF",
        "#DE4968FF",
        "#F76F5CFF",
        "#FE9F6DFF",
        "#FECE91FF",
        "#FCFDBFFF",
        "#FDE725FF",
        "#67CC5CFF",
        "#440154FF"
    )) +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    facet_grid(~fibertype, scales = "free_x", labeller = labeller(fibertype = c("mhc1" = "Type I", "mhc2" = "Type II"))) +
    theme(text = ggplot2::element_text(size = 8),
          strip.text = element_text(size = 8),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.key.size = ggplot2::unit(1, units = "mm"),
          legend.position = "right",
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(0,-5,0,-7),
          legend.title = element_blank()) +
    xlab("Sample") +
    ylab("MYH expression (%)")

ggsave(plot = myh, here::here('figures/figure_1/figure1_myh.pdf'), height = 60, width = 190, units = "mm")
