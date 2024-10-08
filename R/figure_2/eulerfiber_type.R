#FEMALE TYPE I VS TYPE II FIBERS AT BASELINE

#load results of female type I and type II fiber differential expression
female_fibers <- vroom::vroom(here::here("results/differential_expression_female_fibertypes.csv")) %>%
    dplyr::rename("genes" = "...1") |>
    mutate(
        regulated = case_when(
            logFC > 0 & q < 0.05 ~ "type2",
            logFC < 0 & q < 0.05 ~ "type1",
            q > 0.05 ~ "no"))

female_fast <- female_fibers |>
    dplyr::filter(adj.P.Val < 0.05) |>
    dplyr::filter(logFC > 0) |>
    dplyr::pull(genes)

female_slow <- female_fibers |>
    dplyr::filter(adj.P.Val < 0.05) |>
    dplyr::filter(logFC < 0) |>
    dplyr::pull(genes)

#MALE TYPE I VS TYPE II FIBERS AT BASELINE

#load results of male type I and type II fiber differential expression
male_fibers <- vroom::vroom(here::here("results/differential_expression_male_fibertypes.csv")) %>%
    dplyr::rename("genes" = "...1") |>
    mutate(
        regulated = case_when(
            logFC > 0 & q < 0.05 ~ "type2",
            logFC < 0 & q < 0.05 ~ "type1",
            q > 0.05 ~ "no"))

male_fast <- male_fibers |>
    dplyr::filter(adj.P.Val < 0.05) |>
    dplyr::filter(logFC > 0) |>
    dplyr::pull(genes)

male_slow <- male_fibers |>
    dplyr::filter(adj.P.Val < 0.05) |>
    dplyr::filter(logFC < 0) |>
    dplyr::pull(genes)


# Euler plot DE proteins between fiber types ------------------------------

euler_list <- list(
    female_slow,
    female_fast,
    male_slow,
    male_fast
)

names(euler_list) <- c(
    "Female type I",
    "Female type II",
    "Male type I",
    "Male type II"
)

ggplotify::as.ggplot(plot(
    eulerr::euler(
        euler_list
    ),
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
            c("#000000", "#000000", "#FF7518", "#FF7518")
        ),
        lex = 3
    ),
))

ggplot2::ggsave(here::here("figures/figure_2/euler_DE_fibertype_sex.pdf"),
                units = "mm",
                height = 60,
                width = 60)
