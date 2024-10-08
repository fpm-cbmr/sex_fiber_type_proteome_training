source(here::here("R/Library.R"))

# Main effect of sex ------------------------------------------------------

#load normalized log2 data
df <- vroom::vroom(here::here("data/data_log2_normalized.csv")) %>%
    rename_with(~ "gene", colnames(.)[1]) %>%
    column_to_rownames("gene")

#load metadata and create groups
metadata <- vroom::vroom(here::here("data-raw/Metadata_REX_pooled_fibers.csv")) %>%
    tidyr::unite(col = "group",
                 c(sex, fibertype), sep ="_", remove = FALSE, na.rm = FALSE) %>%
    column_to_rownames("sample_id") %>%
    dplyr::filter(subject != "rex09")

#Create SummarizedExperiment
se <- PhosphoExperiment(assay = list(df), colData=metadata)

se_pre <- se[,se$trial == 'pre']

design <- model.matrix(~0+ se_pre$sex)

colnames(design)=c("female", "male")

correlation <- duplicateCorrelation(assay(se_pre), design, block=se_pre$subject)

contrast <- makeContrasts(male - female,
                          levels = design)

fit <- eBayes(lmFit(assay(se_pre), design, block=se_pre$subject, correlation=correlation$consensus))

fit2 <- eBayes(contrasts.fit(fit, contrast))

#main effect of sex
baseline_male_vs_female <- topTable(fit2, coef = 1, number = Inf, sort.by = "logFC") %>%
    dplyr::mutate(xiao=10^-(sqrt(log10(1/(P.Value^logFC))^2))) %>%
    dplyr::mutate(protein = row.names(.),
                  q = qvalue(.$P.Value)$qvalues,
                  qiao = qvalue(.$xiao)$qvalues,
                  "-log10p" = -log10(.$P.Value),
                  regulated = case_when(
                      q < 0.05 ~ "yes",
                      q > 0.05 ~ "no")) %>%
    dplyr::arrange(desc(logFC))

write.csv(baseline_male_vs_female, here::here("results/sex_main_effect.csv"))


# Main effect of training -------------------------------------------------

#load normalized log2 datadf <- vroom::vroom(here::here("data/data_log2_normalized.csv")) %>%
rename_with(~ "gene", colnames(.)[1]) %>%
    column_to_rownames("gene")

#load metadata and create groups
metadata <- vroom::vroom(here::here("data-raw/Metadata_REX_pooled_fibers.csv")) %>%
    tidyr::unite(col = "group",
                 c(sex, fibertype), sep ="_", remove = FALSE, na.rm = FALSE) %>%
    dplyr::filter(sample_id %in% colnames(df)) |>
    column_to_rownames("sample_id")

#Create SummarizedExperiment
se <- PhosphoExperiment(assay = list(df), colData=metadata)

TS_interaction <- paste(se$trial, se$fibertype, sep=".")

TS_interaction <- factor(TS_interaction, levels=c("pre.mhc1","post.mhc1","pre.mhc2","post.mhc2"))

design_interaction <- model.matrix(~0+ TS_interaction)

colnames(design_interaction)=c("pre.mhc1","post.mhc1","pre.mhc2","post.mhc2")

correlation_interaction <- duplicateCorrelation(assay(se), design_interaction, block=se$subject)

fit_interaction <- eBayes(lmFit(assay(se), design_interaction, block=se$subject, correlation=correlation_interaction$consensus))

contrast_interaction <- makeContrasts(
    trial = (post.mhc1 - pre.mhc1 + post.mhc2 - pre.mhc2) / 2,
    mhc1 = post.mhc1 - pre.mhc1,
    mhc2 = post.mhc2 - pre.mhc2,
    interaction=(post.mhc1 - pre.mhc1)-(post.mhc2 - pre.mhc2),
    levels = design_interaction)

fit2_interaction <- eBayes(contrasts.fit(fit_interaction, contrast_interaction))

#Main effect of trial
results_trial <- topTable(fit2_interaction, coef = "trial", number = Inf) %>%
    dplyr::mutate(xiao=10^-(sqrt(log10(1/(P.Value^logFC))^2))) %>%
    dplyr::mutate(protein = row.names(.),
                  q = qvalue(.$P.Value)$qvalues,
                  qiao = qvalue(.$xiao)$qvalues,
                  "-log10p" = -log10(.$P.Value),
                  regulated_xiao = ifelse(xiao < 0.05, "+", ""),
                  regulated_qiao = ifelse(qiao < 0.05, "+", ""),
                  regulated_q = ifelse(q < 0.05, "+", "")) %>%
    dplyr::arrange(desc(logFC))

write.csv(results_trial, here::here("results/trial_main_effect.csv"))
