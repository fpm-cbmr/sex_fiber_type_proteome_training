source(here::here("R/Library.R"))

#load normalized log2 data
data_log2 <- vroom::vroom(here::here("data/data_log2_normalized.csv")) %>%
    rename_with(~ "gene", colnames(.)[1]) %>%
    column_to_rownames("gene")

#load metadata, remove REX09 and create groups
metadata <- vroom::vroom(here::here("data-raw/Metadata_REX_pooled_fibers.csv")) %>%
    filter(!subject == "rex09") %>%
    tidyr::unite(col = "group",
                 c(sex, fibertype), sep ="_", remove = FALSE, na.rm = FALSE) %>%
    column_to_rownames("sample_id")

#create summarized experiment
se <- SummarizedExperiment(
    assays = list(counts = as.matrix(data_log2)),
    colData = DataFrame(metadata)
)

#FEMALE TYPE I FIBERS PRE TO POST

#filter summarized experiment to only include female_mhc1
se_f_mhc1 <- se [, se$group == "female_mhc1"]

#create design matrix for female_mhc1
design_f_mhc1 <- model.matrix(~0 + se_f_mhc1$trial)
colnames(design_f_mhc1) <- c("post", "pre")

#define comparison using contrast
contrast_f_mhc1 <- makeContrasts(post - pre, levels = design_f_mhc1)

#correlation between samples from same subject
correlation_f_mhc1 <- duplicateCorrelation(assay(se_f_mhc1), design_f_mhc1, block = se_f_mhc1$subject)

#Use Bayes statistics to assess differential expression
fit_f_mhc1 <- eBayes(lmFit(assay(se_f_mhc1), design_f_mhc1, block = se_f_mhc1$subject, correlation = correlation_f_mhc1$consensus.correlation))
ebayes_f_mhc1 <- eBayes(contrasts.fit(fit_f_mhc1, contrast_f_mhc1))

#extract results
results_f_mhc1 <- topTable(ebayes_f_mhc1, coef = 1, number = Inf, sort.by = "logFC") %>%
    dplyr::mutate(xiao=10^-(sqrt(log10(1/(P.Value^logFC))^2))) %>%
    dplyr::mutate(protein = row.names(.),
                  q = qvalue(.$P.Value)$qvalues,
                  qiao = qvalue(.$xiao)$qvalues,
                  "-log10p" = -log10(.$P.Value),
                  regulated_xiao = ifelse(xiao < 0.05, "+", ""),
                  regulated_qiao = ifelse(qiao < 0.05, "+", ""),
                  regulated_q = ifelse(q < 0.05, "+", "")
    )%>%
    arrange(desc(logFC))

write.csv(results_f_mhc1, here::here("results/differential_expression_female_mhc1.csv"))


#FEMALE TYPE II FIBERS PRE TO POST

#filter summarized experiment to only include female_mhc2
se_f_mhc2 <- se [, se$group == "female_mhc2"]

#create design matrix for female_mhc2
design_f_mhc2 <- model.matrix(~0 + se_f_mhc2$trial)
colnames(design_f_mhc2) <- c("post", "pre")

#define comparison using contrast
contrast_f_mhc2 <- makeContrasts(post - pre, levels = design_f_mhc2)

#correlation between samples from same subject
correlation_f_mhc2 <- duplicateCorrelation(assay(se_f_mhc2), design_f_mhc2, block = se_f_mhc2$subject)

#Use Bayes statistics to assess differential expression
fit_f_mhc2 <- eBayes(lmFit(assay(se_f_mhc2), design_f_mhc2, block = se_f_mhc2$subject, correlation = correlation_f_mhc2$consensus.correlation))
ebayes_f_mhc2 <- eBayes(contrasts.fit(fit_f_mhc2, contrast_f_mhc2))

#extract results
results_f_mhc2 <- topTable(ebayes_f_mhc2, coef = 1, number = Inf, sort.by = "logFC") %>%
    dplyr::mutate(xiao=10^-(sqrt(log10(1/(P.Value^logFC))^2))) %>%
    dplyr::mutate(protein = row.names(.),
                  q = qvalue(.$P.Value)$qvalues,
                  qiao = qvalue(.$xiao)$qvalues,
                  "-log10p" = -log10(.$P.Value),
                  regulated_xiao = ifelse(xiao < 0.05, "+", ""),
                  regulated_qiao = ifelse(qiao < 0.05, "+", ""),
                  regulated_q = ifelse(q < 0.05, "+", "")
    )%>%
    arrange(desc(logFC))

write.csv(results_f_mhc2, here::here("results/differential_expression_female_mhc2.csv"))


#MALE TYPE I FIBERS PRE TO POST

#filter summarized experiment to only include male_mhc1
se_m_mhc1 <- se [, se$group == "male_mhc1"]

#create design matrix for male_mhc1
design_m_mhc1 <- model.matrix(~0 + se_m_mhc1$trial)
colnames(design_m_mhc1) <- c("post", "pre")

#define comparison using contrast
contrast_m_mhc1 <- makeContrasts(post - pre, levels = design_m_mhc1)

#correlation between samples from same subject
correlation_m_mhc1 <- duplicateCorrelation(assay(se_m_mhc1), design_m_mhc1, block = se_m_mhc1$subject)

#Use Bayes statistics to assess differential expression
fit_m_mhc1 <- eBayes(lmFit(assay(se_m_mhc1), design_m_mhc1, block = se_m_mhc1$subject, correlation = correlation_m_mhc1$consensus.correlation))
ebayes_m_mhc1 <- eBayes(contrasts.fit(fit_m_mhc1, contrast_m_mhc1))

#extract results
results_m_mhc1 <- topTable(ebayes_m_mhc1, coef = 1, number = Inf, sort.by = "logFC") %>%
    dplyr::mutate(xiao=10^-(sqrt(log10(1/(P.Value^logFC))^2))) %>%
    dplyr::mutate(protein = row.names(.),
                  q = qvalue(.$P.Value)$qvalues,
                  qiao = qvalue(.$xiao)$qvalues,
                  "-log10p" = -log10(.$P.Value),
                  regulated_xiao = ifelse(xiao < 0.05, "+", ""),
                  regulated_qiao = ifelse(qiao < 0.05, "+", ""),
                  regulated_q = ifelse(q < 0.05, "+", "")
    )%>%
    arrange(desc(logFC))

write.csv(results_m_mhc1, here::here("results/differential_expression_male_mhc1.csv"))


#MALE TYPE II FIBERS PRE TO POST

#filter summarized experiment to only include male_mhc2
se_m_mhc2 <- se [, se$group == "male_mhc2"]

#create design matrix for male_mhc1
design_m_mhc2 <- model.matrix(~0 + se_m_mhc2$trial)
colnames(design_m_mhc2) <- c("post", "pre")

#define comparison using contrast
contrast_m_mhc2 <- makeContrasts(post - pre, levels = design_m_mhc2)

#correlation between samples from same subject
correlation_m_mhc2 <- duplicateCorrelation(assay(se_m_mhc2), design_m_mhc2, block = se_m_mhc2$subject)

#Use Bayes statistics to assess differential expression
fit_m_mhc2 <- eBayes(lmFit(assay(se_m_mhc2), design_m_mhc2, block = se_m_mhc2$subject, correlation = correlation_m_mhc2$consensus.correlation))
ebayes_m_mhc2 <- eBayes(contrasts.fit(fit_m_mhc2, contrast_m_mhc2))

#extract results
results_m_mhc2 <- topTable(ebayes_m_mhc2, coef = 1, number = Inf, sort.by = "logFC") %>%
    dplyr::mutate(xiao=10^-(sqrt(log10(1/(P.Value^logFC))^2))) %>%
    dplyr::mutate(protein = row.names(.),
                  q = qvalue(.$P.Value)$qvalues,
                  qiao = qvalue(.$xiao)$qvalues,
                  "-log10p" = -log10(.$P.Value),
                  regulated_xiao = ifelse(xiao < 0.05, "+", ""),
                  regulated_qiao = ifelse(qiao < 0.05, "+", ""),
                  regulated_q = ifelse(q < 0.05, "+", "")
    )%>%
    arrange(desc(logFC))

write.csv(results_m_mhc2, here::here("results/differential_expression_male_mhc2.csv"))

#FEMALE TYPE I VS TYPE II FIBERS

#filter summarized experiment to only include females and pre
se_f_pre <- se [, se$sex == "female"]
se_f_pre <- se_f_pre [, se_f_pre$trial == "pre"]

#create design matrix for female fiber types
design_f_pre <- model.matrix(~0 + se_f_pre$fibertype)
colnames(design_f_pre) <- c("typeI", "typeII")

#define comparison using contrast
contrast_f_pre <- makeContrasts(typeII - typeI, levels = design_f_pre) #type II will appear as positive Log2FC, type I as negative Log2FC

#correlation between samples from same subject
correlation_f_pre <- duplicateCorrelation(assay(se_f_pre), design_f_pre, block = se_f_pre$subject)

#Use Bayes statistics to assess differential expression
fit_f_pre <- eBayes(lmFit(assay(se_f_pre), design_f_pre, block = se_f_pre$subject, correlation = correlation_f_pre$consensus.correlation))
ebayes_f_pre <- eBayes(contrasts.fit(fit_f_pre, contrast_f_pre))

#extract results
results_f_pre <- topTable(ebayes_f_pre, coef = 1, number = Inf, sort.by = "logFC") %>%
    dplyr::mutate(xiao=10^-(sqrt(log10(1/(P.Value^logFC))^2))) %>%
    dplyr::mutate(protein = row.names(.),
                  q = qvalue(.$P.Value)$qvalues,
                  qiao = qvalue(.$xiao)$qvalues,
                  "-log10p" = -log10(.$P.Value),
                  regulated_xiao = ifelse(xiao < 0.05, "+", ""),
                  regulated_qiao = ifelse(qiao < 0.05, "+", ""),
                  regulated_q = ifelse(q < 0.05, "+", "")
    )%>%
    arrange(desc(logFC))

write.csv(results_f_pre, here::here("results/differential_expression_female_fibertypes.csv"))

#MALE TYPE I VS TYPE II FIBERS

#filter summarized experiment to only include males and pre
se_m_pre <- se [, se$sex == "male"]
se_m_pre <- se_m_pre [, se_m_pre$trial == "pre"]

#create design matrix for female fiber types
design_m_pre <- model.matrix(~0 + se_m_pre$fibertype)
colnames(design_m_pre) <- c("typeI", "typeII")

#define comparison using contrast
contrast_m_pre <- makeContrasts(typeII - typeI, levels = design_m_pre) #type II will appear as positive Log2FC, type I as negative Log2FC

#correlation between samples from same subject
correlation_m_pre <- duplicateCorrelation(assay(se_m_pre), design_m_pre, block = se_m_pre$subject)

#Use Bayes statistics to assess differential expression
fit_m_pre <- eBayes(lmFit(assay(se_m_pre), design_m_pre, block = se_m_pre$subject, correlation = correlation_m_pre$consensus.correlation))
ebayes_m_pre <- eBayes(contrasts.fit(fit_m_pre, contrast_m_pre))

#extract results
results_m_pre <- topTable(ebayes_m_pre, coef = 1, number = Inf, sort.by = "logFC") %>%
    dplyr::mutate(xiao=10^-(sqrt(log10(1/(P.Value^logFC))^2))) %>%
    dplyr::mutate(protein = row.names(.),
                  q = qvalue(.$P.Value)$qvalues,
                  qiao = qvalue(.$xiao)$qvalues,
                  "-log10p" = -log10(.$P.Value),
                  regulated_xiao = ifelse(xiao < 0.05, "+", ""),
                  regulated_qiao = ifelse(qiao < 0.05, "+", ""),
                  regulated_q = ifelse(q < 0.05, "+", "")
    )%>%
    arrange(desc(logFC))

write.csv(results_m_pre, here::here("results/differential_expression_male_fibertypes.csv"))


#MALE VS FEMALE TYPE I FIBERS

#filter summarized experiment to only include type I fibers and pre
se_mhc1_pre <- se [, se$fibertype == "mhc1"]
se_mhc1_pre <- se_mhc1_pre [, se_mhc1_pre$trial == "pre"]

#create design matrix for type I fiber difference between sexes
design_mhc1_pre <- model.matrix(~0 + se_mhc1_pre$sex)
colnames(design_mhc1_pre) <- c("female", "male")

#define comparison using contrast
contrast_mhc1_pre <- makeContrasts(male - female, levels = design_mhc1_pre) #male will appear as positive Log2FC, female as negative Log2FC

#correlation between samples from same subject should not influence results as no samples are from same subject
correlation_mhc1_pre <- duplicateCorrelation(assay(se_mhc1_pre), design_mhc1_pre, block = se_mhc1_pre$subject)

#Use Bayes statistics to assess differential expression
fit_mhc1_pre <- eBayes(lmFit(assay(se_mhc1_pre), design_mhc1_pre, block = se_mhc1_pre$subject, correlation = correlation_mhc1_pre$consensus.correlation))
ebayes_mhc1_pre <- eBayes(contrasts.fit(fit_mhc1_pre, contrast_mhc1_pre))

#extract results
results_mhc1_pre <- topTable(ebayes_mhc1_pre, coef = 1, number = Inf, sort.by = "logFC") %>%
    dplyr::mutate(xiao=10^-(sqrt(log10(1/(P.Value^logFC))^2))) %>%
    dplyr::mutate(protein = row.names(.),
                  q = qvalue(.$P.Value)$qvalues,
                  qiao = qvalue(.$xiao)$qvalues,
                  "-log10p" = -log10(.$P.Value),
                  regulated_xiao = ifelse(xiao < 0.05, "+", ""),
                  regulated_qiao = ifelse(qiao < 0.05, "+", ""),
                  regulated_q = ifelse(q < 0.05, "+", "")
    )%>%
    arrange(desc(logFC))

write.csv(results_mhc1_pre, here::here("results/differential_expression_mhc1_sex.csv"))


#MALE VS FEMALE TYPE II FIBERS

#filter summarized experiment to only include type II fibers and pre
se_mhc2_pre <- se [, se$fibertype == "mhc2"]
se_mhc2_pre <- se_mhc2_pre [, se_mhc2_pre$trial == "pre"]

#create design matrix for type II fiber difference between sexes
design_mhc2_pre <- model.matrix(~0 + se_mhc2_pre$sex)
colnames(design_mhc2_pre) <- c("female", "male")

#define comparison using contrast
contrast_mhc2_pre <- makeContrasts(male - female, levels = design_mhc2_pre) #male will appear as positive Log2FC, female as negative Log2FC

#correlation between samples from same subject will be set to zero and not influence results as no samples are from same subject
correlation_mhc2_pre <- duplicateCorrelation(assay(se_mhc2_pre), design_mhc2_pre, block = se_mhc2_pre$subject)

#Use Bayes statistics to assess differential expression
fit_mhc2_pre <- eBayes(lmFit(assay(se_mhc2_pre), design_mhc2_pre, block = se_mhc2_pre$subject, correlation = correlation_mhc2_pre$consensus.correlation))
ebayes_mhc2_pre <- eBayes(contrasts.fit(fit_mhc2_pre, contrast_mhc2_pre))

#extract results
results_mhc2_pre <- topTable(ebayes_mhc2_pre, coef = 1, number = Inf, sort.by = "logFC") %>%
    dplyr::mutate(xiao=10^-(sqrt(log10(1/(P.Value^logFC))^2))) %>%
    dplyr::mutate(protein = row.names(.),
                  q = qvalue(.$P.Value)$qvalues,
                  qiao = qvalue(.$xiao)$qvalues,
                  "-log10p" = -log10(.$P.Value),
                  regulated_xiao = ifelse(xiao < 0.05, "+", ""),
                  regulated_qiao = ifelse(qiao < 0.05, "+", ""),
                  regulated_q = ifelse(q < 0.05, "+", "")
    )%>%
    arrange(desc(logFC))

write.csv(results_mhc2_pre, here::here("results/differential_expression_mhc2_sex.csv"))
