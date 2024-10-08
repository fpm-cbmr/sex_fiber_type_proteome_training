source(here::here("R/Library.R"))

#Load datasets of differential expression pre to post
female_mhc1 <- vroom::vroom(here::here("results/differential_expression_female_mhc1.csv"))
female_mhc2 <- vroom::vroom(here::here("results/differential_expression_female_mhc2.csv"))
male_mhc1 <- vroom::vroom(here::here("results/differential_expression_male_mhc1.csv"))
male_mhc2 <- vroom::vroom(here::here("results/differential_expression_male_mhc2.csv"))

##FEMALE TYPE I FIBERS##

#create ranked protein list
gsea_list_f_t1 <- as.numeric(female_mhc1$logFC)
names(gsea_list_f_t1) = as.character(female_mhc1$protein)
gsea_list_f_t1 <- gsea_list_f_t1[!is.na(gsea_list_f_t1)]

#gene set enrichment analysis biological process of female type I fibers
gsea_f_t1_bp <- gseGO(
    geneList = gsea_list_f_t1,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pvalueCutoff = 0.05,
    keyType = "SYMBOL",
    eps = 0,
)


#Filter for redundancy of GO-terms. Terms with more than 50% similarity are omitted and the most significant is used.
gsea_f_t1_bp <- clusterProfiler::simplify(gsea_f_t1_bp, cutoff=0.5, by="qvalue", select_fun=min)

#save as data frame
gsea_f_t1_bp <- as.data.frame(gsea_f_t1_bp)

write.csv(gsea_f_t1_bp, here::here("results/gsea_bp_female_mhc1.csv"))

##FEMALE TYPE II FIBERS##

#create ranked protein list
gsea_list_f_t2 <- as.numeric(female_mhc2$logFC)
names(gsea_list_f_t2) = as.character(female_mhc2$protein)
gsea_list_f_t2 <- gsea_list_f_t2[!is.na(gsea_list_f_t2)]

#gene set enrichment analysis biological process of female type II fibers
gsea_f_t2_bp <- gseGO(
    geneList = gsea_list_f_t2,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pvalueCutoff = 0.05,
    keyType = "SYMBOL",
    eps = 0,
)

#Filter for redundancy of GO-terms. Terms with more than 50% similarity are omitted and the most significant is used.
gsea_f_t2_bp <- clusterProfiler::simplify(gsea_f_t2_bp, cutoff=0.5, by="qvalue", select_fun=min)

#save as data frame
gsea_f_t2_bp <- as.data.frame(gsea_f_t2_bp)

write.csv(gsea_f_t2_bp, here::here("results/gsea_bp_female_mhc2.csv"))

##MALE TYPE I FIBERS##

#create ranked protein list
gsea_list_m_t1 <- as.numeric(male_mhc1$logFC)
names(gsea_list_m_t1) = as.character(male_mhc1$protein)
gsea_list_m_t1 <- gsea_list_m_t1[!is.na(gsea_list_m_t1)]

#gene set enrichment analysis biological process of male type I fibers
gsea_m_t1_bp <- gseGO(
    geneList = gsea_list_m_t1,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pvalueCutoff = 0.05,
    keyType = "SYMBOL",
    eps = 0,
)

#Filter for redundancy of GO-terms. Terms with more than 50% similarity are omitted and the most significant is used.
gsea_m_t1_bp <- clusterProfiler::simplify(gsea_m_t1_bp, cutoff=0.5, by="qvalue", select_fun=min)

#save as data frame
gsea_m_t1_bp <- as.data.frame(gsea_m_t1_bp)

write.csv(gsea_m_t1_bp, here::here("results/gsea_bp_male_mhc1.csv"))

##MALE TYPE II FIBERS##

#create ranked protein list
gsea_list_m_t2 <- as.numeric(male_mhc2$logFC)
names(gsea_list_m_t2) = as.character(male_mhc2$protein)
gsea_list_m_t2 <- gsea_list_m_t2[!is.na(gsea_list_m_t2)]

#gene set enrichment analysis biological process of male type II fibers
gsea_m_t2_bp <- gseGO(
    geneList = gsea_list_m_t2,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pvalueCutoff = 0.05,
    keyType = "SYMBOL",
    eps = 0,
)

#Filter for redundancy of GO-terms. Terms with more than 50% similarity are omitted and the most significant is used.
gsea_m_t2_bp <- clusterProfiler::simplify(gsea_m_t2_bp, cutoff=0.5, by="qvalue", select_fun=min)

#save as data frame
gsea_m_t2_bp <- as.data.frame(gsea_m_t2_bp)

write.csv(gsea_m_t2_bp, here::here("results/gsea_bp_male_mhc2.csv"))

#Load differential expression datasets of sex and fibertype differences at baseline
female_fibertypes <- vroom::vroom(here::here("results/differential_expression_female_fibertypes.csv"))
male_fibertypes <- vroom::vroom(here::here("results/differential_expression_male_fibertypes.csv"))
sex_mhc1 <- vroom::vroom(here::here("results/differential_expression_mhc1_sex.csv"))
sex_mhc2 <- vroom::vroom(here::here("results/differential_expression_mhc2_sex.csv"))

##FEMALE FIBERTYPES##

#create ranked protein list
gsea_list_female <- as.numeric(female_fibertypes$logFC)
names(gsea_list_female) = as.character(female_fibertypes$protein)
gsea_list_female <- gsea_list_female[!is.na(gsea_list_female)]

#gene set enrichment analysis biological process of female fibertypes
gsea_female_bp <- gseGO(
    geneList = gsea_list_female,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pvalueCutoff = 0.05,
    keyType = "SYMBOL",
    eps = 0,
)

#Filter for redundancy of GO-terms. Terms with more than 50% similarity are omitted and the most significant is used.
gsea_female_bp <- clusterProfiler::simplify(gsea_female_bp, cutoff=0.5, by="qvalue", select_fun=min)

#save as data frame
gsea_female_bp <- as.data.frame(gsea_female_bp)

write.csv(gsea_female_bp, here::here("results/gsea_bp_female_fibertypes.csv"))

##MALE FIBERTYPES##

#create ranked protein list
gsea_list_male <- as.numeric(male_fibertypes$logFC)
names(gsea_list_male) = as.character(male_fibertypes$protein)
gsea_list_male <- gsea_list_male[!is.na(gsea_list_male)]

#gene set enrichment analysis biological process of male fibertypes
gsea_male_bp <- gseGO(
    geneList = gsea_list_male,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pvalueCutoff = 0.05,
    keyType = "SYMBOL",
    eps = 0,
)

#Filter for redundancy of GO-terms. Terms with more than 50% similarity are omitted and the most significant is used.
gsea_male_bp <- clusterProfiler::simplify(gsea_male_bp, cutoff=0.5, by="qvalue", select_fun=min)

#save as data frame
gsea_male_bp <- as.data.frame(gsea_male_bp)

write.csv(gsea_male_bp, here::here("results/gsea_bp_male_fibertypes.csv"))

##TYPE I SEX DIFFERENCES##

#create ranked protein list
gsea_list_mhc1 <- as.numeric(sex_mhc1$logFC)
names(gsea_list_mhc1) = as.character(sex_mhc1$protein)
gsea_list_mhc1 <- gsea_list_mhc1[!is.na(gsea_list_mhc1)]

#gene set enrichment analysis biological process of type I difference between sexes
gsea_mhc1_bp <- gseGO(
    geneList = gsea_list_mhc1,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pvalueCutoff = 0.05,
    keyType = "SYMBOL",
    eps = 0,
)

#Filter for redundancy of GO-terms. Terms with more than 50% similarity are omitted and the most significant is used.
gsea_mhc1_bp <- clusterProfiler::simplify(gsea_mhc1_bp, cutoff=0.5, by="qvalue", select_fun=min)

#save as data frame
gsea_mhc1_bp <- as.data.frame(gsea_mhc1_bp)

write.csv(gsea_mhc1_bp, here::here("results/gsea_bp_mhc1_sex.csv"))

##TYPE II SEX DIFFERENCES##

#create ranked protein list
gsea_list_mhc2 <- as.numeric(sex_mhc2$logFC)
names(gsea_list_mhc2) = as.character(sex_mhc2$protein)
gsea_list_mhc2 <- gsea_list_mhc2[!is.na(gsea_list_mhc2)]

#gene set enrichment analysis biological process of type II difference between sexes
gsea_mhc2_bp <- gseGO(
    geneList = gsea_list_mhc2,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pvalueCutoff = 0.05,
    keyType = "SYMBOL",
    eps = 0,
)

#Filter for redundancy of GO-terms. Terms with more than 50% similarity are omitted and the most significant is used.
gsea_mhc2_bp <- clusterProfiler::simplify(gsea_mhc2_bp, cutoff=0.5, by="qvalue", select_fun=min)

#save as data frame
gsea_mhc2_bp <- as.data.frame(gsea_mhc2_bp)

write.csv(gsea_mhc2_bp, here::here("results/gsea_bp_mhc2_sex.csv"))
