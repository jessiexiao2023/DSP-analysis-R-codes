library(ConsensusClusterPlus)
library(survival)
library(survminer)
library(SpatialDecon)
library(immunedeconv)

# 1. model ----------------------------------------------------------------

 
all_data_splitByMS[3] %>% 
  map2(names(.), 
       wgcna_heatmap = wgcna_plots_heatmap,
       function(dat, name, wgcna_heatmap){browser()
         
         # cluster
         wgcna_cluster <- 
           wgcna_heatmap[[name]]$tree_col %>% {
             dat <- .
             cutree(dat, h =
                      ifelse(str_detect(name, "MSI"), 
                             22, 
                             ifelse(str_detect(name, "MSS"), 
                                    18, 
                                    stop("Only MSI and MSS supported.")
                             )
                      )
             )
           } %>% 
           enframe(value = "WGCNA_cluster") %>%
           mutate(WGCNA_cluster = as.character(WGCNA_cluster)) 
         
         # 筛选WGCNA基因
         gene <- wgcna_heatmap[[name]]$tree_row$labels
         
         counts_ana <- dat[['counts_ana']][gene, ]
         meta_ana_gene <- dat[['meta_ana_gene']][gene, ]
         meta_ana_ROI <- dat[['meta_ana_ROI']]
         
         meta_ana_ROI <-
           meta_ana_ROI %>%
           left_join(wgcna_cluster, by = c(ID = 'name'))
         rownames(meta_ana_ROI) <- meta_ana_ROI$ID
         
         
         # name精细划分
         title <- name
         name <- name %>% str_split(" -> ") %>% .[[1]]
         location <- name[1]
         name <- name[2]
         
         # 验证cluster
         pheatmap::pheatmap(
           mat = counts_ana %>% log(),
           color = colorRampPalette(c("green", "black", "red"))(50),
           scale = "row",
           border_color = NA,
           fontsize = 8,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method= "ward.D2",
           annotation_col = meta_ana_ROI[, c("Group5", "WGCNA_cluster"), drop=FALSE],
           annotation_colors = meta_ana_ROI[, c("Group5", "WGCNA_cluster"), drop=FALSE] %>% 
             makeColors(random = TRUE, detect_continuous = FALSE)
         )
         
         
         ###  建模 cycle 1  ----------------------------
         x <- counts_ana %>% t
         
         y <- meta_ana_ROI[["WGCNA_cluster"]] %>%
           fct_recode(cold = "1", "non-cold" = "2", "non-cold" = "3")
         
         # 2  4  6 10 16 17 18 20 22 24 27 29 33 35 38
         # 2  3  6  8 15 17 23 26 29 31 32 36 40 41 43
         set.seed(1234)
         test_idx <- sample(seq.int(length(y)), ceiling(length(y)/3), replace = FALSE)
         test_idx_used <- seq.int(length(y)) %in% test_idx 
         
         useAll = TRUE # FALSE
         manualStandardize = TRUE
         if(useAll){
           x_train <- x
           y_train <- y
           
           x_test <- x
           y_test <- y
           
           test_idx_used <- TRUE
         }else{
           x_train <- x[!test_idx_used, ,drop=F]
           y_train <- y[!test_idx_used]
           
           x_test <- x[test_idx_used, ,drop=F]
           y_test <- y[test_idx_used]
         }
         
         if(manualStandardize){
           x_train <- scale(x_train)
           x_test  <- scale(x_test)
         }
         
         lasso0 <- glmnet(x_train, 
                          y_train, 
                          family = "binomial", 
                          alpha = 1, 
                          standardize = !manualStandardize)
         pdf("lasso 1.pdf", height = 5, width = 6)
         plot(lasso0, xvar = "lambda")
         dev.off()
         
         cv.lasso <- cv.glmnet(x_train, 
                               y_train,
                               family = "binomial", 
                               alpha = 1, 
                               standardize = !manualStandardize,
                               nfolds = length(y_train))
         pdf("lasso 1 cross validation.pdf", height = 5, width = 6)
         plot(cv.lasso)
         dev.off()
         
         lasso <- glmnet(x_train,
                         y_train, 
                         family = "binomial", 
                         alpha = 1, 
                         standardize = !manualStandardize, 
                         lambda = cv.lasso$lambda.1se)
         
         lasso$lambda
         lasso$beta
         lasso$a0
         coef(lasso)
         
         coefs <- coef(cv.lasso, s = cv.lasso$lambda.1se) 
         
         cbind(1, x_train) %*% coefs
         cbind(1, x_train) %*% coefs %>% apply(1, function(x) exp(x)/(1+exp(x)) ) %>% t
         cbind(1, x_train) %*% coefs %>% apply(1, function(x) exp(x)/(1+exp(x)) ) %>% sapply(function(x)ifelse(x < 0.5, 1, 2)) 
         
         pred <- predict(lasso, x_train, type = "link")
         pred <- predict(lasso, x_train, type = "response")
         pred <- predict(lasso, x_train, type = "class")
         
         none_zero <- coefs %>% .[apply(. != 0, 1, any), ]
         
         lasso_gene <- names(none_zero)[-1]
         
         list(coefs, none_zero) %>%
           set_names(c("full", "ultimate")) %>%
           map(~set_colnames(.x, paste0("Cluster", colnames(.x)))) %>%
           map(as.matrix) %>%
           map(as.data.frame) %>%
           map(rownames_to_column, "Coefficient") %>%
           openxlsx::write.xlsx(str_glue("{location} {name} model coef.xlsx"),
                                overwrite = TRUE,
                                asTable = TRUE)
         
         # 测试建模
         pheatmap::pheatmap(
           mat = counts_ana[lasso_gene, ] %>% log,
           scale = "row",
           color = colorRampPalette(c("green", "black", "red"))(50),
           clustering_method= "ward.D2",
           border_color = NA,
           annotation_col = meta_ana_ROI[, c("Group5", "WGCNA_cluster"), drop=FALSE],
           annotation_colors = meta_ana_ROI[, c("Group5", "WGCNA_cluster"), drop=FALSE] %>% 
             makeColors(random = TRUE, detect_continuous = FALSE)
         )
         
         lasso_dat <- predict(lasso, x_test, type = "response") %>%
           as.data.frame() %>%
           rownames_to_column("ID") %>%
           left_join(meta_ana_ROI[test_idx_used, c("WGCNA_cluster", "ID"), drop=FALSE]) %>%
           mutate(Group = ifelse(s0 >0.5, 'non-cold', 'cold'))
         
         lasso_dat %>% 
           ggplot(aes(
             x = Group, 
             y = s0
           )) +
           geom_boxplot(outlier.colour = NA) +
           geom_jitter(aes(color = WGCNA_cluster)) +
           ggpubr::theme_pubclean()
         
         ###  建模 cycle 2  ----------------------------
         # omit 1
         cold_idx <- y == "cold"
         x_2 <- x[!cold_idx, ]
         y_2 <- meta_ana_ROI[["WGCNA_cluster"]] %>%
           .[!cold_idx] %>%
           fct_recode(hot = "2", intermediate = "3")
         
         if(manualStandardize){
           x_2 <- scale(x_2)
         }  
         
         lasso2 <- glmnet(x_2, 
                          y_2, 
                          family = "binomial", 
                          alpha = 1, 
                          standardize = !manualStandardize)
         pdf("lasso 2.pdf", height = 5, width = 6)
         plot(lasso2, xvar = "lambda")
         dev.off()
         
         cv.lasso2 <- cv.glmnet(x_2,
                                y_2, 
                                family = "binomial", 
                                alpha = 1, 
                                standardize = !manualStandardize,
                                nfolds = length(y_2))
         pdf("lasso 2 cross validation.pdf", height = 5, width = 6)
         plot(cv.lasso2)
         dev.off()
         
         lasso2 <- glmnet(x_2,
                          y_2, 
                          family = "binomial", 
                          alpha = 1, 
                          standardize = !manualStandardize, 
                          lambda = cv.lasso2$lambda.1se)
         
         coefs2 <- coef(cv.lasso2, s = cv.lasso2$lambda.1se) 
         
         cbind(1, x_2) %*% coefs2
         cbind(1, x_2) %*% coefs2 %>% apply(1, function(x) exp(x)/(1+exp(x)) ) %>% t
         cbind(1, x_2) %*% coefs2 %>% apply(1, function(x) exp(x)/(1+exp(x)) ) %>% sapply(function(x)ifelse(x < 0.5, 1, 2)) 
         
         pred2 <- predict(lasso2, x_2, type = "link")
         pred2 <- predict(lasso2, x_2, type = "response")
         pred2 <- predict(lasso2, x_2, type = "class")
         
         none_zero2 <- coefs2 %>% .[apply(. != 0, 1, any), ]
         
         lasso_gene2 <- names(none_zero2)[-1]
         
         
         # 测试建模
         pheatmap::pheatmap(
           mat = counts_ana[lasso_gene2, !cold_idx] %>% log,
           scale = "row",
           color = colorRampPalette(c("green", "black", "red"))(50),
           clustering_method= "ward.D2",
           border_color = NA,
           annotation_col = meta_ana_ROI[!cold_idx, c("Group5", "WGCNA_cluster"), drop=FALSE],
           annotation_colors = meta_ana_ROI[!cold_idx, c("Group5", "WGCNA_cluster"), drop=FALSE] %>% 
             makeColors(random = TRUE, detect_continuous = FALSE)
         )
         
         lasso_dat2 <- predict(lasso2, x_2, type = "response") %>%
           as.data.frame() %>%
           rownames_to_column("ID") %>%
           left_join(meta_ana_ROI[, c("WGCNA_cluster", "ID"), drop=FALSE]) %>%
           mutate(Group = ifelse(s0 > 0.5, 'intermediate', 'hot'))
         
         lasso_dat2 %>% 
           ggplot(aes(
             x = Group, 
             y = s0
           )) +
           geom_boxplot(outlier.colour = NA) +
           geom_jitter(aes(color = WGCNA_cluster)) +
           ggpubr::theme_pubclean()
         
         # 合并两个模型
         coefs_total <- cbind(coefs, coefs2)
         
         none_zero_total <- coefs_total %>% .[apply(. != 0, 1, any), ]
         
         lasso_gene_total <- rownames(none_zero_total)[-1]
         
         list(coefs_total, none_zero_total) %>%
           set_names(c("full", "ultimate")) %>%
           map(~set_colnames(.x, paste0("Model", seq_along(colnames(.x))))) %>%
           map(as.matrix) %>%
           map(as.data.frame) %>%
           map(rownames_to_column, "Coefficient") %>%
           openxlsx::write.xlsx(str_glue("{location} {name} two round model coef.xlsx"),
                                overwrite = TRUE,
                                asTable = TRUE)
         
         # 测试建模
         meta_ana_ROI$Group <- meta_ana_ROI$WGCNA_cluster %>% fct_recode(cold = '1', hot = '2', intermediate = '3')
         colnames(counts_ana) %<>% 
           str_remove("_.*$") %>% 
           str_split("-") %>%
           map(
             function(x){
               x[2] %<>% str_pad(width = 3, pad = "0")
               x
             }
           ) %>%
           map_chr(~ .x %>% str_c(collapse = "-"))
         rownames(meta_ana_ROI) <- meta_ana_ROI$ID <- colnames(counts_ana)
         
         p_cluster <-
           pheatmap::pheatmap(
             mat = counts_ana[lasso_gene_total, ] %>% log,
             scale = "row",
             color = colorRampPalette(c("green", "black", "red"))(50),
             clustering_method= "ward.D2",
             border_color = NA,
             silent = TRUE,
             annotation_col = meta_ana_ROI[,"Group",drop=F], # meta_ana_ROI[, c("Group5", "Group", "WGCNA_cluster"), drop=FALSE],
             annotation_colors = list(Group = c(
               cold         = '#47B5C2', # '#3F74B8', # '#47B5C2'
               intermediate = '#FA9346', # '#FBA09F', # '#FA9346'
               hot          = '#E61F22'
             ))
           )
         
         p_nocluster <-
           pheatmap::pheatmap(
             mat = counts_ana[c(setdiff(lasso_gene_total, "DNMT3A"), "DNMT3A"), ] %>% log,
             scale = "row",
             color = colorRampPalette(c("green", "black", "red"))(50),
             clustering_method= "ward.D2",
             cluster_rows = FALSE,
             border_color = NA,
             silent = TRUE,
             annotation_col = meta_ana_ROI[,"Group",drop=F], 
             annotation_colors = list(Group = c(
               cold         = '#47B5C2', 
               intermediate = '#FA9346', 
               hot          = '#E61F22'
             ))
           )
         
         
         pdf("14 gene cluster heatmap add cold-hot group.pdf", height = 8, width = 12)
         print(p_cluster)
         dev.off()
         
         pdf("14 gene nocluster heatmap add cold-hot group.pdf", height = 8, width = 12)
         print(p_nocluster)
         dev.off()

         saveRDS(
           object = list(coefs_total, none_zero_total) %>%
             set_names(c("full", "ultimate")) %>%
             map( ~ set_colnames(.x, paste0( "Model", seq_along(colnames(.x))))),
           file = str_glue("two steps coefs.rds")
         )
         
         # enrich by CTA
         CTA_db <-
           dat[["meta_ana_gene"]][, c('Probe_name','TargetGroup')] %>% 
           deframe() %>%
           map(str_split, pattern = ",\\s+") %>% 
           map(~ .x[[1]]) %>%
           map(~ .x[-1]) %>%
           enframe(name = 'Gene', value = "Pathway") %>%
           unnest(Pathway) %>% 
           dplyr::select(
             term = Pathway, 
             gene = Gene
           ) 
         
         panel14_enrich <- 
           clusterProfiler::enricher(
             gene = lasso_gene_total,
             TERM2GENE = CTA_db,
             pvalueCutoff = 0.5,
             qvalueCutoff = 0.5
           ) 
         
         panel14_enrich@result %>%
           openxlsx::write.xlsx("14 panel CTA enrich.xlsx")
         
         
         panel14_enrich@result %>%
           dplyr::filter(pvalue < 0.1) %>%
           dplyr::select(ID, geneID) %>%
           mutate(geneID = str_split(geneID, "/")) %>%
           unnest(geneID) %>%
           summarise(gene = names(table(geneID)), Counts = as.numeric(table(geneID))) %>%
           ggplot(aes(x = reorder(gene, -Counts), y = Counts)) +
           geom_col(width = 0.55) +
           scale_y_continuous(expand = c(0,0)) +
           ggpubr::theme_pubr() +
           labs(
             x = NULL,
             title = 'Counts gene in signifcant pathway',
             caption = "Signifcant p.value cutoff: 0.1"
           )
         
         panel14_enrich %>%
           enrichplot::dotplot(orderBy = "p.adjust") + 
           scale_y_discrete(limits = function(x) rev(x))
         
         ggsave("gene number in signicant pathway.pdf", width = 8, height = 5)
         
         pdf("70 panel CTA enrich dotplot.pdf", height = 6, width = 7)
         clusterProfiler::enricher(
           gene = wgcna_plots_heatmap$`Tumor -> MSI`$tree_row$labels,
           TERM2GENE = CTA_db,
         ) %>% 
           enrichplot::dotplot(orderBy = "p.adjust") + scale_y_discrete(limits = function(x) rev(x))
         dev.off()
         
         
       })



# 2. validation by TCGA dataset ----------------------------------------
model_coefs <- readRDS("two steps coefs.rds")
model_gene  <- model_coefs$ultimate %>% rownames() %>% .[-1]
model_gene_df <- meta_ana_gene[meta_ana_gene$Probe_name %in% model_gene, c("Probe_name", "GeneID")]

TCGA_UCEC <- readRDS("TCGA/UCEC.rds")


### does all model gene in TCGA dataset?
model_gene %in% TCGA_UCEC$meta_gene$gene_name %>% table
model_gene %>% .[!.%in% TCGA_UCEC$meta_gene$gene_name]

# Yes, all model genes existed in TCGA dataset.
model_gene %in% TCGA_UCEC$meta_gene$gene_name %>% table


### predict TCGA 
tcga_x <- t(TCGA_UCEC$counts[model_gene, ])

model_tcga_cycle1 <- model_coefs$ultimate %>% .[.[,1] != 0, 1, drop=FALSE]

y_tcga_cycle1 <-
  tcga_x[, rownames(model_tcga_cycle1)[-1]]  %>%
  apply(2, scale) %>%
  cbind(1, .) %>% 
  {. %*% model_tcga_cycle1} %>% 
  apply(1, function(x) exp(x)/(1+exp(x)) ) %>%  
  sapply(function(x) ifelse(x < 0.5, "cold", "non-cold")) 

### cycle2 ------
tcga_cycle2_idx <- y_tcga_cycle1 != "cold"
tcga_x_cycle2 <- tcga_x[tcga_cycle2_idx, ]

model_tcag_cycle2 <- model_coefs$ultimate %>% .[.[, 2] != 0, 2, drop=FALSE]


y_tcag_cycle2 <- 
  tcga_x_cycle2[, rownames(model_tcag_cycle2)[-1]]  %>%
  apply(2, scale) %>%
  cbind(1, .) %>% 
  {. %*% model_tcag_cycle2} %>% 
  apply(1, function(x) exp(x)/(1+exp(x)) ) %>%  
  sapply(function(x) ifelse(x < 0.5, "hot", "intermediate")) 

y_tcga_pred <- y_tcga_cycle1
y_tcga_pred[tcga_cycle2_idx] <- y_tcag_cycle2
names(y_tcga_pred) <- rownames(tcga_x)
y_tcga_pred <- y_tcga_pred %>% enframe(name = "barcode", value = "group") 

### DE of immune cell between cold-hot  ------
y_pred_id <- y_tcga_pred$barcode

TCGA_UCEC_cibersort <-
  TCGA_UCEC$res_cibersort %>%
  dplyr::select(1, all_of(y_pred_id)) %>% 
  pivot_longer(
    -cell_type,
    names_to = "barcode",
    values_to = "cell_abundance"
  ) %>% 
  group_by(barcode) %>%
  mutate(cell_proportion = cell_abundance/sum(cell_abundance)*100) %>%
  mutate(sample = str_sub(barcode, 14, 15)) %>%
  dplyr::filter(sample == "01")

TCGA_UCEC_cibersort_withPred <-
  left_join(
    TCGA_UCEC_cibersort, 
    y_tcga_pred,
    by = "barcode"
  ) 
TCGA_UCEC_cibersort_withPred %>%
  ggplot(aes(x = group, y = cell_abundance)) +
  geom_boxplot() +
  ggpubr::stat_compare_means(comparisons = list(
    c("cold", "hot"),
    c("cold", "intermediate"),
    c("intermediate", "hot")
  )) +
  facet_wrap(~cell_type, scales = "free_y")

p_cibersort <-
  TCGA_UCEC_cibersort_withPred %>%
  mutate(group = factor(group, levels = c("cold", "intermediate", "hot"))) %>%
  ggpubr::ggboxplot(
    x = "group",
    y = "cell_abundance",
    color = "group",
    facet.by = "cell_type",
    # nrow = 3,
    scales = "free_y"
  ) +
  ggpubr::stat_compare_means(vjust = -2) +
  ggpubr::stat_compare_means(
    aes(label = ..p.signif..),
    # hide.ns = TRUE,
    step.increase = 0.1,
    comparisons = list(
      c("cold", "hot"),
      c("cold", "intermediate"),
      c("intermediate", "hot")
    )) +
  ggpubr::theme_pubclean() +
  theme(panel.grid.major.y = element_blank()) +
  scale_y_continuous(limits = function(x) c(x[1], 1.3*x[2])) +
  scale_color_manual(values = c(cold = "green", hot = "red", intermediate = "gray5")) +
  labs(x = NULL)

ggsave(plot = p_cibersort,
       filename = "TCGA_UCEC_cibersort_validation.pdf",
       height = 15, width = 20)


p_cibersort_props <-
  TCGA_UCEC_cibersort_withPred %>%
  mutate(group = factor(group, levels = c("cold", "intermediate", "hot"))) %>%
  ggpubr::ggboxplot(
    x = "group",
    y = "cell_proportion",
    color = "group",
    facet.by = "cell_type",
    # nrow = 3,
    scales = "free_y"
  ) +
  ggpubr::stat_compare_means(vjust = -2) +
  ggpubr::stat_compare_means(
    aes(label = ..p.signif..),
    # hide.ns = TRUE,
    step.increase = 0.1,
    comparisons = list(
      c("cold", "hot"),
      c("cold", "intermediate"),
      c("intermediate", "hot")
    )) +
  ggpubr::theme_pubclean() +
  theme(panel.grid.major.y = element_blank()) +
  scale_y_continuous(limits = function(x) c(x[1], 1.3*x[2])) +
  scale_color_manual(values = c(cold = "green", hot = "red", intermediate = "gray5")) +
  labs(x = NULL)

ggsave(plot = p_cibersort_props,
       filename = "TCGA_UCEC_cibersort_validation_cell_proportion.pdf",
       height = 15, width = 20)

### survival analysis -------------------------------------------------------
TCGA_UCEC_survival_raw <-
  TCGA_UCEC$meta_roi %>%
  mutate(barcode2 = rownames(.)) %>%
  # 合并MSI和clinical数据
  mutate(bcr_patient_barcode2 = str_sub(bcr_patient_barcode, end = 12)) %>%
  left_join(TCGA_UCEC$msi, by = c(bcr_patient_barcode2 = 'bcr_patient_barcode')) %>% 
  left_join(TCGA_UCEC$clic[, c("days_to_new_tumor_event_after_initial_treatment",
                               "new_tumor_event_after_initial_treatment",
                               'bcr_patient_barcode'), 
                           drop = FALSE],
            by = c(bcr_patient_barcode2 = 'bcr_patient_barcode')) 

TCGA_UCEC_survival <-
  TCGA_UCEC_survival_raw %>%
  dplyr::select(
    barcode,
    barcode2,
    sample_type, 
    vital_status, 
    days_to_death,
    days_to_last_follow_up,
    days_to_new_tumor_event_after_initial_treatment,
    new_tumor_event_after_initial_treatment,
    mononucleotide_and_dinucleotide_marker_panel_analysis_status 
  ) %>%
  mutate(status = ifelse(vital_status == 'Dead', "1", "0")) %>%
  mutate(
    OS = case_when(
      status == 1 ~ days_to_death,
      status == 0 ~ days_to_last_follow_up,
      TRUE ~ NA_integer_
    )
  ) %>%
  mutate(OS_status = status) %>%
  mutate(PFS_status = case_when(
    new_tumor_event_after_initial_treatment == "YES" ~ "1",
    TRUE ~ OS_status
  )) %>%
  mutate(PFS = case_when(
    new_tumor_event_after_initial_treatment == "YES" ~ days_to_new_tumor_event_after_initial_treatment,
    TRUE ~ OS
  )) %>%
  dplyr::filter(!is.na(OS))

TCGA_UCEC_survival_withPred <-
  TCGA_UCEC_survival %>%
  left_join(y_tcga_pred, by = c(barcode2 = 'barcode')) %>%
  dplyr::filter(sample_type == "Primary Tumor") # only use '01' sample
dplyr::filter(OS <= 3000)

TCGA_UCEC_survival_withPred %>% openxlsx::write.xlsx("TCGA_UCEC.survival.full.xlsx")

TCGA_UCEC_survival_withPred %>%
  mutate(bcr_patient_barcode = str_sub(barcode2, end = 12)) %>%
  distinct(bcr_patient_barcode, .keep_all = TRUE) %>%
  openxlsx::write.xlsx("TCGA_UCEC.survival.full.byPatient.xlsx")

# 生存信息status保证是数值类型
TCGA_UCEC_survival_withPred$status <- as.numeric(TCGA_UCEC_survival_withPred$status)

# 手动构造生存分析的拟合公式
pattern <- str_glue("Surv(OS, status) ~ group") 
formu <- pattern %>% as.formula()

# 拟合
fit <- survfit(formu, data=TCGA_UCEC_survival_withPred)
fit$call$formula <- formu

# 差异分析
diff <- survdiff(formu, data=TCGA_UCEC_survival_withPred)
diff$call$formula <- formu

p.value <- 1 - pchisq(diff$chisq, length(diff$n) -1)

TCGA_UCEC_survival_withPred$group %>%
  unique %>%
  set_names(., .) %>%
  map(function(x){
    dat <- TCGA_UCEC_survival_withPred %>% dplyr::filter(group != x)
    diff <- survdiff(formu, data=dat)
    diff$call$formula <- formu
    
    p.value <- 1 - pchisq(diff$chisq, length(diff$n) -1)
    p.value
  }) %>%
  enframe() %>%
  group_by(name) %>%
  mutate(name = name %>% map(~ str_c(setdiff(name, .x), collapse = " - "))) %>%
  unnest(name)

survminer::ggsurvplot(
  fit = fit,
  data = TCGA_UCEC_survival_withPred,
  pval = TRUE, 
  risk.table = TRUE
) 

fit %>% 
  unclass() %>%
  .[c("time", "surv", "strata", "n.censor")] %>%
  imap(~{
    if(.y == "strata"){
      grp <- names(.x) %>% str_extract("=.*$") %>% str_remove("^=")
      res <- c()
      for(i in seq_along(grp)){
        res <- c(res, rep(grp[i], times = .x[i]))
      }
      res
    }else{
      .x
    }
  }) %>%
  as.data.frame() %>%{df <- .
  ggplot(df, aes(x = time, y = surv, color = strata)) +
    # geom_line() +
    # geom_path()+
    geom_step() +
    theme_pubclean() +
    scale_y_continuous(limits = c(0, 1)) +
    geom_point(data = df %>% dplyr::filter(n.censor == 1), shape = 3)
  }

# 3. 14gene panel生存分析 TCGA UCEC---------
tcga_surv_14gene <-
  model_gene %>%
  set_names(., .) %>%
  map(function(x){
    surv_dat <- TCGA_UCEC_survival
    counts_dat <- TCGA_UCEC$counts
    
    match_patient <- all(colnames(counts_dat) == surv_dat$barcode2)
    if(!match_patient){
      inter_patient <- intersect(colnames(counts_dat), surv_dat$barcode2)
      stopifnot(length(inter_patient) > 0)
      
      surv_dat <- 
        cbind(
          surv_dat[match(inter_patient, surv_dat$barcode2), ],
          counts_dat[x, inter_patient] %>% t
        )
      surv_dat$group <- ifelse(surv_dat[[x]] > median(surv_dat[[x]]), "High", "Low")
    }
    
    # 手动构造生存分析的拟合公式
    p_ls <- c("OS", "PFS") %>%
      set_names(., .) %>%
      map(function(y){
        # 生存信息status保证是数值类型
        surv_dat[[str_glue('{y}_status')]] <- as.numeric(surv_dat[[str_glue('{y}_status')]])
        
        pattern <- str_glue("Surv({y}, {y}_status) ~ group") 
        formu <- pattern %>% as.formula()
        
        # 拟合
        fit <- survfit(formu, data=surv_dat)
        fit$call$formula <- formu
        
        # plot
        p <- survminer::ggsurvplot(
          fit = fit,
          data = surv_dat,
          pval = TRUE, 
          risk.table = TRUE
        )
        p$plot <- p$plot + ggtitle(str_glue("{x} --- {y}"))
        p
      })
    p_ls
    
  })

pdf("14 gene TCGA UCEC survival.xlsx", width = 8, height = 7)
tcga_surv_14gene %>% 
  flatDoubleList(sep = "---") %>%
  walk(print)
dev.off()


