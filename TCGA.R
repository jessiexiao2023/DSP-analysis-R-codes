library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)
library(cBioPortalData)
source('tools.R')
# remotes::install_github("omnideconv/immunedeconv")
library(immunedeconv)

 

# data download ---------------------------------------
allProj <- c(
  "TCGA-UCEC",
  "TCGA-READ",
  "TCGA-COAD"
)

allData <-
  allProj %>%
  set_names(., .) %>%
  map(function(proj){
    
    # gene expression query
    query <- GDCquery(
      project       = proj,
      data.category = "Transcriptome Profiling", 
      data.type     = "Gene Expression Quantification",
      workflow.type = "STAR - Counts"
      # sample.type   = NA # Primary Tumor, ... 
    )
    # msi query
    query_mis <- GDCquery(project = proj, 
                          data.category = "Other",
                          legacy = TRUE,
                          access = "open",
                          data.type = "Auxiliary test"
    ) 
    
    # query clinical
    query_clic <- GDCquery(project = proj, 
                           data.category = "Clinical",
                           data.format = "BCR XML"
    ) 
    
    # save query
    saveRDS(
      list(RNA = query, msi = query_mis, clic = query_clic),
      file = str_glue("{proj}-query.rds")
    )
    
    # download
    GDCdownload(
      query           = query,
      method          = "api",
      files.per.chunk = 5
    )
    
    GDCdownload(query_mis)
    GDCdownload(query_clic)
    
    # prepare
    dat <- GDCprepare(
      query         = query,
      save          = TRUE,
      save.filename = str_glue("{proj}.rda"),
      remove.files.prepared = FALSE
    )
    msi_results <- GDCprepare_clinic(query_mis, "msi")
    clic <- GDCprepare_clinic(query_clic, "follow_up")
    
    list(dat = dat, 
         msi = msi_results,
         clic = clic)
  })


# batch prepare -----------------------------------------------------------
allData_final <-
  allData %>%
  imap(function(tcga, name){
    TCGA_dat <- tcga$dat
    
    dat <- assay(TCGA_dat, "fpkm_uq_unstrand")
    meta_roi  <- colData(TCGA_dat)
    meta_gene <- rowData(TCGA_dat)
    
    stopifnot(all(rownames(dat) == rownames(meta_gene)))
    stopifnot(all(colnames(dat) == rownames(meta_roi)))
    
    # duplicated gene -> mean
    duplicated_gene <- unique(meta_gene$gene_name[duplicated(meta_gene$gene_name)])
    dat_1 <-
      cbind(dat, meta_gene[, "gene_name", drop=FALSE]) %>%
      .[.$gene_name %in% duplicated_gene, ] %>% 
      as.data.frame() %>%
      group_by(gene_name) %>%
      summarise(across(everything(), mean))
    dat_2 <- 
      cbind(dat, meta_gene[, "gene_name", drop=FALSE]) %>%
      .[!.$gene_name %in% duplicated_gene, ] %>%
      as.data.frame()
    
    dat_final <- rbind(
      dat_1,
      dat_2
    ) %>%
      remove_rownames() %>%
      column_to_rownames("gene_name")
    
    
    meta_gene_final <- meta_gene[match(rownames(dat_final), meta_gene$gene_name), ]
    rownames(meta_gene_final) <- meta_gene_final$gene_name
    
    
    # only use MSI data
    # Indeterminate         MSI-H         MSI-L           MSS
    msi <- tcga$msi %>% 
      dplyr::rename(
        Patient_ID = bcr_patient_barcode,
        MS_status = mononucleotide_and_dinucleotide_marker_panel_analysis_status
      ) %>% 
      .[.$MS_status %in% c('MSI-H'), 'Patient_ID']
    
    
    meta_roi_final <- meta_roi[meta_roi$patient %in% msi, ]
    meta_roi_final$barcode2 <- meta_roi_final$barcode %>% make.names()
    rownames(meta_roi_final) <- meta_roi_final$barcode2 
    
    dat_final <- dat_final[, rownames(meta_roi_final)]
    
    stopifnot(all(rownames(dat_final) == rownames(meta_gene_final)))
    stopifnot(all(colnames(dat_final) == rownames(meta_roi_final)))
    
    
    
    # final dat
    all_data_final <-
      list(
        counts = dat_final,
        meta_roi  = meta_roi_final %>% as.data.frame(),
        meta_gene = meta_gene_final %>% as.data.frame(),
        msi  = tcga$msi,
        clic = tcga$clic
      )
    all_data_final
  })


# ImmuneDecon -------------------------------------------------------------
set_cibersort_binary("data/CIBERSORT.R")
set_cibersort_mat("data/LM22.txt")

allData_final_addCibersort <-
  allData_final %>%
  imap(function(dat, name){
    mat  <- dat$counts
    res_cibersort_abs <- deconvolute(data.frame(mat), "cibersort_abs")
    
    dat$res_cibersort <- res_cibersort_abs
    dat
  })