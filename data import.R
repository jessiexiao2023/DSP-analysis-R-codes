library(tidyverse)
library(magrittr)
source("tools.R")


majorPath <- "CTApipeline/XX_exportDat_V6"

pval_cutoff <- 0.05  
fc_cutoff    <- 1.5
y_axis <- "pval"
User_columns <- "Group5"
 
# 构造路径

reportPath <-
  list.files(majorPath) %>%
  str_subset("^Report") %>%
  str_subset(., "total", negate = TRUE) %>%
  str_c(majorPath, ., sep="/") 

name <- reportPath %>%
  basename() %>%
  str_remove("Report-") 

names(reportPath) <- name

# 导入数据
# 以Location拆分的报告
all_data <-
  reportPath %>%
  map(function(path) {
    c('de_list_single',
      'de_list_double',
      'QC_V2_targetCounts',
      'counts_ana',
      'meta_ana_gene',
      'meta_ana_ROI') %>%
      set_names(., .) %>%
      map(function(file) {
        file <- paste0(file, ".rds")
        tryCatch(
          readRDS(file.path(path, file)),
          error = function(e) return(NULL)
        )
      })
  })

# 将Group1 - Group5的双分组设置为单分组数据
all_data <- 
  all_data %>% 
  map(function(data){
    if("de_list_double" %in% names(data) &&
       !is.null(data$de_list_double)){
      double <- data$de_list_double$`Group1 - Group5` %>% flatDoubleList(sep = " -> ")
      data$de_list_single$`MSI-MSS-TILs` <- double
      
      data
    }
  })

# 将MSI和MSS作为分组标识提到第一分组
all_data_splitByMS <- 
  all_data %>%
  map(function(x){
    msi_idx <- x$meta_ana_ROI$Group1 == "MSI"
    mss_idx <- x$meta_ana_ROI$Group1 == "MSS"
    
    list(
      MSI = list(
        counts_ana    = x$counts_ana[, msi_idx,drop=FALSE],
        meta_ana_gene = x$meta_ana_gene,
        meta_ana_ROI  = x$meta_ana_ROI[msi_idx, ,drop=FALSE],
        QC_V2_targetCounts = x$QC_V2_targetCounts[, msi_idx,drop=FALSE]
      ),
      MSS = list(
        counts_ana    = x$counts_ana[, mss_idx,drop=FALSE],
        meta_ana_gene = x$meta_ana_gene,
        meta_ana_ROI  = x$meta_ana_ROI[mss_idx, ,drop=FALSE],
        QC_V2_targetCounts    = x$QC_V2_targetCounts[, mss_idx,drop=FALSE]
      )
    )
  }) %>%
  flatDoubleList(sep=" -> ")


### de_gene
de_gene_list <- 
  all_data %>% 
  map2(names(.), function(dat, name){
    counts_ana <- dat[['counts_ana']] 
    
    meta_ana_gene <- dat[['meta_ana_gene']]
    meta_ana_ROI <- dat[['meta_ana_ROI']]
    
    de_list <- dat[["de_list_single"]]
    
    map(de_list, function(de){
      map(de, function(dat){
        data <- dat %>%
          volcanoTable(
            neg_num = Inf,
            pos_num = Inf,
            y_axis = y_axis,
            fc_cutoff = fc_cutoff,
            pval_cutoff = pval_cutoff
          )
        
        ##获取基因
        data$Gene
        
      })
    })
  })

### de_list

de_alldat_list <- 
  all_data %>% 
  map2(names(.), function(dat, name){
    counts_ana <- dat[['counts_ana']] 
    
    meta_ana_gene <- dat[['meta_ana_gene']]
    meta_ana_ROI <- dat[['meta_ana_ROI']]
    
    de_list <- dat[["de_list_single"]]
    
    map(de_list, function(de){
      map(de, function(dat){
        data <- dat %>%
          volcanoTable(
            full = TRUE,
            # neg_num = Inf,
            # pos_num = Inf,
            y_axis = y_axis,
            fc_cutoff = fc_cutoff,
            pval_cutoff = pval_cutoff
          )
        
        ##获取数据
        data
        
      })
    })
  })
  