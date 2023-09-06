# 0 GeomxTools ------------------------------------------------------------
readInitialData <- function(path, 
                            combine = TRUE,
                            SegmentSheetName = 'SegmentProperties',
                            ProbeSheetname   = 'BioProbeCountMatrix'){
  # test path
  is_exist <- purrr::map_lgl(path, file.exists)
  if(length(which(!is_exist)) > 0){
    message(stringr::str_glue("File(s) {glue::glue_collapse(path[!is_exist], sep = ', ', last = 'and')} \
                              not exist and will be ignored!"))
    path <- path[is_exist]
  }
  
  if(length(path) < 1) stop("None of the file path exist!")
  
  if(length(path) < 2) {
    
    if(combine) message("Only one path exist, file will not combine!")
    combine <- FALSE
    
  }
  
  # read data
  dat <- map(path, function(p){
    SegmentProperties   <- readxl::read_excel(p,
                                              sheet = SegmentSheetName)
    BioProbeCountMatrix <- readxl::read_excel(p,
                                              sheet = ProbeSheetname)
    list(
      SegmentProperties   = SegmentProperties,
      BioProbeCountMatrix = BioProbeCountMatrix
    )
  })
  
  # combine?
  if(!combine) return(dat[[1]])
  
  dat <- purrr::transpose(dat)
  dat$SegmentProperties <- bind_rows(dat$SegmentProperties)
  
  # make sure identical order of probename
  ROI_label_prefix <- getROILabelPrefix(dat$SegmentProperties$SegmentDisplayName)
  dat$BioProbeCountMatrix <- map_at(dat$BioProbeCountMatrix, -1, function(x) dplyr::select(x, "ProbeName", starts_with(ROI_label_prefix)))
  dat$BioProbeCountMatrix  <- purrr::reduce(dat$BioProbeCountMatrix , function(x, y) full_join(x, y, by = 'ProbeName'))
  
  return(dat)
}

# bind

# get ROI label
getROILabelPrefix <- function(x){
  x <- as.character(x)
  x <- substr(x, 1, 5)
  
  x_most <- sort(table(x))[1]
  return(names(x_most))
}

# construct new GeomxDataSet object
newGeomxDataSet <- function(InitialData,
                            Module     = "Standard",
                            annotation = "GeoMx_Hs_CTA",
                            check      = FALSE,
                            dimLabels  = c("RTS_ID", "Sample_ID"),
                            MinArea    = 16000,
                            MinNuclei  = 200,
                            shiftedByOne = TRUE){
  
  SegmentProperties   <- InitialData$SegmentProperties
  BioProbeCountMatrix <- InitialData$BioProbeCountMatrix
  
  pheno <- SegmentProperties %>% 
    dplyr::select(`slide name` = SlideName, 
                  ScanLabel, 
                  ROILabel, 
                  ROIID, 
                  SegmentLabel, 
                  SegmentID,
                  SegmentDisplayName, 
                  area   = AOISurfaceArea, 
                  nuclei = AOINucleiCount,
                  ROICoordinateX, 
                  ROICoordinateY) %>% 
    mutate(area   = as.numeric(area), 
           nuclei = as.numeric(nuclei),
           ROICoordinateX = as.numeric(ROICoordinateX), 
           ROICoordinateY = as.numeric(ROICoordinateY)) %>%
    as.data.frame()
  rownames(pheno) <- pheno$SegmentDisplayName
  
  protocol <- SegmentProperties %>%
    dplyr::select(segmentDisplayName = SegmentDisplayName, 
                  Raw = RawReads, 
                  Aligned = AlignedReads, 
                  DeduplicatedReads,
                  Trimmed = TrimmedReads, 
                  Stitched = StitchedReads,
                  SequencingSaturation, 
                  SequencingSetID,
                  umiQ30 = UMIQ30, 
                  rtsQ30 = RTSQ30) %>%
    mutate(across(-c(segmentDisplayName, SequencingSetID), as.numeric)) %>%
    # mutate(NTC = 3000) %>% # 如果需要出现QCFlags.HighNTC, 则需要有NTC这一列
    as.data.frame()
  rownames(protocol) <- protocol$segmentDisplayName
  
  feature <- BioProbeCountMatrix %>%
    dplyr::select(ProbeName, 
                  ProbeDisplayName, 
                  TargetName, 
                  HUGOSymbol, 
                  Accessions,
                  GenomeBuild,
                  GenomicPosition, 
                  AnalyteType,
                  CodeClass, 
                  ProbePool, 
                  TargetGroup) %>%
    mutate(Module = Module,
           RTS_ID = ProbeName) %>%
    as.data.frame() %>%
    column_to_rownames("ProbeName")
  
  ROILabelPrefix <- getROILabelPrefix(pheno$`slide name`)
  if(nchar(ROILabelPrefix) < 1) stop('slide name is null or not in Segment data of Initial Data!')
  
  assay <- BioProbeCountMatrix %>%
    dplyr::select(ProbeName, starts_with(ROILabelPrefix)) %>%
    as.data.frame() %>%
    column_to_rownames("ProbeName")
  stopifnot(all(colnames(assay) == rownames(pheno)))
  
  experiment <- 
    Biobase::MIAME(name = "", 
                   other = list(AnalyteType = "RNA",
                                MinArea = 16000,
                                MinNuclei = 200,
                                shiftedByOne = shiftedByOne)) 
  
  # new object 
  InitalObject <-
    NanoStringGeoMxSet(assayData = as.matrix(assay),
                       phenoData = pheno %>% AnnotatedDataFrame(),
                       featureData = feature %>% AnnotatedDataFrame(),
                       experimentData = experiment,
                       annotation = annotation,
                       protocolData = protocol %>% AnnotatedDataFrame(),
                       check = check, 
                       dimLabels = dimLabels)
  InitalObject
}

# tools for suffix duplicated items
str_distinct <- function(x, sep = "_"){
  
  if(length(x) < 1) return(x)
  if(sum(duplicated(x)) < 1) return(x)
  
  dat <- data.frame(x = x)
  
  dat2 <-
    dat %>% 
    group_by(x) %>%
    mutate(suffix = as.character(1:n())) %>%
    mutate(suffix = ifelse(suffix == "1", "", str_c(sep, suffix))) %>%
    mutate(x2 = str_c(x, suffix)) 
  
  return(dat2$x2)
}

calculate_loq <- function(counts, 
                          gene_label = NA,
                          dataset = c("CTA", "WTA")){
  # browser()
  dataset <- match.arg(dataset)
  
  if(dataset == "CTA"){
    neg_probe <- "Negative Probe"
  }else{
    neg_probe <- "NegProbe-WTX"
  }
  
  if(!is.na(gene_label)){
    gene <- counts[, gene_label, drop = TRUE]
    counts <- counts[, -match(gene_label, colnames(counts))]
  }else{
    gene <- rownames(counts)
  }
  
  neg_idx <- which(gene %in% neg_probe)
  
  if(length(neg_idx) == 0){
    stop(str_glue("No negative probe found: {neg_probe}"))
  }
  
  if(length(neg_idx) == 1){
    message(str_glue("Only one negative probe found, return negative data..."))
    neg_counts <- t(counts[neg_idx, , drop=FALSE])
    return(data.frame(loq = neg_counts, row.names = rownames(neg_counts)))
  }
  
  loq <- counts[neg_idx, ] %>% 
    apply(2, function(x){
      x <- as.numeric(x)
      
      dat <- log(x)
      MEAN <-  mean(dat) %>% exp()
      SD <- sd(dat) %>% exp()
      
      MEAN * (SD ^ 2)
    })
  
  loq %>% 
    as.data.frame() %>% 
    `colnames<-`("loq")
  
}

# 1 QC --------------------------------------------------------------------
##获取gene的metadata
getMetaGene <- function(raw_data){

  raw_data %>% dplyr::select(!starts_with("MIPCC"))
  
}

##获取ROI的metadata
getMetaROI <- function(raw_data){
  raw_data %>% 
    dplyr::rename(Scan_ID = ScanLabel,
           ROI_ID = ROILabel,
           Segment_ID = SegmentLabel,
           Original_ID = SegmentDisplayName) 
}

## 验证anno的列是否填写内容
isNA <- function(df){
  
  if(is.data.frame(df) || is.matrix(df)){
    if(ncol(df) > 1) message("Data Frame should only have one column!")
  }
  
  ratio <- length(which(unlist(is.na(df))))/length(unlist(df))
  ifelse(ratio > 0.9, TRUE, FALSE)
} 


## 获取客户注释的列名
getUserColumnName <- function(anno_data, singleLevel = F){
  cols_name <- colnames(anno_data)
  cols_name <- setdiff(cols_name, c("No.",  "Scan_ID",  "ROI_ID",  "Segment_ID")) # 固定列，去除
  
  cols_name <- map_chr(cols_name,  # 不为NA的列返回名字，否则标记为"NA"
                       ~ifelse(!isNA(anno_data[, .]), ., "NA"))
  cols_name <- cols_name[!cols_name=="NA"] 
  # 去除只有一个水平的列
  if(singleLevel){
    idx <- map_lgl(cols_name,                       
                   ~ifelse(length(unique(anno_data[, ., drop=T])) == 1, TRUE, FALSE))
    cols_name <- cols_name[!idx]
  }
  cols_name
}

## ROI注释信息
annoMetaROI <- function(meta_ROI, anno_data){
  
  # 验证是否有patientID
  if (! "Patient_ID" %in% colnames(anno_data)) message("Annotation table should have column Patient_ID!")
  
  if(isNA(anno_data[,"Patient_ID", drop=F])){
    slide_levels <- anno_data$Scan_ID %>% { factor(., levels = unique(.)) } %>% as.numeric()
    anno_data <- mutate(anno_data, Patient_ID = paste0("Patient", slide_levels))
  }
  
  col_omit <- c("No.")
  # 筛选合并数据
  # Segment_ID不为空:AOI
  # if(F){
  if("Segment_ID" %in% colnames(anno_data) &&
     !isNA(anno_data[,"Segment_ID", drop=F])){
    
    anno_data <- anno_data[!map_lgl(anno_data, isNA)]
    anno_data <- dplyr::select(anno_data, -any_of(col_omit)) %>%
      mutate(Patient_ID = str_replace_na(Patient_ID, replacement = ""),
             ROI_ID = str_replace_na(ROI_ID, replacement = ""),
             Segment_ID = str_replace_na(Segment_ID, replacement = "")) %>%
      mutate(ID = str_c(Patient_ID, ROI_ID, Segment_ID, sep = "_")) %>%
      mutate(ID = str_replace(ID, "Geometric Segment", "GS") %>%
               str_to_upper() %>% 
               str_remove("\\+$"))
    
    meta_ROI <- left_join(meta_ROI, anno_data, by=c("Scan_ID", "ROI_ID", "Segment_ID"))
    if(any(is.na(meta_ROI$ID))) {
      missing_roi <- meta_ROI[is.na(meta_ROI$ID), c("Scan_ID", "ROI_ID", "Segment_ID")]
      missing_label <- with(missing_roi, str_c(Scan_ID, ROI_ID, Segment_ID, sep = " | "))
      
      stop(str_glue("{length(missing_label)} ROI(s) had no annotation: \n{str_c(missing_label, collapse = ', ')}"))
    }
    meta_ROI <- as.data.frame(meta_ROI)
    meta_ROI$ID <- str_distinct(meta_ROI$ID)
    rownames(meta_ROI) <- meta_ROI$ID
    return(meta_ROI)
  }
  # Segment_ID为空：ROI
  anno_data <- anno_data[!map_lgl(anno_data, isNA)]
  
  anno_data <- dplyr::select(anno_data, -any_of(col_omit)) %>%
    mutate(Patient_ID = str_replace_na(Patient_ID, replacement = ""),
           ROI_ID = str_replace_na(ROI_ID, replacement = ""))%>%
    mutate(ID = str_c(Patient_ID, ROI_ID, sep = "_")) %>%
    mutate(ID = str_replace(ID, "Geometric Segment", "GS") %>%
             str_to_upper() %>% 
             str_remove("\\+$"))
  
  
  # 自定义ROI ID后缀
  if("ID_Suffix" %in% colnames(anno_data)) {
    
    anno_data <- anno_data %>% 
      mutate(ID_Suffix = str_replace_na(ID_Suffix, replacement = "")) %>%
      mutate(ID = str_c(ID, ID_Suffix, sep = "_"))
  }
  
  meta_ROI <- left_join(meta_ROI, anno_data, by=c("Scan_ID", "ROI_ID"))
  if(any(is.na(meta_ROI$ID))) {
    missing_roi <- meta_ROI[is.na(meta_ROI$ID), c("Scan_ID", "ROI_ID", "Segment_ID")]
    missing_label <- with(missing_roi, str_c(Scan_ID, ROI_ID, Segment_ID, sep = " | "))
    
    stop(str_glue("{length(missing_label)} ROI(s) had no annotation: \n{str_c(missing_label, collapse = ', ')}"))
  }
  meta_ROI <- as.data.frame(meta_ROI)
  meta_ROI$ID <- str_distinct(meta_ROI$ID)
  rownames(meta_ROI) <- meta_ROI$ID
  meta_ROI
}

##获取数据count
getCount <- function(raw_data, meta_gene, meta_ROI, version = "1.0"){
  # 搜索ROI信息和Counts之间的分界线
  idx <- which(raw_data[, 1] == "#Probe Group")
  if(raw_data[idx, 2] != "#Analyte type") message(paste0("Is #Analyte type at ", idx-1,"st row?"))
  if(raw_data[idx, 3] != "#CodeClass") message(paste0("Is #CodeClass at ", idx,"st row?"))
  if(raw_data[idx, 4] != "ProbeName (display name)") message(paste0("Is ProbeName (display name) at ", idx-1,"st row?"))
  # 
  if(version == "1.0"){
    raw_data %>%
      dplyr::select(-c(1:4)) %>%
      slice( (idx+1):nrow(.) ) %>% 
      map_dfc(as.numeric) %>%{
        . <- as.data.frame(.)
        colnames(.) <- meta_ROI$ID
        rownames(.) <- meta_gene$Probe_name
        . }
  }else{
    message("Version Not Supported!")
  }
}

# Compute Positive Norm Factor
computerERCC <- function(count_data, meta_gene){
  ERCC_pos_gene <- count_data[meta_gene$Probe_name =="HYB-POS", ]
  ERCC_pos_gene <- t(ERCC_pos_gene)
  ERCC_factor <- mean(unlist(ERCC_pos_gene))/ERCC_pos_gene
  colnames(ERCC_factor) <- "Positive Norm Factor"
  as.data.frame(ERCC_factor)
}

# ERCC tranformation
normERCC <- function(count_data, ERCC_factor){
  rowN <- rownames(count_data)
  colN <- colnames(count_data)
  
  lapply(seq_along(count_data), function(x){
    count_data[,x] * ERCC_factor[x, 1]
  })  %>% do.call(cbind, .) %>% {
    rownames(.) <- rowN
    colnames(.) <- colN
    .
  } -> QC_data
  as.data.frame(QC_data)
}

# Normalizaiton
Normalization <- function(count_data, meta_gene){
  HK_idx <- meta_gene$Probe_name %in% findHKGene(meta_gene) 
  HK_counts <- count_data[HK_idx, ]
  HK_GeoMean_counts <- HK_counts %>% map_dfc(geoMean) %>% unlist %>% mean
  HK_factor <- map_dfc(HK_counts, ~ HK_GeoMean_counts/geoMean(.)) %>% t
  
  rowN <- rownames(count_data)
  colN <- colnames(count_data)
  
  lapply(seq_along(count_data),function(x){
    count_data[,x] * HK_factor[x,1]
  })  %>% do.call(cbind, .) %>% {
    rownames(.) <- rowN
    colnames(.) <- colN
    .
  } -> Norm_data
  
  Norm_data
}


# 1.1 HouseKeeping counts distribution ------------------------------------

##找到HouseKeeping Gene
findHKGene <- function(meta_gene, assay = "RNA"){
  if(assay == "RNA" || assay == "Protein"){
    meta_gene %>% 
      dplyr::filter(Code_class == "Control") %>%
      dplyr::select(Probe_name) %>%
      unlist()
  }else if(assay == "Protein"){
    NULL #待完善
    # 由于RNA和Protein的HK gene取法相同，不再需要特殊指定assay参数。
  }else{
    message("Assay should only be 'RNA' or 'Protein'")
  }

}

##找到Background Gene
findBGGene <- function(meta_gene, assay = "RNA"){
  if(assay == "RNA" || assay == "Protein" ){
    meta_gene %>% 
      dplyr::filter(Analyte_type != "SpikeIn") %>%
      dplyr::filter(Code_class == "Negative") %>%
      dplyr::select(Probe_name) %>%
      unlist()
  }else{
    message("Assay should only be 'RNA' or 'Protein'")
  }
  
}
##找到Negtive Gene
findNegGene <- function(meta_gene, assay = "RNA"){
  if(assay == "RNA" || assay == "Protein" ){
    meta_gene %>% 
      dplyr::filter(Analyte_type == "SpikeIn") %>%
      dplyr::filter(Code_class == "Negative") %>%
      dplyr::select(Probe_name) %>%
      unlist()
  }else{
    message("Assay should only be 'RNA' or 'Protein'")
  }
  
}

### 绘图函数
#1. 直方图
.histPlot <- function(data, geneclass){
  data %>% ggplot(aes(x=Geo_Mean)) +
    geom_histogram(color="black", fill="grey") +
    theme_minimal()+
    theme(panel.grid = element_blank()) +
    xlim(0, max(data$Geo_Mean) * 1.2) +
    labs(x = paste("AOI/ROIs'log2", geneclass, "Geomean (Counts)"),
         y = "Frequency")
}
#2. 分布图
.distbtPlot <- function(data, column, geneclass){
  ggplot(data, aes(x = reorder(1:nrow(data), Geo_Mean),
                y = Geo_Mean, 
                fill = .data[[column]]) )+
    geom_col(color = "#7F7F7F", width = 0.8) +
    theme_minimal()+
    labs(x = "Individual AOI/ROIs",
         y = paste("Log2", geneclass, "Geomean (Counts)"))+
    theme(axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          legend.position = "top",
          panel.grid = element_blank()) +
    scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(n_distinct(data[[column]])))
    # scale_fill_brewer(palette = "Dark2")
}

countDistributionPlot <- function(count_data, meta_gene, meta_ROI, genelist, 
                                  columns, geneclass = "HK"){
  ### 构造数据
  count_df <- count_data %>%
    dplyr::filter(rownames(.) %in% genelist) %>% 
    map_dbl( ~log2(geoMean(.)) ) %>% 
    data.frame(Geo_Mean = .) %>%
    merge(meta_ROI, by = "row.names")
  
  ### 绘图
  .histPlot(count_df, geneclass = geneclass) -> p1
  map(columns,   ~ .distbtPlot(count_df, column = ., geneclass = geneclass)) -> p2
  
  ### 拼图
  require(patchwork)

  plot_list <- vector(length = length(p2)+1, mode = "list")
  plot_list[[1]] <- p1
  for(i in seq_along(p2) ) plot_list[[i+1]] <- p2[[i]]
  patchwork::wrap_plots(plot_list, ncol = 1)

}

# 1.2 HK_bg_cor ------------------------------------------------------
.corScatterPlot <- function(data, column){
  data %>% ggplot(aes(x = HK, y = BG, color = .data[[column]]))+
    geom_point(size = 5, alpha=0.8) +
    theme_bw()+
    # scale_color_brewer(palette = "Dark2")+
    scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(n_distinct(data[[column]]))) +
    theme(legend.position = "top",
          panel.grid = element_blank())+
    labs(x = "Control Geomean (Counts)",
         y = "Background Geomean (Counts)")
}

signalBackgroundPlot <- function(count_data, meta_gene, meta_ROI, columns){
  # HK_count
  HK_gene <- findHKGene(meta_gene)
  HK_count <- count_data %>%
    dplyr::filter(rownames(.) %in% HK_gene) %>% 
    map_dbl( ~log2(geoMean(.)) ) %>% 
    data.frame(HK = .)
  
  # BG_count
  BG_gene <- findBGGene(meta_gene)
  BG_count <- count_data %>%
    dplyr::filter(rownames(.) %in% BG_gene) %>% 
    map_dbl( ~log2(geoMean(.)) ) %>% 
    data.frame(BG = .)
  
  # count
  HK_BG_conunt <- bind_cols(HK_count, BG_count, meta_ROI)
  
  # 绘图
  # .corScatterPlot(HK_BG_conunt, columns)
  plot_list <- map(columns,   ~ .corScatterPlot(HK_BG_conunt, .)) 
  patchwork::wrap_plots(plot_list, ncol = 2)
}

# 1.3 Signal-to-background Ratio ------------------------------------------


SNRPlot <- function(count_data, meta_gene, meta_ROI){
  # SNR为进行绘图的count矩阵
  BG_gene <- findBGGene(meta_gene)
  idx <- match(BG_gene, meta_gene$Probe_name)
  
  SNR <- count_data %>% 
    apply(2, function(x){
      x <- x/geoMean(x[idx])
      log2(x)
    }) %>% t %>% 
    as_tibble()
  
  SNR <- pivot_longer(SNR, cols = everything(),
               names_to = "Gene",
               values_to = "Counts") 
  
  SNR <- SNR %>% { 
    seq_gene <- reorder(.$Gene, .$Counts)  %>% levels() # 按照counts值对Gene排序
    levels_gene <- c(sort(BG_gene), setdiff(seq_gene, BG_gene)) # 将BG基因的顺序排到最前面
    mutate(., Gene = factor(Gene, levels = levels_gene)) # 将Gene调整为特定顺序的factor
  }
  # 绘图
  SNR %>% {
    ggplot(., aes(x=Gene, y=Counts)) +
      geom_boxplot(outlier.shape = NA)+
      geom_jitter(size = 1, color = "#8B0202", width = 0.2, alpha = 0.6) +
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            panel.grid = element_blank()) +
      labs(x=NULL, y="Log2 Signal-to-Background Ratio") +
      geom_hline(yintercept = 0)+
      geom_vline(xintercept = length(BG_gene)+0.5, linetype = 2 )
  }
}

# 1.4 correlation of HK or BG ------------------------------------------
corrMaxtrixPlot <- function(count_data, meta_gene, meta_ROI, genelist, column, 
                            assay = "RNA", geneclass = "HK"){
  ## 构造数据
  rowN <- rownames(count_data)
  count_data <- count_data %>% map_dfc(log2) %>% 
    as.data.frame() %>%
    {  rownames(.) <- rowN; . } %>%
    t %>% as.data.frame() %>%
    merge(meta_ROI, by= "row.names") %>%
    column_to_rownames("Row.names")
  
  ## 去除只出现一次的水平值
  # ggpairs绘图时，调用cor.test，要求每个水平的case数大于2
  bad_levels <- count_data[[column]] %>% table %>% .[. < 3] %>% names()
  if(length(bad_levels) > 0)   count_data <- count_data %>% dplyr::filter(! .data[[column]] %in% bad_levels)
  
  ## 绘图
  GGally::ggpairs(count_data, columns = genelist,
                  mapping = aes(color = .data[[column]]),
                  upper = list(continuous = GGally::wrap("cor", align_percent = 0.5)),
                  diag = list(continuous=GGally::wrap("densityDiag", alpha=0.6))
                  # diag = "blank",
                  # axisLabels = "internal"
                  ) +
    ggtitle(paste(assay, "counts\ncolored by", column)) +
    theme_bw() +
    theme(title = element_text(size = 15))
}

# 1.5 total correlation ------------------------------------------
corrMaxtrixTotalPlot <- function(count_data, meta_gene, meta_ROI, column, assay = "RNA", Quantile = 0.75,
                                 HK_except = NULL, Neg_except = NULL){
  # HK_count
  HK_gene <- findHKGene(meta_gene) %>% setdiff(HK_except)
  HK_count <- count_data %>%
    dplyr::filter(rownames(.) %in% HK_gene) %>% 
    map_dbl( ~geoMean(.) ) %>% log2 %>%
    data.frame(HK = .)
  
  # Neg_count
  
  Neg_gene <- findNegGene(meta_gene) %>% setdiff(Neg_except)
  Neg_count <- count_data %>%
    dplyr::filter(rownames(.) %in% Neg_gene) %>% 
    map_dbl( ~geoMean(.) ) %>% log2 %>%
    data.frame(NegProbe = .)
  
  # Q3 data
  topGene <- count_data %>%
    apply(1, geoMean) %>% log2 %>%{
      Q3 <- quantile(., probs = Quantile)
      .[. >= Q3]
    } %>% names()
  
  topCounts <- count_data %>%
    # dplyr::filter(rownames(.) %in% topGene) %>%
    # map_dbl( ~geoMean(.) ) %>% log2 %>%
    # data.frame(`Q3` = ., check.names = FALSE)
    apply(2, function(x)quantile(x, Quantile))%>% 
    log2 %>%
    data.frame(Q3 = ., check.names = FALSE)

  # count
  data_count <- merge(HK_count, Neg_count, by = "row.names") %>% column_to_rownames("Row.names")
  data_count <- merge(data_count, topCounts, by = "row.names") %>% column_to_rownames("Row.names")
  data_count <- merge(data_count, meta_ROI, by = "row.names")
  data_count <- data_count %>% mutate(Area=as.numeric(AOISurfaceArea)) %>%
    mutate(Nuclei_count = as.numeric(AOINucleiCount))
  
  ## 去除只出现一次的水平值
  bad_levels <- data_count[[column]] %>% table %>% .[. < 3] %>% names()
  if(length(bad_levels) > 0)   data_count <- data_count %>% dplyr::filter(! .data[[column]] %in% bad_levels)
  
  
  # 绘图
  cor_columns <- c("HK", "NegProbe","Q3",
                   "Area", "Nuclei_count")
  GGally::ggpairs(data_count, columns = cor_columns,
                  mapping = aes(color = .data[[column]]),
                  upper = list(continuous = GGally::wrap("cor", align_percent = 0.5)),
                  diag = list(continuous = GGally::wrap("densityDiag", alpha=0.6)) ) +
    ggtitle(paste(assay, "counts\ncolored by", column)) +
    theme_bw() +
    theme(title = element_text(size = 15))
}

# 1.6 ROI filter ------------------------------------------
roiFilter <- function(meta_ROI, assay = "RNA", results = "df"){
  if (!all(c("Row.names", "Positive Norm Factor") %in% colnames(meta_ROI))) 
    message("Has ERCC factor computation done?")
  
  if (assay == "RNA" ){
    Area_limit = 16000
    Nuclei_limit = 200
  }else if(assay == "Protein"){
    Area_limit = 1600
    Nuclei_limit = 20
  }else{
    message("Wrong assay name provided!")
  }
  
  ROI_filter <-
    meta_ROI %>% dplyr::select(Row.names, Area, Nuclei_count, `Positive Norm Factor`) %>%
    mutate(Postive_factor = 
             case_when(`Positive Norm Factor` > 4  ~ "Fail",
                       `Positive Norm Factor` < 0.2 ~ "Fail",
                       `Positive Norm Factor` > 3  ~ "Caution",
                       `Positive Norm Factor` < 0.3 ~ "Caution",
                       TRUE ~ "ERCCPass")) %>%
    mutate(Group = ifelse(Postive_factor != "Fail" &
                            Area > .env[["Area_limit"]] &
                            Nuclei_count > .env[["Nuclei_limit"]], "Pass", "Fail"))
  
  p <- ROI_filter %>% ggplot(aes(x = Area, y = Nuclei_count, 
                            size = Postive_factor, fill = Group)) + 
    geom_point(alpha = 0.6, shape = 21) +
    theme_bw() +
    geom_vline(xintercept = Area_limit, linetype = 2, alpha = 0.5) +
    geom_hline(yintercept = Nuclei_limit, linetype = 2, alpha = 0.5) +
    annotate("text", x = Inf, y = Nuclei_limit, label = glue::glue("Nuclei counts = {Nuclei_limit}  "), hjust = 1, vjust = -1) +   
    annotate("text", x = Area_limit, y = Inf, label = glue::glue("  Area = {Area_limit}"), hjust = 0, vjust = 1.5) + 
    theme(axis.text = element_blank(),
          panel.grid = element_blank(),
          legend.position = "top",
          legend.box = "vertical",
          legend.spacing = unit(0.1, 'cm')) +
    scale_fill_manual(values = c("Pass" = "blue", "Fail" = "red"),
                      breaks = c("Pass", "Fail")) +
    scale_size_manual(values = c("ERCCPass" = 6, "Caution" = 4, "Fail" = 2), 
                      breaks = c("ERCCPass", "Caution", "Fail")) +
    labs(y = "Nuclei counts") +
    guides( fill = guide_legend(override.aes = list(size = 4))) +
    ggrepel::geom_text_repel(data = dplyr::filter(ROI_filter, Group == "Fail"),
                             aes(label = Row.names), size = 3) +
    scale_y_log10()+
    scale_x_log10()
  
  if(results == "df"){
    return(ROI_filter)
  }else if(results == "plot"){
    return(p)
  }else{
    message("Parameter results should only be 'df' or 'plot'!")
  }
  
}
# 2 Analysis --------------------------------------------------------------
##去除Control样本
removeControl <- function(count_data, meta_gene){
  # 筛选数据,只要Endogenous数据
  idx <- which(meta_gene$Code_class == "Endogenous")
  df <- count_data[idx, ]
  
  rownames(df) <- meta_gene[idx, "Probe_name", drop=T]
}

# 2.1 volcano ---

## 计算几何平均值函数
geoMean <- function(x){
  # prod <- max(cumprod(x))
  # prod^(1/length(x))
  x <- log(x)
  exp(mean(x))
}


## 单分组差异分析
singleVolcano <- function(df, meta_gene, meta_ROI, variableA, p.adj.method = "BH", method = "t"){
  ## method: t wilcox DEseq2 edgeR limma
  ## 目前不支持edgeR limma
  
  ## 约定标准：不需要进行分析的levels在annotation中标记为Others
  #  分析时会去掉Others的ROI
  stopifnot(variableA %in% colnames(meta_ROI))
  
  idx <- meta_ROI[, variableA, drop=TRUE] != "Others"
  df <- df[, idx, drop=FALSE]
  meta_ROI <- meta_ROI[idx, ,drop=FALSE]

  ## 差异分析的总水平
  subGroup_levels <- if(is.factor(meta_ROI[,variableA])){
    levels(meta_ROI[,variableA])
  }else{
    unique(meta_ROI[,variableA])
  }

  ## 差异分析变量转为因子变量
  subGroup <- if(is.factor(meta_ROI[,variableA])){
    meta_ROI[, variableA, drop=T]
  }else{
    as.factor(meta_ROI[, variableA, drop=T])
  }

  ## 水平小于1则不进行统计检验
  if( length(subGroup_levels) <= 1){
    de_list <- NULL
  }else{
    subGroup_combn <- combn(1:length(subGroup_levels), 2) # 每列是所有的因子组合
    
    ##差异分析
    de_list <- vector(mode = "list") # 差异表达结果列表
    
    ###### <DESeq2方法> ######
    if(tolower(method) == "deseq2"){
      # original data
      cts <- round(df)
      colData <- meta_ROI[, variableA, drop = FALSE]
      cts <- cts[, rownames(colData)]
      
      # 防止colData里面有非法R命名，先命名为X1 X2...
      if(ncol(colData) > 1) message("Multi group deseq2: old name would not be changed...")
      lev_names <- colData %>% unlist(use.names = FALSE) %>% as.character()
      old2new <- unique(lev_names)
      names(old2new) <- paste0("X", 1:length(old2new))
      new_names <- fct_recode(lev_names, !!!old2new)
      colData[[1]] <- new_names
      
      # 构造object
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts, colData = colData,
                                    design = as.formula(paste("~", variableA)))
      # filter & fit
      keep <- which(rowMeans(DESeq2::counts(dds)) >= 1 )
      dds <- dds[keep,]
      
      # fit: 需要样本数大于拟合变量的水平数
      if(ncol(DESeq2::counts(dds)) <= length(subGroup_levels)) return(NULL)
      dds <- DESeq2::DESeq(dds)
      
      # 数据
      for (i in seq_along(as.data.frame(subGroup_combn))){
        ## 差异比较的两组
        de_group_levels <- subGroup_combn[ , i, drop=T]
        de_group_levels <- subGroup_levels[de_group_levels]
        de_name <- paste(de_group_levels[2], de_group_levels[1], sep = "_vs_") #差异比较组别命名
        
        ## 进行比较的两组数据，如果数量小于2则不进行统计检验
        de_group_sub <- subGroup[subGroup == de_group_levels[1] |
                                   subGroup == de_group_levels[2] ]
        de_group_sub <- as.character(de_group_sub) # 因子变量的table函数会保留所有的因子水平的统计
        if( all(table(de_group_sub) <= 1) ) next()
        
        # 将label修成成DESeq2中使用的名字
        de_group_levels <- de_group_levels %>% fct_recode(!!!old2new) %>%as.character()
        
        
        ## results函数获取差异数据
        res <- DESeq2::results(dds, contrast = c(as.character(variableA), 
                                                 de_group_levels[2], 
                                                 de_group_levels[1]) )
        res <- as.data.frame(res)
        
        # Gene,	Pvalue,	Log2FC,	Adjusted.pvalue,	-log10 pvalue,	-log10 adjusted pvalue
        de_df <- res %>% rownames_to_column(var = "Gene") %>% 
          dplyr::rename(Log2FC = log2FoldChange,  
                 Pvalue = pvalue,
                 Adjusted.pvalue = padj) %>%
          mutate(`-log10 pvalue` = -log10(Pvalue)) %>%
          mutate(`-log10 adjusted pvalue` = -log10(Adjusted.pvalue))
        
        de_list[[de_name]] <- de_df
      } #END for "for"
      return(de_list)
    } #END for "if"
      
    
    
    ###### <t.test/wilcox.test方法> ######
    
    for (i in seq_along(as.data.frame(subGroup_combn))){
      ## 差异比较的两组
      de_group_levels <- subGroup_combn[ , i, drop=T]
      de_group_levels <- subGroup_levels[de_group_levels]
      de_name <- paste(de_group_levels[2], de_group_levels[1], sep = "_vs_") #差异比较组别命名
      
      ## 挑选出差异比较的cases
      idx_sub <- c(which(subGroup == de_group_levels[1] |
                           subGroup == de_group_levels[2]) )
      de_group_sub <- subGroup[idx_sub]
      de_group_sub <- factor(de_group_sub, levels = de_group_levels)
      
      ## 进行比较的两组数据，如果数量小于2则不进行统计检验
      if( any(table(de_group_sub) <= 1) ) next()
      
      ## t.test计算p值
      df_sub <- df[, idx_sub]
      pval <- apply(df_sub, 1, function(x){
        pval <- t.test(x~de_group_sub)$p.value
      }) %>% cbind() %>%
        as.data.frame() %>%
        { rownames(.) <- rownames(df_sub); . } %>%
        rownames_to_column(var = "Gene") %>%
        dplyr::rename(Pvalue = 2)
      
      ## 计算p.adjust
      p.adj <- p.adjust(pval$Pvalue, method = p.adj.method)
      
      ## 计算log2FC
      log2df <- apply(df, 1, function(x){ # df,不是def_sub
        # subGroup_level <- levels(meta_ROI[,variableA])
        Test <- x[meta_ROI[,variableA] == de_group_levels[2] ]
        Control <- x[meta_ROI[,variableA] == de_group_levels[1] ]
        log2df <- mean(Test)/mean(Control)
        log2(log2df)
      })
      
      ## 合并p,p_adj,log2FC，并计算-log10 p等数据
      # Gene,	Pvalue,	Log2FC,	Adjusted.pvalue,	-log10 pvalue,	-log10 adjusted pvalue
      de_df <- transform(pval, Log2FC = log2df, Adjusted.pvalue = p.adj) %>%
        mutate(`-log10 pvalue` = -log10(Pvalue)) %>%
        mutate(`-log10 adjusted pvalue` = -log10(Adjusted.pvalue))
      
      de_list[[de_name]] <- de_df
    } # End for "for"
    
  }#End for "else"
  

  return(de_list)
}

## 双分组差异分析
doubleVolcano <- function(df, meta_gene, meta_ROI, variableA, 
                              variableB, ...){
  
  majorGroup_levels <- if(is.factor(meta_ROI[,variableB, drop=T])){
    levels(meta_ROI[,variableB, drop=T])
  }else{
    unique(meta_ROI[,variableB, drop=T]) 
  }
  
  df_list <- map(majorGroup_levels, function(x){
    # 筛选数据
    idx <- which(meta_ROI[,variableB, drop=T] == x)
    df_n <- df[,idx,drop=FALSE]
    
    meta_ROI_n <- meta_ROI[idx, ,drop=FALSE]
    # 分析
    singleVolcano(df_n, meta_gene, meta_ROI_n, variableA, ...)
  })
  names(df_list) <- majorGroup_levels
  df_list
}

## 差异分析主函数
volcanoAnalysis <- function(count_data, meta_gene, meta_ROI, 
                            variableA, variableB = NULL, ...){
  # 筛选数据,只要Endogenous数据
  idx <- which(meta_gene$Code_class == "Endogenous")
  df <- count_data[idx, ]
  
  df <- as.data.frame(df)
  rownames(df) <- meta_gene[idx, "Probe_name", drop=T]
  
  # 差异分析 
  if(is.null(variableB)) {
    singleVolcano(df = df, meta_gene = meta_gene, meta_ROI = meta_ROI, variableA = variableA, ...)
  }else{
    doubleVolcano(df, meta_gene, meta_ROI, variableA, variableB, ...)
  }
}

## volcano绘图函数
volcanoPlot <-
  function(data,
           title = NULL,
           subtitle = NULL,
           y_axis = "padj",
           label_pos_num = 5,
           label_neg_num = 5,
           pval_cutoff = 0.05,
           fc_cutoff = 2,
           pointsize = 5,
           pointalpha = 0.6,
           color = NULL,
           ...) {
      
    if(is.null(data)) return(NULL)
    
    if(is.null(color)){
      color <- c("Negative" = "blue",
                 'NegOnlyP' = "#B3B3FF",
                 "Non-sig"  = "gray",
                 'PosOnlyP' = "#FFB3B3",
                 "Positive" = "red")
    } 
    
    pval_cutoff = -log10(pval_cutoff)
    fc_cutoff   = log2(fc_cutoff)
  
    
  ## y_axis控制使用p值还是调整p值的-log10进行绘图
  y_axis <- case_when(y_axis == "padj" ~ "-log10 adjusted pvalue" ,
                      y_axis == "pval" ~ "-log10 pvalue",
                      TRUE ~ "illegal")
  if(y_axis == "illegal") message("parameter y_axis should be pval or padj!")
  
  ## x轴
  x_axis <- c('Log2', 'Log2FC')
  if(!any(x_axis %in% colnames(data)) ) stop('Please check log2fc`name in DiffExpr data!') 
  x_axis <- x_axis[x_axis %in% colnames(data)][1]
  
  
  ##数据变换，添加分组信息
  data <- data %>%
    mutate(Group = case_when( .data[[y_axis]] > pval_cutoff & .data[[x_axis]] >= fc_cutoff ~ "Positive",
                              .data[[y_axis]] > pval_cutoff & .data[[x_axis]] <  fc_cutoff & .data[[x_axis]] > 0 ~ "PosOnlyP",
                              .data[[y_axis]] > pval_cutoff & .data[[x_axis]] < -fc_cutoff ~ "Negative",
                              .data[[y_axis]] > pval_cutoff & .data[[x_axis]] > -fc_cutoff & .data[[x_axis]] < 0 ~ "NegOnlyP",
                              TRUE ~ "Non-sig")) %>%
    mutate(Group = factor(Group, levels = c("Negative", "NegOnlyP", "Non-sig", "PosOnlyP", "Positive")))
  
  #########绘图:p1 p2，p1和p2的代码基本上一样：logFC取反，subtitle取反，左右label取反
  # 类似于：
  # p1: Tumor_vs_Stroma
  # p2: Stroma_vs_Tumor
  
  ##p1
  data %>% {
    # 需要标注的点：Posive和Negative各前5个点
    data_pos <- dplyr::filter(., Group == "Positive") %>% 
      top_n(label_pos_num,  -Adjusted.pvalue) 
    data_neg <- dplyr::filter(., Group == "Negative") %>%  
      top_n(label_neg_num,  -Adjusted.pvalue) 
    data_lab <- bind_rows(data_pos, data_neg)
    
    # 火山图左右注释
    if(is.null(subtitle)){
      left_text <- NULL
      right_text <- NULL
    }else{
      subtitle %>% str_split(pattern = "_vs_", simplify = T) %>%{
        # 结果为一行二列矩阵
        # 第一个元素标记在火山图右侧
        # 第二个元素标记在火山图左侧
        left_text <<- paste0("  ",  .[1,2])
        right_text <<- paste0(.[1,1], "  ")
      }
    }
    
    
    # 绘图
    ggplot(., aes(x = .data[[x_axis]], y = .data[[y_axis]])) +
      geom_point(
        aes(color = Group, size = Group),
        # size = pointsize,
        alpha = pointalpha,
        show.legend = FALSE
      ) +
      theme_bw() +
      theme(...) +
      theme(legend.position = "null") +
      geom_vline(xintercept = c(-fc_cutoff, fc_cutoff),linetype=2)+
      geom_hline(yintercept = pval_cutoff, linetype=2) +
      ggrepel::geom_text_repel(
        data = data_lab,
        mapping = aes(label=Gene) ) +
      scale_color_manual(values = color)  +
      scale_size_manual(values = c(
        Positive = pointsize,
        Negative = pointsize,
        PosOnlyP = max(1, round(pointsize/2, 1)),
        NegOnlyP = max(1, round(pointsize/2, 1)),
        "Non-sig" = 1
      ))+
      labs(title = title, subtitle = subtitle,
           x = "log2FC") +
      annotate("text", label = right_text, 
               x = Inf, y = 0, 
               size =5, colour = "black", 
               hjust="right") +
      annotate("text", label = left_text, 
               x = -Inf, y = 0, 
               size = 5, colour = "black", 
               hjust="left") -> p
    
    ##调整坐标轴范围
    y_lim <- max(data[[y_axis]], na.rm = TRUE)
    if(y_lim < 5) p <- p + ylim(NA, 5)
    
    x_lim <- range(data[[x_axis]], na.rm = TRUE)
    if(x_lim[1] > -2) p <- p + xlim(-2, NA)
    if(x_lim[2] < 2) p <- p + xlim(NA, 2)
    if(x_lim[1] > -2 && x_lim[2] < 2) p <- p + xlim(-2, 2)
    
    p1 <<- p
  }   
  
  ##p2
  temp <- color['Negative']
  color['Negative'] <- color['Positive']
  color['Positive'] <- temp
  
  temp <- color['NegOnlyP']
  color['NegOnlyP'] <- color['PosOnlyP']
  color['PosOnlyP'] <- temp
  
  
  data %>% {
    # 需要标注的点：Posive和Negative各前5个点
    data_pos <- dplyr::filter(., Group == "Positive") %>% 
      top_n(label_pos_num,  -Adjusted.pvalue) 
    data_neg <- dplyr::filter(., Group == "Negative") %>%  
      top_n(label_neg_num,  -Adjusted.pvalue) 
    data_lab <- bind_rows(data_pos, data_neg)
    
    # 火山图左右注释
    if(is.null(subtitle)){
      left_text <- NULL
      right_text <- NULL
    }else{
      subtitle %>% str_split(pattern = "_vs_", simplify = T) %>%{ # 与p1图顺序相反
        # 结果为一行二列矩阵
        # 第一个元素标记在火山图左侧left
        # 第二个元素标记在火山图右侧right
        left_text <<- paste0("  ",  .[1,1])
        right_text <<- paste0(.[1,2], "  ")
        subtitle <<- paste0(.[1,2], "_vs_", .[1,1])
      }
    }
    
    
    # 绘图
    ggplot(., aes(x = - .data[[x_axis]], y = .data[[y_axis]])) +  # Log2FC取负号
      geom_point(aes(color = Group, size = Group),
                 # size = pointsize,
                 alpha = pointalpha,
                 show.legend = FALSE) +
      theme_bw() +
      theme(...) +
      theme(legend.position = "null") +
      geom_vline(xintercept = c(-fc_cutoff,fc_cutoff),linetype=2)+
      geom_hline(yintercept = pval_cutoff, linetype=2) +
      ggrepel::geom_text_repel(
        data = data_lab,
        mapping = aes(label=Gene) ) +
      scale_color_manual(values = color)  +
      scale_size_manual(values = c(
        Positive = pointsize,
        Negative = pointsize,
        PosOnlyP = max(1, round(pointsize/2, 1)),
        NegOnlyP = max(1, round(pointsize/2, 1)),
        "Non-sig" = 1
      ))+
      labs(title = title, subtitle = subtitle,
           x = "log2FC") +
      annotate("text", label = right_text, 
               x = Inf, y = 0, 
               size =5, colour = "black", 
               hjust="right") +
      annotate("text", label = left_text, 
               x = -Inf, y = 0, 
               size = 5, colour = "black", 
               hjust="left") -> p
    
    ##调整坐标轴范围
    y_lim <- max(data[[y_axis]], na.rm = TRUE)
    if(y_lim < 5) p <- p + ylim(NA, 5)
    
    x_lim <- range(- data[[x_axis]], na.rm = TRUE) # Log2FC取负号
    if(x_lim[1] > -2) p <- p + xlim(-2, NA)
    if(x_lim[2] < 2) p <- p + xlim(NA, 2)
    if(x_lim[1] > -2 && x_lim[2] < 2) p <- p + xlim(-2, 2)
    
    p2 <<- p
  }
  
  return(patchwork::wrap_plots(p1, p2, nrow = 1) )
}


## volcano绘图函数
multiVolcanoPlot <- function(plot_data, User_columns){
  names(plot_data) <- User_columns
  map(seq_along(plot_data), function(data){
    enframe <- plot_data[[data]]
    col_name <- colnames(enframe)
    name <- col_name[1]
    value <- col_name[2]
    
    map2(enframe[[name]], enframe[[value]], function(x,y){
      volcanoPlot(data = y, title = names(plot_data)[data], subtitle = x)
    })
  }) -> plot_list 

  names(plot_list) <- User_columns

  map(plot_list, function(x){
      width = 8
      height = 6 * length(x)
      pdf(paste0(names(x),".pdf"), height = height, width = width)
      patchwork::wrap_plots(x, ncol = 1)
      dev.off()
  })
}

## 差异表达显著基因表
volcanoTable <- function(data,
                         pos_num = 5,
                         neg_num = 5,
                         full = FALSE,
                         y_axis = "padj",
                         pval_cutoff = 0.05,
                         fc_cutoff = 2,
                         DE_groups = NULL,
                         Major_Group = NULL
){
  
    if(is.null(data)) return(NULL)
    
    pval_cutoff = -log10(pval_cutoff)
    fc_cutoff   = log2(fc_cutoff)
  
    ## y_axis控制使用p值还是调整p值返回差异基因
    y_axis <- case_when(y_axis == "padj" ~ "-log10 adjusted pvalue" ,
                        y_axis == "pval" ~ "-log10 pvalue",
                        TRUE ~ "illegal")
    if(y_axis == "illegal") message("parameter y_axis should be pval or padj!")
    
    
    ## x轴
    x_axis <- c('Log2', 'Log2FC')
    if(!any(x_axis %in% colnames(data)) ) stop('Please check log2fc`name in DiffExpr data!') 
    x_axis <- x_axis[x_axis %in% colnames(data)][1]
    
    
    ##数据变换，添加分组信息：上调、下调及不显著
    data <- data %>%
      mutate(Group = case_when( .data[[y_axis]] > pval_cutoff & .data[[x_axis]] > fc_cutoff ~ "Positive",
                                .data[[y_axis]] > pval_cutoff & .data[[x_axis]] < -fc_cutoff ~ "Negative",
                                TRUE ~ "Non-sig")) %>%
      mutate(Group = factor(Group, levels = c("Negative", "Non-sig", "Positive")))
    
    ##绘图
    # 需要标注的点：Posive和Negative各前5个点
    # 如果full为真，则返回全部数据
    stopifnot(is.logical(full))
    
    if(!full){
      data_pos <- dplyr::filter(data, Group == "Positive") %>% 
        top_n(pos_num,  -Adjusted.pvalue) 
      
      data_neg <- dplyr::filter(data, Group == "Negative") %>%  
        top_n(neg_num,  -Adjusted.pvalue) 
      
      data_lab <- bind_rows(data_pos, data_neg)
      # 添加差异比较信息
      data_lab <- data_lab %>% mutate(Major_Group = Major_Group,
                                      DE_groups = DE_groups)
      if(nrow(data_lab) == 0){
        return(NULL) # 无显著基因返回NULL
      }else{
        return(data_lab)
      }
      
    }else{
      # 添加差异比较信息
      data <- data %>% mutate(Major_Group = Major_Group,
                              DE_groups = DE_groups)
      return(data)
    }
    
}

##差异表达boxPlot
boxPlot <- function(count_data, meta_gene = NULL, meta_ROI = NULL, 
                    genelist, User_columns){
  
  if(!is.data.frame(count_data)) count_data <- as.data.frame(count_data)
  # 设置count_data列名
  if(!is.null(meta_gene) && is.data.frame(meta_gene)){
    stopifnot(ncol(count_data) == nrow(meta_ROI))
    colnames(count_data) <- meta_ROI$ID
  }
  # 设置count_data行名
  if(!is.null(meta_gene) && is.data.frame(meta_gene)){
    stopifnot(nrow(count_data) == nrow(meta_gene))
    rownames(count_data) <- meta_gene$Probe_name
  }
  
  # Patient_ID作为第一个boxplot图
  if(!"Patient_ID" %in% User_columns &&
     "Patient_ID" %in% colnames(meta_ROI)){
    User_columns <- c("Patient_ID", User_columns)
  }
  # 数据框变换
  df <- count_data %>% t %>% as.data.frame() %>%
    bind_cols(dplyr::select(meta_ROI, any_of(User_columns))) %>%
    pivot_longer(cols = -any_of(User_columns),
                 names_to = c("Gene"),
                 values_to = c("Counts"))
  # 绘图
  envir <- new.env()
  assign('first_plot', value = TRUE, envir = envir)
  
  # 散点注释
  jitter_anno <- User_columns[1]
  
  map(User_columns, function(column){
    df %>% dplyr::filter(Gene %in% genelist) %>%
      dplyr::filter(.data[[column]] != 'Others') %>% # 去除Others的ROI
      ggplot(aes(x = .data[[column]], y = Counts, fill=.data[[column]])) +
      geom_boxplot(outlier.color = NA, alpha = 0.7)+
      theme_bw()+
      labs(x = column,
           y = "Counts") +
      theme(legend.position = "top",
            legend.box = 'vertical') -> p_single
    
    # 除第一张图，均不显示散点图的legend
    first_plot <- get('first_plot', envir = envir)
    if(first_plot){ 
      assign('first_plot', value = FALSE, envir = envir)
      return(p_single + 
               geom_jitter(aes(color = .data[[jitter_anno]]), width = 0.25) +
               theme(axis.text.x = element_text(angle = 90)) )
    }
    p_single + geom_jitter(aes(color = .data[[jitter_anno]]), width = 0.25, show.legend = FALSE)
    
  }) -> p_list
  
  ## 拼图
  if(length(unique(meta_ROI[[jitter_anno]])) <= 10 ||
     length(p_list) < 1){
    # 第一张图的水平数少于10
    # 或者只有一张图
    p <- patchwork::wrap_plots(p_list, nrow = 1)
  }else{
    p_label <- LETTERS[1:length(p_list)]
    p_label <- c(p_label[1], p_label)
    
    p <- patchwork::wrap_plots(p_list, 
                               nrow = 1, 
                               design = paste0(p_label, collapse = ""))
    
  }
  
  p
  
}

# 2.2 Heatmap ---
heatmapPlot <- function(data, meta_gene, meta_ROI, select_cols, 
                        color = NULL,
                        fontsize_row = 6,
                        fontsize_col = 6,
                        cluster_rows = TRUE,
                        cluster_cols = TRUE,
                        clustering_method = "ward.D2",
                        scale = "row", 
                        title = NA,
                        border_color = NA,
                        sampleCluster = FALSE, # 样本聚类图：也就是对ROI取相关性或者欧氏距离所做的矩阵热图
                        sort = FALSE, # 将数据的ROI进行排序，以使得不同的组并在一起，cluster_cols为FALSE时生效
                        transform = 'log',
                        # show_rownames = TRUE,
                        # show_colnames = TRUE,
                        ...){
  # 需要注释的列
  # select_cols
  
  # 单一注释时，去掉Others的ROI
  if(length(select_cols) == 1 ) {
    stopifnot(select_cols %in% colnames(meta_ROI))
    
    idx <- meta_ROI[, select_cols, drop=TRUE] != "Others"
    
    data <- data[, idx, drop=FALSE]
    meta_ROI <- meta_ROI[idx, ,drop=FALSE]
  }
  
  # 列注释
  anno_col <- meta_ROI %>%
    as_tibble() %>%
    dplyr::select(any_of(select_cols), ID) %>%
    column_to_rownames("ID")
  
  # 行注释，未使用
  if(!is.null(meta_gene)){
    anno_row <- meta_gene %>% 
      dplyr::select(Code_class, Probe_name) %>%
      remove_rownames()%>%
      column_to_rownames(var = "Probe_name")
  }
  
  # 数量小于2不进行cluster
  if(nrow(data) <= 2) cluster_rows <- FALSE
  if(ncol(data) <= 2) cluster_cols <- FALSE
  if(any(dim(data) == 1)){
    oneDim <- TRUE
  }else{
    oneDim <- FALSE
  }
  
  # 若行或列的值一样，则无法进行scale, 人工增加一些扰动
  has_same_rows <- apply(data, 1, function(x) all(near(x[1],x)))
  has_same_cols <- apply(data, 2, function(x) all(near(x[1],x)))
  
  if(any(has_same_rows) && cluster_rows && !oneDim) {
    data <- data[!has_same_rows, ,drop=FALSE] 
    if(!is.null(meta_gene)){
      meta_gene <- meta_gene[!has_same_rows, ]
    }
  }
  if(any(has_same_cols) && cluster_cols && !oneDim) {
    data <- data[, !has_same_cols,drop=FALSE]
    meta_ROI <- meta_ROI[!has_same_cols, ]
  }
  
  # 样本聚类图
  if(sampleCluster){ # 欧式距离聚类图
    colnames(data) <- meta_ROI$ID
    
    sampleDist <- dist(t(data))
    
    if(is.null(color)){
      color <- colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(100)
    }
    pheatmap::pheatmap(as.matrix(sampleDist), 
                       color = color,
                       fontsize_row = fontsize_row,
                       fontsize_col = fontsize_col,
                       cluster_rows = cluster_rows,
                       cluster_cols = cluster_cols,
                       clustering_method = clustering_method,
                       scale = "none",
                       main = title,
                       annotation_col = anno_col,
                       border_color = border_color,
                       clustering_distance_rows = sampleDist,
                       clustering_distance_cols = sampleDist,
                       # show_rownames = show_rownames,
                       # show_colnames = show_colnames,
                       ...)
    
    return()
  }
  
  # 数据变换
  if(!is.null(transform) ){
    data <- eval(expr(apply(data, 2, !!as.symbol(transform))))
  }
  
  # 绘图 正常的热图log处理
  colnames(data) <- meta_ROI$ID
  if(is.numeric(sort) || sort ){
    ## 如果为TRUE，则以anno_ROI的第一列为顺序调整ROI顺序
    if(is.logical(sort)) sort <- 1
    
    sortROI <- anno_col %>% arrange(.data[[ colnames(anno_col)[1] ]]) %>% rownames()
    data <- data[, sortROI]
  }
  
  if(is.null(color)){
    color <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
  }
  
  pheatmap::pheatmap(data, 
                     color = color,
                     fontsize_row = fontsize_row,
                     fontsize_col = fontsize_col,
                     cluster_rows = cluster_rows,
                     cluster_cols = cluster_cols,
                     clustering_method = clustering_method,
                     scale = scale,
                     main = title,
                     annotation_col = anno_col, 
                     border_color = border_color, 
                     # show_rownames = show_rownames,
                     # show_colnames = show_colnames,
                     ...)
  
}


# 拆分单一注释热图

splitHeatmapPlot <- function(data, meta_gene, meta_ROI, 
                             variableA, variableB, 
                             anno_cols = variableA, title = TRUE, 
                             fontsize_col = 8, fontsize_row = 8, ...){
  # VariableA：ROI注释
  # variableB：筛选数据的变量
  
  # 所有水平，代表画图数量
  vars <-  meta_ROI[ ,variableB, drop = T] %>% unique()
  
  # 去除只有一个水平的ROI
  bad_levels <- meta_ROI[[variableB]] %>% table %>% .[. == 1] %>% names()
  vars <- vars[! vars %in% bad_levels]
  
  # 绘图
  walk(vars, function(var){
    # 数据筛选
    idx <- which(meta_ROI[ ,variableB, drop = T] == var)
    df <- data[,idx]
    meta_ROI_n <- meta_ROI[idx, ]
    
    # 构造标题
    if(title){
      title <- paste0("Heatmap of ", variableA, 
                      " in only ", var, " (one level of ", variableB, ") data")
    }else{
      title <- NA
    }
    
    # 列注释
    if(anno_cols != variableA) {
      anno_cols <- c(variableA, anno_cols)
    }
    
    # 绘图
    heatmapPlot(df, meta_gene, meta_ROI_n, anno_cols, title = title,
                fontsize_col = fontsize_col, fontsize_row = fontsize_row, ...)
  })
}
  
# 2.3 PCA ---
pcaPlot <- function(count_data, meta_ROI, select_cols, 
                    addEllipses = T, 
                    repel = T,
                    highVariableRatio = NULL,
                    palette = "Dark2", 
                    geom = c("point", "text"), ...){
  # 分组,如果有多个列，只取第一个
  select_cols <- match.arg(select_cols, select_cols)
  Group <- factor(meta_ROI[,select_cols, drop=T])
  
  # 去除Others水平的ROI
  idx <- meta_ROI[, select_cols, drop=TRUE] != "Others"
  
  if(length(which(idx)) < length(Group)){ # 有Others
    count_data <- count_data[, idx, drop=FALSE]
    meta_ROI <- meta_ROI[idx, , drop=FALSE]
    Group <- factor(meta_ROI[,select_cols, drop=T])
  }
  
  # highVariableRatio
  if(!is.null(highVariableRatio)){
    if(is.character(highVariableRatio)) highVariableRatio <- as.numeric(highVariableRatio)
    if(is.na(as.numeric)){
      # as.numeric, 无法解析会产生NA
      highVariableRatio <- 1
    }
    
    # 获得sd的排名，降序排列
    variableGeneOrder <- count_data %>% apply(1, sd) %>% order(decreasing = TRUE)
    geneNumSelected <- round(nrow(count_data) * highVariableRatio)
    
    # 选择前<highVariableRatio*nrow>个排名的基因
    highVariableIdx <- variableGeneOrder < geneNumSelected
    
    count_data <- count_data[highVariableIdx, ,drop=FALSE]
  }
  
  
  # 绘图
  FactoMineR::PCA(t(count_data), graph = F) %>%
    factoextra::fviz_pca_ind(habillage = Group,
                             addEllipses = addEllipses, 
                             repel = repel,
                             palette = palette,
                             geom = geom, ...)  

}

# 2.4 --
getDiffGene <- function(de, grp,
                        y_axis = 'p_val',
                        log10pval_cutoff = 1.3,
                        log2fc_cutoff    = 1){
  
  ## y轴
  y_axis <- case_when(y_axis == "padj" ~ "-log10 adjusted pvalue" ,
                      y_axis == "pval" ~ "-log10 pvalue",
                      TRUE             ~ "illegal")
  if(y_axis == 'illegal') stop('"y_axis" can only be "pval" or "padj"!')
  
  ## x轴
  x_axis <- c('Log2', 'Log2FC')
  if(!any(x_axis %in% colnames(de)) ) stop('Please check log2fc`name in DiffExpr data!') 
  x_axis <- x_axis[x_axis %in% colnames(de)][1]
  
  ##数据变换，添加分组信息：上调、下调及不显著
  de <- de %>%
    mutate(Group = case_when( .data[[y_axis]] > {log10pval_cutoff} & .data[[x_axis]] > {log2fc_cutoff} ~ "Positive",
                              .data[[y_axis]] > {log10pval_cutoff} & .data[[x_axis]] < -{log2fc_cutoff} ~ "Negative",
                              TRUE ~ "Non-sig")) %>%
    mutate(Group = factor(Group, levels = c("Negative", "Non-sig", "Positive")))
  
  ##获取基因
  de %>% dplyr::filter(Group != 'Non-sig') %>% .$Gene
  
}

# zzz ---------------------------
## comoute entropy
entropy <- function(x){
  stopifnot(length(x) > 0 )
  stopifnot(all(x <= 1 && x >= 0))
  res <- 0
  for(i in x){
    res <- res + -i*log2(i)
  }
  return(res)
  
}

#  define func.
flatDoubleList <- function(x, sep = "---", ...){
  res <- list()
  
  name1 <- names(x)
  if(is.null(name1)) name1 <- seq_along(x)
  name2 <- x %>% 
    map(function(x){
      name <- names(x)
      if(is.null(name)) name <- seq_along(x)
      name
    }) 
  
  for(i in seq_along(x)){
    for(j in seq_along(x[[i]])){
      val <- x[[i]][[j]]
      name <- paste(name1[i], name2[[i]][j], sep = sep)
      
      res[[name]] <- val
    }
  }
  res
}

## define surg -> ggplot
ggsurv2gg <- function(ggsurv){
  stopifnot(class(ggsurv) %in% c("ggsurvplot", "ggsurv", "list" ))
  
  gg_idx <- ggsurv %>% map_lgl(is.ggplot)
  
  if(sum(gg_idx) < 1) stop("No ggplot in ggsurvplot object.")
  
  gg_list <- ggsurv[gg_idx]
  
  patchwork::wrap_plots(gg_list,
                        ncol = 1,
                        heights = c(8, 2))
  
}

## define print.pheatmap
print.pheatmap <- function(p){
  grid::grid.newpage()
  pheatmap:::print.pheatmap(p)
}

 
## function
cor.multiTest <- function(x, y){
  row_num <- ncol(x)
  col_num <- ncol(y)
  
  stopifnot(row_num > 1)
  stopifnot(col_num > 1)
  
  cor.mat <- matrix(NA, nrow = row_num, ncol = col_num)
  
  for( i in 1:row_num){
    for(j in 1:col_num){
      cor.mat[i, j] <- cor.test(x[,i,drop=TRUE], y[,j,drop=TRUE])$p.value
    }
  }
  
  rownames(cor.mat) <- colnames(x)
  colnames(cor.mat) <- colnames(y)
  
  cor.mat
}