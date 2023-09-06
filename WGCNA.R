library(tidyverse)
library(magrittr)
library(WGCNA)
library(pheatmap)
source('tools.R') 

wgcna_plots_heatmap <-
  all_data_splitByMS[3:4] %>% 
  map2(names(.), function(dat, name){
    
    counts_ana <- dat[['counts_ana']] 
    
    meta_ana_gene <- dat[['meta_ana_gene']]
    meta_ana_ROI <- dat[['meta_ana_ROI']]
    # name
    title <- name
    name <- name %>% str_split(" -> ") %>% .[[1]]
    location <- name[1]
    name <- name[2]
    
    # omit all Others
    others_idx <- (
      meta_ana_ROI$Group1 == "Others" &
        meta_ana_ROI$Group3 == "Others" &
        meta_ana_ROI$Group4 == "Others"
    )
    counts_ana <- counts_ana[, !others_idx, drop=FALSE]
    meta_ana_ROI <- meta_ana_ROI[!others_idx, ,drop=FALSE]
    
    
    datExpr0 <- t(counts_ana)
    
    ## data input: expression dat ##
    # test na
    gsg = goodSamplesGenes(datExpr0, verbose = 3);
    gsg$allOK
    
    # clust for outlier samples
    sampleTree = hclust(dist(datExpr0), method = "average");
    plot(sampleTree, 
         main = "Sample clustering to detect outliers", 
         sub="",
         xlab="", 
         cex.lab = 1.5,
         cex.axis = 1.5, 
         cex.main = 2)
    
    if(name == "MSI"){
      height = 22500
    }else if(name == "MSS"){
      height = 30000
    }else{
      message("Others group name...")
      browser()
    }
    
    abline(h = height, col = "red");
    
    clust = cutreeStatic(sampleTree, cutHeight = height, minSize = 10)
    table(clust)
    
    keepSamples = (clust==1)
    datExpr = datExpr0[keepSamples, ]
    nGenes = ncol(datExpr)
    nSamples = nrow(datExpr)
    
    # plot
    sampleTree2 = hclust(dist(datExpr), method = "average");
    sampleTree2$labels %<>% str_remove("_.*$")
    
    par(mar = c(5, 5, 4, 2))
    plot(sampleTree2, 
         main = "Sample clustering", 
         sub="",
         xlab="", 
         cex.lab = 1.5,
         cex.axis = 1.5, 
         cex.main = 2)
    ## data input: clinal trait ##
    datTraits <- meta_ana_ROI[, "Group5", drop=FALSE] #[,c("Group1", "Group3", "Group4")]
    
    # make sure that order was identity betwent expression mat and clinical trait
    Samples = rownames(datExpr);
    traitRows = match(Samples, rownames(datTraits));
    if(sum(is.na(traitRows)) > 0){
      stop("Some ROIs have no datTraits...")
    }
    datTraits = datTraits[traitRows, ,drop=FALSE];
    
    
    # final test for expr_mat and clinal_trati
    sampleTree2 = hclust(dist(datExpr), method = "average")
    traitColors = labels2colors(datTraits);
    plotDendroAndColors(sampleTree2, traitColors,
                        groupLabels = names(datTraits),
                        main = "Sample dendrogram and trait heatmap")
    
    ## pick soft factor ##
    powers = c(c(1:10), seq(12, 30, by=2))
    
    sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
    
    local({
      par(mfrow = c(1,2));
      cex1 = 0.8;
      # Scale-free topology fit index as a function of the soft-thresholding power
      plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
           xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
           main = paste("Scale independence"));
      text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
           labels=powers,cex=cex1,col="red");
      # this line corresponds to using an R^2 cut-off of h
      abline(h=0.85,col="red")
      # Mean connectivity as a function of the soft-thresholding power
      plot(sft$fitIndices[,1], sft$fitIndices[,5],
           xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
           main = paste("Mean connectivity"))
      text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
      
    })
    
    net = blockwiseModules(datExpr, 
                           maxBlockSize = 5000,
                           power = ifelse(is.na(sft$powerEstimate), 10,
                                          sft$powerEstimate),
                           TOMType = "unsigned", 
                           minModuleSize = 30,
                           reassignThreshold = 0,
                           numericLabels = TRUE,
                           pamRespectsDendro = FALSE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = "wgcnaTOM",
                           verbose = 3)
    nGenes = ncol(datExpr);
    nSamples = nrow(datExpr);
    
    moduleColors = labels2colors(net$colors)
    
    plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                        "Module colors",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
    
    
    design <- str_glue("~0+{colnames(datTraits)}") %>%
      map(as.formula) %>%
      map(~model.matrix(.x, data = datTraits)) %>%
      map(as.data.frame) %>%
      bind_cols()
    
    colnames(design) <- colnames(design) %>% str_remove("^Group[0-9]")

    design <- datTraits
    design$Group5 <- factor(design$Group5, levels = c("TILs-low", "TILs-high")) %>% as.numeric()
    
    
    
    MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
    MEs = orderMEs(MEs0)
    
    plotEigengeneNetworks(cbind(MEs, design), #%>% dplyr::rename(MEbrown = 2), 
                          setLabels = paste0(
                            name, 
                            " Eigengene adjacency plot"
                          ),
                          marDendro = c(1,10,1,2),
                          marHeatmap = c(0.1,10,2,2), 
                          plotDendrograms = T, 
                          xLabelsAngle = 90)
    
    
    moduleTraitCor = cor(MEs, design, use = "p");
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
    
    saveRDS(
      list(
        moduleTraitCor = moduleTraitCor,
        moduleTraitPvalue = moduleTraitPvalue
      ),
      file = str_glue("{name} moduleTraitCor and moduleTraitPvalue.rds")
    )
    
    sizeGrWindow(10,6)
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                       signif(moduleTraitPvalue, 1), ")", sep = "");
    dim(textMatrix) = dim(moduleTraitCor)

    par(mar = c(6, 10, 3, 3));
    labeledHeatmap(Matrix = moduleTraitCor,
                   xLabels = colnames(design) %>% str_remove("^.*Clinical"),
                   yLabels = names(MEs),
                   ySymbols = names(MEs),
                   colorLabels = FALSE,
                   colors = greenWhiteRed(50),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 0.8,
                   xLabelsAngle = 0,
                   xLabelsAdj = 0.5,
                   zlim = c(-1,1),
                   main = paste("Module-trait relationships"))
    
    ## gene and clinical trait ##
    
    modNames = substring(names(MEs), 3)
    geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
    MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
    names(geneModuleMembership) = paste("MM", modNames, sep="");
    names(MMPvalue) = paste("p.MM", modNames, sep="");
    
    geneTraitSignificance = as.data.frame(cor(datExpr, design, use = "p"));
    GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
    names(geneTraitSignificance) = paste("GS.", colnames(design), sep="");
    names(GSPvalue) = paste("p.GS.", colnames(design), sep="");
    
    
    # filter gene by high GS and MM
    
    if(name == "MSI" && location == "Tumor"){
      module <-  c("brown", "green" ) 
    }else if(name == "MSS" && location == "Tumor"){
      module <-  c("blue" ) 
    }else{
      message("Others group name...")
      browser()
    }
    
    column = match(module, modNames);
    
    filterHubGene <- function(MM,
                              GS,
                              clinicalTrait = "Clincal",
                              module = "black",
                              moduleColors = NULL){
      module = module
      modNames = substring(names(MM), 3)
      
      column = match(module, modNames);
      moduleGenes = moduleColors==module;
      
      par(mfrow = c(1,1));
      verboseScatterplot(abs(MM[moduleGenes, column]),
                         abs(GS[moduleGenes, str_glue("GS.{clinicalTrait}")]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = str_glue("Gene significance for {clinicalTrait}"),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
    }
    
    filterHubGene(MM = geneModuleMembership,
                  GS = geneTraitSignificance,
                  clinicalTrait = "Group5",
                  module = module[1],
                  moduleColors = moduleColors)
    
    abline(h = 0.1, col = "red")
    abline(v = 0.7, col = "red")
    
    # hub gene for T_Cell_Infiltration 
    gene_list <- 
      module %>%
      set_names(., .) %>%
      map(function(x){ 
        moduleGenes <- moduleColors==x
        x <- paste0("MM", x)
        idx <- abs(geneModuleMembership[[x]]) > 0.7 & abs(geneTraitSignificance$GS.Group5) > 0.1
        hub_gene <- colnames(datExpr)[idx & moduleGenes]
        
        hub_gene
      })
    
    select_col <- colnames(datTraits)
    title <- name 
    
    color <- c(
      colorRampPalette(c("green", "#008200"))(20),
      colorRampPalette(c("#007700", "#000500"))(10),
      "black",
      colorRampPalette(c("#050000", "#820000"))(10),
      colorRampPalette(c("#8C0000", "red"))(20)
    )
    
    annotation_colors <- list(
      Module = module %>% set_names(., .)
    )
    
    anno_row  <- gene_list %>%
      enframe(name = "Module", value = "Gene") %>%
      unnest(Gene) %>%
      column_to_rownames("Gene")
    
    # filter
    roi_idx <- meta_ana_ROI[["Group1"]] == title
    target_list <- rownames(anno_row)
    
    # plot
    p <- tryCatch(
      heatmapPlot(
        data = counts_ana[target_list, roi_idx, drop = FALSE],
        meta_gene = meta_ana_gene[target_list, , drop = FALSE],
        meta_ROI = meta_ana_ROI[roi_idx, , drop = FALSE],
        select_cols = select_col %>% union(c("block", "Group1", "Group2", "Group3", "Group4", "Group5")),
        color = color,
        annotation_colors = annotation_colors,
        annotation_row = anno_row,
        clustering_distance_rows = "euclidean", # "correlation", #'euclidean',
        clustering_distance_cols = "euclidean", # "correlation", #'euclidean',
        fontsize_col = 10,
        fontsize_row =  6, # ifelse(length(target_list) > 25, 7, 10),
        show_rownames = TRUE,
        show_colnames = T,
        title =  title,
        silent = TRUE
      ),
      error = function(x)
        return(NULL)
    )
    
    return(p)
  })