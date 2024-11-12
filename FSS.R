# ----------------------------------------------------------------------------------
# Pancreatic Disease Center, Ruijin Hospital, SHSMU, Shanghai, China.
# Author: Dongjie Chen
# Date: 2024-11-03
# Email: chen_dj@sjtu.edu.cn
# ----------------------------------------------------------------------------------

rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

# ----------------------------------------------------------------------------------
# Step 1: Calculate FSx
# ----------------------------------------------------------------------------------

# Load necessary data
cellular_senescence_gs <- readRDS("cellular_senescence_gene.rds")
scRNA_list <- readRDS("40_pancancer_scRNA_seq_dat.rds")

# Calculate FSx
FSx_list <- pbapply::pblapply(
  1:length(scRNA_list),
  FUN = function(x) {
    sce <- scRNA_list[[x]]
    Idents(sce) <- "Celltype"
    sce <- subset(sce, idents = "Malignant cells")
    counts <- sce@assays$RNA@data
    
    S_CAF_score <- gsva(
      expr = counts,
      gset.idx.list = cellular_senescence_gs,
      kcdf = "Gaussian",
      parallel.sz = 60
    )
    S_CAF_score <- as.data.frame(t(getAUC(S_CAF_score)))
    
    lmGenes <- data.frame(
      gene = rownames(counts),
      coef = NA,
      p = NA
    )
    
    for (i in 1:nrow(counts)) {
      cor <- cor.test(counts[i, ], S_CAF_score$cellular_senescence_gs, method = "spearman")
      lmGenes$coef[i] <- cor$estimate
      lmGenes$p[i] <- cor$p.value
    }
    
    lmGenes$p.adjust <- p.adjust(lmGenes$p, method = "BH")
    return(lmGenes)
  }
)
names(FSx_list) <- names(scRNA_list)

# ----------------------------------------------------------------------------------
# Step 2: Calculate FSy
# ----------------------------------------------------------------------------------

FSy_list <- pbapply::pblapply(
  1:length(scRNA_list),
  FUN = function(x) {
    sce <- scRNA_list[[x]]
    sce$group <- ifelse(
      sce$Celltype == "CAF",
      "CAF",
      "control"
    )
    Idents(sce) <- "group"
    future::plan("multisession", workers = 60)
    
    DE <- FindMarkers(
      sce,
      ident.1 = "CAF",
      group.by = "group",
      logfc.threshold = 0.25,
      min.pct = 0.1,
      base = exp(1)
    )
    DE <- DE[DE$p_val_adj < 1e-05, ]
    return(DE)
  }
)
names(FSy_list) <- names(scRNA_list)

# ----------------------------------------------------------------------------------
# Step 3: LM.SIG Construction
# ----------------------------------------------------------------------------------

# Verify FSx and FSy have identical structure
identical(names(FSx_list), names(FSy_list))

# Construct LM.SIG
ls_FSn <- pbapply::pblapply(
  1:length(FSx_list),
  FUN = function(x) {
    FSx <- FSx_list[[x]]
    FSy <- FSy_list[[x]]
    
    FSx <- FSx[FSx$coef > 0 & FSx$p.adjust < 1e-05, ]
    FSy <- rownames(FSy[FSy$avg_logFC > 0.25, ])
    FSy <- FSy[!grepl("^RP[SL]", FSy, ignore.case = FALSE)]  # Exclude ribosomal proteins
    
    FSn <- FSx[FSx$gene %in% FSy, ]
    return(FSn)
  }
)

# Combine all gene data
allGenes <- Reduce(rbind, ls_FSn)
allGenes <- unique(allGenes$gene)
allGenesDf <- data.frame(gene = allGenes)

# Merge coefficients from each dataset
for (i in 1:length(ls_FSn)) {
  allGenesDf <- left_join(
    allGenesDf,
    ls_FSn[[i]][, c("gene", "coef")],
    by = "gene"
  )
}

# Prepare final signature
allGenesDf <- allGenesDf[!is.na(allGenesDf$gene), ]
rownames(allGenesDf) <- allGenesDf$gene
allGenesDf <- allGenesDf[, -1]
colnames(allGenesDf) <- names(FSx_list)

genelist <- allGenesDf
genelist$all_gmean <- compositions::geometricmeanRow(genelist[, 1:length(ls_FSn)])  # Geometric mean

# Filter genes by geometric mean threshold
sig <- genelist[genelist$all_gmean > 0.25, ]
sig <- sig[order(sig$all_gmean, decreasing = TRUE), ]
FSS <- rownames(sig)

# Save results
saveRDS(FSS, file = "FSS.rds")
