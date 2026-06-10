
suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
  library(ggplot2)
  library(ggrepel)
  library(matrixStats)
  library(RColorBrewer)
})

# ----------------------------
# Helpers
# ----------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

canonicalize_sample_id <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\.bam$", "", x, ignore.case = TRUE)
  x <- gsub("\\.sorted.*$", "", x, ignore.case = TRUE)
  x
}

fix_metadata_sample_ids <- function(metadata, counts_colnames) {
  metadata$SampleID <- canonicalize_sample_id(metadata$SampleID)
  counts_colnames <- canonicalize_sample_id(counts_colnames)

  # Manual repairs observed in the uploaded files
  repair_map <- c(
    "HCC00_017_PRE" = "HCC002_017_PRE",
    "HCC006_002"    = "HCC006_002_PRE"
  )

  idx <- metadata$SampleID %in% names(repair_map)
  metadata$SampleID[idx] <- unname(repair_map[metadata$SampleID[idx]])

  missing_in_counts <- setdiff(metadata$SampleID, counts_colnames)
  extra_in_counts   <- setdiff(counts_colnames, metadata$SampleID)

  list(
    metadata = metadata,
    missing_in_counts = missing_in_counts,
    extra_in_counts = extra_in_counts
  )
}

strip_ensembl_version <- function(x) sub("\\..*$", "", x)

lncrna_biotypes <- c(
  "lncRNA",
  "3prime_overlapping_ncRNA",
  "antisense",
  "antisense_RNA",
  "bidirectional_promoter_lncRNA",
  "lincRNA",
  "macro_lncRNA",
  "non_coding",
  "processed_transcript",
  "sense_intronic",
  "sense_overlapping"
)


volcano_plot <- function(res_df, file, title = "Volcano plot",
                         padj_cutoff = 0.05, lfc_cutoff = 1) {
  plot_df <- res_df
  plot_df$neglog10_padj <- -log10(plot_df$padj)
  plot_df$neglog10_padj[is.infinite(plot_df$neglog10_padj)] <- NA_real_

  plot_df$status <- "NS"
  plot_df$status[plot_df$padj < padj_cutoff & plot_df$log2FoldChange >= lfc_cutoff] <- "Up"
  plot_df$status[plot_df$padj < padj_cutoff & plot_df$log2FoldChange <= -lfc_cutoff] <- "Down"

  label_df <- subset(plot_df, status != "NS")
  label_df <- label_df[order(label_df$padj, -abs(label_df$log2FoldChange)), ]
  label_df <- head(label_df, 20)

  p <- ggplot(plot_df, aes(x = log2FoldChange, y = neglog10_padj, color = status)) +
    geom_point(alpha = 0.7, size = 1.2, na.rm = TRUE) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    
    scale_color_manual(values = c(Down = "#2C7BB6", NS = "grey70", Up = "#D7191C")) +
    labs(title = title, x = "log2 fold change", y = "-log10 adjusted p-value") +
    theme_bw(base_size = 12)

  print(p)
}

plot_pca_custom <- function(vsd, metadata, intgroup, file, title = "PCA") {
  pca_data <- plotPCA(vsd, intgroup = intgroup, returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))

  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = .data[[intgroup[1]]])) +
    geom_point(size = 3) +
    labs(
      title = title,
      x = paste0("PC1: ", percent_var[1], "% variance"),
      y = paste0("PC2: ", percent_var[2], "% variance"),
      color = intgroup[1]
    ) +
    theme_bw(base_size = 12)

  if (length(intgroup) > 1) {
    p <- p + aes(shape = .data[[intgroup[2]]])
  }
  print(p)
}

plot_sample_distance_heatmap <- function(vsd, metadata, ann_cols, file, title = "Sample distance heatmap") {
  sample_dists <- dist(t(assay(vsd)))
  mat <- as.matrix(sample_dists)
  rownames(mat) <- colnames(vsd)
  colnames(mat) <- colnames(vsd)

  annotation_col <- metadata[, ann_cols, drop = FALSE]
  annotation_col <- as.data.frame(annotation_col)
  rownames(annotation_col) <- rownames(metadata)

  pheatmap(
    mat,
    annotation_col = annotation_col,
    annotation_row = annotation_col,
    clustering_distance_rows = sample_dists,
    clustering_distance_cols = sample_dists,
    main = title
  )
  
}

plot_ma_to_file <- function(res, file, title = "MA plot") {
 
  plotMA(res, ylim = c(-6, 6), main = title)
}

plot_dispersion_to_file <- function(dds, file) {
  
  plotDispEsts(dds)
  
}

plot_heatmap_block <- function(mat, metadata, ann_cols, title, file) {
  annotation_col <- metadata[, ann_cols, drop = FALSE]
  annotation_col <- as.data.frame(annotation_col)
  rownames(annotation_col) <- rownames(metadata)

  
  pheatmap(
    mat,
    scale = "row",
    annotation_col = annotation_col,
    show_rownames = FALSE,
    cluster_cols = F,
    cluster_rows = F,
    fontsize_col = 9,
    main = title
  )
  
}

# ----------------------------
# Main pipeline
# ----------------------------

run_deseq_subset <- function(counts_file,
                             metadata_file,
                             subset_filters = list(),
                             contrast_col,
                             contrast_levels,
                             sample_id_col = "SampleID",
                             gene_id_col = NULL,
                             design_covariates = NULL,
                             min_count = 10,
                             min_samples = 3,
                             alpha = 0.05,
                             lfc_shrink = TRUE,
                             shrink_type = "apeglm",
                             annotation_df = NULL,
                             annotation_cache = "gene_annotation_cache.csv",
                             outdir = "deseq2_run",
                             pca_intgroup = NULL) {

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  qc_dir <- file.path(outdir, "qc")
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

  counts <- read.csv(counts_file, check.names = FALSE, stringsAsFactors = FALSE, row.names = 'geneId')
  metadata <- readxl::read_xlsx(metadata_file)

  if (is.null(gene_id_col)) {
    gene_id_col <- colnames(counts)[1]
  }

  colnames(counts)[1] <- gene_id_col
  metadata[[sample_id_col]] <- canonicalize_sample_id(metadata[[sample_id_col]])
  count_sample_names <- canonicalize_sample_id(colnames(counts)[-1])

  fixed <- fix_metadata_sample_ids(metadata, count_sample_names)
  metadata <- fixed$metadata

  if (length(fixed$missing_in_counts) > 0) {
    message("Metadata samples not found in counts and dropped: ",
            paste(fixed$missing_in_counts, collapse = ", "))
  }

  if (length(subset_filters) > 0) {
    for (nm in names(subset_filters)) {
      if (!nm %in% colnames(metadata)) {
        stop("Subset variable not found in metadata: ", nm)
      }
      keep_values <- subset_filters[[nm]]
      metadata <- metadata[metadata[[nm]] %in% keep_values, , drop = FALSE]
    }
  }

  metadata <- metadata[metadata[[sample_id_col]] %in% count_sample_names, , drop = FALSE]

  if (!contrast_col %in% colnames(metadata)) {
    stop("contrast_col not present in metadata: ", contrast_col)
  }

  if (length(contrast_levels) != 2) {
    stop("contrast_levels must be length 2: c(test, reference)")
  }

  metadata <- metadata[metadata[[contrast_col]] %in% contrast_levels, , drop = FALSE]
  metadata[[contrast_col]] <- factor(metadata[[contrast_col]], levels = contrast_levels)

  if (nrow(metadata) < 4) {
    stop("Too few samples after subsetting: ", nrow(metadata))
  }

  group_table <- table(metadata[[contrast_col]])
  if (any(group_table < 2)) {
    stop("Each contrast level must have at least 2 samples. Current counts: ",
         paste(names(group_table), group_table, sep = "=", collapse = "; "))
  }

  #print(counts)
  print(metadata[[sample_id_col]])
  #print(colnames(counts))
  
  intersect_cols <- intersect(x = colnames(counts), y = metadata[[sample_id_col]])
  counts_sub <- counts[, c( metadata[[sample_id_col]]), drop = FALSE]
  rownames(counts_sub) <- rownames(counts_sub)
  counts_sub[[gene_id_col]] <- NULL

  counts_mat <- as.matrix(counts_sub)
  storage.mode(counts_mat) <- "integer"

  metadata <- as.data.frame(metadata)
  rownames(metadata) <- metadata[[sample_id_col]]

  keep <- rowSums(counts_mat >= min_count) >= min_samples
  counts_mat <- counts_mat[keep, , drop = FALSE]

  design_vars <- c(design_covariates %||% character(0), contrast_col)
  design_vars <- unique(design_vars)

  missing_design_vars <- setdiff(design_vars, colnames(metadata))
  if (length(missing_design_vars) > 0) {
    stop("These design variables are missing from metadata: ",
         paste(missing_design_vars, collapse = ", "))
  }

  metadata_design <- metadata[, design_vars, drop = FALSE]
  metadata_design[] <- lapply(metadata_design, function(x) {
    if (is.character(x) || is.logical(x)) factor(x) else x
  })

  design_formula <- as.formula(paste("~", paste(design_vars, collapse = " + ")))

  dds <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData = metadata,
    design = design_formula, tidy = F
  )

  dds <- DESeq(dds)
  
  ###### For CORR
  vsd <- vst(dds)
  expr <- assay(vsd)
  ######

  if (lfc_shrink) {
    if (!requireNamespace("apeglm", quietly = TRUE) && shrink_type == "apeglm") {
      message("apeglm not installed; falling back to unshrunk log2FC.")
      lfc_shrink <- FALSE
    }
  }

  contrast_name <- resultsNames(dds)[grep(paste0("^", contrast_col, "_"), resultsNames(dds))][1]
  res <- results(dds, contrast = c(contrast_col, contrast_levels[1], contrast_levels[2]), alpha = alpha)

  if (lfc_shrink) {
    res <- lfcShrink(
      dds,
      contrast = c(contrast_col, contrast_levels[1], contrast_levels[2]),
      res = res,
      type = shrink_type
    )
  }

  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  #res_df$gene_id_clean <- strip_ensembl_version(res_df$gene_id)
  full_degs <- res_df

  
  lnc_annotation <- read.delim('mart_export_all_lnc.txt', header = T)
  lnc_annotation <- lnc_annotation[
    !duplicated(lnc_annotation$Gene_stable_ID_version),
  ]
  #print(head(lnc_annotation))
  
  #print(head(res_df))
  
  

  res_df <- merge(
    res_df,
    lnc_annotation,
    by.x = "gene_id",
    by.y = "Gene_stable_ID_version",
    all.x = FALSE,
    sort = FALSE
  )
  print(head(res_df))

  if (!"gene_name" %in% colnames(res_df)) {
    res_df$gene_name <- res_df$gene_id_clean
  }
  res_df$gene_name[is.na(res_df$gene_name) | res_df$gene_name == ""] <- res_df$gene_id_clean[is.na(res_df$gene_name) | res_df$gene_name == ""]

  res_df <- res_df[order(res_df$padj, -abs(res_df$log2FoldChange)), ]

  norm_counts <- counts(dds, normalized = TRUE)
  

  vsd <- vst(dds, blind = FALSE)
  vst_mat <- assay(vsd)

  ann_cols <- unique(c(contrast_col, names(subset_filters), pca_intgroup %||% contrast_col))
  ann_cols <- ann_cols[ann_cols %in% colnames(metadata)]

  plot_dispersion_to_file(dds, file.path(qc_dir, "dispersion_plot.png"))
  plot_ma_to_file(res, file.path(qc_dir, "MA_plot.png"),
                  title = paste0(contrast_levels[1], " vs ", contrast_levels[2]))
  plot_sample_distance_heatmap(
    vsd, metadata, ann_cols = ann_cols,
    file = file.path(qc_dir, "sample_distance_heatmap.png")
  )

  pca_groups <- pca_intgroup %||% c(contrast_col, names(subset_filters))
  pca_groups <- unique(pca_groups[pca_groups %in% colnames(metadata)])
  if (length(pca_groups) == 0) pca_groups <- contrast_col
  plot_pca_custom(
    vsd, metadata, intgroup = pca_groups,
    file = file.path(qc_dir, "PCA.png"),
    title = paste0("PCA: ", contrast_levels[1], " vs ", contrast_levels[2])
  )

  volcano_plot(
    res_df,
    file = file.path(qc_dir, "volcano.png"),
    title = paste0(contrast_levels[1], " vs ", contrast_levels[2], " (", contrast_col, ")")
  )
  #print(head(res_df))

  ordered_res <- res_df[order(res_df$log2FoldChange, decreasing = TRUE), ]
  top_up_ids <- head(ordered_res$gene_id[!is.na(ordered_res$log2FoldChange)], 50)
  top_down_ids <- head(rev(ordered_res$gene_id[!is.na(ordered_res$log2FoldChange)]), 50)
  
  #actuals
  actual_up_id <- ordered_res$gene_id[
    !is.na(ordered_res$log2FoldChange) &
      !is.na(ordered_res$padj) &
      ordered_res$log2FoldChange > 1 &
      ordered_res$padj < 0.05
  ]

  actual_down_id <- ordered_res$gene_id[
    !is.na(ordered_res$log2FoldChange) &
      !is.na(ordered_res$padj) &
      ordered_res$log2FoldChange < -1 &
      ordered_res$padj < 0.05
  ]
  
  

  top_up_ids <- intersect(top_up_ids, rownames(vst_mat))
  top_down_ids <- intersect(top_down_ids, rownames(vst_mat))
  
  print(head(top_up_ids))

  if (length(top_up_ids) >= 2) {
    plot_heatmap_block(
      vst_mat[top_up_ids, , drop = FALSE],
      metadata = metadata,
      ann_cols = ann_cols,
      title = "Top 50 upregulated genes",
      file = file.path(qc_dir, "heatmap_top50_up.png")
    )
  }

  if (length(top_down_ids) >= 2) {
    plot_heatmap_block(
      vst_mat[top_down_ids, , drop = FALSE],
      metadata = metadata,
      ann_cols = ann_cols,
      title = "Top 50 downregulated genes",
      file = file.path(qc_dir, "heatmap_bottom50_down.png")
    )
  }
  print(top_down_ids)
  
  up_ <- subset(res_df, gene_id %in% top_up_ids)
  down_ <- subset(res_df, gene_id %in% top_down_ids)
  
  #lncTops
  upLNC <- subset(res_df, (gene_id %in% top_up_ids) & !is.na(Gene_stable_ID) )
  downLNC <- subset(res_df, (gene_id %in% top_down_ids) & !is.na(Gene_stable_ID) )

  return(list('UP' = up_, 'DOWN' = down_, 'ALL' = ordered_res,
              'UP_LIST' = top_up_ids, 'DOWN_LIST' = top_down_ids,
              'UP_ACT' = actual_up_id, 'DOWN_ACT' = actual_down_id,
              'FULL_DEGS' = full_degs, FOR_CORR = expr,
              'upLNC'= upLNC, 'downLNC' = downLNC))
}


##for correlation

makeCorrelation <- function(direction, for_corr, cutoff = 0.8, n_samp_min = 10) {
  
  
  #corr
  lncRNA_list <- subset(direction, !is.na(Gene_stable_ID))$gene_id #Gene_stable_ID
  
  FOR_CORR <- for_corr
  lncEXP <- FOR_CORR[lncRNA_list,]
  
  #remove low overall count
  keep <- rowSums(FOR_CORR >= 10) >= n_samp_min #Keep genes with at least 10 counts, in at least 5 samples
  FOR_CORR <- FOR_CORR[keep, ]
  
  #remove low variance
  gene_var <- apply(FOR_CORR, 1, var)
  FOR_CORR <- FOR_CORR[gene_var > 0.75, ]
  dim(FOR_CORR)
  FOR_CORR <- rbind(FOR_CORR, lncEXP)
  
  
  dim(FOR_CORR)
  final_CORR <- cor(t(FOR_CORR))
  final_CORR_UT <- upper.tri(final_CORR)
  
  corr_table <- data.frame(
    gene1 = rownames(final_CORR)[row(final_CORR)[final_CORR_UT]],
    gene2 = colnames(final_CORR)[col(final_CORR)[final_CORR_UT]],
    correlation = final_CORR[final_CORR_UT]
  )
  
  dim(corr_table)
  head(corr_table)
  
  
  lnc_corr <- subset(
    corr_table,
    abs(correlation) > cutoff &
      (
        (gene1 %in% lncRNA_list & !(gene2 %in% lncRNA_list)) |
          (gene2 %in% lncRNA_list & !(gene1 %in% lncRNA_list))
      )
  )
  
  return(lnc_corr)
}
