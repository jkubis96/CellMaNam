#' @import openxlsx dplyr ggplot2 pheatmap reshape2



NULL

#' Load marker data from a file
#' Marker data structure available on GitHub (https://github.com/jkubis96/CellMaNam)
#'
#' @param data_path Path to the input file (CSV, TSV, TXT, or XLSX)
#' @return A dataframe containing the marker data
#' @export
load_markers_data <- function(data_path) {


  ext <- tools::file_ext(data_path)

  if (".csv" %in% ext) {
    data <- fread(data_path)
  } else if (".tsv" %in% ext || ".txt" %in% ext) {
    data <- fread(data_path, sep = "\t")
  } else if (ext == "xlsx") {
    data <- read.xlsx(data_path)
  } else {
    stop("Unsupported file extension")
  }

  return(data)

}



#' Load marker hierarchy from a file
#' Marker hierarchy data structure available on GitHub (https://github.com/jkubis96/CellMaNam)
#'
#' @param data_path Path to the input file (CSV, TSV, TXT, or XLSX)
#' @return A dataframe containing the marker hierarchy
#' @export
load_markers_hierarchy <- function(data_path) {


  ext <- tools::file_ext(data_path)

  if (".csv" %in% ext) {
    data <- fread(data_path)
  } else if (".tsv" %in% ext || ".txt" %in% ext) {
    data <- fread(data_path, sep = "\t")
  } else if (ext == "xlsx") {
    data <- read.xlsx(data_path)
  } else {
    stop("Unsupported file extension")
  }

  return(data)

}



#' Create a data format for marker occurrences, which will be used to identify and remove the most common markers across the dataset.
#'
#' @param markers_data Dataframe containing marker data
#' @param features_col Name of the column containing features
#' @param cell_column Name of the column containing cell annotations
#' @return Dataframe with feature occurrence counts
#' @export
calc_occurrence <- function(markers_data, features_col = NULL, cell_column = NULL) {

  if (!is.null(features_col) && !(features_col %in% colnames(markers_data))) {
    stop(paste("Error: Column", features_col, "not found in the dataset."))
  }

  if (!is.null(cell_column) && !(cell_column %in% colnames(markers_data))) {
    stop(paste("Error: Column", cell_column, "not found in the dataset."))
  }

  if (!is.null(cell_column) && !is.null(features_col)) {
    markers_data <- markers_data %>%
      select(all_of(c(cell_column, features_col)))
  }


  colnames(markers_data)[colnames(markers_data) == cell_column] <- 'cell_annotation'
  colnames(markers_data)[colnames(markers_data) == features_col] <- 'features'


  markers_data <- distinct(markers_data)

  markers_data <- markers_data %>%
    group_by(features) %>%
    mutate("features_occurrence" := n()) %>%
    ungroup()

  return(markers_data)

}



#' Identify the most frequently occurring features for each cell annotation.
#' The input data must be derived from the calc_occurrence() function.
#'
#' @param markers_occ Dataframe containing feature occurrences
#' @param top_n Number of thresholds for selecting the top features
#' @return Filtered dataframe with top occurring features
#' @export
select_top_occ <- function(markers_occ, top_n = 2) {

  final_df <- data.frame()

  for (l in unique(markers_occ$cell_annotation)) {
    tmp <- markers_occ[markers_occ$cell_annotation %in% l,]
    n <- sort(as.numeric(unique(tmp$features_occurrence)))

    if (top_n < length(n)) {
      n <- n[[top_n]]
    } else {
      n <- n[[length(n)]]

    }

    tmp <- tmp[tmp$features_occurrence <= n,]

    final_df <- rbind(final_df, tmp)
  }

  return(final_df)
}


#' Identify markers from RNA sequencing (scRNAseq) data.
#' The input dataframe should be tailored to the type of RNA sequencing used:
#' For single-cell RNA sequencing (scRNA-seq) where each cell is sorted and sequenced like RNAseq each column should represent an individual cell.
#' In the case of Drop-seq data with cluster information, the dataframe should include the average expression per cluster, with each cluster represented by its name.
#'
#' @param rna_seq_df Dataframe with RNA sequencing (scRNAseq) data
#' @param exclude_genes List of genes to exclude
#' @param exclude_mt Boolean, whether to exclude mitochondrial genes. Default: TRUE
#' @param include_markers List of markers to include. Recommended all markers inside in cell markers data used for naming.
#' @param stat Statistical method for filtering ('mean', 'median', 'q1', 'q3'). Default: 'mean'
#' @return Dataframe of identified markers
#' @export
find_markers_rna_seq <- function(rna_seq_df, exclude_genes = NULL, exclude_mt = TRUE, include_markers = NULL, stat = 'mean') {

  if (is.atomic(exclude_genes)) {

    rna_seq_df <- rna_seq_df[!toupper(rownames(rna_seq_df)) %in% toupper(exclude_genes),]

  }


  if (is.atomic(include_markers)) {

    rna_seq_df <- rna_seq_df[toupper(rownames(rna_seq_df)) %in% toupper(include_markers),]

  }


  if (exclude_mt == TRUE) {

    rna_seq_df <- rna_seq_df[!grepl('MT.|MT-',toupper(rownames(rna_seq_df))),]

  }




  colnames(rna_seq_df) <- make.unique(colnames(rna_seq_df))


  cell_annotation <- c()
  markers <- c()

  for (c in colnames(rna_seq_df)) {

    if (stat == 'mean') {

      tmp_markers <- as.vector(rownames(rna_seq_df)[rna_seq_df[,c] >= mean(rna_seq_df[,c][rna_seq_df[,c] > 0])])

    } else if (stat == 'median') {

      tmp_markers <- as.vector(rownames(rna_seq_df)[rna_seq_df[,c] >= median(rna_seq_df[,c][rna_seq_df[,c] > 0])])

    } else if (stat == 'q1') {

      tmp_markers <- as.vector(rownames(rna_seq_df)[rna_seq_df[,c] >= quantile(rna_seq_df[,c][rna_seq_df[,c] > 0] ,0.25)])

    } else if (stat == 'q3') {

      tmp_markers <- as.vector(rownames(rna_seq_df)[rna_seq_df[,c] >= quantile(rna_seq_df[,c][rna_seq_df[,c] > 0] ,0.75)])

    }

    markers <- c(markers, tmp_markers)
    cell_annotation <- c(cell_annotation, as.vector(rep(c, length(tmp_markers))))


  }


  results <- data.frame(cell_annotation,markers)

}



#' Annotate cell markers based on statistical significance.
#' This function annotates each cell using the provided markers for cell naming, based on the markers returned by `find_markers_rna_seq()`.
#' Alternatively, you can use a different set of markers (e.g., obtained from Seurat) by loading them into the function.
#' For custom cell marker sets, you must modify the 'cell_markers' dataframe, ensuring that the cell names are in the 'cell_annotation' column and the genes are listed under a separate 'features' column.
#'
#' @param cell_markers Dataframe containing marker information for each cell
#' @param markers_occ Dataframe with marker occurrences in different cell types
#' @param max_genes Integer, the maximum number of genes to consider
#' @return A dataframe with cell annotations, associated cell types, p-values, and gene percentages
#' @details Uses Fisher's exact test to assess marker significance.
#' @export
get_annotation <- function(cell_markers, markers_occ, max_genes) {


  cell_annotation <- c()
  cell_name <- c()
  p_val <- c()
  genes_perc <- c()



  for (cell in unique(cell_markers$cell_annotation)) {

    p_value_fish_tmp <- c()
    p_value_bin_tmp <- c()
    genes_perc_tmp <- c()

    for (cell_type in unique(markers_occ$cell_annotation)) {

      cell_name <- c(cell_name, cell_type)

      tmp_case <- cell_markers$markers[cell_markers$cell_annotation %in% cell]
      case_genes <-  length(tmp_case[toupper(tmp_case) %in% toupper(markers_occ$features[markers_occ$cell_annotation %in% cell_type])])
      db_genes <- length(cell_markers$markers[cell_markers$cell_annotation %in% cell])


      case_group <- length(markers_occ$features[markers_occ$cell_annotation %in% cell_type])
      db_group <- max_genes



      fisher_matrix <- matrix(c(case_genes, db_genes - case_genes,
                                case_group, db_group - case_group),
                              nrow = 2, byrow = TRUE)


      fisher_result <- fisher.test(fisher_matrix, alternative = "greater")

      p_value_fish_tmp <- c(p_value_fish_tmp, fisher_result$p.value)

      genes_perc_tmp <- c(genes_perc_tmp, round(case_genes / case_group,2))





    }


    cell_annotation <- c(cell_annotation, rep(cell, length(genes_perc_tmp)))
    p_val <- c(p_val, p_value_fish_tmp)
    genes_perc <- c(genes_perc, genes_perc_tmp)



  }



  annotation_df <- data.frame(cell_annotation, cell_name, p_val, genes_perc)


  return(annotation_df)

}





#' Perform cell type annotation based on annotation (get_annotation()).data.
#' If necessary, you can use hierarchical data (via `load_markers_hierarchy()`) to find and construct hierarchical relationships in cell names.
#'
#' @param annotation_data Dataframe with cell annotations and p-values (get_annotation())
#' @param hierarchy_data Dataframe representing the hierarchy of cell types (load_markers_hierarchy() or NULL if not needed)
#' @param p_val Threshold p-value for selecting significant annotations. Default: 0.01
#' @param level Integer, specifying the level in the hierarchy to use. Default: 1
#' @param hierarchy Boolean, whether to use hierarchical annotation. If you want to use hierarchical data for name creating. If not the names will based on 'level'. Default: TRUE
#' @return A dataframe with assigned cell types, weighted p-values, and log-transformed p-values
#' @details Uses hierarchical relationships and statistical thresholds to classify cell types.
#' @export
cell_typing <- function(annotation_data, hierarchy_data, p_val = 0.01, level = 1, hierarchy = TRUE) {

  input_annotations <- annotation_data$cell_annotation

  annotation_data <- annotation_data[annotation_data$p_val <= p_val,]

  annotation_data$p_val[annotation_data$p_val == 0] = min(annotation_data$p_val)*0.9

  full_names <- c()
  annotation <- c()
  weighted_p_val <- c()
  weighted_pct <- c()



  if (is.data.frame(hierarchy_data)) {

    if (hierarchy == FALSE) {


      annotation_data = annotation_data[annotation_data$cell_name %in% hierarchy_data[,level], ,drop = FALSE]

      full_data <- data.frame(annotation = annotation_data$cell_annotation,
                              full_names = annotation_data$cell_name,
                              weighted_p_val = annotation_data$p_val,
                              weighted_pct = annotation_data$genes_perc,
                              log2_p_val = -log2(annotation_data$p_val),
                              completed = -log2(annotation_data$p_val) * annotation_data$genes_perc)

      for (miss in input_annotations) {

        if (!miss %in% full_data$annotation) {
          full_data <- rbind(full_data, data.frame(annotation = c(miss),
                                                   full_names = c('Undefined'),
                                                   weighted_p_val = c(-1),
                                                   weighted_pct = c(0),
                                                   log2_p_val = c(0),
                                                   completed = c(0))
          )

        }
      }

      return(full_data)

    }

    if (level > 1) {

      hierarchy_data <- hierarchy_data[,level:length(hierarchy_data), drop = F]
    }

    for (ar in unique(annotation_data$cell_annotation)) {

      anno_tmp <- annotation_data[annotation_data$cell_annotation %in% ar,]

      cells <- unique(hierarchy_data[,1])

      for (c in cells) {

        tmp_h <- hierarchy_data[hierarchy_data[,1] %in% c,]


        for (r in 1:nrow(tmp_h)) {

          tmp_full_names <- c()
          tmp_pval <- c()
          tmp_perc <- c()


          for (cs in colnames(tmp_h)) {

            tmp_markers <- anno_tmp[anno_tmp$cell_name %in%  tmp_h[r,cs],]

            if (length(tmp_markers) > 0) {

              if (!tmp_h[r,cs] %in% tmp_full_names) {


                tmp_full_names <- c(tmp_full_names, tmp_h[r,cs])
                tmp_pval <- c(tmp_pval, tmp_markers$p_val)
                tmp_perc <- c(tmp_perc, tmp_markers$genes_perc)




              }


            }



          }

          if (length(tmp_full_names) > 0) {

            full_names <- c(full_names, paste(tmp_full_names, collapse = " -> "))
            annotation <- c(annotation, ar)
            weighted_p_val <- c(weighted_p_val, prod(tmp_pval))
            weighted_pct <- c(weighted_pct, mean(tmp_perc))


          }

        }


      }


    }


    full_data <- data.frame(annotation,
                            full_names,
                            weighted_p_val,
                            weighted_pct,
                            log2_p_val = -log2(weighted_p_val),
                            completed = -log2(weighted_p_val) * weighted_pct)

    for (miss in input_annotations) {

      if (!miss %in% full_data$annotation) {
        full_data <- rbind(full_data, data.frame(annotation = c(miss),
                                                 full_names = c('Undefined'),
                                                 weighted_p_val = c(-1),
                                                 weighted_pct = c(0),
                                                 log2_p_val = c(0),
                                                 completed = c(0))
        )

      }
    }

    return(full_data)


  } else {


    full_data <- data.frame(annotation = annotation_data$cell_annotation,
                            full_names = annotation_data$cell_name,
                            weighted_p_val = annotation_data$p_val,
                            weighted_pct = annotation_data$genes_perc,
                            log2_p_val = -log2(annotation_data$p_val),
                            completed = -log2(annotation_data$p_val) * annotation_data$genes_perc)


    for (miss in input_annotations) {

      if (!miss %in% full_data$annotation) {
        full_data <- rbind(full_data, data.frame(annotation = c(miss),
                                                 full_names = c('Undefined'),
                                                 weighted_p_val = c(-1),
                                                 weighted_pct = c(0),
                                                 log2_p_val = c(0),
                                                 completed = c(0))
        )

      }
    }

    return(full_data)


  }


}




#' Generate a heatmap of cell type annotations
#'
#' @param cell_types_data Dataframe containing cell type data (cell_typing())
#' @param value Column name used for heatmap values. Default: 'completed'
#' @return A heatmap object displaying the relationship between full names and annotations
#' @export
heat_map_names <- function(cell_types_data, value = 'completed') {


  df_wide <- dcast(cell_types_data, full_names ~ annotation, value.var = value)

  rownames(df_wide) <- df_wide$full_names
  df_wide$full_names <- NULL
  df_wide[is.na(df_wide)] = 0

  pheat <- pheatmap(as.matrix(df_wide),
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    scale = "column",
                    color = colorRampPalette(c("blue", "white", "red"))(50))


  return(pheat)

}


