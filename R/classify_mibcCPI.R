#' Single-sample classifier of CPI-response related classes
#'
#' Classifies RNA-seq from muscle-invasive bladder cancer samples normalized
#' with log2(FPKM+1) into one of 5 classes
#' proposed by Robertson et al (2021, doi: tba), related to response to
#' neoadjuvant treatment with checkpoint inhibitors (CPI).
#'
#' @param gexp gene expression matrix, normalized by log2(FPKM+1)
#' @param gene_id type of gene identifier. One of: 'hgnc_symbol',
#' 'ensembl_gene_id', 'entrezid'. Default: 'hgnc_symbol
#'
#' @return A data.frame with classification predictions and probabilities.
#' @importFrom stringr str_replace_all
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate select rename relocate rename_all case_when pull
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom stats predict
#' @import glmnet caret
#' @examples
#' data(tcga_blca_ex)
#' classify_mibcCPI(tcga_blca_ex)
#' @export
classify_mibcCPI <- function(gexp,
                        gene_id = c("hgnc_symbol", "ensembl_gene_id", "entrezid")[1]) {
    #-- Transform into symbol
    if(gene_id != "hgnc_symbol") {
        ids <- gene_annot %>% pull(hgnc_symbol, {{gene_id}})
        ids_present <- intersect(rownames(gexp), names(ids))
        gexp <- gexp[ids_present,]
        rownames(gexp) <- ids[rownames(gexp)]
    }

    #-- Fix names
    rownames(gexp) <- str_replace_all(rownames(gexp), "-", "_")

    #-- Address missing
    missing <- feats[! feats %in% rownames(gexp)]
    val <- mean(as.matrix(gexp[feats[feats %in% rownames(gexp)],]), na.rm = TRUE)

    if(length(missing) > 0) {
        msg <- paste("There are markers missing from this dataset.
Missing genes:", paste(missing, collapse = ", "))
        warning(msg)
    }

    mss_mat <- matrix(val, nrow = length(missing), ncol = ncol(gexp),
                      dimnames = list(missing, colnames(gexp)))
    gexp <- rbind(gexp, mss_mat)
    #-- Rearrange
    data <- gexp[feats,] %>% t() %>% as.data.frame()
    #-- Make prediction
    pred_cls <- predict(glm_model, data)
    probs_cls <- predict(glm_model, data, type = "prob") %>%
        rename_all(~ paste0("probability_s", .))
    #-- Make return
    class_res <- cbind(prediction = pred_cls, probs_cls) %>%
        as.data.frame() %>%
        rownames_to_column("id") %>%
        as_tibble()

    class_res <- class_res %>%
        mutate(prediction = factor(paste0("S", prediction), levels = c("S1", "S2", "S4", "S5", "S3")),
               class_name = case_when(
                   prediction == "S1" ~ "LumE",
                   prediction == "S2" ~ "LumR",
                   prediction == "S4" ~ "MycU",
                   prediction == "S5" ~ "StroR",
                   prediction == "S3" ~ "ImmBas"
               ) %>% factor(levels = c("LumE", "LumR", "MycU", "StroR", "ImmBas"))) %>%
        relocate(class_name, .after = prediction)

    return(class_res)
}
#' Get mibcCPI classification features
#'
#' Returns a vector of 100 gene symbols with the classification features
#'
#' @param gene_id type of gene identifier. One of: 'hgnc_symbol',
#' 'ensembl_gene_id', 'entrezid'. Default: 'hgnc_symbol
#'
#' @examples
#' classification_features()
#' @export
classification_features <- function(gene_id = c("hgnc_symbol", "ensembl_gene_id", "entrezid")[1]) {
    if(!gene_id %in% c("hgnc_symbol", "ensembl_gene_id", "entrezid")) {
        stop('`gene_id` must be one of: "hgnc_symbol", "ensembl_gene_id", "entrezid"')
    }
    class_feats <- as.data.frame(gene_annot)[,gene_id]

    return(class_feats)
}

#' A 100-gene, 81 sample subset of the TCGA-BLCA dataset
#'
#' @format A matrix with 100 rows (genes) and 81 columns (samples).
#' Downloaded using TCGAbiolinks package on 2021-07-09, and normalized as:
#' log2(FPKM+1), with gene IDs mapped to HGNC symbols using Ensembl version 101.
#' @references Robertson et al., Cell, 2017 doi: 10.1016/j.cell.2017.09.007
"tcga_blca_ex"
