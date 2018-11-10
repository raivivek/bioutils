# disable strings as factors, but re-enable upon exit
old <- options(stringsAsFactors = FALSE)
on.exit(options(old), add = TRUE)

#' Ensembl ids 2 hgnc
#'
#' \code{ensembl2hgnc} converts Ensembl gene ids to hgnc symbols. If no
#' hgnc symbol then uses external_gene_name. If no external_gene_name then
#' uses Ensembl gene id.
#'
#' @param ensembl_gene_ids Character vector.
#'     List of Ensembl gene ids to get hgnc symbols for.
#' @param host Character.
#'     Ensembl biomaRt host.
#' @param drop_dot_ensembl_id Logical.
#'     Drop "." from ENSG00000072310.12
#'
#' @return Character.
#'     The corresponding hgnc symbols to the Ensembl ids.
#'
#' @examples
#' ensembl2hgnc(c("ENSG00000109501", "ENSG00000058453", "ENSG00000030066"))
#'
#' @importFrom biomaRt useMart
#' @importFrom biomaRt getBM
#' @export
ensembl2hgnc <- function(
        ensembl_gene_ids,
        host = "grch37.ensembl.org",
        drop_dot_ensembl_id = TRUE
    ) {

    ensembl_gene_ids <- as.character(ensembl_gene_ids)
    ensembl_gene_ids_original <- ensembl_gene_ids
    if (drop_dot_ensembl_id) {
        ensembl_gene_ids <- unlist(lapply(ensembl_gene_ids,
            FUN = function(x) { return(strsplit(x, '\\.')[[1]][1]) }
        ))
    }
    # put this in a tryCatch block in case biomaRt is down!
    gene_info_tmp <- tryCatch({
        ensembl <- biomaRt::useMart(
            biomart = "ENSEMBL_MART_ENSEMBL",
            host = host,
            path = "/biomart/martservice",
            dataset = "hsapiens_gene_ensembl"
        )
        #biomaRt::listAttributes(ensembl)
        #cols <- c("ensembl_gene_id", "hgnc_symbol", "gene_biotype",
        #   "description")
        cols <- c("ensembl_gene_id", "hgnc_symbol", "external_gene_name")
        gene_info_tmp <- biomaRt::getBM(
            attributes = cols,
            mart = ensembl,
            filters = c("ensembl_gene_id"),
            values = unique(ensembl_gene_ids)
        )
    }, error = function(cond) {
        warning(paste0("Problems with biomaRt, probably due to a down",
                       " website. Returning ensembl_ids rather than hgnc_id."))
        return(data.frame())
    })

    # case where no results
    if (nrow(gene_info_tmp) == 0) {
        return_vector <- ensembl_gene_ids_original
        names(return_vector) <- ensembl_gene_ids_original
    } else {
        # make the gene_info_tmp the length of ensembl_gene_ids... so that
        # the return vector is the same length
        # this command also expands cases when there are duplicates.
        gene_info_tmp <- merge(
            data.frame(
                "ensembl_gene_id" = ensembl_gene_ids,
                "ensembl_gene_ids_original" = ensembl_gene_ids_original
            ),
            gene_info_tmp,
            by = c("ensembl_gene_id"),
            all.x = T,
            all.y = T
        )
        # if hgnc_symbol is empty, replace it with external_gene_name
        filt <- is.na(gene_info_tmp$hgnc_symbol) |
            gene_info_tmp$hgnc_symbol == ""
        if (any(filt)) {
            gene_info_tmp$hgnc_symbol[filt] <- as.character(
                gene_info_tmp$external_gene_name[filt]
            )
        }
        gene_info_tmp$external_gene_name <- NULL  # del external_gene_name
        # if hgnc_symbol is *still* empty, replace it with ensembl_gene_id
        filt <- is.na(gene_info_tmp$hgnc_symbol) |
            gene_info_tmp$hgnc_symbol == ""
        if (any(filt)) {
            gene_info_tmp$hgnc_symbol[filt] <- as.character(
                gene_info_tmp$ensembl_gene_id[filt]
            )
        }
        # get the return vector
        return_vector <- gene_info_tmp$hgnc_symbol
        names(return_vector) <- gene_info_tmp$ensembl_gene_ids_original
    }

    return(return_vector[ensembl_gene_ids_original])
}

