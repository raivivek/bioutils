# disable strings as factors, but re-enable upon exit
old <- options(stringsAsFactors = FALSE)
on.exit(options(old), add = TRUE)

#' rsid to variant info.
#'
#' \code{rsid2varinfo} converts rs ids to chromosome position, the reference
#' allele, and the alternate allele.
#'
#' @param rs_ids Character vector.
#'     List of rs gene ids to get information for.
#' @param host Character.
#'     Ensembl biomaRt host.
#' @param var_attributes Character vector.
#'     List of attributes to select from Ensembl.
#'
#' @return Data.frame.
#'     Data.frame of the attributes requested.
#'
#' @examples
#' rsid2varinfo(c("rs8050136", "rs780094", "rs16942"))
#'
#' @importFrom biomaRt useMart
#' @importFrom biomaRt getBM
#' @importFrom plyr rename
#' @export
rsid2varinfo <- function(
        rs_ids,
        host = "grch37.ensembl.org",
        var_attributes = c("refsnp_id", "chr_name", "chrom_start", "allele")
    ) {

    rs_ids <- as.character(rs_ids)

    # put this in a tryCatch block in case biomaRt is down!
    snp_info_tmp <- tryCatch({
        ensembl <- biomaRt::useMart(
            biomart = "ENSEMBL_MART_SNP",
            host = host,
            dataset = "hsapiens_snp"
        )
        #biomaRt::listAttributes(ensembl)
        #allele_1 is the ancestral allele
        #cols <- c("refsnp_id", "chr_name", "chrom_start", "allele", "allele_1")
        snp_info_tmp <- biomaRt::getBM(
            attributes = var_attributes,
            filters = "snp_filter",
            values = as.character(rs_ids),
            mart = ensembl
        )
        # remove funky instances
        snp_info_tmp <- subset(snp_info_tmp, !(grepl("HSCHR", chr_name)))
        snp_info_tmp <- subset(snp_info_tmp, !(grepl("NOVEL_TEST", chr_name)))
    }, error = function(cond) {
        warning(paste0("Problems with biomaRt, probably due to a down",
                       " website. Returning data.frame of rs_ids."))
        return_df <- data.frame("rs_id" = rs_ids)
        row.names(return_df) <- rs_ids
        return(return_df)
    })

    # case where no results
    if (nrow(snp_info_tmp) == 0) {
        return_df <- data.frame("rs_id" = rs_ids)
        row.names(return_df) <- rs_ids
    } else {
        snp_info_tmp <- merge(
            data.frame("refsnp_id" = rs_ids),
            snp_info_tmp,
            by = c("refsnp_id"),
            all.x = T,
            all.y = T
        )
        # the "allele" column from ensembl corresponds to the reference /
        # alternate alleles according to the forward strand
        if ("allele" %in% colnames(snp_info_tmp)) {
            allele_df <- stringr::str_split_fixed(
                snp_info_tmp[["allele"]], '\\/', 2
            )
            colnames(allele_df) <- c("ref_allele", "alt_allele")
            snp_info_tmp$allele <- NULL
            snp_info_tmp <- cbind(snp_info_tmp, allele_df)
        }
        row.names(snp_info_tmp) <- as.character(snp_info_tmp$refsnp_id)
        return_df <- plyr::rename(snp_info_tmp, c("refsnp_id" = "rs_id"))
    }

    return(return_df[rs_ids, ])
}

