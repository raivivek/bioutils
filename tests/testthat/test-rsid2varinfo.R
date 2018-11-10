context('rsid2varinfo')

df_test1 <- data.frame(
    "rs_id" = c("rs8050136", "rs780094", "rs1694", "rs66698963"),
    "chr_name" = c(16, 2, 18, 11),
    "chrom_start" = c(53816275, 27741237, 25398056, 61602537),
    "ref_allele" = c("C", "T", "A", "ACTTCTCCCTGCCTCCCCAGGG"),
    "alt_allele" = c("A", "C", "T", "-")
)
row.names(df_test1) <- df_test1$rs_id


test_that('rsid2varinfo returns correct data.fram', {
    expect_equal(
        rsid2varinfo(c("rs8050136", "rs780094", "rs1694", "rs66698963")),
        df_test1
    )
})
