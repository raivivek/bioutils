context('ensembl2hgnc')
library(bioutils)

test_that('ensembl2hgnc returns correct hgnc', {
    expect_equal(
        ensembl2hgnc('ENSG00000109501'),
        c('ENSG00000109501' = 'WFS1')
    ) # ignore.case = TRUE
})

test_that('ensembl2hgnc returns correct hgnc', {
    expect_equal(
        ensembl2hgnc(c('ENSG00000058453.12', 'ENSG00000030066.9')),
        c('ENSG00000058453.12' = 'CROCC', 'ENSG00000030066.9' = 'NUP160')
    )
})
