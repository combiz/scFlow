context("annotation")
test_that("ensembl mapping works", {
  expect_equal(map_ensembl_gene_id(ensembl_ids = "ENSG00000130707",
                                   mappings = "external_gene_name")$external_gene_name,
               "ASS1")
  expect_equal(map_ensembl_gene_id(ensembl_ids = "ENSG00000130707")$ensembl_gene_id,
               "ENSG00000130707")
})
