context("annotation")
test_that("ensembl mapping works", {
  expect_equal(map_ensembl_gene_id(ensembl_ids = "ENSG00000130707",
                                   mappings = "external_gene_name",
                                   species = "human")$external_gene_name,
               "ASS1")
  expect_equal(map_ensembl_gene_id(ensembl_ids = "ENSG00000130707",
                                   species = "human")$ensembl_gene_id,
               "ENSG00000130707")
})
