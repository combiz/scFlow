#!/usr/bin/env Rscript
# Perform impacted pathway analysis on the differential expression result
#  Nurun Fancy <n.fancy@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

options(mc.cores = parallel::detectCores())

##  ............................................................................
##  Load packages                                                           ####
library(argparse)
library(scFlow)
library(cli)

##  ............................................................................
##  Parse command-line arguments                                            ####

# create parser object
parser <- ArgumentParser()

# specify options
required <- parser$add_argument_group("Required", "required arguments")
optional <- parser$add_argument_group("Optional", "required arguments")

required$add_argument(
  "--gene_file",
  help = "full path to the gene file",
  metavar = ".tsv",
  required = TRUE,
  default = "NULL"
)

required$add_argument(
  "--enrichment_tool",
  help = "one or more enrichment tools",
  metavar = "WebGestaltR",
  required = TRUE,
  default = "WebGestaltR"
)

required$add_argument(
  "--enrichment_method",
  help = "name of the enrichment method used for webgestaltr",
  metavar = "ORA,GSEA",
  required = TRUE,
  default = "ORA"
)

required$add_argument(
  "--enrichment_database",
  help = "name of the enrichment databases",
  metavar = "GO_Biological_Process,GO_Cellular_Component,GO_Molecular_Function",
  required = TRUE,
  default = "KEGG"
)



### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####

args <- parser$parse_args()

args$enrichment_method <- strsplit(args$enrichment_method, ",")[[1]]
args$enrichment_tool <- strsplit(args$enrichment_tool, ",")[[1]]
args$enrichment_database <- strsplit(args$enrichment_database, ",")[[1]]
args$gene_file <- strsplit(args$gene_file, ",")[[1]]
args <- purrr::map(args, function(x) {
  if (length(x) == 1) {
    if (toupper(x) == "TRUE") {
      return(TRUE)
    }
    if (toupper(x) == "FALSE") {
      return(FALSE)
    }
    if (toupper(x) == "NULL") {
      return(NULL)
    }
  }
  return(x)
})

##  ............................................................................
##  Start impacted pathway analysis(IPA)                                    ####

output_dir <- file.path(getwd(), "ipa")
report_dir <- file.path(getwd())

dir.create(output_dir)
dir.create(report_dir)

for (gene_file in args$gene_file) {
  enrichment_result <- find_impacted_pathways(
    gene_file = gene_file,
    enrichment_tool = args$enrichment_tool,
    enrichment_method = args$enrichment_method,
    enrichment_database = args$enrichment_database,
    is_output = TRUE,
    output_dir = output_dir
  )
  report_name <-  tools::file_path_sans_ext(gene_file)
  report_fp <- paste0(report_name, "_scflow_ipa_report")
  report_impacted_pathway(
    res = enrichment_result,
    report_folder_path = report_dir,
    report_file = report_fp
  )
  cli::cli_text(c(
    "{cli::col_green(symbol$tick)} Analysis complete, output is found at: ",
    "{.file {output_dir}}"
  ))
}
