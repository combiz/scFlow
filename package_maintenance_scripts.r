#install.packages(c("devtools", "roxygen2", "usethis", "spelling"))
library(devtools)
library(roxygen2)
library(usethis)
library(spelling)


## Initial package creation steps
#create_package("~/Documents/scflow")
#options(usethis.full_name = "Combiz Khozoie")
#use_git()
#use_github()
#use_mit_license()
#use_readme_rmd()
#use_travis()
#use_testthat()
#use_spell_check()
#use_news_md()
#use_pkgdown()

## Ongoing package dev
use_r("function_name")
use_test("function_name")

use_package("dplyr") # add to DESCRIPTION
#use_vignette("How to do this cool analysis") # later

### BEFORE EVERY COMMIT
#Restart R Session Cmd+Shift+F10 (Ctrl+Shift+F10 for Windows)
#Document Package Cmd+Shift+D (Ctrl+Shift+D for Windows)
#Check Package Cmd+Shift+E (Ctrl+Shift+E for Windows)

## before every release
# knit the readme.Rmd
use_version()
pkgdown::build_site()
