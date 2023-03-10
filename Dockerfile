#LABEL maintainer="Combiz Khozoie, Ph.D. c.khozoie@imperial.ac.uk, Alan Murphy, a.murphy@imperial.ac.uk"

## Use rstudio installs binaries from RStudio's RSPM service by default,
## Uses the latest stable ubuntu, R and Bioconductor versions. Created on unbuntu 20.04, R 4.0 and BiocManager 3.12
FROM rocker/rstudio:4.2


## Add packages dependencies
RUN apt-get update \
	&& apt-get install -y --no-install-recommends apt-utils \
	&& apt-get install -y --no-install-recommends \
	## Basic deps
	gdb \
	libxml2-dev \
	python3-pip \
	libz-dev \
	liblzma-dev \
	libbz2-dev \
	libpng-dev \
	libgit2-dev \
	## sys deps from bioc_full
	pkg-config \
	fortran77-compiler \
	byacc \
	automake \
	curl \
	## This section installs libraries
	libpcre2-dev \
	libnetcdf-dev \
	libhdf5-serial-dev \
	libfftw3-dev \
	libopenbabel-dev \
	libopenmpi-dev \
	libxt-dev \
	libudunits2-dev \
	libgeos-dev \
	libproj-dev \
	libcairo2-dev \
	libtiff5-dev \
	libreadline-dev \
	libgsl0-dev \
	libgslcblas0 \
	libgtk2.0-dev \
	libgl1-mesa-dev \
	libglu1-mesa-dev \
	libgmp3-dev \
	libhdf5-dev \
	libncurses-dev \
	libbz2-dev \
	libxpm-dev \
	liblapack-dev \
	libv8-dev \
	libgtkmm-2.4-dev \
	libmpfr-dev \
	libmodule-build-perl \
	libapparmor-dev \
	libprotoc-dev \
	librdf0-dev \
	libmagick++-dev \
	libsasl2-dev \
	libpoppler-cpp-dev \
	libprotobuf-dev \
	libpq-dev \
	libperl-dev \
	## software - perl extentions and modules
	libarchive-extract-perl \
	libfile-copy-recursive-perl \
	libcgi-pm-perl \
	libdbi-perl \
	libdbd-mysql-perl \
	libxml-simple-perl \
	libmysqlclient-dev \
	default-libmysqlclient-dev \
	libgdal-dev \
	## new libs
	libglpk-dev \
	## Databases and other software
	sqlite \
	openmpi-bin \
	mpi-default-bin \
	openmpi-common \
	openmpi-doc \
	tcl8.6-dev \
	tk-dev \
	default-jdk \
	imagemagick \
	tabix \
	ggobi \
	graphviz \
	protobuf-compiler \
	jags \
	## Additional resources
	xfonts-100dpi \
	xfonts-75dpi \
	biber \
	libsbml5-dev \
	## qpdf needed to stop R CMD Check warning
	qpdf \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*


RUN install2.r -e \
argparse \
assertthat \
BiocManager \
cli \
covr \
cowplot \
data.table \
devtools \
DirichletReg \
dplyr \
DT \
english \
enrichR \
forcats \
formattable \
future \
future.apply \
gdtools \
ggdendro \
ggplot2 \
ggpubr \
ggrepel \
ggridges \
Hmisc \
httr \
ids \
igraph \
knitr \
leaflet \
lme4 \
magrittr \
Matrix \
paletteer \
patchwork \
plyr \
prettydoc \
purrr \
qs \
R.utils \
RANN \
rcmdcheck \
Rcpp \
RcppArmadillo \
RcppEigen \
RcppParallel \
RcppProgress \
remotes \
rlang \
rliger \
rmarkdown \
Rtsne \
scales \
sctransform \
Seurat \
snow \
spelling \
stringr \
testthat \
threejs \
tibble \
tidyr \
tidyselect \
tidyverse \
UpSetR \
utils \
vroom \
WebGestaltR

## Install Bioconductor packages
COPY ./misc/requirements-bioc.R .
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
   libfftw3-dev \
   gcc && apt-get clean \
 && rm -rf /var/lib/apt/lists/*
RUN Rscript -e 'requireNamespace("BiocManager"); BiocManager::install(ask=F);' \
&& Rscript requirements-bioc.R

## Install from GH the following
RUN installGithub.r chris-mcginnis-ucsf/DoubletFinder \
ropensci/plotly \
cole-trapnell-lab/monocle3 \
theislab/kBET \
jlmelville/uwot \
NathanSkene/One2One \
hhoeflin/hdf5r \
mojaveazure/loomR \
ropensci/bib2df \
cvarrichio/Matrix.utils \
neurogenomics/scFlowExamples \
neurogenomics/scFlowData

## Install scFlow package
# Copy description
WORKDIR scFlow
ADD . .

# Run R CMD check - will fail with any errors or warnings
RUN Rscript -e "devtools::check()"
# Install R package from source
RUN Rscript -e "remotes::install_local()"
RUN rm -rf *
