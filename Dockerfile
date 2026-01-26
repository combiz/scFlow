#LABEL maintainer="Combiz Khozoie, Ph.D. c.khozoie@imperial.ac.uk, Alan Murphy, a.murphy@imperial.ac.uk"

## Use rstudio installs binaries from RStudio's RSPM service by default,
## Uses the latest stable ubuntu, R and Bioconductor versions. Created on unbuntu 20.04, R 4.3 and BiocManager 3.18
FROM rocker/rstudio:4.3


## Add packages dependencies
RUN apt-get update \
	&& apt-get install -y --no-install-recommends apt-utils \
	&& apt-get install -y --no-install-recommends \
	## Basic deps
	build-essential \
	gfortran \
	libcurl4-openssl-dev \
    libssl-dev \
    libzstd-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
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
	curl \
	## This section installs libraries
	libpcre2-dev \
	libfftw3-dev \
	libopenmpi-dev \
	libxt-dev \
	libudunits2-dev \
	libgeos-dev \
	libproj-dev \
	libcairo2-dev \
	libtiff5-dev \
	libreadline-dev \
	libgsl-dev \
	libgl1-mesa-dev \
	libglu1-mesa-dev \
	libhdf5-dev \
	libncurses-dev \
	liblapack-dev \
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
	gcc \
	nodejs \
	npm \
	## qpdf needed to stop R CMD Check warning
	qpdf \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

RUN pip install stratocumulus \
&& curl https://sdk.cloud.google.com > install.sh \
&& bash install.sh --disable-prompts \
&& curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" \
-o "awscliv2.zip" \
&& unzip awscliv2.zip \
&& ./aws/install \
&& rm -rf awscliv2.zip \
&& rm -rf /tmp/*

ENV MAKEFLAGS="-j2"
ENV OMP_NUM_THREADS=1

RUN install2.r -e source \
Matrix \
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
knitr \
leaflet \
magrittr \
lme4 \
igraph \
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
SeuratObject \
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
WebGestaltR \
apcluster \
&& rm -rf /tmp/downloaded_packages

## Install Bioconductor packages
COPY ./misc/requirements-bioc.R .
RUN Rscript -e 'requireNamespace("BiocManager"); BiocManager::install(ask=F);' \
&& Rscript requirements-bioc.R \
&& rm -rf /tmp/downloaded_packages

## Install from GH the following
RUN installGithub.r NathanSkene/EWCE \
chris-mcginnis-ucsf/DoubletFinder \
ropensci/plotly \
cole-trapnell-lab/monocle3 \
theislab/kBET \
jlmelville/uwot \
hhoeflin/hdf5r \
ropensci/bib2df \
cvarrichio/Matrix.utils \
&& rm -rf /tmp/downloaded_packages

## Install scFlow package
# Copy description
WORKDIR scFlow
ADD . .

# Run R CMD check - will fail with any errors or warnings
RUN Rscript -e "devtools::check(vignettes = FALSE)"
# Install R package from source
RUN Rscript -e "remotes::install_local()"
RUN rm -rf *
