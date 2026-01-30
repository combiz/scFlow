#LABEL maintainer="Nurun Fancy, n.fancy@imperial.ac.uk Combiz Khozoie, Ph.D. c.khozoie@imperial.ac.uk, Alan Murphy, a.murphy@imperial.ac.uk"

## Use rstudio installs binaries from RStudio's RSPM service by default,
## Uses the latest stable ubuntu, R and Bioconductor versions. Created on unbuntu 20.04, R 4.3 and BiocManager 3.18
FROM rocker/rstudio:4.5


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
	pkg-config \
	curl \
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
	libhwloc-dev \
	libopenblas-dev \
	libncurses-dev \
	liblapack-dev \
	libmagick++-dev \
	libsasl2-dev \
	libpoppler-cpp-dev \
	libprotobuf-dev \
	libpq-dev \
	libperl-dev \
	libarchive-extract-perl \
	libfile-copy-recursive-perl \
	libcgi-pm-perl \
	libdbi-perl \
	libdbd-mysql-perl \
	libxml-simple-perl \
	libmysqlclient-dev \
	default-libmysqlclient-dev \
	libgdal-dev \
	libglpk-dev \
	sqlite3 \
    libsqlite3-dev \
	openmpi-bin \
	tcl8.6-dev \
	tk-dev \
	default-jdk \
	imagemagick \
	tabix \
	graphviz \
	protobuf-compiler \
	jags \
	xfonts-100dpi \
	xfonts-75dpi \
	biber \
	libsbml5-dev \
	gcc \
	nodejs \
	npm \
	qpdf \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

RUN pip install --break-system-packages stratocumulus \
&& curl https://sdk.cloud.google.com > install.sh \
&& bash install.sh --disable-prompts \
&& curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" \
-o "awscliv2.zip" \
&& unzip awscliv2.zip \
&& ./aws/install \
&& rm -rf awscliv2.zip \
&& rm -rf /tmp/*

ENV MAKEFLAGS="-j1"
ENV OMP_NUM_THREADS=1
ARG GITHUB_PAT
ENV GITHUB_PAT=${GITHUB_PAT}

## Install CRAN, Bioconductor and github packages 
COPY ./misc/requirements-bioc.R .

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
HighFive \
ids \
knitr \
leaflet \
magrittr \
lme4 \
paletteer \
patchwork \
plyr \
prettydoc \
purrr \
R.utils \
RANN \
rcmdcheck \
rmarkdown \
Rtsne \
scales \
sctransform \
Seurat \
SeuratObject \
snow \
spelling \
speedglm \
stringr \
testthat \
threejs \
tidyr \
tidyverse \
UpSetR \
vroom \
WebGestaltR \
apcluster \
&& Rscript -e "install.packages('RcppPlanc', repos = c( \
  linux = 'https://welch-lab.r-universe.dev/bin/linux/noble/4.5/', \
  sources = 'https://welch-lab.r-universe.dev', \
  cran = 'https://cloud.r-project.org' \
))" \
&& Rscript requirements-bioc.R \
&& installGithub.r NathanSkene/EWCE \
qsbase/qs \
chris-mcginnis-ucsf/DoubletFinder \
ropensci/plotly \
cran/grr \
bnprks/BPCells/r \
cole-trapnell-lab/monocle3 \
theislab/kBET \
jlmelville/uwot \
hhoeflin/hdf5r \
ropensci/bib2df \
cvarrichio/Matrix.utils \
welch-lab/liger \
&& rm -f requirements-bioc.R \
&& rm -rf /tmp/downloaded_packages

## Install scFlow package
WORKDIR /home/rstudio/scFlow
ADD . .

# Run R CMD check - will fail with any errors or warnings
RUN Rscript -e "devtools::check(vignettes = FALSE)" \
&& Rscript -e "remotes::install_local()" \
&& apt-get purge -y build-essential gfortran gcc \
&& apt-get autoremove -y \
&& apt-get clean \
&& rm -rf *