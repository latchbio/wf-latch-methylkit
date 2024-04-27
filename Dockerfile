# DO NOT CHANGE
from 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:fe0b-main

workdir /tmp/docker-build/work/

shell [ \
    "/usr/bin/env", "bash", \
    "-o", "errexit", \
    "-o", "pipefail", \
    "-o", "nounset", \
    "-o", "verbose", \
    "-o", "errtrace", \
    "-O", "inherit_errexit", \
    "-O", "shift_verbose", \
    "-c" \
]
env TZ='Etc/UTC'
env LANG='en_US.UTF-8'

arg DEBIAN_FRONTEND=noninteractive

# Latch SDK
# DO NOT REMOVE
run pip install latch==2.38.8
run mkdir /opt/latch

# Install rig the R installation manager
run \
    curl \
        --location \
        --fail \
        --remote-name \
        https://github.com/r-lib/rig/releases/download/latest/rig-linux-latest.tar.gz && \
    tar \
        --extract \
        --gunzip \
        --file rig-linux-latest.tar.gz \
        --directory /usr/local/ && \
    rm rig-linux-latest.tar.gz

# Install R
run rig add release # Change to any R version

# Install R dependencies
# copy environment.R /opt/latch/environment.R
# run Rscript /opt/latch/environment.R

run pip install pandas
run apt-get install -y libxml2-dev

run R -e 'install.packages("jsonlite")'
run R -e 'library(jsonlite)'
run R -e "install.packages('BiocManager', dependencies=TRUE, repos='http://cran.rstudio.com/')"
run R -e 'BiocManager::install("methylKit", ask=FALSE)'
# run R -e "install.packages('devtools')"
# library(devtools)
# install_github("al2na/methylKit", build_vignettes=FALSE,
#   repos=BiocManager::repositories(),
#   dependencies=TRUE)

# Copy workflow data (use .dockerignore to skip files)
copy . /root/


# Latch workflow registration metadata
# DO NOT CHANGE
arg tag
# DO NOT CHANGE
env FLYTE_INTERNAL_IMAGE $tag

workdir /root
