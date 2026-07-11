# syntax=docker/dockerfile:1

FROM bioconductor/bioconductor_docker:devel
LABEL Name=spiDE Version=0.99.0
RUN mkdir /root/rpkg
COPY . /root/rpkg
RUN R -e "devtools::install('/root/rpkg', dependencies=T)"
