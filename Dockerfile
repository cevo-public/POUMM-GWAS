# Base image https://hub.docker.com/u/rocker/
FROM rocker/r-base:latest

## install R-packages
COPY scripts/R/install_packages.R .
RUN Rscript install_packages.R

## create output directory
RUN mkdir -p output

## copy analysis scripts
COPY scripts/ scripts

## run script
CMD Rscript scripts/R/simulate_for_estimator_accuracy.R