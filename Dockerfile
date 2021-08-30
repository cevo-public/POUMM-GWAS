# Base image https://hub.docker.com/u/rocker/
FROM rocker/r-base:latest

## install R-packages
COPY install_packages.R .
RUN Rscript install_packages.R

## create directories
RUN mkdir -p functions
RUN mkdir -p analysis
RUN mkdir -p output

## copy files for analysis
COPY functions/* functions
COPY analysis/simulate_for_estimator_accuracy.R analysis/simulate_for_estimator_accuracy.R

## run script
CMD Rscript analysis/simulate_for_estimator_accuracy.R