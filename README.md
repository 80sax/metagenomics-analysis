# METAGENOMIC SEQUENCING DATA ANALYSIS

This project aims to analyse metagenomic sequencing data, specifically for the microbiome in coastal water samples

## Setup - Package installation
The project uses renv to create a virtual R environment, there are two methods to setup the environment:
- In the terminal run:
```bash
Rscript setup.R --prepare
```
- Open an R terminal and run:
```R
renv::restore()
```

## Run unit tests
- Unit tests are implemented with the [testthat](https://testthat.r-lib.org/) framework.
- Test are located in the folder `tests`
- To execute the tests run:
```bash
./run_tests.sh
```