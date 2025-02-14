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

## Workflow
Maybe not all samples are avilable and the anlysis is expensive, so it makes sense to track and control which stages of the analysis have been done to avoid repeating already performed steps in the same samples.

The workflow is tracked by default in the `data/state` directory. This location can be change in the `state_path` in the file `config.yaml`, under the section `state`.

To start an analysis from scratch, it is needed to run the following in a R console:
```R
source("workflow/utils/status.R")
internal_create_state_file()
```

The function `internal_create_state_file` creates a new empty json file with the initial metadata and increase the version of the file. So previous tracking files can be kept.