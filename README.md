# NEON-NIWO-misclass

Here we develop a joint classification-occupancy model for multispecies survey datasets. We illustrate the model using carabid pitfall trap data from the National Ecological Observatory Network (NEON) collected at Niwot Ridge, CO, USA.

## Publication

https://doi.org/10.1111/2041-210X.13858

The manuscript describing this work is published in *Methods in Ecology & Evolution*

**Estimating species misclassification with occupancy dynamics and encounter rates: A semi-supervised, individual-level approach**
Anna I. Spiers, J. Andrew Royle, Christa L. Torrens, Maxwell B. Joseph

Code is archived on Zenodo at https://doi.org/10.5281/zenodo.6394878

## Prerequisites

To run this compendium, you'll need R installed.

## Quickstart

This compendium produces the article's figures based off of the carabid case study results. 

It does not produce the paper nor simulation output. For simulation output, run ```simulations.R```. This will take a long time depending on your CPU bandwidth. For simulation figures, run ```sim_figs.R```

### Navigate to repository

Download the repository. In terminal, navigate to the repository.

### Execute the compendium

To run the analysis run GNU Make from the terminal:

```bash
make
```

## Details

Here's what's inside: 

### `data/`

A place for local storage of publicly available NEON data. Data are ignored by version control.

### `figures/` 

A location for generated figures. 

### `output/` 

JAGS model output

### `source/` 

custom functions

### `Makefile` 

This file contains instructions for how to execute the compendium. 

