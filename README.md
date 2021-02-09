# NEON-NIWO-misclass

Here we develop a joint classification-occupancy model for multispecies datasets. We illustrate the model using carabid pitfall trap data from the National Ecological Observatory Network (NEON).

## Prerequisites

To run this compendium, you'll need R installed

## Quickstart

This compendium produces the article's main figures. It does not produce the paper.

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

A place for `raw/` and `clean/` data. Cleaned data are ignored by version control.

### `figures/` 

A location for generated figures. 

### `output/` 

JAGS model output

### `source/` 

custom functions

### `Makefile` 

This file contains instructions for how to execute the compendium. 

