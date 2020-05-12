# README for sdhdf\_tools

## About

This repository contains a collection of tools for interrogating
Single Dish Hierarchical Data Format (SDHDF) files

## Prerequisites

The C code requires the following dependencies to be installed:
* ERFA (https://github.com/liberfa/erfa) or the SOFA library (http://www.iausofa.org/)
* HDF5 library (https://www.hdfgroup.org/downloads/hdf5) 
* PGPLOT library (http://www.astro.caltech.edu/~tjp/pgplot/)
* calceph library (https://www.imcce.fr/recherche/equipes/asd/calceph)

The python modules require Python 3

## Branches

The master branch requires the ERFA library (SOFA replacement) by default
Checkout the 'sofa\_build' branch if you need to use the SOFA library

## Test data

Accompanying this repository is the sdhdf\_test\_data repository:
* https://bitbucket.csiro.au/scm/cpda/sdhdf\_test\_data.git

## More information

If you have any comments/queries, contact lawrence.toomey@csiro.au or george.hobbs@csiro.au

