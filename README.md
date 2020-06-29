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

The python modules require the following dependencies to be installed:
* Python 3

The code has been tested with the following library versions installed:
* liberfa 1.3.0-2
* libhdf5 1.10.0
* pgplot5 5.2.2-19.3
* calceph 3.0.0

## Branches

The master branch requires the ERFA library (SOFA replacement) by default
Checkout the 'sofa\_build' branch if you need to use the SOFA library

## Test data

Accompanying this repository is the sdhdf\_test\_data repository:
* https://bitbucket.csiro.au/scm/cpda/sdhdf\_test\_data.git

## More information

If you have any comments/queries, contact lawrence.toomey@csiro.au or george.hobbs@csiro.au

