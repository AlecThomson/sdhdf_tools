README for INSPECTA
===================


# About
---
**INSPECTA** is an 'INtegrated SDHDF Processing Engine in C for Telescope data Analysis'. Formerly known as 'sdhdfProc', INSPECTA is a software package to read, manipulate and process radio astronomy data in Spectral-Domain Hierarchical Data Format (SDHDF).

**Author**    George.Hobbs@csiro.au
**Copyright** CSIRO 2020


## Prerequisites
---
* PGPLOT (http://astro.caltech.edu/~tjp/pgplot/)
* HDF5 library (http://hdfgroup.org/HDF5/)
* ERFA library (https://github.com/liberfa/erfa)
* calceph library (https://www.imcce.fr/recherche/equipes/asd/calceph/)


## Installing
---
```
./bootstrap
./configure
make
make install
```


## Routines
---
### sdhdf_describe
### sdhdf_extract
### sdhdf_quickdump
### sdhdf_onoff
### sdhdf_sum
### sdhdf_modify
### sdhdf_identify
### sdhdf_cal
### sdhdf_listlines
### sdhdf_primaryFluxCal
### sdhdf_gainCurve
### sdhdf_applyCal
### sdhdf_calProc
### sdhdf_pointing
### sdhdf_autoFlag
### sdhdf_flag
### sdhdf_plotScan
### sdhdf_plotWide
### sdhdf_plotMultiSpec
### sdhdf_plotSpectrum
