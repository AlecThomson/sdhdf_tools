README for sdhdfProc
=======================

# About
---
**sdhdfProc** is a software package to read, manipulate and process radio astronomy data in Single Dish Hierarchical Data Format (SDHDF) 

**Author**    George.Hobbs@csiro.au  
**Copyright** CSIRO 2020 

## Prerequisites
---
* PGPLOT
* HDF5 library
* SOFA library
* calceph library

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
