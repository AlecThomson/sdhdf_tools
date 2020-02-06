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
### Note that if any libraries are not present in default paths then use, e.g.:
### ./configure CFLAGS="-I/u/hob044/software/new_c/sdhdfProc/sofa/20190722/c/src/ -I/pulsar/psr/software/current/stretch/include/ -I/u/hob044/software/new_c/hdf5-1.10.4/src/" LDFLAGS="-L/u/hob044/software/new_c//hdf5-1.10.4/src/ -L/u/hob044/software/new_c/sdhdfProc/sofa/20190722/c/src/"
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
