#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Read and write SDHDF tables"""

__author__ = ["Danny Price", "Alec Thomson"]

from astropy.table import QTable
from astropy.time import Time
from .exceptions import VerificationError
import h5py

# Ignore astropy warnings
import warnings
warnings.filterwarnings('ignore', category=Warning, append=True)

class SDHDFTable(QTable):
    """ Class for sdhdf_table """
    attrs: dict = {}

    def __init__(self, sdhdf_dataset: h5py.Dataset, *args, **kwargs):
        """ Initialise SDHDFTable from h5py dataset """
        if hasattr(sdhdf_dataset, 'attrs'):
            self.attrs = dict(sdhdf_dataset.attrs)
        super().__init__(sdhdf_dataset[:], *args, **kwargs)
        self._apply_units()

    def _apply_units(self):
        """ Read attributes and convert columns into astropy Quanitites where possible """
        for col in self.colnames:
            if col + '_UNIT' in self.attrs:
                try:
                    if col == 'MJD':
                        # NOTE: MJD currently expressed as Quantity (days), should be Time (mjd)
                        # TODO: Fix in subsequent SDHDF version
                        self[col] = Time(self[col], format='mjd')
                    else:
                        self[col].unit = self.attrs[col + '_UNIT']
                        self[col] = self[col].astype('float64')
                except:
                    # Handling of time. This is a bit messy.
                    # First, if something like ELAPSED_TIME has a unit,
                    # e.g. seconds, then we treat it as a quantity
                    # Second, we currently have to guess based on column names
                    # TODO: Fix in susbsequent versions of SDHDF!
                    try:
                        if 'TIME' in col or 'DATE' in col:
                            self[col] = Time(self[col], scale='utc')
                        elif 'UTC' in col or 'AEST' in col:
                            pass  # This is a HH:MM:SS, can't convert into astropy Time()
                        if 'MJD' in col:
                            self[col] = Time(self[col], format='mjd')
                    except:
                        raise VerificationError(f"Cannot convert {self.name} {col} into astropy time")