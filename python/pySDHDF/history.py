#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SDHDF history utilities"""

from dataclasses import dataclass
import datetime
import inspect
from pathlib import Path
import pkg_resources
import socket
from typing import List, Optional, Tuple, Union
import warnings

import h5py
import matplotlib.pyplot as plt
from IPython import embed
import numpy as np
import pandas as pd
from astropy.stats import sigma_clip, mad_std
from astropy.table import Table
from tqdm.auto import tqdm
from xarray import DataArray, Variable, Dataset

def generate_history_row(
) -> pd.DataFrame:
    """Generate a history row.
    Returns:
        pd.DataFrame: History row.
    """

    # Get the calling function from inspect
    process_name = inspect.stack()[1][3]
    # Get the calling function's docstring
    process_description = inspect.stack()[1][0].f_locals["self"].__doc__
    # Get the calling function's arguments
    process_arguments = str(inspect.stack()[1][0].f_locals["self"].__dict__)

    history_row = pd.DataFrame(
        {
            "DATE": datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S"),
            "PROC": process_name,
            "PROC_DESCR": process_description,
            "PROC_ARGS": process_arguments,
            "PROC_HOST": socket.getfqdn(),
        },
        index=[0]
    )
    return history_row