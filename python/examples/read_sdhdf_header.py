#!/usr/bin/env python
import argparse

import h5py
import pandas as pd
from astropy.table import QTable

__version__ = "2.2"
__author__ = "Lawrence Toomey"


def print_header(tb):
    """
    Format SDHDF header metadata output

    :param astropy.QTable tb: astropy.QTable metadata object
    :return: None
    """
    params = []
    for col in tb.colnames:
        params.append([col, tb[col][0]])
    df = pd.DataFrame(params, columns=("-- Key --", "-- Value --"))
    print(df)


def read_sdhdf_header(f):
    """
    Read SDHDF primary header metadata

    :param string f: Path to SDHDF file
    :return: None
    """
    try:
        with h5py.File(f, "r") as h5:
            tb = QTable.read(h5, path="/metadata/primary_header")
            print("Displaying primary header for file:\n%s\n" % f)
            print_header(tb)
    except Exception as e:
        print("ERROR: failed to read file %s" % f, e)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--filename", help="Path to SDHDF file to read [required]", required=True
    )
    args = ap.parse_args()

    read_sdhdf_header(args.filename)
