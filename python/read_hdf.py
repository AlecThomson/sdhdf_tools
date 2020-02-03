# list data and data types of an HDF DataSet

import h5py
import pandas as pd
from subprocess import Popen, PIPE


def read_hdf(f, dset):
    h5 = h5py.File(f, 'r')
    ds = h5[dset]
    #attr = ds.attrs
    df = pd.DataFrame(pd.read_hdf(f, key=dset))
    process = Popen(['h5dump', '-H', '-d', dset, f],
                    stdout=PIPE,
                    stderr=PIPE,
                    universal_newlines=True)
    stdout, stderr = process.communicate()
    print("------------------------------------------------------------------")
    print("Data:")
    print("------------------------------------------------------------------")
    print(df)
    print("------------------------------------------------------------------")
    print("Python data types:")
    print("------------------------------------------------------------------")
    print(df.dtypes)
    print("------------------------------------------------------------------")
    print("HDF data types:")
    print(stdout)
    print("------------------------------------------------------------------")
    print("HDF Attributes:")
    print("------------------------------------------------------------------")


    h5.close()


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--filename",
                    help="Path to HDF file to read",
                    required=True)
    ap.add_argument("--dset",
                    help="HDF dataset to read",
                    required=True)
    args = ap.parse_args()
    read_hdf(args.filename, args.dset)
