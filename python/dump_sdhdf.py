import h5py
from astropy.io import ascii
from astropy.table import Table
import numpy as np
import re

__version__ = '1.7'
__author__ = 'Lawrence Toomey'


def dump_sdhdf(f, out_format):

    with h5py.File(f, 'r') as h5:
        for dset in h5.keys():
            print(h5[dset].keys())
            if re.search('data', str(dset)):
                data = np.array(h5[dset])
                n_pol = data.shape[2]
                n_chan = data.shape[3]
                sl_arr = np.empty((n_pol, n_chan))
                for j in range(len(data[2])):
                    #sl = np.array((1, 2))
                    print(data[0, 0, j, :].shape)
                    sl_arr[j] = data[0, 0, j, :]
                    #np.mean(sl, axis=0)
                    #bb = np.ones((1, 2))
                    #bb = np.reshape((-1))
                    #np.savetxt(str(dset) + '.' + out_format, bb, fmt='%f')
                
                sl_arr.reshape((n_pol, n_chan))
                np.savetxt(str(dset) + '.' + out_format, sl_arr, fmt='%f')

                #print(data[0, 0, 0, :])
                #for sl in data:
                #    print(sl.shape)
                    #np.savetxt(str(dset) + '.' + out_format, sl, fmt='%f')
                #np.savetxt(str(dset) + '.' + out_format, np.vstack([data[0,0,0,:], data[0,0,1,:]]), delimiter='\n', fmt='%f')
                #data.tofile(str(dset) + '.' + out_format, sep=' ', format='%f')
            if re.search('frequency_SB', str(dset)):
                print(dset)
                data = np.array(h5[dset])
                np.savetxt(str(dset) + '.' + out_format, data[0], delimiter='\n', fmt='%f')


if __name__ == '__main__':
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument('--filename', help='Path to SDHDF file to read', required=True)
    ap.add_argument('--format', help='Output file format [Default ascii]', default='ascii')
    args = ap.parse_args()

    dump_sdhdf(args.filename, args.format)
