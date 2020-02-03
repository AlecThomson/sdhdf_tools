import h5py
import argparse
import glob
import os
from astropy.table import QTable
from tabulate import tabulate

__version__ = '1.9'
__author__ = 'Lawrence Toomey'


def quick_list(dir):
	row_data = []
	hdr = 'File', 'Project_ID', 'Beam', 'N_bands', 'Source'

	for fpath in glob.glob(dir + '/*.hdf'):
		try:
			with h5py.File(fpath, 'r') as h5:
				ph = QTable.read(h5, path='/metadata/primary_header')
				bp = QTable.read(h5, path='/metadata/beam_params')
				f = os.path.basename(fpath)
				nbeams = len(bp)
				pid = ph['PID'][0]

				for row in range(0, nbeams):
					beam = bp['LABEL'][row]
					nbands = bp['N_BANDS'][row]
					source = bp['SOURCE'][row]
					row = f, pid, beam, nbands, source
					row_data.append(row)

		except Exception as e:
			print('ERROR: failed to read file %s' % f, e)

	print(tabulate(row_data, headers=hdr))


if __name__ == '__main__':
	ap = argparse.ArgumentParser()
	ap.add_argument('--dir', help='Path to directory containing SDHDF files [required]', required=True)
	args = ap.parse_args()

	quick_list(args.dir)
