import numpy as np
import pylab as plt
import pandas as pd
import h5py
import os
import matplotlib as mpl
from astropy.table import QTable

__version__ = '1.9'
__author__ = 'Lawrence Toomey'

# increase matplotlib chunk size above default -
# this improves speed and prevents Agg rendering failure
# when plotting large data sets
mpl.rcParams['agg.path.chunksize'] = 20000


def print_hdr(tb):
    """
        Print header metadata

        :param astropy.QTable tb: astropy.QTable metadata object
        :return None
    """
    params = []
    for col in tb.colnames:
        params.append([col, tb[col][0]])
    df = pd.DataFrame(params, columns=('-- Key --', '-- Value --'))
    print(df)


def read_sdhdf_header(f_pth, dset_pth):
    """
        Read SDHDF header metadata

        :param string f_pth: Path to SDHDF file
        :param HDF dataset dset_pth: Path to HDF metadata dataset
        :return None
    """
    try:
        with h5py.File(f_pth, 'r') as h5:
            tb = QTable.read(h5, path=dset_pth)
            print('Displaying primary header for file:\n%s\n' % f_pth)
            print_hdr(tb)
    except Exception as e:
        print('ERROR: failed to read file %s' % f_pth, e)


def save_plot(plt, f_pth):
    """
        Save a plot to disk

        :param matplotlib.pyplot plt: matplotlib.pyplot object
        :param string f_pth: Path to SDHDF file to plot
        :return None
    """
    plot_prefix = os.path.basename(f_pth)
    plt.savefig(plot_prefix + '.png')

    return None


def setup_figure():
    """
        Setup the parameters of the plot figure, such as spacing,
        viewport size, number of rows/columns etc.

        :return None
    """
    fig_n_rows = 2
    fig_n_cols = 26
    fig_height = 1.5
    fig_width = 12
    fig_wspace = 0
    fig_hspace = 0

    fig, axs = plt.subplots(fig_n_rows, fig_n_cols)
    fig.set_size_inches(fig_width, fig_height)
    plt.subplots_adjust(wspace=fig_wspace, hspace=fig_hspace)

    return fig, axs


def verify_sdhdf(f):
    """
        Plot the spectra (uncalibrated flux vs frequency),
        and waterfall (time vs frequency) for each polarisation product,
        for a specified sub-band from an SDHDF format file,
        and zoom if specified by the user (Default is no zoom)

        :param string f: Name of SDHDF file to read
        :return None
    """
    h5 = h5py.File(f, 'r')
    bp = QTable.read(h5, path='/metadata/beam_params')

    # get primary header metadata
    read_sdhdf_header(f, '/metadata/primary_header')

    # create quick-look plots for web monitor
    # loop over the beams
    for beam in range(0, len(bp)):
        beam_label = 'beam_' + str(beam)
        sb_avail = QTable.read(h5, path=beam_label + '/metadata/band_params')

        # set up the figure
        fig, axs = setup_figure()

        # loop over the sub-bands
        for sb_id in range(0, 26):
            sb_label = 'band_SB' + str(sb_id)
            if sb_label in sb_avail['LABEL']:
                sb_data = beam_label + '/' + sb_label + '/astronomy_data/data'
                sb_freq = beam_label + '/' + sb_label + '/astronomy_data/frequency'
                op = QTable.read(h5, path=beam_label + '/' + sb_label + '/metadata/obs_params')

                # plot spectra
                # plot_spectra()
                lc = 'b'
                np.mean(h5[sb_data], axis=0)
                plt.yscale('log')
                axs[0, sb_id].plot(h5[sb_freq][:], h5[sb_data][0, 0, 0, :],
                            linewidth=0.5, color=lc)

                axs[0, sb_id].set_title(sb_id, fontsize=6)
                axs[0, sb_id].axis('off')

                # plot waterfall data (time/frequency)
                # plot_waterfall()
                if len(op['MJD']) > 1:
                    av = np.mean(h5[sb_data], axis=4)
                    axs[1, sb_id].imshow(av[:, 0, 0, :],
                       aspect='auto',
                       extent=(h5[sb_freq][0],
                               h5[sb_freq][-1],
                               op['MJD'][0],
                               op['MJD'][-1]),
                       interpolation='nearest')
                    axs[1, sb_id].axis('off')
            else:
                axs[0, sb_id].set_title(sb_id, fontsize=6)
                axs[0, sb_id].axis('off')
                axs[1, sb_id].axis('off')

        # write plot to disk
        save_plot(plt, f)

    # check ATOA ingest validation

    # other checks


if __name__ == '__main__':
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument('--filename', help='Path to SDHDF file to read',
                    required=True)
    args = ap.parse_args()

    verify_sdhdf(args.filename)
