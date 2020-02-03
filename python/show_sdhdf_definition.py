import h5py
import argparse
import csv
import os
from datetime import datetime
from tabulate import tabulate
from astropy.table import QTable

__version__ = '1.9'
__author__ = 'Lawrence Toomey'

hr = '#' * 79
dte = datetime.now().strftime('%Y')


def add_rows_to_csv(f_name, rows):
    f_csv = open(f_name, 'a')
    writer = csv.writer(f_csv)
    writer.writerows(rows)
    f_csv.close()

    return None


def add_line(f_name, line):
    f = open(f_name, 'a')
    f.write(line + '\n')

    return None


def read_csv(f_name):
    row_data = []
    f_csv = open(f_name, 'r')
    reader = csv.reader(f_csv, delimiter=',')
    for row in reader:
        row_data.append(row)
    print(tabulate(row_data))

    return None


def show_sdhdf_definition(f):
    defn_list = []
    try:
        with h5py.File(f, 'r') as h5:
            # Definition overview
            sdhdf_ver = bytes(h5['/metadata/primary_header']['HDR_DEFN_VERSION']).decode()
            f_name = os.path.basename(f)
            defn_csv = 'SDHDF_definition_' + sdhdf_ver + '.csv'

            add_line(defn_csv, hr)
            add_line(defn_csv, 'SDHDF Definition Overview')
            add_line(defn_csv, hr)
            add_line(defn_csv, 'SDHDF Definition Version: %s' % sdhdf_ver)
            add_line(defn_csv, 'Author: %s' % __author__)
            add_line(defn_csv, 'Copyright: CSIRO %s' % dte)
            add_line(defn_csv, hr)

            # File object meta-data overview
            add_line(defn_csv, '\n' + hr)
            add_line(defn_csv, 'SDHDF File Overview')
            add_line(defn_csv, 'HDF_Object_Name, HDF_Object_Type, Description')
            add_line(defn_csv, hr)
            defn_list.append([f_name, 'File', h5.attrs['DESCRIPTION']])
            add_rows_to_csv(defn_csv, defn_list)
            add_line(defn_csv, hr)

            # SDHDF structure overview
            add_line(defn_csv, '\n' + hr)
            add_line(defn_csv, 'SDHDF Structure Overview')
            add_line(defn_csv, 'HDF_Object_Name, HDF_Object_Type, Description')
            add_line(defn_csv, hr)

            for k in h5:
                defn_list = []
                if 'Group' in str(type(h5[k])):
                    defn_list.append(['/' + k, 'Group', h5[k].attrs['DESCRIPTION']])
                    for a in list(h5[k].keys()):
                        if 'Group' in str(type(h5[k][a])):
                            defn_list.append(['/' + k + '/' + a, 'Group',
                                              h5[k][a].attrs['DESCRIPTION']])
                            for b in list(h5[k][a].keys()):
                                if 'Group' in str(type(h5[k][a][b])):
                                    defn_list.append(['/' + k + '/' + a + '/' + b, 'Group',
                                                      h5[k][a][b].attrs['DESCRIPTION']])
                                if 'Dataset' in str(type(h5[k][a][b])):
                                    tb = QTable.read(h5, path='/' + k + '/' + a + '/' + b)
                                    defn_list.append(['/' + k + '/' + a + '/' + b, 'Dataset',
                                                      tb.meta['DESCRIPTION']])
                        elif 'Dataset' in str(type(h5[k][a])):
                            tb = QTable.read(h5, path='/' + k + '/' + a)
                            defn_list.append(['/' + k + '/' + a, 'Dataset',
                            tb.meta['DESCRIPTION']])

                    add_rows_to_csv(defn_csv, defn_list)
                    add_line(defn_csv, hr)

            defn_list = []
            for k in h5:
                if 'Dataset' in str(type(h5[k])):
                    tb = QTable.read(h5, path=k)
                    defn_list.append(['/' + k, 'Dataset', tb.meta['DESCRIPTION']])

            add_rows_to_csv(defn_csv, defn_list)

#            # SDHDF structure detail
#            add_line(defn_csv, '\n' + hr)
#            add_line(defn_csv, 'SDHDF Structure Detail')
#            add_line(defn_csv, 'HDF_Object_Name, HDF_Object_Type, Description')
#            add_line(defn_csv, hr)

#            for k in h5:
#                defn_list = []
#                if 'Group' in str(type(h5[k])):
#                    defn_list.append(['/' + k, 'Group', h5[k].attrs['DESCRIPTION']])
#                    for a in h5[k].attrs.keys():
#                        defn_list.append([a, 'Attribute', h5[k].attrs[a]])#
#
#                    for a in list(h5[k].keys()):
#                        defn_list.append(['/' + k + '/' + a, 'Dataset', h5[k][a].attrs['DESCRIPTION']])
#                        for attr in h5[k][a].attrs.keys():
#                            if 'KEY' in str(h5[k][a].attrs[attr]):
#                                attr_dict_val = h5[k][a].attrs[attr]
#                            else:
#                                if attr == 'CLASS':
#                                    attr_dict_val = bytes(h5[k][a].attrs[attr]).decode()
#                                elif attr == 'REFERENCE_LIST' or attr == 'DIMENSION_LIST':
#                                    attr_dict_val = attr
#                                else:
#                                    attr_dict_val = h5[k][a].attrs[attr]#
#
#                            defn_list.append([attr, 'Attribute', attr_dict_val])#
#
#                    add_rows_to_csv(defn_csv, defn_list)
#                    add_line(defn_csv, hr)

            add_line(defn_csv, '\n' + hr)
            add_line(defn_csv, 'File %s \nconforms to the SDHDF definition %s' % (f, sdhdf_ver))
            add_line(defn_csv, hr)
            add_line(defn_csv, 'Output written to %s' % defn_csv)
            add_line(defn_csv, hr)

            # read back csv file
            read_csv(defn_csv)

    except Exception as e:
        print('ERROR: failed to read file %s - '
              'contents does not match SDHDF definition.' % f, e)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--filename',
                    help='Path to SDHDF file to read',
                    required=True)
    args = ap.parse_args()

    show_sdhdf_definition(args.filename)
