import h5py
import argparse
import csv
import os
from datetime import datetime
from tabulate import tabulate

__version__ = '1.9'
__author__ = 'Lawrence Toomey'

hr = '#' * 79
dte = datetime.now().strftime('%Y')


def add_rows_to_csv(f_name, rows):
    f_csv = open(f_name, 'a')
    writer = csv.writer(f_csv, lineterminator="\n")
    writer.writerows(rows)
    f_csv.close()

    return None


def add_line(f_name, line):
    f = open(f_name, 'a')
    f.write(line + "\n")

    return None


def read_csv(f_name):
    row_data = []
    f_csv = open(f_name, 'r')
    reader = csv.reader(f_csv, delimiter=',')
    for row in reader:
        row_data.append(row)
    print(tabulate(row_data))

    return None


def append_attributes(obj_attrs, defn_list):
    for attr in obj_attrs:
        defn_list.append([attr, 'Attribute', obj_attrs[attr]])

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

            # File object metadata overview
            add_line(defn_csv, "\n" + hr)
            add_line(defn_csv, 'SDHDF File Overview')
            add_line(defn_csv, 'HDF_Object_Name, HDF_Object_Type, Value')
            add_line(defn_csv, hr)
            defn_list.append([f_name, 'File', h5.attrs['DESCRIPTION']])
            add_rows_to_csv(defn_csv, defn_list)
            add_line(defn_csv, hr)

            # SDHDF structure overview
            add_line(defn_csv, "\n" + hr)
            add_line(defn_csv, 'SDHDF Structure Overview')
            add_line(defn_csv, 'HDF_Object_Name, HDF_Object_Type, Value')
            add_line(defn_csv, hr)

            # Loop over the groups and datasets and retrieve the attributes
            for k in h5:
                defn_list = []
                # root group
                if 'Group' in str(type(h5[k])):
                    pth_str = '/' + k
                    defn_list.append([pth_str, 'Group', ''])
                    append_attributes(h5[k].attrs, defn_list)
                    # beam, config, metadata groups
                    for a in list(h5[k].keys()):
                        if 'Group' in str(type(h5[k][a])):
                            pth_str = '/' + k + '/' + a
                            defn_list.append([pth_str, 'Group', ''])
                            append_attributes(h5[k][a].attrs, defn_list)
                            # band group
                            for b in list(h5[k][a].keys()):
                                if 'Group' in str(type(h5[k][a][b])):
                                    pth_str = '/' + k + '/' + a + '/' + b
                                    defn_list.append([pth_str, 'Group', ''])
                                    append_attributes(h5[k][a][b].attrs, defn_list)
                                    # data group
                                    for c in list(h5[k][a][b].keys()):
                                        if 'Dataset' in str(type(h5[k][a][b][c])):
                                            pth_str = '/' + k + '/' + a + '/' + b + '/' + c
                                            defn_list.append([pth_str, 'Dataset', ''])
                                            append_attributes(h5[k][a][b][c].attrs, defn_list)
                                elif 'Dataset' in str(type(h5[k][a][b])):
                                    pth_str = '/' + k + '/' + a + '/' + b
                                    defn_list.append([pth_str, 'Dataset', ''])
                                    append_attributes(h5[k][a][b].attrs, defn_list)
                        elif 'Dataset' in str(type(h5[k][a])):
                            pth_str = '/' + k + '/' + a
                            defn_list.append([pth_str, 'Dataset', ''])
                            append_attributes(h5[k][a].attrs, defn_list)

                    add_rows_to_csv(defn_csv, defn_list)
                    add_line(defn_csv, hr)

            add_rows_to_csv(defn_csv, defn_list)
            add_line(defn_csv, "\n" + hr)
            add_line(defn_csv, "File %s \nconforms to the SDHDF definition %s" % (f, sdhdf_ver))
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
