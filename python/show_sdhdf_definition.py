#!/usr/bin/env python
import h5py
import argparse
import csv
import os
from datetime import datetime
from tabulate import tabulate

__version__ = '2.2'
__author__ = 'Lawrence Toomey'

hr = '-' * 79
dte = datetime.now().strftime('%Y')


def add_rows_to_csv(f_name, rows):
    """
    Bulk write text to a CSV file

    :param string f_name: Path to output CSV file
    :param list rows: List of rows to write to CSV
    :return: None
    """
    f_csv = open(f_name, 'a')
    writer = csv.writer(f_csv, lineterminator="\n")
    writer.writerows(rows)
    f_csv.close()

    return None


def add_line(f_name, line):
    """
    Append to a CSV file

    :param string f_name: Path to output CSV file
    :param string line: Line to append to CSV
    :return: None
    """
    f = open(f_name, 'a')
    f.write(line + "\n")

    return None


def read_csv(f_name):
    """
    Read a CSV file and print to stdout

    :param string f_name: Path to CSV file to read
    :return: None
    """
    row_data = []
    f_csv = open(f_name, 'r')
    reader = csv.reader(f_csv, delimiter=',')
    for row in reader:
        row_data.append(row)
    print(tabulate(row_data))

    return None


def append_attributes(obj_attrs, defn_list):
    """
    Sanitise the attribute if required and append to a list

    :param list obj_attrs: List of attributes to append
    :param defn_list: List to append to
    :return: None
    """
    for attr in obj_attrs:
        if attr in ('DIMENSION_LIST', 'REFERENCE_LIST'):
            print('Ignoring attribute %s' % attr)
        elif attr == 'DIMENSION_LABELS':
            vl = []
            for v in obj_attrs[attr]:
                vl.append(v.decode("utf-8"))
            defn_list.append([attr, 'Attribute', vl])
        elif attr == 'CLASS':
            vl = obj_attrs[attr].decode("utf-8")
            defn_list.append([attr, 'Attribute', vl])
        else:
            defn_list.append([attr, 'Attribute', obj_attrs[attr]])

    return None


def show_sdhdf_definition(f, sb, output):
    """
    Display the formal SDHDF definition contained in an SDHDF file

    :param string f: Path to SDHDF file to read
    :param string sb: Sub-band to read
    :param bool output: Print output to stdout [True|False]
    :return: None
    """
    defn_list = []
    try:
        with h5py.File(f, 'r') as h5:
            # Definition overview
            sdhdf_ver = bytes(h5['/metadata/primary_header']['HDR_DEFN_VERSION']).decode()
            f_name = os.path.basename(f)
            defn_csv = 'SDHDF_definition_' + sdhdf_ver + '.csv'

            add_line(defn_csv, 'SDHDF Definition Overview,-,-')
            add_line(defn_csv, 'SDHDF Definition Version,%s,-' % sdhdf_ver)
            add_line(defn_csv, 'Author,%s,-' % __author__)
            add_line(defn_csv, 'Copyright,CSIRO %s,-' % dte)

            # File object metadata overview
            add_line(defn_csv, "-, -, -")
            add_line(defn_csv, 'SDHDF File Overview,-,-')
            add_line(defn_csv, 'HDF_Object_Name,HDF_Object_Type,Value')
            append_attributes(h5.attrs, defn_list)
            defn_list.append(['-', '-', '-'])
            add_rows_to_csv(defn_csv, defn_list)

            # SDHDF structure overview
            add_line(defn_csv, "-, -, -")
            add_line(defn_csv, 'SDHDF Structure Overview,-,-')
            add_line(defn_csv, 'HDF_Object_Name,HDF_Object_Type,Value')
            add_line(defn_csv, "-, -, -")

            # check existence of requested sub-band
            root_keys = list(h5.keys())
            beam = root_keys[0]
            if sb != 'all':
                if sb not in list(h5[beam].keys()):
                    raise Exception('ERROR: no such sub-band %s' % sb)

            # Loop over the groups and datasets and retrieve the attributes
            for k in h5:
                print('Processing root group %s' % k)
                defn_list = []
                # root group
                if 'Group' in str(type(h5[k])):
                    pth_str = '/' + k
                    defn_list.append([pth_str, 'Group', '-'])
                    append_attributes(h5[k].attrs, defn_list)
                    defn_list.append(['-', '-', '-'])
                    # beam, config, metadata groups
                    for a in list(h5[k].keys()):
                        if a == sb or 'SB' not in a or sb != 'all':
                            print('-> Processing group %s' % a)
                            if 'Group' in str(type(h5[k][a])):
                                pth_str = '/' + k + '/' + a
                                defn_list.append([pth_str, 'Group', '-'])
                                append_attributes(h5[k][a].attrs, defn_list)
                                defn_list.append(['-', '-', '-'])
                                # band group
                                for b in list(h5[k][a].keys()):
                                    if 'Group' in str(type(h5[k][a][b])):
                                        pth_str = '/' + k + '/' + a + '/' + b
                                        defn_list.append([pth_str, 'Group', '-'])
                                        append_attributes(h5[k][a][b].attrs, defn_list)
                                        defn_list.append(['-', '-', '-'])
                                        # data group
                                        for c in list(h5[k][a][b].keys()):
                                            if 'Dataset' in str(type(h5[k][a][b][c])):
                                                pth_str = '/' + k + '/' + a + '/' + b + '/' + c
                                                defn_list.append([pth_str, 'Dataset', '-'])
                                                append_attributes(h5[k][a][b][c].attrs, defn_list)
                                                defn_list.append(['-', '-', '-'])
                                    elif 'Dataset' in str(type(h5[k][a][b])):
                                        pth_str = '/' + k + '/' + a + '/' + b
                                        defn_list.append([pth_str, 'Dataset', '-'])
                                        append_attributes(h5[k][a][b].attrs, defn_list)
                                        defn_list.append(['-', '-', '-'])
                            elif 'Dataset' in str(type(h5[k][a])):
                                pth_str = '/' + k + '/' + a
                                defn_list.append([pth_str, 'Dataset', '-'])
                                append_attributes(h5[k][a].attrs, defn_list)
                                defn_list.append(['-', '-', '-'])
                        else:
                            if 'SB' in a:
                                print('Skipping group %s' % a)

                    add_rows_to_csv(defn_csv, defn_list)

            add_line(defn_csv, "\n" + hr)
            add_line(defn_csv, "PASS: File %s conforms to SDHDF definition v%s" % (f_name, sdhdf_ver))
            add_line(defn_csv, 'Output written to %s' % defn_csv)

            # read back csv file
            if bool(output) is True:
                read_csv(defn_csv)
            else:
                print("PASS: File %s conforms to SDHDF definition v%s" % (f_name, sdhdf_ver))

    except Exception as e:
        print('ERROR: failed to read %s - '
              'file does not conform to the SDHDF definition' % f, e)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--filename',
                    help='Path to SDHDF file to read',
                    required=True)
    ap.add_argument('--subband',
                    help='Sub-band to read [Default=all]',
                    required=False,
                    default='all')
    ap.add_argument('--output',
                    help='Print output to stdout [True|False]',
                    required=True,
                    default=False)
    args = ap.parse_args()

    show_sdhdf_definition(args.filename, args.subband, args.output)
