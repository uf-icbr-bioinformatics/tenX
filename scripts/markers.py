#!/usr/bin/env python

import sys
import csv
import os.path
import xlsxwriter

HEADER = ["Gene", "avg_log2FC", "p_val_adj", "pct.1", "pct.2"]

def add_header(ws, bold):
    idx = 0
    for hdr in HEADER:
        ws.write(0, idx, hdr, bold)
        idx += 1

# ['TNNT3', '2.71669193463137e-242', '1.9504026200521', '0.99', '0.446', '4.76996769882575e-238', '0', 'TNNT3']

def store_row(ws, wr, row):
    ws.write(wr, 0, row[0])
    ws.write(wr, 1, row[2])
    ws.write(wr, 2, row[5])
    ws.write(wr, 3, row[3])
    ws.write(wr, 4, row[4])

def main(matrixfile, xlsxfile):
    if not os.path.isfile(matrixfile):
        sys.stderr.write("Error: input file {} not found.\n".format(matrixfile))
        sys.exit(1)

    sys.stderr.write("{} => {}\n".format(matrixfile, xlsxfile))
    workbook = xlsxwriter.Workbook(xlsxfile, {'strings_to_numbers': True})
    try:
        workbook.set_properties({'author': 'A. Riva, ariva@ufl.edu', 'company': 'DiBiG - ICBR Bioinformatics'}) # these should be read from conf or command-line
        bold = workbook.add_format({'bold': 1})

        curr_clust = ""
        with open(matrixfile, "r") as f:
            c = csv.reader(f, delimiter='\t')
            next(c)             # skip header
            for row in c:
                clust = row[6]
                if clust != curr_clust:
                    ws = workbook.add_worksheet("Cluster " + clust)
                    curr_clust = clust
                    add_header(ws, bold)
                    wr = 1
                store_row(ws, wr, row)
                wr += 1

    finally:
        workbook.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
