#!/usr/bin/env python

import sys
import csv

WANTED = ['Estimated Number of Cells', 'Mean Reads per Cell', 'Median Genes per Cell', 'Fraction Reads in Cells', 'Number of Reads']

def getMetrics(filename):
    values = []
    with open(filename, "r") as f:
        c = csv.reader(f)
        hdr = c.__next__()
        data = c.__next__()
        for w in WANTED:
            p = hdr.index(w)
            values.append(data[p])
    for d in values:
        sys.stdout.write("<TD align='right'>{}</TD>".format(d))
    sys.stdout.write("\n")

def getAllMetrics(filenames):
    first = True
    for fn in filenames:
        smp = fn.split("/")[0] # we know first component of path is sample name
        with open(fn, "r") as f:
            hdr = f.readline()
            if first:
                sys.stdout.write("Sample," + hdr)
                first = False
            sys.stdout.write(smp + "," + f.readline())

if __name__ == "__main__":
    if sys.argv[1] == "-s":
        getAllMetrics(sys.argv[2:])
    else:
        getMetrics(sys.argv[1])
