#!/usr/bin/env python

import sys
import csv

WANTED = ['Estimated Number of Cells', 'Mean Reads per Cell', 'Median Genes per Cell', 'Fraction Reads in Cells', 'Number of Reads']
WANTED_ATAC = ['cells_detected', 'median_fragments_per_cell', 'frac_fragments_overlapping_targets', 'frac_cut_fragments_in_peaks', 'num_fragments']

class Metrics(object):
    wanted = []
    mode = "r"                  # r = scRNA-Seq, a = scATAC-Seq, s = Visium
    filenames = []
    what = "m"                  # a = all metrics

    def __init__(self, args):
        self.wanted = WANTED
        self.filenames = []

        for a in args:
            if a == "-s":
                self.what = "a"
            elif a == "-a":
                self.mode = "a"
                self.wanted = WANTED_ATAC
            elif a == "-v":
                self.mode = "v"
            else:
                self.filenames.append(a)

    def run(self):
        if self.what == "a":
            self.getAllMetrics()
        else:
            self.getMetrics()

    def getMetrics(self):
        values = []
        with open(self.filenames[0], "r") as f:
            c = csv.reader(f)
            hdr = c.__next__()
            data = c.__next__()
            for w in self.wanted:
                p = hdr.index(w)
                values.append(data[p])
        for d in values:
            sys.stdout.write("<TD align='right'>{}</TD>".format(d))
        sys.stdout.write("\n")

    def getAllMetrics(self):
        first = True
        for fn in self.filenames:
            smp = fn.split("/")[0] # we know first component of path is sample name
            with open(fn, "r") as f:
                hdr = f.readline()
                if first:
                    sys.stdout.write("Sample," + hdr)
                    first = False
                sys.stdout.write(smp + "," + f.readline())

if __name__ == "__main__":
    M = Metrics(sys.argv[1:])
    M.run()
