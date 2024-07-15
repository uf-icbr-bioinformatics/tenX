#!/usr/bin/env python

import sys
import csv

class F(object):
    name = ""
    fmt = "r"

    def __init__(self, spec):
        if ":" in spec:
            parts = spec.split(":")
            self.name = parts[0]
            self.fmt = parts[1]
        else:
            self.name = spec

    def write(self, value):
        if self.fmt == "r":
            return value
        elif self.fmt == "d":
            return "{:,}".format(int(value))
        elif self.fmt == "f":
            return "{:.1f}".format(float(value))
        elif self.fmt == "p":
            return "{:.1f}%".format(100.0*float(value))

WANTED = ['Estimated Number of Cells', 'Mean Reads per Cell', 'Median Genes per Cell', 'Fraction Reads in Cells', 'Number of Reads']
WANTED_ATAC = ['cells_detected:d', 'median_fragments_per_cell:f', 'frac_fragments_overlapping_targets:p', 'frac_cut_fragments_in_peaks:p', 'num_fragments:d']
WANTED_VISIUM = ['Number of Spots Under Tissue:d', 'Mean Reads per Spot:f', 'Median Genes per Spot:f', 'Fraction of Spots Under Tissue:p', 'Fraction Reads in Spots Under Tissue:p', 'Number of Reads:d']
WANTED_HD = ['Reads Mapped to Probe Set:f', 'Fraction Reads in Squares Under Tissue:p', 'Genes Detected:d', 'Number of Reads:d']

class Metrics(object):
    wanted = []
    mode = "r"                  # r = scRNA-Seq, a = scATAC-Seq, v = Visium, d = Visium HD
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
                self.wanted = WANTED_VISIUM
            elif a == "-d":
                self.mode = "d"
                self.wanted = WANTED_HD
            else:
                self.filenames.append(a)
        self.wanted = [ F(spec) for spec in self.wanted ]

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
                p = hdr.index(w.name)
                values.append(w.write(data[p]))
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
