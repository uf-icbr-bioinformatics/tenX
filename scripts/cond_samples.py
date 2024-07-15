#!/usr/bin/env python

import sys
import csv

def main(samplesheet, condition):
    with open(samplesheet, "r") as f:
        c = csv.reader(f)
        hdr = next(c)
        samplecol = hdr.index("Sample")
        if "Condition" in hdr:
            condcol = hdr.index("Condition")
            samples = []
            for row in c:
                if row[condcol] == condition:
                    samples.append(row[samplecol])
            sys.stdout.write(",".join(samples) + "\n")
        else:
            for row in c:
                if row[samplecol] == condition:
                    sys.stdout.write(condition + "\n")

main(sys.argv[1], sys.argv[2])
