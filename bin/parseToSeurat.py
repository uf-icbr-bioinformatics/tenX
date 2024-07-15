#!/usr/bin/env python

import sys
import csv

def main():
    for line in sys.stdin:
        if line[0] == '%':
            sys.stdout.write(line)
        else:
            row = line.split(" ")
            sys.stdout.write("{} {} {}".format(row[1], row[0], row[2]))

main()

        
