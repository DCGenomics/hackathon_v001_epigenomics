#!/bin/python

import sys
import subprocess


def main():
    print "Starting Pipeline"
    print
    print "Extracting gene regions"

    refgene_fname = "refGene.txt"
    output_fname = "regions.tab"
    args = ["python","generateRegions.py",refgene_fname]

    proc = subprocess.Popen(args,stdout=subprocess.PIPE)

    out = proc.communicate()[0]

    fout = open(output_fname,"w")
    for line in out:
        fout.write(line)
    fout.close()
    print "Done extracting regions"

    print "Extracting features in regions"

    print "training model"

    print "testing model"

    print "..."







if __name__ == "__main__":
    main()