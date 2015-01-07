#!/bin/python

import sys
import os
import subprocess


def validateFeatureParamFiles(fname):
    f = open(fname)
    #header
    f.readline()
    errors = []
    for line in f:
        line = line.strip()
        items = line.split("\t")
        print items
        if not os.path.isfile(items[2]):
            errors.append("%s does not exist"%(items[2]))

        print items
    print fname
    return len(errors) == 0,errors
    

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
    
    
    

    feature_param_fname = "/home/cemmeydan/inputTest3.txt"
    print "Feature Extraction"
    isFeaturesValid,errorMsgs = validateFeatureParamFiles(feature_param_fname)
    
    if not isFeaturesValid:
        for val in errorMsgs:
            print errorMsgs
        print
        sys.exit(1)
    

    print "training model"

    print "testing model"

    print "..."







if __name__ == "__main__":
    main()