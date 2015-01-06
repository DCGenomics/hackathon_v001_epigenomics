#!/bin/python

import sys
import subprocess


def validateFeatureParamFiles(fname):
    print fname
    

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