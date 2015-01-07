#!/bin/python

import sys
import os
import subprocess


def validateFeatureParamFiles(fname):
    required = set("RNA")
    
    patient_map = {}
    feature_warnings = []
    
    f = open(fname)
    #header
    f.readline()
    errors = []
    for line in f:
        line = line.strip()
        items = line.split("\t")
        #print items
        
        patient_id = items[0]
        
        if patient_id not in patient_map:
            patient_map[patient_id] = []
        
        data_type = items[1]
        
        patient_map[patient_id].append(data_type)
        
        
        bw_fname = items[2]
        
        if not os.path.isfile(bw_file):
            errors.append("%s does not exist"%(bw_file))            
        elif not bw_fname[:-2] == "bw":
            errors.append("%s does not appear to be a bigwig file"%(bw_file))
    
    
    for p in patient_map:
        hasRna = False
        pkeys = patient_map[p]
        
        for val in pkeys:
            if val in required and val == "RNA":
                hasRna == True
            
            if not hasRna:
                errors.append("data set %s has no RNA"%(p))
            
            
            if len(pkeys) == 2 and hasRna:
                feature_warnings.append("warning: data set %s has only two features, three or more recomended."%(p))

        #print items
    #print fname
    return len(errors) == 0,errors,feature_warnings
    

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
    
    print "Feature Extraction"
    feature_param_fname = "/home/cemmeydan/inputTest3.txt"
    isFeaturesValid,errorMsgs,feature_warnings = validateFeatureParamFiles(feature_param_fname)
    
    if len(feature_warnings) > 0:
        for warn in feature_warnings:
            print warn
    
    if not isFeaturesValid:
        for val in errorMsgs:
            print errorMsgs
        print
        sys.exit(1)
        
    args = []
    #args = ["R","generateRegions.py",refgene_fname]

    #proc = subprocess.Popen(args,stdout=subprocess.PIPE)
    #out = proc.communicate()[0]

    print "training model"

    print "testing model"

    print "..."







if __name__ == "__main__":
    main()