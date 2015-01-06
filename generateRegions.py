#!/bin/python

import sys

def parseRefGene(fname):
    f = open(fname)
    for line in f:
        line = line.strip()
        items = line.split("\t")
        yield items
    f.close()
    
def gene_id(rec):
    return rec[1]

def chrom(rec):
    return rec[2]


def start(rec):
    return int(rec[4])


def end(rec):
    return int(rec[5])


def strand(rec):
    return rec[3]


def upstreamRegion(rec,distance):
    outrec = []
    if strand(rec) == "+":
        s = start(rec) - distance
        outrec.append(s)
        outrec.append(start(rec))
    else:
        s = end(rec)
        outrec.append(s)
        outrec.append(s+distance)
    return outrec

def downstreamRegion(rec,distance):
    outrec = []
    if strand(rec) == "-":
        s = start(rec) - distance
        outrec.append(s)
        outrec.append(start(rec))
    else:
        s = end(rec)
        outrec.append(s)
        outrec.append(s+distance)
    return outrec

def genOutRecs(rec,regions):
    new_recs = []
    for akey in regions:
        s,e = regions[akey]
        new_rec = []
        new_rec.append(gene_id(rec))
        new_rec.append(akey)
        new_rec.append(chrom(rec))
        new_rec.append(str(s))
        new_rec.append(str(e))
        new_recs.append(new_rec)
    return new_recs

def transformGene(rec,upstream,downstream):
    upReg = upstreamRegion(rec,upstream)
    downReg = downstreamRegion(rec,downstream)
    geneBody = [start(rec),end(rec)]
    
    #out_recs = genOutRecs(rec, [upReg,geneBody,downReg])
    out_map = {}
    out_map["up"] = upReg
    out_map["body"] = geneBody
    out_map["down"] = downReg
    #out_recs = genOutRecs(rec, [upReg,geneBody,downReg])
    out_recs = genOutRecs(rec, out_map)
            
    return out_recs
        
    

def main():
    refgene_fname = sys.argv[1]
    
    up   = 5000
    down = 5000
    
    count = 0
    for rec in parseRefGene(refgene_fname):
        regions = transformGene(rec,up,down)
        
        for vals in regions:
            outs = []
            outs.append(str(count))
            outs.extend(vals)
            
            count += 1
            print "\t".join(outs)
            

if __name__ == "__main__":
    main()