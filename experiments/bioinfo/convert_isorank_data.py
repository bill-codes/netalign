#!/usr/bin/env python

"""
Convert isorank data into smat data 

The isorank data from Bonnie Berger comes as a set of "tab" files
for the graphs and a "pair_match.txt" file for the bipartite graph.

## History
 2009-05-28: Initial coding
"""

__author__ = 'David F. Gleich'

import pprint

# TODO Implement better parsing
filesuffix = '_normalized_blast'
g1 = 'dmela'
g2 = 'scere'
l = g1+'-'+g2+filesuffix+'.evals'

# Outline
# 1.  parse pair information and build initial vertex hashes for g1, g2
# 2.  parse g1 and finish the vertex hash, and output edges
# 2.  parse g2 and finish the vertex hash, and output edges


def count_lines(filename):
    """ Count the lines in a file. """
    lineno = 0
    file = open(filename,'rt')
    for line in file:
        lineno += 1
    return lineno
    
def read_tab(filename,curmap=None):
    if curmap is None:
        nverts = 0
        name2id = {}
    else:
        nverts = len(curmap)
        name2id = curmap
    edges = []
    file = open(filename,'rt')
    file.readline() # skip one line
    for line in file:
        line = line.rstrip()
        parts = line.split()
        v1 = name2id.get(parts[0],None)
        if v1 is None:
            v1 = nverts
            name2id[parts[0]] = nverts
            nverts += 1
        v2 = name2id.get(parts[1],None)
        if v2 is None:
            v2 = nverts
            name2id[parts[1]] = nverts
            nverts += 1
        edges.append((v1,v2))
    file.close()
    return nverts, edges, name2id
    
def write_smat(filename, nverts, edges):
    file = open(filename, 'wt')
    file.write('%i %i %i\n'%(nverts,nverts,len(edges)))
    for e in edges:
        file.write('%i %i 1\n'%(e[0],e[1]))
    file.close()
    
def file_field_to_map(filename,curmap=None,field=0,skiplines=0):
    """ Build a map from a particular field in a file. """
    file = open(filename, 'rt')
    if skiplines>0:
        for i in xrange(skiplines):
            file.readline()
    nelem = 0     
    if curmap is None:
        map = {}
    else:
        map = curmap
        nelem = max(map.values())
    for line in file:
        line = line.rstrip()
        parts = line.split()
        if parts[field] not in map:
            map[parts[field]] = nelem
            nelem += 1
    return map
            
def main():
    g1filename = g1+filesuffix+'.tab'
    g2filename = g2+filesuffix+'.tab'
    lfilename = l
    # build maps for g1 and g2, this is stupidly inefficient, sigh.
    g1map = file_field_to_map(g1filename, field=0, skiplines=1)
    #pprint.pprint(g1map)
    g1map = file_field_to_map(g1filename, field=1, skiplines=1, curmap=g1map)
    g2map = file_field_to_map(g2filename, field=0, skiplines=1)
    g2map = file_field_to_map(g2filename, field=1, skiplines=1, curmap=g2map)
    
    g1map = file_field_to_map(lfilename, field=0, curmap=g1map)
    g2map = file_field_to_map(lfilename, field=1, curmap=g2map)
    
    #print g1map['dm4971']
    
    g1n, g1e, g1map = read_tab(g1+filesuffix+'.tab',g1map)
    g2n, g2e, g2map = read_tab(g2+filesuffix+'.tab',g2map)
    write_smat(g1+'.smat', g1n, g1e)
    write_smat(g2+'.smat', g2n, g2e)

    nmatches = count_lines(l)
    matchfile = open(l,'rt')
    matchsmatfile = open(g1+'-'+g2+'.smat','wt')
    matchsmatfile.write('%i %i %i\n'%(g1n,g2n,nmatches))
    for line in matchfile:
        line = line.rstrip()
        parts = line.split()
        v1 = g1map[parts[0]]
        v2 = g2map[parts[1]]
        matchsmatfile.write('%i %i %s\n'%(v1,v2,parts[2]))
    matchsmatfile.close()
            
if __name__=='__main__':
    main()

