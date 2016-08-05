#!/usr/bin/env python
# we will get bed file of positions and report the homologs
# for all the species from maf file

import argparse
from collections import defaultdict
from maf_extractor import intersect
from synteny_blocks.model import MAF_Entry
from synteny_blocks.model import parse_bed
import math

def process(bed_file, maf_file):
    regions = parse_bed(bed_file)
    e_homologs = defaultdict(list)
    species = set()
    with open(maf_file) as maf:
        maf_entries = []
        for line in maf:
            line = line.strip()
            if not line or '#' in line:
                continue
            if line[0] == 'a':
                if maf_entries:
                    intersected_regions = intersect(maf_entries, regions)
                    if intersected_regions:
                        for e in intersected_regions:
                            e_homologs[e].append(maf_entries)
                    maf_entries = []
            else:
                line = line.split()
                line = line[1:]
                genome = line[0].split('.')[0]
                species.add(genome)
                chrom = '.'.join(line[0].split('.')[1:])
                maf_entries.append(MAF_Entry(genome, chrom, int(line[1]), int(line[2]),\
                                             line[3], int(line[4]), line[5]))
    return e_homologs, list(species)



def print_out_matrix(e_homologs, species):
    header = '\t'.join(['regions'] + species)
    print header
    for e in e_homologs:
        l = e.to_string() + '\t'
        maf_entries = e_homologs[e]
        for s in species:
            aligned_length = 0
            for maf in maf_entries:
                s_entries = filter(lambda x: x.genome == s, maf)
                for entry in s_entries:
                    aligned_length += entry.length
            #tunable:
            if aligned_length > math.fabs(e.end - e.start)/2:
                l += '+\t'
            else:
                l += '-\t'
        print l


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('bed', help='bed file')
    parser.add_argument('in_maf', help='input maf file')

    args = parser.parse_args()
    homologs, species = process(args.bed, args.in_maf)
    print_out_matrix(homologs, species)
