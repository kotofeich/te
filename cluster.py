#!/usr/bin/env python

from argparse import ArgumentParser
import os
from collections import defaultdict,Counter
from time import gmtime, strftime

MIN_SIZE_OF_COMMON_FEATURES = 1
SPECIES_NAMES = ['CanFam','CroCro','PanLeo','PanPar','PanOnc','PanTig','CarCar','LynPar','AciJub','PriBen','PriViv','FelCat']

#path - path to the clusters for the given feature
#candidates_info = defaultdict
def parse_candidate_ids2clusters(path, candidates_info, feature_name):
    with open(path) as f:
        cluster_id = '-2'
        for line in f:
            line = line.strip()
            if '#' == line[0]:
                continue
            if line[-1] == ':':
                cluster_id = line[:-1]
            else:
                candidate_id = line
                candidates_info[candidate_id].add((feature_name, cluster_id))
    return candidates_info
        
# parse file in the input folder
# and create list of the features 
# with dictionary of candidate_ids to cluster_id
def build_candidate_ids2features(folder):
    files = os.listdir(folder)
    #clustered_features = {}
    candidates_info = defaultdict(set)
    for f in files:
        feature_name = f.split('.')[0]
        path = os.path.join(folder,f)
        candidates_info = parse_candidate_ids2clusters(path, candidates_info, feature_name)
    #for k in candidates_info.keys():
    #    print k, candidates_info[k]
    return candidates_info


def are_compatible(cluster_set1, cluster_set2):
    names1 = map(lambda x: x[0], cluster_set1)
    names2 = map(lambda x: x[0], cluster_set2)
    common_names = set(names1) & set(names2)
    if len(common_names) < MIN_SIZE_OF_COMMON_FEATURES:
        return False
    #print 'common_names:', len(common_names)
    for name in common_names:
        cluster1 = filter(lambda x: x[0] == name, cluster_set1)
        cluster2 = filter(lambda x: x[0] == name, cluster_set2)
        if cluster1[0][1] != cluster2[0][1]:
            return False
    return True

#for each id find list of compatible ids
def check_compatibility_for_all_pairs(candidates_info):
    compatibility_dict = defaultdict(set)
    cand_ids = candidates_info.keys()    
    for i in cand_ids:
        compatibility_dict[i]
        for j in cand_ids:
            if i == j:
                continue
            if are_compatible(candidates_info[i], candidates_info[j]):
                compatibility_dict[i].add(j)
                compatibility_dict[j].add(i)
    return compatibility_dict

def do_not_overlap(cluster_set1, cluster_set2):
    names1 = map(lambda x: x[0], cluster_set1)
    names2 = map(lambda x: x[0], cluster_set2)
    common_names = set(names1) & set(names2)
    if len(common_names) == 0:
        return True
    return False

def check_candidates_with_nonoverlapping_features(candidates_info):
    cand_ids = candidates_info.keys()
    nonoverlapping = defaultdict(set)
    for i in cand_ids:
        compatibility_dict[i]
        for j in cand_ids:
            if i == j:
                continue
            if do_not_overlap(candidates_info[i], candidates_info[j]):
                nonoverlapping[i].add(j)
                nonoverlapping[j].add(i)
    return nonoverlapping

def check_connections(v, c, compatibility_dict):
    adj_vs = compatibility_dict[v]
    for v_compatible_set in c:
        if not v_compatible_set in adj_vs:
            return False
    return True

def check_connections_keeping_nonoverlapping_features(v, c, compatibility_dict, nonoverlapping_dict):
    adj_vs = compatibility_dict[v]
    nonoverlapping_vertices = nonoverlapping_dict[v]
    compatible = 0
    for v_compatible_set in c:
        if v_compatible_set in nonoverlapping_vertices:
            continue
        if not v_compatible_set in adj_vs:
            return False
        else:
            compatible += 1
    if compatible > 0:
        return True
    return False

def get_full_adjacent_subgraphs_from_vertex(start_vertex, compatibility_dict, nonoverlaping_features_dict):
    compatible_sets_for_vertex = []
    v_compatibility_list = compatibility_dict[start_vertex]
    if len(v_compatibility_list) == 0:
        return [[start_vertex]]
    #for each vertex adjacent to the start_vertex...
    for v in v_compatibility_list:
        added = False
        for c in compatible_sets_for_vertex:
            #..check if v is adjacent to all vertices in the accumulated compatible_set c
            # in order to comply the paper we will also allow the case when once vertex is not connected to
            # some others because it lacks the corresponding features. it still goes into the compatible set
            # if it's connected at least to one vertex in that set
            ##if check_connections(v, c, compatibility_dict):
            if check_connections_keeping_nonoverlapping_features(v, c, compatibility_dict, nonoverlapping_features_dict):
                #if so append v to compatible_set c
                c.append(v)
                added = True
        if not added:
            compatible_sets_for_vertex.append([start_vertex, v])
    return map(sorted,compatible_sets_for_vertex) 

#we don't call them clicks because our compatible sets are not necessarily disjoint
def find_compatible_sets(compatibility_dict, nonoverlaping_features_dict):
    vertices = compatibility_dict.keys()
    compatible_set_candidates = []
    print strftime("%Y-%m-%d %H:%M:%S", gmtime())
    print '    finding full adjacent subgraphs for each vertex...'
    for v in vertices:
        print 'start_vertex:', v
        v_candidates = get_full_adjacent_subgraphs_from_vertex(v, compatibility_dict, nonoverlaping_features_dict)
        for c in v_candidates:
            if not c in compatible_set_candidates:
                compatible_set_candidates.append(c)
    #print 'unsorted candidates:', compatible_set_candidates
    print strftime("%Y-%m-%d %H:%M:%S", gmtime())
    print '    sorting candidates...'
    sorted_candidates = sorted(compatible_set_candidates, key=len, reverse=True)
    used = set()
    print strftime("%Y-%m-%d %H:%M:%S", gmtime())
    #print 'found clusters:'
    result = []
    for c in sorted_candidates:
        intersected = set(c) & used
        if len(intersected) < len(c):
            for e in c:
                print e
            print 
            result.append(c)
            used.update(set(c))
    return result

def filter_ambigious(compatible_sets):
    ambigious = set()
    for c in compatible_sets:
        for c_other in compatible_sets:
            if c == c_other:
                continue
            overlap = set(c) & set(c_other)
            if overlap:
                ambigious |= overlap
    unambigious = []
    for c in compatible_sets:
        upd_c = set()
        for e in c:
            if not e in ambigious:
                upd_c.add(e)
        unambigious.append(sorted(upd_c))
    return unambigious, ambigious

def count_species(s, names):
    counts = []
    for n in names:
        cnt = 0
        for x in s:
            if n in x:
                cnt += 1
        counts.append(cnt)
    return counts

#check if cluster contains only one species
def find_uniq_clusters(clusters, names):
    #results = defaultdict(set)
    for c in clusters:
        counts = count_species(c, names)
        notnull = filter(lambda x: x!=0, counts)
        if len(notnull) == 1:
            print names[counts.index(notnull[0])], ':', c

def find_occurences(c, ordered_names):
    oc = [0] * len(ordered_names)
    i = 0
    for name in ordered_names:
        for x in c:
            if name in x:
                oc[i] = 1
        i += 1
    return oc

def count_occurences(c, ordered_names):
    oc = [0] * len(ordered_names)
    i = 0
    for name in ordered_names:
        for x in c:
            if name in x:
                oc[i] += 1
        i += 1
    return oc

def process_proportions(c, ordered_names):
    oc = [0] * len(ordered_names)
    i = 0
    #print 'ordered_names:', ordered_names
    for name in ordered_names:
        for x in c:
            if name in x:
                oc[i] += 1
        i += 1
    oc = map(float,oc)
    value = sum(oc)
    #print 'oc float', oc
    #print 'value',value
    oc = map(lambda x:x/value, oc)
    return oc


#if there is 0 between 1s
def not_regular(ar):
    if ar[0] == 0 and ar[1] == 1:
        return False
    if ar[-1] == 0 and ar[-2] == 1:
        return False
    for i in range(1, len(ar)-2):
        if ar[i] == 1 and ar[i+1] == 0 and ar[i+2] == 1:
            return True
    return False
        

def filter_suspicious_clusters(ordered_names, clusters):
    for c in clusters:
        oc = find_occurences(c, ordered_names)
        if not_regular(oc):
            print 'suspicious:', oc
            for e in c:
                print e

def process_compatible_sets_as_you_like(sorted_compatible_sets):
    for c_1 in sorted_compatible_sets:
        for c_2 in sorted_compatible_sets:
            pass

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('clustered_features',help='folder to clustered features')
    args = parser.parse_args()
    print strftime("%Y-%m-%d %H:%M:%S", gmtime())
    print 'building index of candidates to features...'
    candidates_info = build_candidate_ids2features(args.clustered_features)
    #for e in candidates_info:
    #    print 'number_of_features:', len(candidates_info[e])
    print strftime("%Y-%m-%d %H:%M:%S", gmtime())
    print 'building compatibility relation...'
    compatibility_dict = check_compatibility_for_all_pairs(candidates_info)
    nonoverlapping_features_dict = check_candidates_with_nonoverlapping_features(candidates_info)
    #for e in compatibility_dict.keys():
    #    print e, ':', compatibility_dict[e]
    print strftime("%Y-%m-%d %H:%M:%S", gmtime())
    print 'finding compatible sets...'
    compatible_sets = find_compatible_sets(compatibility_dict, nonoverlapping_features_dict)
    print 'number of compatible sets', len(compatible_sets)
    print strftime("%Y-%m-%d %H:%M:%S", gmtime())
    print 'filtering ambigious sets...'
    unambigious_compatile_sets, ambigious_vertices = filter_ambigious(compatible_sets)
    print strftime("%Y-%m-%d %H:%M:%S", gmtime())
    print 'unambigious sets:'
    for c in unambigious_compatile_sets:
        for e in c:
            print e,
        print
    print 'ambigious:'
    for v in ambigious_vertices:
        print v
    #filter_suspicious_clusters(SPECIES_NAMES, compatible_sets) 
    #species_counts = defaultdict(list) 
    #for c in compatible_sets:
    #    counts = count_species(c, SPECIES_NAMES)
    #    for i in range(len(counts)):
    #        species_counts[SPECIES_NAMES[i]].append(counts[i])

    #for s in species_counts:
    #    print s, ':', Counter(species_counts[s])

    #find_uniq_clusters(compatible_sets, SPECIES_NAMES)

    #for c in compatible_sets:
    #    oc = count_occurences(c, SPECIES_NAMES)
    #    if sum(oc) == 0:
    #        raise Exception("Cluster is empty!")
    #    z = zip(SPECIES_NAMES, oc)
    #    for e in sorted(z, key = lambda x:x[0]):
    #        print e[0],':',e[1], 
    #    print

    
    #for c in compatible_sets:
    #    prop = process_proportions(c, SPECIES_NAMES)
    #    z = zip(SPECIES_NAMES, prop)
    #    for e in sorted(z, key = lambda x:x[0]):
    #        print e[0],':', '%1.3f'%e[1], 
    #    print
