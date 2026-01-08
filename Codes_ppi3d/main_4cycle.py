from data_processing import *
from stats_graph import *

import random
import math
import numpy as np
import time

########## BASIC (OUTDATED) VERSION, ONLY COUNTING 4-CYCLES

# my_file_path = "../Networks/ppi3d_c95_human_net.txt"
# #my_file_path = "../Networks/ppi3d_c95_human_bsi_net.txt"
# #my_file_path = "../Networks/toy.txt"

# my_ppi_el = read_data(my_file_path)
# my_ppi_ael = create_adjacency_list(my_ppi_el)

# def identify_l3_structures(ael):
#     path3  = 0
#     cycle4 = 0 
#     for x in ael:
#         for u in ael[x]:
#             if x < u:
#                 for v in ael[u]:
#                     if x < v and v != u:
#                         for y in ael[v]:
#                             if x < y and u < y and y != v:
#                                 path3 += 1
#                                 if x in ael[y]:
#                                     cycle4 += 1
#                                     print(x,u,v,y)
#     return path3 , cycle4

# print(identify_l3_structures(my_ppi_ael))
# print(avg_degree(dd(my_ppi_ael)))
# print(density(my_ppi_ael))
# print(num_edges(my_ppi_ael))


########## ALTERNATIVE (OUTDATED) VERSION, COUNTING 4-CYCLES AND REGISTERING THEIR PDBID

# my_file_path = "../Networks/ppi3d_c95_human_net_pdbid.txt"

# my_ppi_el = read_data(my_file_path)
# my_ppi_el_pdbid = read_data_pdbid(my_file_path)
# my_ppi_ael = create_adjacency_list(my_ppi_el)

# def identify_l3_structures(ael):
#     path3  = 0
#     cycle4 = 0
#     store_cycles = set()
#     for x in ael:
#         for u in ael[x]:
#             if x < u:
#                 for v in ael[u]:
#                     if x < v and v != u:
#                         for y in ael[v]:
#                             if x < y and u < y and y != v:
#                                 path3 += 1
#                                 if x in ael[y]:
#                                     cycle4 += 1
#                                     store_cycles.add((x,u,v,y))
#     return path3 , cycle4 , store_cycles 

# my_path3 , my_cycle4 , my_store_cycles = identify_l3_structures(my_ppi_ael)

# find the pdbids of the 4 edges in the cycle

def pdb_ids(dict_edges_pdbid , store_cycles):
    store_cycles_pdbids = dict()
    for (x,u,v,y) in store_cycles:
        store_cycles_pdbids[(x,u,v,y)] = (dict_edges_pdbid[x,u],dict_edges_pdbid[u,v],dict_edges_pdbid[v,y],dict_edges_pdbid[y,x])
    return store_cycles_pdbids

# analyze the types of cycles according to their pdbid diversity

def analyze_pdb_ids(store_cycles_pdbids):
    cycle_types=[0,0,0,0]
    for cycle in store_cycles_pdbids:
        a, b, c, d = store_cycles_pdbids[cycle]
        print(a,b,c,d)
        cycle_types[len({a,b,c,d})-1] +=1
    return cycle_types

#my_store_cycles_pdbids = pdb_ids(my_ppi_el_pdbid , my_store_cycles)
#my_cycle_types = analyze_pdb_ids(my_store_cycles_pdbids)
#print(my_cycle_types)


########## ALTERNATIVE VERSION, COUNTING 4-CYCLES AND REGISTERING THEIR PDBID, AND MINIMIZING PDBID DIVERSITY OVER THE CYCLES

# loading data

my_file_path = "../Networks/ppi3d_c95_human_net_pdbid.txt"
my_dico_path = "../Networks/ppi3d_c95_human_net_pdbid_dico.txt"

my_ppi_el = read_data(my_file_path)
my_dico_all_pdbid = read_data_all_pdbid(my_dico_path)

my_ppi_ael = create_adjacency_list(my_ppi_el)

# counting 3-edges paths, 4-edges cycles and storing 4-cycles in a set of quadruplets

def identify_l3_structures(ael):
    path3  = 0
    cycle4 = 0
    store_cycles = set()
    for x in ael:
        for u in ael[x]:
            if x < u:
                for v in ael[u]:
                    if x < v and v != u:
                        for y in ael[v]:
                            if x < y and u < y and y != v:
                                path3 += 1
                                if x in ael[y]:
                                    cycle4 += 1
                                    store_cycles.add((x,u,v,y))
    return path3 , cycle4 , store_cycles 

# function to test all pdbids combinations to minimize their diversity in one cycle

def test_all (L1, L2, L3, L4):
    minsize_tuple = (L1[0], L2[0], L3[0], L4[0])
    minsize = 4
    for elt1 in L1:
        for elt2 in L2:
            for elt3 in L3:
                for elt4 in L4:
                    current_tuple = (elt1, elt2, elt3, elt4)
                    if len({elt1, elt2, elt3, elt4}) < minsize:
                        minsize_tuple = current_tuple
                        minsize = len({elt1, elt2, elt3, elt4})
    return minsize_tuple

# function to apply the previous function to all cycles in the storage

def pdb_ids(D_all , store_cycles):
    store_cycles_pdbids = dict()
    for (x,u,v,y) in store_cycles:
        min_xu , max_xu = min(x,u) , max(x,u)
        min_uv , max_uv = min(u,v) , max(u,v)
        min_vy , max_vy = min(v,y) , max(v,y)
        min_yx , max_yx = min(y,x) , max(y,x)
        store_cycles_pdbids[(x,u,v,y)] = test_all (D_all[(min_xu , max_xu)],\
                                                   D_all[(min_uv , max_uv)],\
                                                   D_all[(min_vy , max_vy)],\
                                                   D_all[(min_yx , max_yx)])
    return store_cycles_pdbids

# analyze the diversity of cycle types accordig to their diversity of pdbids

def analyze_pdb_ids(store_cycles_pdbids):
    cycle_types=[0,0,0,0]
    for cycle in store_cycles_pdbids:
        a, b, c, d = store_cycles_pdbids[cycle]
        print(a,b,c,d)
        cycle_types[len({a,b,c,d})-1] +=1
    return cycle_types

my_path3 , my_cycle4 , my_store_cycles = identify_l3_structures(my_ppi_ael)
my_store_cycles_pdbids = pdb_ids(my_dico_all_pdbid , my_store_cycles)
my_cycle_types = analyze_pdb_ids(my_store_cycles_pdbids)
print(my_cycle_types)
