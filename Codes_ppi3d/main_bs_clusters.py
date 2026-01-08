from data_processing import *
from stats_graph import *

import random
import math
import numpy as np
import time

my_file_path = "../Networks/ppi3d_c95_human_net.txt"

my_file_with_bsc_path = "../Networks/ppi3d_c95_human_net_with_bsc40.txt"
#my_file_with_bsc_path = "../Networks/ppi3d_c95_human_net_with_bsc70.txt"
#my_file_with_bsc_path = "../Networks/ppi3d_c95_human_net_with_bsc95.txt"

### reading EDGE LIST

my_ppi_el = read_data(my_file_path)

### transform EDGE LIST format into ADJACENCY EDGE LIST format

my_ppi_ael = create_adjacency_list(my_ppi_el)

### reading EDGE LIST with BSC

my_dict_edge_bsc = read_data_bsc(my_file_with_bsc_path)

### find and register 4-nodes cycles

def identify_l3_structures (ael):
    memo = set()
    path3  = 0
    cycle4 = 0 
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
                                    memo.add((x,u,v,y))
    return path3 , cycle4 , memo

### for each cycle, categorize it regarding the bsc compatibilities

def categorize_cycle_bsc (memo, D_bsc):
    cat_1 = set() # cat_1: binding sites of individual nodes consistent
    cat_1_1 = set() # cat_1_1: bsc(u)-bsc(y) and bsc(x)-bsc(v) consistent (full jigsaw)
    cat_1_2 = set() # cat_1_2: bsc(u)-bsc(y) xor bsc(x)-bsc(v) consistent (partial jigsaw)
    cat_1_3 = set() # cat_1_3: bsc(u)-bsc(y) and bsc(x)-bsc(v) not consistent (no jigsaw)
    cat_2 = set() # cat_2: binding sites of individual nodes inconsistent
    
    for (x,u,v,y) in memo:
        
        # consistence in x binding sites
        x_bsc_xu = {bsc for bsc , _ in D_bsc[(x,u)]}
        x_bsc_yx = {bsc for _ , bsc in D_bsc[(y,x)]}
        x_bsc = x_bsc_xu & x_bsc_yx

        # consistence in u binding sites
        u_bsc_xu = {bsc for _ , bsc in D_bsc[(x,u)]}
        u_bsc_uv = {bsc for bsc , _ in D_bsc[(u,v)]}
        u_bsc = u_bsc_xu & u_bsc_uv

        # consistence in v binding sites
        v_bsc_uv = {bsc for _ , bsc in D_bsc[(u,v)]}
        v_bsc_vy = {bsc for bsc , _ in D_bsc[(v,y)]}
        v_bsc = v_bsc_uv & v_bsc_vy

        # consistence in y binding sites
        y_bsc_vy = {bsc for _ , bsc in D_bsc[(v,y)]}
        y_bsc_yx = {bsc for bsc , _ in D_bsc[(y,x)]}
        y_bsc = y_bsc_vy & y_bsc_yx        

        if (len(x_bsc) > 0 and len(u_bsc) > 0 and len(v_bsc) > 0 and len(y_bsc) > 0):
            cat_1.add((x,u,v,y))    
            
            print(x,u,v,y,x_bsc & v_bsc,u_bsc & y_bsc)
            
            if (len(x_bsc & v_bsc) > 0) and (len(u_bsc & y_bsc) > 0): 
                cat_1_1.add((x,u,v,y))
                #print (x_bsc & v_bsc , u_bsc & y_bsc)
            elif (len(x_bsc & v_bsc) > 0) or (len(u_bsc & y_bsc) > 0):
                cat_1_2.add((x,u,v,y))
            else :
                cat_1_3.add((x,u,v,y))
            
        else :
            cat_2.add((x,u,v,y))
            
    return cat_1 , cat_2 , cat_1_1 , cat_1_2 , cat_1_3

my_path3 , my_cycle4 , my_memo = identify_l3_structures(my_ppi_ael)
my_cat_1 , my_cat_2 , my_cat_1_1 , my_cat_1_2 , my_cat_1_3 = categorize_cycle_bsc(my_memo, my_dict_edge_bsc)
print(len(my_cat_1) , len(my_cat_2) , "-" , len(my_cat_1_1) , len(my_cat_1_2) , len(my_cat_1_3))





