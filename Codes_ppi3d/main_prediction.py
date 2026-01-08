from data_processing import *
from pair_scores import *
from evaluation import *
from stats_graph import *

import random
import math
import numpy as np
import time
import matplotlib.pyplot as plt

########## LOADING DATA

### FILE PATHS

#my_file_path = "../Networks/ppi3d_c95_human_net.txt"
my_file_path = "../Networks/ppi3d_c40_human_bsi_net.txt"

### LOADING AS GRAPHS

my_ppi_el = read_data(my_file_path)
my_ppi_ael = create_adjacency_list(my_ppi_el)

### SPLITING DATA

my_train_edges, my_test_edges = random_split(my_ppi_el, test_size=0.10)
my_num_edges = len(my_test_edges)

### TRANSFORM TO AEL (Adjacency Edge List) FORMAT

my_train_ael = create_adjacency_list(my_train_edges)
my_test_ael = create_adjacency_list(my_test_edges)
my_potential_ael = create_potential_adj_list(my_train_ael)

### ADAMIC ADAR SCORING

my_aa_scores = sorting_scores(adamic_adar_scores(my_potential_ael, my_train_ael))
my_tp_fp_aa = tp_fp(my_aa_scores, my_test_ael, 100000, 20)
my_pr_rc_aa = pr_rc(my_tp_fp_aa[0], my_tp_fp_aa[1], my_num_edges)

### L3 SCORING

my_l3_scores = sorting_scores(l3_compute_scores(my_train_ael))
my_tp_fp_l3 = tp_fp(my_l3_scores, my_test_ael, 100000, 20)
my_pr_rc_l3 = pr_rc(my_tp_fp_l3[0], my_tp_fp_l3[1], my_num_edges)

### L3Nf1 SCORING

my_l3Nf1_scores = sorting_scores(compute_all_L3Nf1(my_train_ael))
my_tp_fp_l3Nf1 = tp_fp(my_l3Nf1_scores, my_test_ael, 100000, 20)
my_pr_rc_l3Nf1 = pr_rc(my_tp_fp_l3Nf1[0], my_tp_fp_l3Nf1[1], my_num_edges)

### L3Nf2 SCORING

my_l3Nf2_scores = sorting_scores(compute_all_L3Nf2(my_train_ael))
my_tp_fp_l3Nf2 = tp_fp(my_l3Nf2_scores, my_test_ael, 100000, 20)
my_pr_rc_l3Nf2 = pr_rc(my_tp_fp_l3Nf2[0], my_tp_fp_l3Nf2[1], my_num_edges)

### FIGURES GENERATION

# scoring tables

tab1_aa = my_pr_rc_aa[0]
tab2_aa = my_pr_rc_aa[1]

tab1_l3 = my_pr_rc_l3[0]
tab2_l3 = my_pr_rc_l3[1]
#register_image ("l3-" + "ppi3d_c95_human_bsi_homo_net_90-10.png" , tab1_l3, tab2_l3)

tab1_l3Nf1 = my_pr_rc_l3Nf1[0]
tab2_l3Nf1 = my_pr_rc_l3Nf1[1]
#register_image ("l3Nf1-" + "ppi3d_c95_human_net_75-25.png" , tab1_l3Nf1, tab2_l3Nf1)

tab1_l3Nf2 = my_pr_rc_l3Nf2[0]
tab2_l3Nf2 = my_pr_rc_l3Nf2[1]
#register_image ("l3Nf2-" + "ppi3d_c95_human_net_75-25.png" , tab1_l3Nf2, tab2_l3Nf2)

# options

plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve')
plt.grid()
plt.xticks(np.arange(0, 1, 0.05))
plt.yticks(np.arange(0, 1, 0.25))
plt.xlim(0, 0.1)
plt.ylim(0, 1.0)
plt.plot(tab1_aa , tab2_aa , marker = '.' , label="AA")
plt.plot(tab1_l3 , tab2_l3 , marker = '.' , label="L3")
plt.plot(tab1_l3Nf1 , tab2_l3Nf1 , marker = '.' , label="L3Nf1")
plt.plot(tab1_l3Nf2 , tab2_l3Nf2 , marker = '.' ,  label="L3Nf2")
plt.legend(loc="upper right")
plt.savefig("ppi3d_c40_human_bsi_net_90-10.pdf", bbox_inches = "tight")
plt.show()

