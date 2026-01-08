import math
import numpy as np

############### FUNCTIONS TO COMPUTE PREDICTION SCORES

def adamic_adar_scores(potential_adj_list, adj_list):
    """ compute Adamic-Adar scores of all potential links in potential_adj_list using adj_list structure with set intersections """    
    aa_scores = {}
    for k, potential_links in potential_adj_list.items():
        for node_i, node_j in potential_links:
            if (node_i, node_j) not in aa_scores and (node_j, node_i) not in aa_scores:
                neighbors_i = set(adj_list[node_i])
                neighbors_j = set(adj_list[node_j])
                common_neighbors = neighbors_i.intersection(neighbors_j)
                aa_score = 0
                for neighbor in common_neighbors:
                    degree = len(adj_list[neighbor])
                    if degree > 1:
                        aa_score += 1 / math.log(degree)
                aa_scores[(node_i, node_j)] = aa_score
    return aa_scores

def l3_compute_scores(adj_list):
    """ compute L3 scores of all potential links of unconnected nodes at distance 3 using adj_list structure with set intersections """    
    l3_scores = {}
    for u in adj_list:
        ku = len(adj_list[u])
        for v in adj_list[u]:
            if v > u:
                kv = len(adj_list[v])                
                for i in adj_list[u]:
                    if i != v:
                        for j in adj_list[v]:
                            if j !=u and j != i and i not in adj_list[j]:
                                if ((min(i,j),max(i,j))) not in l3_scores:
                                    l3_scores[(min(i,j),max(i,j))] = 1/np.sqrt(ku*kv)
                                else : 
                                    l3_scores[(min(i,j),max(i,j))] += 1/np.sqrt(ku*kv)
    return l3_scores

def f1(a,b):
    """ f1 in the sense of Yuen et al.: symmetric simple ratio"""
    if len(a) == 0:
        return 0
    else:
        return len(a.intersection(b))/len(a)

def f2(a,b):
    """ f2 in the sense of Yuen et al.: jaccard"""
    if len(a.union(b)) == 0:
        return 0
    else:
        return len(a.intersection(b))/len(a.union(b))

def compute_all_L3Nf1(adj_list):
    """ compute L3N scores with f1 of all potential links of unconnected nodes at distance 3 using adj_list structure with set intersections """
    neighbors = {node: set(adj_list[node]) for node in adj_list}
    snd_neighbors = {node: set() for node in adj_list}
    L3N_scores = {}
    for u in adj_list:
        for v in adj_list[u]:
                for i in adj_list[u]:
                    if i != v:
                        snd_neighbors[i].add(v)
    for i in adj_list:
        for v in snd_neighbors[i]: 
            for j in adj_list[v]:
                if j > i and j not in neighbors[i] and (i,j) not in L3N_scores:        
                    pair = (i,j)
                    L3N_score = compute_L3Nf1(i, j, neighbors, snd_neighbors)
                    L3N_scores[pair] = L3N_score                
    return L3N_scores

def compute_L3Nf1(x, y, neighbors, snd_neighbors):
    """ subroutine of compute_all_L3Nf1, which computes L3N with f1 score for one pair only """
    U = neighbors[x].intersection(snd_neighbors[y]) 
    V = neighbors[y].intersection(snd_neighbors[x])
    if not U or not V:
        return 0
    term1 = len(U) / len(neighbors[x])
    term2 = len(V) / len(neighbors[y])
    term3 = 0
    for u in U:
        for v in neighbors[u].intersection(V):
            Nu_V = len(neighbors[u].intersection(V) - {x}) / len(neighbors[u] - {x})
            f_Nu_V = Nu_V if Nu_V != 0 else 1
            Nv_U = len(neighbors[v].intersection(U) - {y}) / len(neighbors[v] - {y})
            f_Nv_U = Nv_U if Nv_U != 0 else 1
            Nx_Nv = len(neighbors[x].intersection(neighbors[v]) - {y}) / len(neighbors[v] - {y})
            f_Nx_Nv = Nx_Nv if Nx_Nv != 0 else 1            
            Ny_Nu = len(neighbors[y].intersection(neighbors[u]) - {x}) / len(neighbors[u] - {x})
            f_Ny_Nu = Ny_Nu if Ny_Nu != 0 else 1
            term3 += f_Nu_V * f_Nv_U * f_Nx_Nv * f_Ny_Nu
    return term1 * term2 * term3

def compute_all_L3Nf2(adj_list):
    """ compute L3N scores with f2 of all potential links of unconnected nodes at distance 3 using adj_list structure with set intersections \
    identical to compute_all_L3Nf1 except for the subroutin call """
    neighbors = {node: set(adj_list[node]) for node in adj_list}
    snd_neighbors = {node: set() for node in adj_list}
    L3N_scores = {}
    for u in adj_list:
        for v in adj_list[u]:
                for i in adj_list[u]:
                    if i != v:
                        snd_neighbors[i].add(v)
    for i in adj_list:
        for v in snd_neighbors[i]: 
            for j in adj_list[v]:
                if j > i and j not in neighbors[i] and (i,j) not in L3N_scores:        
                    pair = (i,j)
                    L3N_score = compute_L3Nf2(i, j, neighbors, snd_neighbors)
                    L3N_scores[pair] = L3N_score                
    return L3N_scores

def compute_L3Nf2(x, y, neighbors, snd_neighbors):
    """ subroutine of compute_all_L3Nf2, which computes L3N with f2 score for one pair only """
    U = neighbors[x].intersection(snd_neighbors[y]) 
    V = neighbors[y].intersection(snd_neighbors[x])
    if not U or not V:
        return 0
    term1 = len(U) / len(neighbors[x]) # similar to L3Nf1
    term2 = len(V) / len(neighbors[y]) # similar to L3Nf1
    term3 = 0
    for u in U:
        for v in neighbors[u].intersection(V): # TO OPTIMIZE
            Nu_V = len(neighbors[u].intersection(V) - {x}) / len(V.union(neighbors[u] - {x}))
            f_Nu_V = Nu_V if Nu_V != 0 else 1
            Nv_U = len(neighbors[v].intersection(U) - {y}) / len(U.union(neighbors[v] - {y}))
            f_Nv_U = Nv_U if Nv_U != 0 else 1
            Nx_Nv = len(neighbors[x].intersection(neighbors[v]) - {y}) / len(neighbors[x].union(neighbors[v] - {y}))
            f_Nx_Nv = Nx_Nv if Nx_Nv != 0 else 1            
            Ny_Nu = len(neighbors[y].intersection(neighbors[u]) - {x}) / len(neighbors[y].union(neighbors[u] - {x}))
            f_Ny_Nu = Ny_Nu if Ny_Nu != 0 else 1
            term3 += f_Nu_V * f_Nv_U * f_Nx_Nv * f_Ny_Nu
    return term1 * term2 * term3
