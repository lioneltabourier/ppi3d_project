from sklearn.model_selection import train_test_split
import random
random_seed = 42

############## GENERAL GRAPH FUNCTIONS

def read_data(file_path):
    """ load graph data (edge list format) as an el (edge list) graph """
    edges = []
    with open(file_path, "r") as file:
        for line in file.readlines():
            line = line.replace('_','0')
            nodes = line.strip().split()
            if nodes[0] != nodes[1]:
                edge = (int(nodes[0]), int(nodes[1]))
                edges.append(edge)
    return edges

def el_to_ael(el):
    """ transform el graph to ael graph """
    ael = dict()
    for node,neigh in el:
        if node in ael:
            if neigh not in ael[node]:
                ael[node].append(neigh)
        else:
            ael[node] = [neigh]
        if neigh in ael:
            if node not in ael[neigh]:
                ael[neigh].append(node)
        else:
            ael[neigh] = [node]
    return ael

def random_split(edges, test_size=0.25, random_state=random_seed):
    """split Edge List data into train set / test set follwing the test_size parameter""" 
    train_edges, test_edges = train_test_split(edges, test_size=test_size, random_state=random_state)
    return train_edges, test_edges

def create_adjacency_list(edge_list):
    """ create an el graph into a sorted ael graph"""
    adj_list = {}
    for edge in edge_list:
        src, dest = edge
        if src not in adj_list:
            adj_list[src] = [dest]
        else:
            if dest not in adj_list[src]:
                adj_list[src].append(dest)
        if dest not in adj_list:
            adj_list[dest] = [src]
        else:
            if src not in adj_list[dest]:
                adj_list[dest].append(src)
    # guarantee that adj_list is sorted
    for src in adj_list:
        adj_list[src].sort()
    return adj_list

def create_potential_adj_list(adj_list):
    """ create an Adjacency List containing for all node k its pairs of unconnected neighbors (= distance-2 pairs) """
    potential_adj_list = {}
    for k, neighbors in adj_list.items():
        potential_links = []
        for i in range(len(neighbors)):
            for j in range(i+1, len(neighbors)):
                node_i, node_j = min(neighbors[i], neighbors[j]) , max(neighbors[i], neighbors[j])  
                if node_j not in adj_list[node_i]: # Check that i and j are not connected
                    potential_links.append((node_i, node_j))
        potential_adj_list[k] = potential_links
    return potential_adj_list

def write_rankingfile (filename , ranking):
    """write a ranking into a file"""
    f = open(filename, "w")
    for elt in ranking:
        f.write(str(elt))
        f.write("\n")
    f.close()    

############## FUNCTIONS SPECIFIC TO PDBID DATA

def read_data_pdbid(file_path):
    """ read pdbid data to make a dictionary of edges of the form (node1,node2):pdbid  ('_' in pdbid are replaced by 0 for compatibility with other codes)
    note that the pdbid is randomly selected among all possible ones for this edge (last one read selected) """
    dict_edges_pdbid = dict()
    with open(file_path, "r") as file:
        for line in file.readlines():
            line = line.replace('_','0')
            nodes = line.strip().split()
            if nodes[0] != nodes[1]:
                edge = (int(nodes[0]), int(nodes[1]))
                edge_reciprocal = (int(nodes[1]), int(nodes[0]))
                dict_edges_pdbid[edge]=nodes[2]
                dict_edges_pdbid[edge_reciprocal]=nodes[2]
    return dict_edges_pdbid
    
def read_data_all_pdbid(file_path):
    """ read pdbid data to make a dictionary of edges of the form (node1,node2):list_pdbids  ('_' in pdbid are replaced by 0 for compatibility with other codes) 
    note that all possible pdbids for an edge are in the list od pdbids"""
    dict_edges_pdbid = dict()
    with open(file_path, "r") as file:
        for line in file.readlines():
            line = line.replace('_','0')
            nodes = line.strip().split()
            if nodes[0] != nodes[1]:
                mini = min(int(nodes[0]) , int(nodes[1]))
                maxi = max(int(nodes[0]) , int(nodes[1]))
                edge = (mini,maxi)
                if edge in dict_edges_pdbid:
                    dict_edges_pdbid[edge].append(nodes[2])
                else:
                    dict_edges_pdbid[edge]=[nodes[2]]
    return dict_edges_pdbid

def read_data_bsc(file_path):
    """ load data as dictionary edge:bsc """
    dict_edges_pdbsc = dict()
    with open(file_path, "r") as file:    
        for line in file.readlines():
            line = line.replace('_','0')
            nodes = line.strip().split()
            if nodes[0] != nodes[1]:
                edge = (int(nodes[0]),int(nodes[1]))
                edge_reciprocal = (int(nodes[1]), int(nodes[0]))
                if edge in dict_edges_pdbsc:
                    dict_edges_pdbsc[edge].append((nodes[2],nodes[3]))
                else:
                    dict_edges_pdbsc[edge]=[(nodes[2],nodes[3])]

                if edge_reciprocal in dict_edges_pdbsc:
                    dict_edges_pdbsc[edge_reciprocal].append((nodes[3],nodes[2]))
                else:
                    dict_edges_pdbsc[edge_reciprocal]=[(nodes[3],nodes[2])]
    return dict_edges_pdbsc
                    
