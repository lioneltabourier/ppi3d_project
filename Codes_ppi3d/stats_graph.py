############## GENERAL FUNCTIONS TO PERFORM STATS ON GRAPHS

def count_nodes(ael):
    """ count the number of nodes of an ael graph """
    counter = 0
    for node in ael:
        if ael[node] != 0:
            counter += 1
    return counter

def num_edges (ael):
    """ count the number of edges of an ael graph """
    m = 0
    for node in ael:
        m += len(ael[node])
    return m/2

def dd (ael):
    """ compute the degree distribution of an ael graph """
    tab = dict()
    for node in ael:
        degree = len(ael[node])
        if degree in tab:
            tab[degree] += 1
        else:
            tab[degree] = 1
    return tab

def avg_degree (tab):
    """ compute the average degree of a graph from its degree table """
    total = 0
    num_nodes = 0
    for i in tab:
        total += i*tab[i]
        num_nodes += tab[i]
    return total/num_nodes , num_nodes

def density (ael):
    """ compute the density of an ael graph """ 
    n = len(ael)
    m = 0
    for node in ael:
        degree = len(ael[node])
        m += degree
    return m/(n*(n-1))

