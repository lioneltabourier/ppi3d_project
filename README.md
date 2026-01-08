# ppi3d_project

## code directory

### usage
all mains usable as follows (parameters are hard coded)

    python3 main_[complete].py

### data_processing.py
general functions to load data and manage formats read from the files in Networks

### evaluation.py
auxiliary functions of main_prediction.py to evaluate prediction quality

### main_4cycle.py
counting 4-cycles in the graph, and categorizing them according to their pdbid diversity

### main_bs_clusters.py
counting 4-cycles in the graph, and categorizing them according to their binding sites compatibility

### main_prediction.py
predicting edges according to tunable split with various metrics (AA, L3 and variants)

### pair_scores.py
auxiliary functions of main_prediction.py to compute rankings of the metrics

### stats_graph.py
shared auxiliary functions to compute general graph properties (n, m, density, ...)

