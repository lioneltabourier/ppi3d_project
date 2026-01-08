import numpy as np
import matplotlib.pyplot as plt
import random

############### PREDICTION EVALUATION METRICS COMPUTATION AND PLOTTING

def sorting_scores(score_tab):
    """ generate a list of pairs by decreasing score value from a dictionary of scores of the form d[(i,j)]=score """
    pre_sorting =  sorted(score_tab.items(), key=lambda v: random.random())
    sorted_scores = sorted(pre_sorting, key=lambda x: x[1], reverse=True)  
    #sorted_scores = sorted(score_tab.items(), key=lambda x: x[1], reverse=True)  
    return list(sorted_scores)

def tp_fp(score_lst, test_ael, num_pred_max, step):
    """ generate a triplet of lists with the number of predictions, number of TP predictions , number of FP predictions every step until num_pred_max """
    lst_tp = []
    lst_fp = []
    lst_pred = []
    num_pred = 0
    num_tp = 0
    num_fp = 0
    for pred in score_lst[:num_pred_max+1] :
        (i,j), score = pred
        if i in test_ael and j in test_ael[i]:
            num_tp += 1
        else:
            num_fp += 1
        num_pred +=1
        if (num_pred-1) % step == 0:
            lst_tp.append(num_tp)
            lst_fp.append(num_fp)
            lst_pred.append(num_pred)
    return lst_pred, lst_tp, lst_fp

def pr_rc(lst_pred, lst_tp, num_edges):
    """ generate a couple of lists with the precision and recall for every number of prediction"""
    lst_pr = []
    lst_rc = []
    for i in range(len(lst_pred)):
        lst_pr.append(lst_tp[i]/lst_pred[i])
        lst_rc.append(lst_tp[i]/num_edges)
    return lst_rc, lst_pr  

### PLOTS

def show_pr (tab1, tab2):
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall Curve')
    plt.grid()
    plt.xticks(np.arange(0, 1, 0.05))
    plt.yticks(np.arange(0, 1, 0.25))
    plt.xlim(0, 0.5)
    plt.ylim(0, 1)
    plt.plot (tab1 , tab2 , marker = '.' , linewidth=1)
    plt.show()

def register_image (filename, tab1, tab2):
    plt.figure(figsize=(12,12))
    plt.plot (tab1 , tab2 , marker = '.' , linewidth=1)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall Curve')
    plt.grid()
    plt.xlim(0, 0.5)
    plt.ylim(0, 1)  
    plt.xticks(np.arange(0, 1.01, 0.05))
    plt.yticks(np.arange(0, 1.01, 0.25))
    plt.savefig(filename, bbox_inches = "tight", marker = '.' , linewidth=1)
    plt.show()
