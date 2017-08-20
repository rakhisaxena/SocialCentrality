import time
#from pylab import *
import sys
from os import listdir
from os.path import isfile, join
from igraph import *
import csv
import pandas as pd
import math
from scipy import asarray as ar,exp
from scipy import stats
import numpy as np
import numpy.random as nprnd
from WeightedCore import *
from mytruss import *
from scipy.stats import beta
import matplotlib.pyplot as plt
from operator import mul


from utils import *

outputdirectory = '../CentralityOutput/'

def getcommcentrality(graph):
    n = graph.vcount()
    comcen=[0] * n
    community = [0]*n
    kin = [0] * n
    kout = [0] * n
    
    
    degree = graph.degree()
    communities = Graph.community_multilevel(graph) 
    level = 0
    for com in communities:
        for v in com:
            community[v] = level
        level = level + 1

    for v in graph.vs:
        vindex = v.index
        for w in graph.neighbors(v):
            if (community[vindex] == community[w]):
               kin[vindex] = kin[vindex] + 1
            else:
               kout[vindex] = kout[vindex] + 1
               
    mu = [0] * len(communities)
    maxkin = [0] * len(communities)
    maxkout = [0] * len(communities)
    R = [0] * len(communities)
    level = -1
    for com in communities:
        level = level + 1
        maxin = 0
        maxout = 0
        for v in com:
            mu[level] = mu[level] + (kout[v]*1.0/degree[v])
            if (maxin < kin[v]):
               maxin = kin[v]
            if (maxout < kout[v]):
               maxout = kout[v]
        mu[level]= mu[level] / len(com)
        R[level] = maxin
        maxkin[level] = maxin
        maxkout[level] = maxout

    for v in graph.vs:
        vindex = v.index
        c = community[vindex]
        if not(maxkout[c] == 0) :
           firstterm = (1+mu[c]) * kin[vindex] 
           secondterm = (1-mu[c]) * kout[vindex]/maxkout[c] * R[c]  
           comcen[vindex] = firstterm + (secondterm*secondterm)
    
    #print comcen 
    return comcen

    

def SortedRankedToFile(graphname, vnames,centrality, measurename):
       directory = outputdirectory+ graphname + "/"
       if not os.path.exists(directory):
           os.makedirs(directory)
       n = len(vnames)
       data = zip(vnames,centrality)
       pddata =  pd.DataFrame(data, index=range(0,n), columns=["Name",measurename])
       #pddata.to_csv(directory + graphname +"."+measurename+".txt",index=False) 
       sortedcol = pddata.sort_values(by=measurename, ascending=False)
       sortedcol[measurename+'DenseRank'] = sortedcol[measurename].rank(method='dense',ascending=False)
       sortedcol[measurename+'MinRank'] = sortedcol[measurename].rank(method='min',ascending=False)
       outputfile = directory + graphname +".sortedranked."+measurename+".txt"
       print outputfile
       sortedcol.to_csv(outputfile,index=False) 
       return




if __name__=="__main__":
    graph = Graph.Read_Ncol(sys.argv[1], directed=False,weights=True)
    graph = graph.simplify()
    comc = getcommcentrality(graph)
    
    graphname = sys.argv[1].split('/')[-1].split('.')[0]
    vnames = [v["name"] for v in graph.vs]
    #print vnames
    SortedRankedToFile(graphname,vnames,comc,'COMC' )
