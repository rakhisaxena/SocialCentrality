import time
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
from mytruss import *
from scipy.stats import beta
import matplotlib.pyplot as plt
from operator import mul


from utils import *

# import local config: Set your local paths in dev_settings.py
DATA_URL=""
SAVE_URL=""
try:
    from dev_settings import *
except ImportError:
    pass

outputdirectory = '../CentralityOutput/'

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

def getCentralities(g) :
       graphname = g.split('/')[-1].split('.')[0]
       directory = outputdirectory+ graphname + "/"
       if not os.path.exists(directory):
           os.makedirs(directory)
       outputfile = directory + graphname +".Centralities"
       timingoutputfile = directory + graphname +".CentralitiesTiming"
       os.remove(outputfile) if os.path.exists(outputfile) else None 
       os.remove(timingoutputfile) if os.path.exists(timingoutputfile) else None 
 
       graph = Graph.Read_Ncol(g, directed=False,weights=True)
       strength = graph.strength(weights = graph.es['weight'])
       if (sum(strength) == 0.0):
          graph.es['weight'] = [1.0] * graph.ecount()     #assign weight 1 to each edge
          strength = graph.strength(weights = graph.es['weight'])


       n = graph.vcount()
       print 'Num Vertices:',n
       graph.vs.select(_degree = 0).delete()
       n = graph.vcount()
       print 'After dropping isolates - Num Vertices:',n
       print 'graph n,m: ', len(graph.vs), len(graph.es)
       vnames = [v["name"] for v in graph.vs]
       vindices = [v.index for v in graph.vs]
       
       Timings = {}   
 
       print 'Computing DegreeCentrality'
       start_time = time.time()
       degreecentrality = graph.degree()
       Timings["DC"] =  time.time() - start_time
       SortedRankedToFile(graphname,vnames,degreecentrality,'DC' )
       
       print 'Computing Laplacian Centrality'
       start_time = time.time()
       lccentrality =  LaplacianCentrality(graph)   
       Timings["LC"] = time.time() - start_time 
       SortedRankedToFile(graphname,vnames,lccentrality,'LC' )
       
       print 'Computing EigenVector Centrality'
       start_time = time.time()
       eccentrality = graph.eigenvector_centrality(directed=False,scale=True,return_eigenvalue=False)
       #eccentrality = graph.eigenvector_centrality(directed=False, weights=graph.es['weight'] , scale=True,return_eigenvalue=False)
       #cost = [1/x for x in graph.es['weight']]
       #eccentrality = graph.eigenvector_centrality(directed=False, weights=cost , scale=True,return_eigenvalue=False)
       Timings["EC"] =  time.time() - start_time
       SortedRankedToFile(graphname,vnames,eccentrality,'EC' )
      print 'Computing NC Centrality'
       start_time = time.time()
       nc  = graph.constraint(weights=graph.es['weight'] )
       ncn = [1.0 if math.isnan(x) or x == 0 else x for x in nc]	

       ncccentrality = [1.0 / x for x in ncn]
       Timings["NC"]   = time.time() - start_time
       SortedRankedToFile(graphname,vnames,ncccentrality,'NC' )
            
       print 'Computing k-truss SC Centrality'
       start_time = time.time()
       tsccentrality, bonding, bridging =  trussSCT(graph)     #,graphname)  
       Timings["TSC"] = time.time() - start_time 
       SortedRankedToFile(graphname,vnames,tsccentrality,'TSC' )
       SortedRankedToFile(graphname,vnames,bonding,'bonding' )
       SortedRankedToFile(graphname,vnames,bridging,'bridging' )
         
       saveDict(Timings, timingoutputfile)
       return          

def LaplacianCentrality(graph, vs=None):
    if vs is None:
        vs = xrange(graph.vcount())
    degrees = graph.degree(mode="all")
    result = []
    for v in vs:
        neis = graph.neighbors(v, mode="all")
        result.append(degrees[v]**2 + degrees[v] + 2 * sum(degrees[i] for i in neis))
    return result    

def getnodeedgetrussness(graph):
    n = graph.vcount()
    m = graph.ecount()
    ktrussdict = ktruss(graph)
    nodetruss = [0] * n
    edgetruss = [0] * m
    for edge in graph.es:
      source = edge.source		
      target = edge.target
      if not (source == target) :
         t = ktrussdict[(source,target)]
         edgetruss[edge.index] = t
      else:
         t = 0		  
      nodetruss[source] = max(nodetruss[source], t)
      nodetruss[target] = max(nodetruss[target], t)
    
    return nodetruss, edgetruss
    
def trussSCT(graph):
    n = graph.vcount()
    degree = graph.degree()
    socialcapital=[0] * n
    
    strength = graph.strength(weights = graph.es['weight'])
           
    if (sum(strength) == 0.0) :
       graph.es['weight'] = [1.0] * graph.ecount()
       strength = graph.strength(weights = graph.es['weight'])
       
    vnames = [v["name"] for v in graph.vs]

    nodetruss,edgetruss = getnodeedgetrussness(graph)
    bonding = [0] * n
    bridging = [0] * n
    netrusses = [set() for i in range(n)] 
    
    for edge in graph.es:
      source = edge.source		
      target = edge.target		
      thisedgeweight = edge['weight']
      
      sourcetruss = nodetruss[source]
      targettruss = nodetruss[target]
      
      thisedgetruss = edgetruss[edge.index]
      
      if (sourcetruss == targettruss == thisedgetruss) :
         bonding[source] += strength[target]  *  targettruss
         bonding[target] += strength[source] * sourcetruss
         
      else :
         bridging[source] += thisedgeweight * targettruss
         bridging[target] += thisedgeweight  * sourcetruss
      
      #netrusses[source].add(targettruss)
      #netrusses[target].add(sourcetruss)
      '''
      m = max(sourcetruss,targettruss) 
      theta = thisedgetruss / (m * 1.0)
      bonding[source] += theta * strength[target] *  targettruss
      bonding[target] += theta * strength[source] * sourcetruss
      if (thisedgetruss == 2):
         bridging[source] += degree[target]
         bridging[target] += degree[source]
      '''          
    #print 'strength', strength
    #print 'nodetruss',nodetruss
    #print 'bonding',bonding
    #print 'bridging', bridging
    #bonding = [a*b for a,b in zip(bonding,nodetruss)]  
    #bridging = [a/(b*1.0) for a,b in zip(bridging,degree)]
    #diversity = [len(x) for x in  netrusses]   
    socialcapital = [a*(b+1)*(c+1) for a,b,c in zip(strength,bonding,bridging)]    #best result : strength,diversity...not sumtrussscore
    #socialcapital = [a*(b+c) for a,b,c in zip(strength,bridging,bonding)] 
    return socialcapital,bonding,bridging


dir_path = join(DATA_URL, sys.argv[1])
if __name__=="__main__":

    graph = Graph.Read_Ncol(sys.argv[1], directed=False,weights=True)
    #graph = graph.simplify()
    #graph.es['weight'] = [1.0] * graph.ecount()     #assign weight 1 to each edge
    #graphname = sys.argv[1].split('/')[-1].split('.')[0]
    #vnames = [v["name"] for v in graph.vs]
    
    #tsc,bo,br = trussSCT(graph)
    #SortedRankedToFile(graphname,vnames,tsc,'TSC' )
    #SortedRankedToFile(graphname,vnames,inf,'IF' )
    #SortedRankedToFile(graphname,vnames,ent,'ENTROPY' )
    #SortedRankedToFile(graphname,vnames,brif,'BRIF' )

    
    getCentralities(sys.argv[1])



