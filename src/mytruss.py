import sys
import time
from os import listdir
from os.path import isfile, join
from igraph import *
import matplotlib.pyplot as plt
import pandas as pd
from utils import *


colorlist = {0:"red", 1:"orange", 2:"green",3:"yellow", 4:"pink",5:"blue", 6:"azure",7:"cyan",8:"magenta",9:"purple",10:"white",\
    11:"peru",12:"sienna",13:"navy",14:"tomato",15:"violet",16:"plum",17:"thistle",18:"orchid",19:"beige",20:"tan",\
    21:"orchid",22:"lavender",23:"goldenrod",24:"khaki",25:"grey",26:"ivory",27:"salmon",28:"royalblue",29:"limegreen",30:"seagreen",\
    31:"burlywood",32:"coral",33:"black",34:"sandybrown",35:"firebrick"\
    }
    
def triangles(G,nodes=None):
    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs= ( (n,G[n]) for n in G.nbunch_iter(nodes) )
    for v,v_nbrs in nodes_nbrs:
        vs=set(v_nbrs) -set([v])
        ntriangles=0
        for w in vs:
            ws=set(G[w])-set([w])
            ntriangles+=len(vs.intersection(ws))
        yield (v,len(vs),ntriangles)

def edge_support(G):
    neighbors=G.neighborhood()    #neighbors_iter
    nbrs=dict((v.index,set(neighbors[v.index])) for v in G.vs)
    support = {}
    for e in G.es:
        nod1,nod2 = e.source, e.target	
        nod1_nbrs = set(nbrs[nod1])-set([nod1])
        nod2_nbrs = set(nbrs[nod2])-set([nod2])
        sup = len(nod1_nbrs.intersection(nod2_nbrs))
        #G[nod1][nod2]['support'] = sup
        support[(nod1,nod2)] = sup
    #print 'support :', support
    return support
        
def ktruss(G):   
    #G = G.simplify()    #assume graph is simple
    #G.to_undirected(mode=False)
    support = edge_support(G)
    edges=sorted(support,key=support.get)
    bin_boundaries=[0]
    curr_support=0
    for i,e in enumerate(edges):
        if support[e]>curr_support:
            bin_boundaries.extend([i]*(support[e]-curr_support))
            curr_support=support[e]
            
    edge_pos = dict((e,pos) for pos,e in enumerate(edges))
    #print 'edge-pos:',edge_pos
    truss={}         ## initial guesses for truss is support
    neighbors=G.neighborhood()    #neighbors_iter
    #print 'neighbors:', neighbors
    nbrs=dict((v.index,(set(neighbors[v.index])-set([v.index]))) for v in G.vs)
    #nbrs=dict((v.index,set(neighbors[v.index])) for v in G.vs)
    #print 'nbrs:', nbrs
    for e in edges:
      #print 'processing edge : ', e, 'support :', support[e], 'pos:', edge_pos[e]
      u,v =e[0], e[1]
      if not(u == v) :
        common_nbrs = set(nbrs[u]).intersection(nbrs[v])
        #print u,v,'common_nbrs',common_nbrs
        for w in common_nbrs:
            if (u,w) in support :             
               e1 = (u,w)
            else :
               e1 = (w,u)
            if (v,w) in support :
               e2 = (v,w)
            else:
               e2 = (w,v)
            pos=edge_pos[e1]
            if support[e1] > support[e] :
               bin_start=bin_boundaries[support[e1]]
               edge_pos[e1]=bin_start
               edge_pos[edges[bin_start]]=pos
               edges[bin_start],edges[pos]=edges[pos],edges[bin_start]
               bin_boundaries[support[e1]]+=1
            #print 'e1',e1,'support:',support[e1], 'pos:', pos, 'new pos:', edge_pos[e1]
            
            pos=edge_pos[e2]
            if support[e2] > support[e] :
               bin_start=bin_boundaries[support[e2]]
               edge_pos[e2]=bin_start
               edge_pos[edges[bin_start]]=pos
               edges[bin_start],edges[pos]=edges[pos],edges[bin_start]
               bin_boundaries[support[e2]]+=1
            #print 'e2',e2,'support:',support[e2], 'pos:', pos, 'new pos:', edge_pos[e2]
              
            support[e1] =  max(support[e], support[e1]-1)     
            support[e2] =  max(support[e], support[e2]-1)   
            
        truss[e] = support[e] + 2 
        nbrs[u].remove(v)
        nbrs[v].remove(u)
    #print 'Truss: ', truss
    #print 'Sorted Truss: ', sorted(truss,key=truss.get)
    return truss


def get_ktrussProbs(g, name):
    print ("Extracting KTruss features: %s" % name)
    trussness = ktruss(g).values()
    n = len(trussness)
    d = {n:trussness.count(n) for n in range(2,max(trussness)+1)} 
    ktrussprobability = [d[key] / (n * 1.0) for key in sorted(d)]
    return (ktrussprobability)

def TrussDistribution(graphs):
    KTrussSignatures = {}
    KTrussTimings = {}   
    t1 = 0
    for g in graphs:
       G = graphs[g]
       print 'n = ',len(G.nodes())
       print 'm = ',len(G.edges())
       start_time = time.time()
       KTrussSignature = get_ktrussProbs(G, g)
       KTrussSignatures[g] = KTrussSignature
       t1 = time.time() - start_time
       KTrussTimings[g] = [t1]
       saveFeature(g,KTrussSignatures[g], sys.argv[1]+"_TrussDistn.txt")
    print (t1, "seconds")
    saveDict(KTrussTimings, sys.argv[1]+"_TrussDistnTimings.txt")
    return

def distance(filename):
    outputfile = os.path.splitext(filename)[0] +"_JSDistanceMatrix.txt" 
    print "outputfle ", outputfile
    os.remove(outputfile) if os.path.exists(outputfile) else None
    
    sigs = readDict(filename, ',')
    print "Read Distributions", filename
    sigkeys= sigs.keys()
    l = len(sigkeys) 
    JSDistanceMatrix = pd.DataFrame(np.ones((l,l)), index=sigkeys, columns=sigkeys)
    for i in range(l):
       g1name = sigkeys[i]
       for j in range(i,l):
            g2name = sigkeys[j]
            JSDistanceMatrix[g1name][g2name] = JensenShannon_dist(sigs[g1name], sigs[g2name]) 
            JSDistanceMatrix[g2name][g1name] = JSDistanceMatrix[g1name][g2name] 
    JSDistanceMatrix.to_csv(outputfile)   
    return JSDistanceMatrix    
 
def getnodetrussness(graph):
    n = graph.vcount()
    ktrussdict = ktruss(graph)
    nodetruss = [0] * n
    for edge in graph.es:
      source = edge.source		
      target = edge.target
      if not (source == target) :
         t = ktrussdict[(source,target)]
      else:
         t = 0		  
      nodetruss[source] = max(nodetruss[source], t)
      nodetruss[target] = max(nodetruss[target], t)
    
    return nodetruss 
       
def runDir():
       dir_path = sys.argv[1]
       print(dir_path)
       graph_files = [f for f in listdir(dir_path) if \
                            isfile(join(dir_path,f)) ]
       graphs = {f: nx.read_edgelist(join(dir_path, f)) for f in graph_files} 
       print("Read Graphs:",graphs.keys())
       TrussDistribution(graphs)
    
def test():
    if (len(sys.argv) > 1) :
      G = Graph.Read_Ncol(sys.argv[1],directed=False)
      graphname = sys.argv[1].split('/')[-1].split('.')[0]
    G = G.simplify()
    print 'n = ',len(G.vs)
    print 'm = ',len(G.es)
    start_time = time.time()
    coreness = GraphBase.coreness(G)
    print 'Time taken to compute Coreness: ', time.time() - start_time
    n = len(coreness)
    d = {n:coreness.count(n) for n in range(2,max(coreness)+1)}
    print 'd: ', d
    corenessdistn = [d[key] / (n * 1.0) for key in sorted(d)]
    print 'corenessdistn:', corenessdistn
    start_time = time.time()
    trussness = ktruss(G).values()
    print 'Time taken to compute Trussness: ', time.time() - start_time
    print 'trussness:', trussness
    n = len(trussness)
    d = {n:trussness.count(n) for n in range(2,max(trussness)+1)}
    print 'd: ', d
    trussdistn = [d[key] / (n * 1.0) for key in sorted(d)]
    print 'trussdistn:', trussdistn
    
    #nx.draw(G,with_labels = True)
    #plt.show()

def gettrussscore(g):    
    graphname = g.split('/')[-1].split('.')[0]
    outfilename = outputdirectory + graphname + ".trussscore"
    out = open(outfilename,'w')
    
    print '\tReading Graph'
    graph = Graph.Read_Ncol(g, directed=False,weights = True)
    n = graph.vcount()
    m = graph.ecount()
    out.write ( 'Num Vertices:' + str(n) + '  Num Edges:' + str(m) +'\n')
    out.write ('Num zero degree vertices:' + str(len(graph.vs.select(_degree = 0))))
    zerodegreevertices = graph.vs.select(_degree = 0)
    graph.vs.select(_degree = 0).delete()
    n = graph.vcount()
    out.write (' After dropping isolates - Num Vertices:' + str(n)+'\n')
    print '\n\nComputing Trussness'
    out.write('\n\Trussness')
    start_time = time.time()
    ktrussdict = ktruss(graph)
    trussness = ktrussdict.values()
    out.write ('Time to compute trussness: '+ str(time.time() - start_time)+'\n')
    out.write ( 'max trussness: ' + str(max(trussness))+'\n')
    out.write ( 'min trussness: '+ str(min(trussness))+'\n')
    l = len(trussness)
    d = {n:trussness.count(n) for n in range(2,max(trussness)+1)}
    out.write( 'trusscount:'+ str(d)+'\n')
    trussdistn = [d[key] / (l * 1.0) for key in sorted(d)]
    out.write( 'trussdistn:'+ str(trussdistn)+'\n')
    
    nodetruss = [0] * n
    
    print '\tComputing Node Trussness'  
    for edge in graph.es:
      source = edge.source		
      target = edge.target
      t = ktrussdict[(source,target)]		
      
      nodetruss[source] = max(nodetruss[source], t)
      nodetruss[target] = max(nodetruss[target], t)
    
    dn = {n:nodetruss.count(n) for n in range(2,max(nodetruss)+1)}
    out.write('nodetrusscount:'+ str({k: v for k, v in dn.items() if v!=0})+'\n')  #don't print zero values  
    out.close()   
    
    #plot(graph,vertex_color=[colorlist[x] for x in nodetruss],vertex_label=[v.index for v in graph.vs ],edge_label=[ktrussdict[(e.source,e.target)]-2 for e in graph.es],layout=graph.layout("kk"))
    
    
outputdirectory = '../CentralityOutput/TrussScores/'   
if __name__=="__main__":
   gettrussscore(sys.argv[1])
