from heapq import *
from igraph import *
import itertools
import time
import sys
from scipy import asarray as ar,exp
import pandas as pd
import csv
import profile

outputdirectory = '../SoCapOutput/'

start_time = time.time()

heap = []                         # list of entries arranged in a heap
entry_finder = {}                 # mapping of vertexs to entries
REMOVED = '<removed-vertex>'      # placeholder for a removed vertex
stack = []
INFINITY = 1000000

def initializeheap():
    heap = []
    entry_finder = {}
    stack = []
    
def add_vertex(vertex, priority=0):
    'Add a new vertex or update the priority of an existing vertex'
    if vertex in entry_finder:
        remove_vertex(vertex)
    #count = next(counter)
    #entry = [priority, count, vertex]
    entry = [priority, vertex]
    entry_finder[vertex] = entry
    heappush(heap, entry)

def remove_vertex(vertex):
    'Mark an existing vertex as REMOVED.  Raise KeyError if not found.'
    entry = entry_finder.pop(vertex)
    entry[-1] = REMOVED

def pop_vertex():
    'Remove and return the lowest priority vertex. Raise KeyError if empty.'
    while heap:
        #priority, count, vertex = heappop(heap)
        priority, vertex = heappop(heap)
        if vertex is not REMOVED:
            del entry_finder[vertex]
            return priority,vertex
    return -1,Vertex('-1')
    #raise KeyError('pop from an empty priority queue')

    
class Vertex:
   'Vertex class'
   def __init__(self,v):
        self.name = v
   	self.status = 'notfound'
   	self.dist = INFINITY
   	self.fpmsgs = {}		#dictionary will have key(hops),value (pkts)
   	self.bpmsgs = {}		#dictionary will have key(dist),value ({hops: pkts})
   	self.activeICEdges = []
   	self.scv = 0.0
    	
   def __str__(self):
        string = str(self.name)+ ','+self.status + ',' + str(self.dist) + ',' +str(self.fpmsgs)+','+str(self.bpmsgs)+','+str(self.activeICEdges) +','+ str(self.scv)
        return string
          
 
def SortedRankedToFile(graphname, vnames,centrality, measurename):
       'Write sorted ranked centrality measure to file'
       directory = outputdirectory+ graphname + "/"
       if not os.path.exists(directory):
           os.makedirs(directory)
       n = len(vnames)
       data = zip(vnames,centrality)
       pddata =  pd.DataFrame(data, index=range(0,n), columns=["Name",measurename])
       pddata.to_csv(directory + graphname +"."+measurename+".txt",index=False) 
       sortedcol = pddata.sort_values(by=measurename, ascending=False)
       sortedcol[measurename+'DenseRank'] = sortedcol[measurename].rank(method='dense',ascending=False)
       sortedcol[measurename+'MinRank'] = sortedcol[measurename].rank(method='min',ascending=False)
       outputfile = directory + graphname +".sortedranked."+measurename+".txt"
       print outputfile
       sortedcol.to_csv(outputfile,index=False) 
       return 

def b(l):					#benefits function  e^ (-lambda * l)
    'benefits function'
    lambd = 1
    return exp(-1 * lambd * l)

 
def computeSCV(v): 
    'social capital of vertex in V for the current source vertex v'
    currentSCV = 0
    for fphops in v.fpmsgs.keys() :
       fppkts = v.fpmsgs[fphops]
       for bpdist in v.bpmsgs.keys():
          dist = bpdist
          msg = v.bpmsgs[bpdist]
          for bphops in msg:
              bppkts = msg[bphops]
              totalPkts = fppkts * bppkts
              totalHops = fphops + bphops;
              currentSCV = currentSCV + (totalPkts * b(dist) / (totalHops + 1))
    return currentSCV
    
        
def SoCap(graphfile):
  Vertices = {}
  G = defaultdict(list)
#  with open(graphfile, 'r') as fi:
#      for line in fi:
  fi = open(graphfile, 'r') 
  reader = list(csv.reader(fi,delimiter='\t'))
  fields = len( reader[0] )
  print 'fields:',fields
  for line in reader:
    if fields == 2:
        l,r = line
																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																							
        c = 1.0
    elif fields == 3:
        l,r,c = line
    G[l].append([r,float(c)])
    G[r].append([l,float(c)])
  
  
  for vname in G.keys():
      vertex = Vertex(vname)
      Vertices[vname] = vertex
  
  n = len(G.keys())
  L = 17  #rough hack ::i have to make it diameter (The diameter in the worst case is |V| - 1)
          #assuming diameter of coauthor networks is around 25 
          #graphs_pub diameter = 17  #Astro-Ph = 14 #Cond-Mat, HepTh,GrQc = 17 # HepPh = 14
  print 'n = ', n, 'L = ',L
  count = 0
  for vname in G.keys():
      v = Vertices[vname]
      #print 'Processing Vertex: ',vname
      initializeheap()
      v.status = 'found'
      v.dist = 0
      v.fpmsgs = {}
      v.bpmsgs = {}
      v.activeICEdges = [];
      v.fpmsgs[0] = 1		#hops <-- O,pkts <-- 1
      add_vertex(v.name, priority = v.dist)  # Q.add(v)
      while (heap):
         priority,uname = pop_vertex()      # u <-- Q.poll()
         if (priority == -1):
            continue
         #print '\t popped ',uname
         u = Vertices[uname]
         u.status = 'closed'
         stack.append(u.name)		# S.push(u)
         if (min(u.fpmsgs) == L):
              continue;
         for wname,edgeweight in G[uname]:      # for each w adjacent to u 
              #print '\t    its neighbor is vertex:',wname
              w = Vertices[wname]
              if (w.status == 'closed') :
                 continue
              if (w.status == 'found'):
                 if (w.dist < (u.dist +  1.0/edgeweight) ):
                      continue
                 if (w.dist > (u.dist+ 1.0/edgeweight)) :
                      w.fpmsgs = {}
                      w.activeICEdges = [];
                      w.dist = u.dist + 1.0/edgeweight;
                      add_vertex(wname, priority = w.dist)   #Q.decreasePriority(w)
              else:                
                 w.dist = u.dist + 1.0/edgeweight
                 w.status = 'found'
                 add_vertex(wname,priority=w.dist)
              for hops in u.fpmsgs.keys():              #for each msg in u.fpmsgs
                 pkts = u.fpmsgs[hops]
                 if (hops == L):
		   continue
		 newhops = hops + 1
		 if not newhops in w.fpmsgs:		#w.fpmsgs.add( hops <- msg.hops() + 1,pkts <- msg.pkts
                    w.fpmsgs[newhops] = pkts
                 else:
                    w.fpmsgs[newhops] += pkts
                 w.activeICEdges.append([u,w]);
                 
      while(stack):
         uname = stack.pop()
         u = Vertices[uname]
         u.bpmsgs[u.dist] = {0:1 }                           #u.bpmsgs.add(dist <- u.dist, hops <- O,pkts <- 1
         u.scv = u.scv + computeSCV(u)
         for w,u in u.activeICEdges :		                   # for each e(w,u) in u.activeICEdge
             for dist in u.bpmsgs.keys() :  
                 msg = u.bpmsgs[dist]                    # for each msg in u.bpmsgs
                 for hops in msg.keys():
                     pkts = msg[hops]
                     newhops = hops + 1
                     if not dist in w.bpmsgs:
                        w.bpmsgs[dist] = {}
		     if not newhops in w.bpmsgs[dist]:
                        w.bpmsgs[dist] = {newhops : pkts}
                     else:
                        w.bpmsgs[dist][newhops] += pkts           # w.bpmsgs.add(dist <- msg.dist, hops <-msg.hops + 1, pkts <- msg.pkts)
                             
         u.status = 'notfound'
         u.dist = INFINITY
         u.fpmsgs = {}
         u.bpmsgs = {}
         u.activeICEdges = []
      #print 'vertex:',v.name,'SCV:',v.scv
      count = count + 1
      if (count % 5 == 0) :
         print 'processed ', count, ' of ', n, ' vertcies in Time:', str(time.time() - start_time)

  vnames = [Vertices[v].name for v in Vertices.keys()]
  socap = [Vertices[v].scv for v in Vertices.keys()]

  graphname = graphfile.split('/')[-1].split('.')[0]
  SortedRankedToFile(graphname,vnames,socap,'SoCap' )
  



SoCap(sys.argv[1])  
print 'Time to compute SoCap: '+ str(time.time() - start_time)
#profile.run('print SoCap(sys.argv[1])')

  

  
  
  
  
  
  
  
  
    
    
    
    
    
    
    
    

  
