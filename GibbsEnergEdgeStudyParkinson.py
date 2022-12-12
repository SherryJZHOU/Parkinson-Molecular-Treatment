import sys
import scipy
import numpy
import networkx
import matplotlib.pyplot as plt

## adjacency lists 
# Compute and print Gibbs free energy on a node-by-node basis. 
# Produce a "dictionary" of node ID and Gibbs free energy.

# Load the network
fh = open('/Users/sherryzhou/CompBio/BIOGRID-ORGANISM-Homo_sapiens-4.4.197.adjlist.txt','rb')

# Convert adjacency list to graph object.
G = networkx.read_adjlist(fh) # networkx function to convert.

fh.close()

genelist = [] # Read in new Keys list of gene id

kh = open('/Users/sherryzhou/CompBio/GibbsEnergy/LungLUSCexpression/tcgagenelist.txt','r')

for line in kh:
    line = line.strip()
    genelist.append(line)

# Print genelist
kh.close()

express = {} # Dictionary


# ************* compute neighbors *************
# this builds a list of neighbors for each node in G
# that data object is not dict. but looks like a list of list.
print("edge count1", "  ", networkx.number_of_edges(G))
# print "node list in G"
# print G.nodes()

Glist = []
Glist = list(G.nodes())
print("Glist_check  ", Glist)
print("lenglist ",len(Glist))
# sys.exit()

glistsize = numpy.real(numpy.shape(Glist)[0])
# print glistsize
# print "size of Glist", " ", len(Glist) # identical function --> numpy.real(numpy.shape(Glist)[0])

# print "copylistsG", " ", Glist

from collections import defaultdict
# it is a list of lists Ndict keyed to a list of genes Glist
Ndict = defaultdict(list) # use Glist as the key set. it's NOT a dectionary

# for k in range(numpy.real(numpy.shape(Glist)[0])):
for (k, node) in enumerate(Glist):
    print("klist", k, "  ", Glist[k])
    # s = Glist[k]
    # G.neighbors(s)
    Ndict[k].append(G.neighbors(node))
print("Glist99  ",Glist)
print("Gneighbors ", Ndict)

numbedges = G.number_of_edges()

# print("Glist ", Glist)
print("shape glist ", numpy.shape(Glist))

# ***** load large expression set *****
f = open('/Users/sherryzhou/CompBio/GibbsEnergy/LungLUSCexpression/expressiondata.txt','r')

table = []
for line in f:
    line = line.strip()
    line = line.split( )
    table.append(line)
f.close()

col0 = [row[0] for row in table] # start counting at 0
stotalarray = []
slnsarray = []
alnsarray = []
gibbsarray = [[]]
nodelist = []
bettilist = []
## for j in range(len(Glist)): 
    ## print("iterate glist:  ",Glist[j])

for e in range(len(table[0])):
    # first build new expression dictionary i.e. one for each expression set.
    print("chkpt e:   ", e)
    elist = [row[e] for row in table] # one column
    # *** normalize the expression vector prior to zipping into dictionary ****
    print("max1 of elist e:  ", e, "  ", max([float(i) for i in elist]), "  ",min([float(i) for i in elist]))
    emax = max([float(i) for i in elist])
    emin = min([float(i) for i in elist])
    for j in range(len(Glist)): 
        print("iterate glist:  ", e, " ", j, " ", Glist[j], "  ", emax, "  ", emin)
    for i in range(int(len(elist))):
        elist[i] = (float(elist[i])-emin)/(emax-emin)
        print("eliste", e, i,  elist[i])
    print("number of rows elist ", len([float(i) for i in elist]))
    print("number of columns elist  ", len(table[0]))
    print("elist888", genelist)
    express = dict(zip(genelist,elist))
    # print express.items()[1]
    print("express dictionary    ", express)

    # *************COMPUTE EXPRESSION-GRAPH PRODUCT *************
  
    Stotal = 0
    alnstotal = 0
    slnstotal = 0
    slns = 0
    alns = 0
    BigS = 0
    gibbsarray = [[]]
    nbunch = []
    threshold = 65 

    #####################################################################################################################
 
    print("eachcolelist ", elist)
    for (k, kth_node) in enumerate(Glist):
        # Glist is the entire list of genes, length is 7000+
        templist = []
        templist = list(Ndict[k][0])
        print("ktemplist ",k, len(Glist))
        listlen = len(templist)
        print("checkpt2  Glist length and templist  ", len(Glist), "   ", templist, "   ", Glist[1], "  ", len(express))
        S = 0
        for i in range(len(templist)):
            if templist[i] not in express:
                continue
            expresskey = templist[i]
            if float(express[expresskey]) == 0:
                express[expresskey] = 0.001
            S = S + express[expresskey] # sum of all expressions
            print("k ,i, express", k, i,express[expresskey])
        
        if  kth_node in express and express[kth_node] != 0 :
            BigS = express[kth_node]/(express[kth_node] + S)
            slns = BigS*numpy.log(BigS)
            alns = express[kth_node]*numpy.log(BigS)
            Stotal = Stotal + BigS
            slnstotal = slns + slnstotal
            alnstotal = alns + alnstotal # total network gibbs energy
        else:
            continue
        print("templist99  ",templist)
        print("degree99  ", networkx.degree(G,kth_node))
        
        # print "Ndict[0]  ", Ndict[k][0], "Glist[k]  ", Glist[k], "  ", alns
        print("expanded chkpt k  ", k, "   single node ", kth_node, "  ", alns, "expression     ", express[kth_node], "   ", numbedges)
        # gibbsarray.append(Glist[k] alns)        
        # alns is devided by degree
        # gibbsarray.append([alns/networkx.degree(G,Glist[k]), Glist[k]])
        gibbsarray.append([alns, kth_node])
        print("stuff to look at:", k, "  ", alns, "   ", kth_node)

    print("Sample  at end of k loop ", e, "Stotal  ", Stotal, "total slns  ", slnstotal, "total alns  ", alnstotal)
    print("lenth of gibbsarray  ", len(gibbsarray))
    print("unsorted gibbsarray  ", gibbsarray)

    # gibbsarray = [(float(x),y) for (x,y) in gibbsarray]
    # gibbsarray.sort() # sorts on first field
    gibbsarray = sorted(gibbsarray) # this is a real command in python to sort arrays

    print("cancer ID,  sorted gibbsarray  ", e, "  ", gibbsarray)

    # select nodes above threshold for plotting and study... the persistent homology
    pathwaysum = 0
    for x in range(1,threshold+1):
        nbunch.append(gibbsarray[x][1])  
        '''
        print("chkpt1 bunch  ", gibbsarray[1][0], "  ", gibbsarray[1][1])
        print("chkpt1 bunch  ", gibbsarray[2][0], "  ", gibbsarray[2][1])
        print("chkpt1 bunch  ", e)
        
        energy1 = gibbsarray[1][0]
        node1 = gibbsarray[1][1]
        energy2 = gibbsarray[2][0]
        node2 = gibbsarray[2][1]
        '''
        pathwaysum = pathwaysum + gibbsarray[x][0]

    print("pathwaygibbs ", e, " pathway gibbs   ", pathwaysum)
    print("chkpt1 bunch  ", gibbsarray[1][0], "  ", gibbsarray[1][1])
    print("chkpt1 bunch  ", gibbsarray[2][0], "  ", gibbsarray[2][1])
    print("chkpt1 bunch  ", e)
        
    energy1 = gibbsarray[1][0]
    node1 = gibbsarray[1][1]
    energy2 = gibbsarray[2][0]
    node2 = gibbsarray[2][1]
    print("example dictionary use  ",express[gibbsarray[1][1]])

    #print("get edges of TP53  ",G.edges("TP53"))

    #fileNameTemplate = r'Network'
    # graphics this is the graphics for the homology subnetwork
    print("nbunch:   ", nbunch)

    # ****** GRAPHICS CODE ******  
    
    #plt.figure(e)
    #networkx.draw(G.subgraph(nbunch))
    #plt.savefig('Network' + str(e))
    
    H = G.subgraph(nbunch)
    print("get edges of TP53  ",H.edges(gibbsarray[1][1]))
    print("len of edge list ", len(H.edges(gibbsarray[1][1])))
    edgelist = H.edges(gibbsarray[1][1])
    print("chunk of list ",edgelist[0][0], edgelist[0][1])
    print("second chunk of list ",edgelist[1][0], edgelist[1][1])
    print("check express ",express[edgelist[1][0]],express[edgelist[1][1]])
    large_exp = 0
    for zz in range(0,len(edgelist)):
        print("express tuple", e, edgelist[zz][1], express[edgelist[zz][1]])
        # find the largest express[edgelist[zz[1]] && in which edgelist[zz][1] is also NOT edgelist[0][1]
        if express[edgelist[zz][1]] > large_exp and edgelist[zz][1] != node1:
            large_exp = express[edgelist[zz][1]]
            large_node = edgelist[zz][1]
    print("large_exp  ", e, "  " ,large_exp, large_node)

    Hc = H.copy()
    print
    print("cancer number:   ", e)
    print("adjacency list  ")
    # save adjacency lists    
    networkx.write_edgelist(H, 'Network'+str(e)+'.txt')

    Hlist = []
    Hlist = list(H.nodes())
    edgelist = [[]]
    
    # cycle basis and degree analysis
    bettin = 0
    betti = 0
    for inode in range(0,threshold): # for each node in homology network
        print(" unedited network  ")
        D = networkx.triangles(H)
        print("number of triangles  ", sum(D.values()))
        print("triangles  ", D)
        print("edge count before    ", networkx.number_of_edges(H))
        cyclesbe4 = len(networkx.cycle_basis(H))
        print("cycles before    ", len(networkx.cycle_basis(H))) #length of dictionary
        Dg = dict(H.degree()) # degree dictionary
        print("max degree ", max(Dg.values()))
        print("  e  ", e, "  max degree node  ", max(Dg, key=Dg.get),  "  max degree  ", max(Dg.values()))
        print(Hlist)
        bettin = 0  # betti nominal
        for knode in range(len(networkx.cycle_basis(H))):
            if len(networkx.cycle_basis(H)[knode]) > 3:
                bettin = bettin + 1 # len(networkx.cycle_basis(H)[x])
    print("cancer ID  betti ", e, "  betti unedited  ", bettin, "  cycles  ", networkx.cycle_basis(H))

    ### KNOCKOUT AND COMPUTATION OF BETTI CHANGE
    betti = 0
    for inode in range(0,threshold):
        print("edited network ")
        nrem = Hlist[inode] # node to remove
        print("node removed   ", Hlist[inode])
        del Hlist[inode] # now it is removed
        print(Hlist) # edited node list
        H = G.subgraph(Hlist) # edited network
        print("edge count after    ", networkx.number_of_edges(H))
        cyclesafter = len(networkx.cycle_basis(H))
        print("cycles after    ", len(networkx.cycle_basis(H)))
        D = networkx.triangles(H)
        # compute betti number
        betti = 0
        for knode in range(0,len(networkx.cycle_basis(H))):
            print("check_betti   ", betti)
            if len(networkx.cycle_basis(H)[knode]) > 3:
                print("cycle length", len(networkx.cycle_basis(H)), "  cancer ID  ", e, "  inode  ", inode, "  knode  ", knode, "first betti  ", betti)
                betti = betti + 1 # len(networkx.cycle_basis(H)[x])
        print("main_results  inode  ", inode, "  cancer ID", e, "  node removed  ", nrem, "  betti edited   ", betti, "  betti unedited  ", bettin) #, "   ", networkx.cycle_basis(H)
        # node removed, nrem; and betti could now be zipped into a dictionary for later searching        

        nodelist.append(nrem)
        bettilist.append(betti)
        
        # reset network graph to original
        H = G.subgraph(nbunch) # makes a new copy to start the process over
        print(Hlist)
        Hlist = list(H.nodes())        
    
    l = []
    J = []
    # do something with the nodelist and bettilist
    Dbetti = dict(zip(nodelist,bettilist)) # a dictionary of nodes and the resulting betti from removing that node
    print("cancer ID  ", e, "  nominal betti  ", bettin)
    
    minv = min(Dbetti.values())
    print(minv)
    print("cancer ID  ", e, "  nominal betti  ", bettin, "  betti targets  ", [(k,v) for (k,v) in Dbetti.items() if v == minv])
    L = [(k) for (k,v) in Dbetti.items() if v == minv]
    print("ell0  ", L[0]) # the node to remove from the pathway network, it is the results in the greatest drop in betti
    
    print("nodelistpriortoremoval  ", nodelist)
    Dnodelist = Dbetti.keys()
    #Dnodelist.remove(str(L[0]))
    Hlist.remove(L[0]) # the new node list to construct a new pathway
    print("Dnodelist  ", Hlist, "  L[0]  ", L[0]) 
    # empty nodelist and bettilist 
    del nodelist[:]
    del bettilist[:]
    
    # compute statistic on little network
    degreehist = networkx.degree_histogram(H) # list
    
    cc = networkx.closeness_centrality(H) # dictionary
    trans = networkx.transitivity(H) #s calar number
    cclique = networkx.graph_clique_number(H) # scalar number

    ### COMPUTE ENTROPY OF PATHWAY NETWORK
    entropy = 0.0
    rangevalue = numpy.shape(degreehist)[0]
    for jj in range(rangevalue):
        if degreehist[jj] > 0:
            intermediate = float(degreehist[jj])/sum(degreehist)
            print("checkpt", intermediate)
            entropy = entropy + intermediate*numpy.log(intermediate)
    entropy = -entropy        
    print("cancer ID  ", e, "  entropy  ", entropy)
    
    ### Now with the edited Hlist construct a new pathway and compute Gibbs
    print("Dnodelist  ", Hlist, "  L[0]  ", L[0]) 

    ########################################################################################################################
    Stotals = 0
    alnstotals = 0
    slnstotals = 0
    slnss = 0
    alnss = 0
    BigsSs= 0
    for kk in range(len(Hlist)):
        # Hlisst is the small list of genes
        templist = []
        templist = list(Ndict[kk][0])
        listlen = len(templist)
        #gibbsarray = [[]]
        print("checkpt2.2  Glist length and templist  ", len(Hlist), "   ", templist, "   ", Hlist[1], "  ", len(express))
        Ss = 0
        for ii in range(len(templist)):
            if express.has_key(templist[ii]) != 1:
                continue
            expresskey = templist[ii]
            if float(express[expresskey]) == 0:
                express[expresskey] = 0.001
            Ss = Ss + express[expresskey] # sum of all expressions
        
        if Hlist[kk] in express and express[Hlist[kk]] != 0 :
            BigSs = express[Hlist[kk]]/(express[Hlist[kk]] + Ss)
            slnss = BigSs*numpy.log(BigSs)
            alnss = express[Hlist[kk]]*numpy.log(BigSs)
            
            Stotals = Stotals + BigSs
            slnstotals = slnss + slnstotals
            alnstotals = alnss + alnstotals # this is the pathway gibbs after knockout of node computed by betti
        else:
            continue
    print("chcekpt5 e  ", e, "  alns  ", alnstotals, "  pathwaygibbs   ", pathwaysum)       

    ######################################################################################################################################################

    ##print("large_exp  ", e, "  " ,large_exp, large_node)

    '''
    adjspectra = networkx.adjacency_spectrum(H) # list of complex numbers
    print " cancer ID  ", e , "  adjacency eigenvalues  ", numpy.real(adjspectra[0]), "  ", numpy.real(adjspectra[1]), "  ", numpy.real(adjspectra[2])

    print("chkpt1 bunch  ", gibbsarray[1][0], "  ", gibbsarray[1][1])
    print("chkpt1 bunch  ", gibbsarray[2][0], "  ", gibbsarray[2][1])
 
    lapspectra = networkx.laplacian_spectrum(H) # list of complex numbers
    print " cancer ID  ", e, "  lapspectra eigenvalues  ", numpy.real(lapspectra[0]), "   ", numpy.real(lapspectra[1]), "   ", numpy.real(lapspectra[2])
    '''

    print("The Output Results  ", "edges   ", energy1,  "  ", node1, " ",energy2, "  ",node2 ," ",e, "large_node", large_node, "large_exp ",\
     large_exp, "expression", express[gibbsarray[1][1]], "  ", express[gibbsarray[2][1]], "  ", "  entropy  ", entropy, " nominal betti  ", bettin,\
     "  pathway gibbs  ", pathwaysum,"total gibbs ", alnstotal, "  number of edges  ", H.number_of_edges(), " transitivity  ", trans, "  clique number  ",\
    cclique, "  delta gibbs  ", alnstotals, "  betti targets  ", [(k,v) for (k,v) in Dbetti.items() if v == minv] )

    sys.stdout.flush()
    '''
    print "Main Results for: ", e, "  nominal betti  ", bettin,  "  adjacency eigenvalues  ", numpy.real(adjspectra[0]), "  ", numpy.real(adjspectra[1]), "  ", numpy.real(adjspectra[2]),\
    "  lapspectra eigenvalues  ", numpy.real(lapspectra[0]), "   ", numpy.real(lapspectra[1]), "   ", numpy.real(lapspectra[2]),\
    "  entropy  ", entropy, " transitivity  ", trans, "  clique number  ", cclique, "  number of edges  ", H.number_of_edges(), "  total gibbs ", alnstotal, "  pathway gibbs   ", pathwaysum
    '''
    print("e, nodes  ", e, "  ", Hlist)
    print(express)

    sys.stdout.flush()
