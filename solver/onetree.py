#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: augustinparjadis
"""

import networkx as nx
import sys, math
import random as rd
import numpy as np
import matplotlib.pyplot as plt



def loadGraph(file):
    pts = []
    with open(file) as f:
        for line in f.readlines():
            if line[0].isnumeric():
                l = line.split()
                pts.append((float(l[1]),float(l[2])))
    N = len(pts)
    edges = []
    for i in range(N):
        for j in range(N):
            if i != j:
                d = math.sqrt( (pts[i][0]-pts[j][0])**2 + (pts[i][1]-pts[j][1])**2 )
                e = (i+1,j+1,d)
                edges.append(e)
        
    
    G = nx.Graph()
    # G.add_nodes_from([i+1 for i in range(N)])
    G.add_weighted_edges_from(edges)
    
    return G

nodes_rem = 0
#returns a minimum 1-tree and its weight over G. If theta penalties given, they are added to weight
def one_tree(G,theta=None,remW=True,gnopen=None):
    nodes_rem = 0
    if theta != None:
        assert(len(theta) == len(G))
        pen = [t for t in theta]
        pen[0] = 0.0
    else:
        pen = [0 for v in range(len(G))]
    
        
    #build min spanning tree
    G1 = G.copy() #graph with modified costs on edges
    for e in G1.edges:
        G1[e[0]][e[1]]['weight'] += pen[e[0]]+pen[e[1]]
    G2 = G1.copy() #graph with modified costs on edges, without node 1
    G2.remove_node(nodes_rem)
    T = nx.minimum_spanning_tree(G2,algorithm="prim")
    
    
    #sort edges
    edges = list(G.edges) #[(u,v),(u,v),...]
    w = []
    for e in edges:
        w.append( G.get_edge_data(*e)['weight']+pen[e[0]]+pen[e[1]] )
    w_edges = [(edges[i][0],edges[i][1],w[i]) for i in range(len(w))]
    w_edges = sorted(w_edges, key=lambda tup: tup[2])
    
    #add two min cost edges adj to {nodes_rem}
    compt = 0
    for we in w_edges:
        if compt < 2 and (we[0] == nodes_rem or we[1] == nodes_rem):
            # print(we)
            T.add_edge(we[0],we[1])
            compt += 1
    
    if gnopen != None:
        G1 = gnopen.copy() #graph with modified costs on edges
        for e in G1.edges:
            G1[e[0]][e[1]]['weight'] += pen[e[0]]+pen[e[1]]
        
    tw = 0
    for e in list(T.edges):
        tw += G1.get_edge_data(*e)['weight']#+pen[e[0]-1]+pen[e[1]-1]
    if remW:
        for p in pen:
            tw -= 2*p
        
    return T,tw



def stepTheta(g,thetas, n):
    
    if n <= 0:
        return thetas
    
    C = 0.1
    t, tw = one_tree(g,theta=thetas)
    deg = [0 for i in thetas]
    for e in t.edges():
        deg[e[0]] += 1
        deg[e[1]] += 1
    
    thetas_ = [t-C*(2-deg[i]) for i,t in enumerate(thetas)]
    
    return stepTheta(g,thetas_, n-1)


def improve_until(g,thetas,dfw):
    
    C = 1.0
    tw,it = 0,0
    while tw < dfw and it<1000:
        it += 1
        t, tw = one_tree(g,theta=thetas)
        deg = [0 for i in thetas]
        for e in t.edges():
            deg[e[0]] += 1
            deg[e[1]] += 1
        thetas_ = [t-C*(2-deg[i]) for i,t in enumerate(thetas)]
        thetas = thetas_
    return tw,it

import torch.nn.functional as F
import dgl.nn.pytorch as dglnnpy
import dgl.nn as dglnn
import torch.nn as nn

class model(nn.Module):
    def __init__(self, in_feats, hid_feats, out_feats):
        super().__init__()
        self.conv1 = dglnn.SAGEConv(
            in_feats=in_feats, out_feats=hid_feats, aggregator_type='mean')
        self.conv2 = dglnn.SAGEConv(
            in_feats=hid_feats, out_feats=hid_feats, aggregator_type='mean')
        self.conv3 = dglnn.SAGEConv(
            in_feats=hid_feats, out_feats=hid_feats, aggregator_type='mean')
        self.lin1 = nn.Linear(hid_feats, hid_feats)
        # self.lin2 = nn.Linear(hid_feats, hid_feats)
        self.lin2 = nn.Linear(hid_feats, out_feats)

    def forward(self, graph, inputs, ew):
        H = F.relu(self.conv1(graph, inputs, edge_weight = ew))
        H = F.relu(self.conv2(graph, H, edge_weight = ew))
        H = F.relu(self.conv3(graph, H, edge_weight = ew))
        H = F.relu(self.lin1(H))
        # H = F.relu(self.lin2(H))
        H = self.lin2(H).squeeze()
        return H



#returns the different edge between the 2 trees, other than e
def edge_diff(G,T1,T2,e):
    c = 0
    eps = ()
    for temp_e in G.edges:
        if temp_e != e and (T1.has_edge(*temp_e) != T2.has_edge(*temp_e)):
            c += 1
            eps = temp_e
    return eps



M = 1000000
#forall edge e in G, return {e:(R,I,epsilon)} with replacement/insertion costs and associated swapped edge
def one_tree_costs(G,theta=None):
    
    T,w1tree = one_tree(G,theta)
    
    e_data = {} #forall e in E, return list {e:(R,I,epsilon)}
    
    for e in G.edges:
        if T.has_edge(*e):
            #compute replacement cost
            print(f'Replacing {e}')
            g = G.copy()
            g.remove_edge(*e)
            t,wt = one_tree(g)
            eps = edge_diff(G,T,t,e)
            e_data[e] = (wt-w1tree,0,eps)
            print(f'Replacing {e} with {eps} costs {wt-w1tree}\n')
        else:
            #compute insertion cost
            print(f'Inserting {e}')
            g = G.copy()
            g[e[0]][e[1]]['weight'] -= M
            t,wt = one_tree(g)
            eps = edge_diff(G,T,t,e)
            e_data[e] = (0,wt-w1tree+M,eps)
            print(f'Inserting {e} with {eps} costs {wt-w1tree+M}\n')
    
    # replacement cost
    # forall edges in T\{1}, set cost to infinity (or remove ?) and get new edge as well as C(T\e)-C(T)
    
    # insertion cost
    # forall edges in G\T, set cost to -infinity and get removed edge as well as C(TU{e})-C(T)
    
    return e_data


PATHS = {}#save paths
def find_path(G,T,edge):#paths do not go through 1
    #find path in T from e1 to e2
    e1,e2 = edge
    
    if edge in T.edges():
        return -1
    
    if (G,T,edge) in PATHS:
        return PATHS[(G,T,edge)]
    
    Tedges = list(T.edges())
    todo = [e1] #nodes
    done = [] #nodes and edges
    paths = [[e1]] #node lists
    
    for a in todo:
        temp = [] #nodes to expand to
        
        for i in Tedges:#look for edges with a
            #don't check edges already used
            if (a in i) and (i not in done) and (nodes_rem not in i):
                i2 = i[0] if i[1] == a else i[1]
                temp.append(i2)
                todo.append(i2)
                done.append(i)
                
        for p in paths:#complete the paths in memory
            if p[-1] == a:
                for t in temp:
                    paths.append(p.copy()+[t])
        done.append(a)
        
        for p in paths:
            if p[-1] == e2:
                PATHS[(G,T,edge)] = p
                return p
        
    return []


def insertion_costs(G,T,theta=None):
    #G has modified costs
    I = {e: 0 for e in G.edges()}
    Ieps = {e: (-1,-1) for e in G.edges()}
    
    #pour tout e dans G\T, trouver chemin Pe qui relie e1 et e2 dans T
    for e in G.edges():
        if e not in T.edges():
            
            if nodes_rem in e:
                max_c1 = -100000000
                max_c1_eps = (-1,-1)
                for te in list(T.edges()):
                    if nodes_rem in te:
                        # max_c1 = max(max_c1,G[te[0]][te[1]]['weight'])
                        if G[te[0]][te[1]]['weight'] > max_c1:
                            max_c1 = G[te[0]][te[1]]['weight']
                            max_c1_eps = te
                I[e] = G[e[0]][e[1]]['weight'] - max_c1#Ce-max(cost 1arc in 1Tree)
                Ieps[e] = max_c1_eps
                
            else:
                P = find_path(G,T,e)
                assert(nodes_rem not in P)
                assert(P != [])
                assert(P != -1)
                
                max_ca = -100000000
                max_ca_eps = (-1,-1)
                for i in range(len(P)-1):
                    a = (P[i],P[i+1])
                    if a not in I:
                        a = a[::-1] #flipped
                    # max_ca = max(max_ca,G[a[0]][a[1]]['weight'])
                    if G[a[0]][a[1]]['weight'] > max_ca:
                        max_ca = G[a[0]][a[1]]['weight']
                        max_ca_eps = a
                    
                #Ie = Ce - max(Ca|a in Pe)
                I[e] = G[e[0]][e[1]]['weight'] - max_ca
                Ieps[e] = max_ca_eps
    
    return I,Ieps
    

def replacement_costs(G,T,theta=None):
    
    #set Re = inf pour tout e dans T
    R = {e: M if e in T.edges() else 0 for e in G.edges()}
    Reps = {e: (-1,-1) for e in G.edges()}
    
    
    #pour tout e dans G\T, trouver chemin Pe qui relie e1 et e2 dans T
    for e in G.edges():
        if e not in T.edges():
            
            if nodes_rem in e:
                for te in list(T.edges()):
                    if nodes_rem in te:
                        if te not in R:
                            te = te[::-1] #flipped
                        if G[e[0]][e[1]]['weight']-G[te[0]][te[1]]['weight'] < R[te]:
                            Reps[te] = e
                        R[te] = min(R[te],G[e[0]][e[1]]['weight']-G[te[0]][te[1]]['weight'])
                        
            else:
                P = find_path(G,T,e)
                assert(nodes_rem not in P)
                assert(P != [])
                assert(P != -1)
                #pour tout a dans Pe, Ra = min(Ra,Ce-Ca)
                for i in range(len(P)-1):
                    a = (P[i],P[i+1])
                    if a not in R:
                        a = a[::-1] #flipped
                    if G[e[0]][e[1]]['weight']-G[a[0]][a[1]]['weight'] < R[a]:
                        Reps[a] = e
                    R[a] = min(R[a],G[e[0]][e[1]]['weight']-G[a[0]][a[1]]['weight'])
        
    return R,Reps
    

def one_tree_costs_fast(G,theta_tree=None,theta_RI_costs=None):
    
    T,w = one_tree(G,theta_tree)
    
    if theta_RI_costs != None:
        assert(len(theta_RI_costs) == len(G))
        pen = [t for t in theta_RI_costs]
    else:
        pen = [0 for v in range(len(G))]
        
    G1 = G.copy() #graph with modified costs on edges
    for e in G1.edges:
        G1[e[0]][e[1]]['weight'] += pen[e[0]]+pen[e[1]]
       
    I,Ieps = insertion_costs(G1,T)
    R,Reps = replacement_costs(G1,T)
    
    e_data = {} #forall e in E, return list {e:(R,I,epsilon)}
    
    for e in G.edges:
        if T.has_edge(*e):
            e_data[e] = (R[e],0,Reps[e])
        else:
            e_data[e] = (0,I[e],Ieps[e])
    
    return e_data
    


def costs_derivative(G,th_true,th_est):
    
    # Gtrue with theta_true
    # R,I -> calculés à partir de f( x(th^),th ) ; réels
    # edata_true = one_tree_costs_fast(G,theta_tree=th_est,theta_RI_costs=th_true)
    edata_true = one_tree_costs_fast(G,theta_tree=th_est,theta_RI_costs=None)
    
    # Gest with theta_est
    # R^,I^ -> calculés à partir de f( x(th^),th^ ) ; >0
    edata_est = one_tree_costs_fast(G,theta_tree=th_est,theta_RI_costs=th_est)
    
    return edata_true,edata_est


#1tree computed thanks to costs/penalties
#get derivatives dx/dtheta -> x is 0/1 in the 1tree, theta are the panlties at nodes
def derivative(G,T,e_data):
    
    #nodes degrees
    deg = [0 for n in T.nodes]
    for e in T.edges:
        deg[e[0]-1] += 1
        deg[e[1]-1] += 1
    # print(deg)
    
    D = [[0 for x in G.edges] for t in G.nodes]
    for e in G.edges:
        if T.has_edge(*e):
            R = max(1/M,e_data[e][0])
            eps = e_data[e][2]
            i1 = e[0]-1
            i2 = e[1]-1
            D[i1][list(G.edges).index(e)] = -(2-deg[i1])/R
            D[i2][list(G.edges).index(e)] = -(2-deg[i2])/R
            D[i1][list(G.edges).index(eps)] = (2-deg[i1])/R
            D[i2][list(G.edges).index(eps)] = (2-deg[i2])/R
        else:
            I = max(1/M,e_data[e][1])
            eps = e_data[e][2]
            i1 = e[0]-1
            i2 = e[1]-1
            D[i1][list(G.edges).index(e)] = -(2-deg[i1])/I
            D[i2][list(G.edges).index(e)] = -(2-deg[i2])/I
            D[i1][list(G.edges).index(eps)] = (2-deg[i1])/I
            D[i2][list(G.edges).index(eps)] = (2-deg[i2])/I
    return D
