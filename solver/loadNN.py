#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: augustinparjadis
"""

import random as rd
import numpy as np
import sys,os,math,time,pickle
import torch.nn.functional as F
import dgl.nn.pytorch as dglnnpy
import dgl.nn as dglnn
import torch.nn as nn
import torch
import dgl
import matplotlib.pyplot as plt
import networkx as nx
from torch.autograd import Function
from torch import autograd
from statistics import mean
import onetree


def getLogits(edges):
    
    k = edges[0]
    edges = edges[1:]
    
    graph = dglTSPGraph(edges,k)
    feats = graph.ndata['x']
    e_weights = graph.edata['weight']
    logits = gnnmodel(graph, feats, e_weights)

    return logits.tolist()


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
    G.add_weighted_edges_from(edges)
    
    return G


class model(nn.Module):
    def __init__(self, in_feats, hid_feats, out_feats, edges_feats, num_heads):
        super().__init__()
        self.conv1 = dglnnpy.conv.EGATConv(in_node_feats=in_feats,
                in_edge_feats=edges_feats,
                out_node_feats=int(hid_feats/num_heads),
                out_edge_feats=edges_feats,
                num_heads=num_heads,bias=True)
        self.conv2 = dglnnpy.conv.EGATConv(in_node_feats=hid_feats,
                in_edge_feats=edges_feats,
                out_node_feats=int(hid_feats/num_heads),
                out_edge_feats=edges_feats,
                num_heads=num_heads)
        self.conv3 = dglnnpy.conv.EGATConv(in_node_feats=hid_feats,
                in_edge_feats=edges_feats,
                out_node_feats=int(hid_feats/num_heads),
                out_edge_feats=edges_feats,
                num_heads=num_heads)
        self.lin1 = nn.Linear(hid_feats, hid_feats)
        self.lin2 = nn.Linear(hid_feats, out_feats)

    def forward(self, graph, inputs, ew):
        H = F.relu(self.conv1(graph, inputs, ew)[0])
        H = torch.reshape(H,(H.shape[0], -1))
        H = F.relu(self.conv2(graph, H, ew)[0])
        H = torch.reshape(H,(H.shape[0], -1))
        H = F.relu(self.conv3(graph, H, ew)[0])
        H = torch.reshape(H,(H.shape[0], -1))
        H = F.relu(self.lin1(H))
        H = self.lin2(H).squeeze()
        return H

loaded_graph = False
first_enc = True
gId = rd.randint(10e7,10e8)

def dglTSPGraph(edges,k):
    
    global loaded_graph
    global first_enc
    file_path = f'dglGraph{gId}.pkl'
    
    
    if not loaded_graph:  
        
        print("Graph init",flush=True)
        E1 = [edges[4*i] for i in range(len(edges)//4)] 
        E2 = [edges[4*i+1] for i in range(len(edges)//4)]
        G = dgl.graph((E1, E2))
        w = [edges[4*i+2] for i in range(len(edges)//4)]

        
        
        #node feature
        n_feat = torch.ones(len(G.nodes()), 32)
        n_feat[0][0] = 0.0
        mean_dist = [0 for i in G.nodes()]
        min_dist = [1000000 for i in G.nodes()]
        max_dist = 0
        for i,e1 in enumerate(G.edges()[0]):
            e2 = G.edges()[1][i]
            dist = w[i]
            if dist > max_dist:
                max_dist = dist
            if dist < min_dist[e1]:
                min_dist[e1] = dist
            if dist < min_dist[e2]:
                min_dist[e2] = dist
            mean_dist[e1] += dist
            mean_dist[e2] += dist
            
        #edge feature
        e_feat = torch.ones(len(G.edges()[0]), 32)
        for i,e1 in enumerate(G.edges()[0]):
            e2 = G.edges()[1][i]
            e_feat[i][0] = w[i] /max_dist
            e_feat[i][1] = 0 #isPresent
            e_feat[i][2] = 0 #isForced
    
        G.edata['weight'] = e_feat
        
        for i in range(len(mean_dist)):
            mean_dist[i] = mean_dist[i]/((len(G.edges()[0])-1)*max_dist)
        for i in range(len(mean_dist)):
            min_dist[i] /= max_dist
        deg = [0 for i in G.nodes()]#deg = [0 for i in Gnx.nodes()]
        Gnx_ = dgl.to_networkx(G).to_undirected()#, edge_attrs=['weight']).to_undirected()
        Gnx = nx.Graph(Gnx_)

        for i,e1 in enumerate(G.edges()[0]):
            e2 = G.edges()[1][i]
            Gnx[e1.item()][e2.item()]['weight'] = G.edata['weight'][i][0].item()

        tree, tw = onetree.one_tree(Gnx)
        for e in tree.edges():
            deg[e[0]] += 1
            deg[e[1]] += 1

        for i,_ in enumerate(n_feat):
            idx, rep = 2,9
            for j in range(rep):
                n_feat[i][idx] = mean_dist[i]
                n_feat[i][idx+1] = min_dist[i]
                n_feat[i][idx+2] = deg[i]/len(G.nodes())
                idx += 3
        
        
        G.ndata['x'] = n_feat
        G = dgl.add_reverse_edges(G,copy_edata=True)
        
        loaded_graph = True
        with open(file_path, 'wb') as file:
            pickle.dump(G, file)
        #with open(f'subgraphsk10n{size_graphs}/dglGraph{gId}.pkl', 'wb') as file:
        #    pickle.dump(G, file)
        

    with open(file_path, 'rb') as file:
        G = pickle.load(file)
    
    
    savesubgraphs = False
    
    ledges = len(edges)
    isForced = [edges[4*i+3] for i in range(ledges//4)]
    if (not savesubgraphs) or first_enc:
        for i in range(ledges//4):
            for _id in [G.edge_ids(edges[4*i],edges[4*i+1]),G.edge_ids(edges[4*i+1],edges[4*i])]:
                dist_ = G.edata['weight'][_id][0]
                for i in range(5):
                    G.edata['weight'][_id][i*4] = dist_
                    G.edata['weight'][_id][1+i*4] = 1.0
                    G.edata['weight'][_id][2+i*4] = isForced[i]

    
    #if savesubgraphs:
        # if first_enc:
        #     with open(f'subgraphsk10/dglGraph{gId}.pkl', 'wb') as file:
        #         pickle.dump(G, file)
        #         first_enc = False
            #with open(f'subgraphsk10n{size_graphs}/dglGraph{gId}_{rd.randint(10e7,10e8)}.pkl', 'wb') as file:
                #pickle.dump(edges, file)
    
    return G


n_graphs = 100
max_coord = 100
from onetree import model
gnnmodel = model(32,32,1)
opt = torch.optim.Adam(gnnmodel.parameters(),lr=0.01,maximize=True)

#size_graphs = 100
#k = 10
#typeg = "sym_geo"
#gnnmodel.load_state_dict(torch.load(f"./models/ULmodeln{size_graphs}_{typeg}_k{k}.pt"))
gnnmodel.load_state_dict(torch.load(f"./model.pt"))

gnnmodel.eval()


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
    G2 = G1.copy() 
    G2.remove_node(nodes_rem)
    T = nx.minimum_spanning_tree(G2,algorithm="prim")
    
    edges = list(G.edges)
    w = []
    for e in edges:
        w.append( G.get_edge_data(*e)['weight']+pen[e[0]]+pen[e[1]] )
    w_edges = [(edges[i][0],edges[i][1],w[i]) for i in range(len(w))]
    w_edges = sorted(w_edges, key=lambda tup: tup[2])
    
    #add two min cost edges adj to {nodes_rem}
    compt = 0
    for we in w_edges:
        if compt < 2 and (we[0] == nodes_rem or we[1] == nodes_rem):
            T.add_edge(we[0],we[1])
            compt += 1
    
    if gnopen != None:
        G1 = gnopen.copy() #graph with modified costs on edges
        for e in G1.edges:
            G1[e[0]][e[1]]['weight'] += pen[e[0]]+pen[e[1]]
        
    tw = 0
    for e in list(T.edges):
        tw += G1.get_edge_data(*e)['weight']
    if remW:
        for p in pen:
            tw -= 2*p
        
    return T,tw
