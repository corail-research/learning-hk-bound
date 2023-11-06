#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: augustinparjadis
"""

import random as rd
import numpy as np
import sys,os,math,pickle,re,copy
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
from HKbound import thetas_C
from onetree import one_tree
from onetree import one_tree_costs_fast
from onetree import costs_derivative
from onetree import stepTheta

#training graphs
gtype = "sym_geo"
size_graphs = 100

def generate_dataset(n_graphs,n_nodes,max_coord,validation_set_p=0.1,generate_graphs=True):
    
    training_set = []
    graph_points = []
    N = n_nodes
    M = n_graphs
    train_dir = "/home/val/Documents/TSP_HeldKarp_Filtering/HK_bound/graphs_training"
    train_graphs_path = ""
    
    if not generate_graphs:
        print('Loading graphs')
        with open(f'{train_dir}/TSPdataset_s{n_graphs}_n{n_nodes}_c{max_coord}.pkl', 'rb') as f:
            training,validation = pickle.load(f)
        if len(training) > 1:
            return training,validation
        else:
            print('loading failed')
            return -1,-1

        for item in os.listdir(train_dir):
            if item.split('_')[0] == f"N{N}":
                print(item)
                c_ = 0
                for _ in os.listdir(train_dir+"/"+item):
                    c_ += 1
                if c_ == M:
                    train_graphs_path = train_dir+"/"+item
                    print(f"Training from {train_graphs_path}")
                    break
        if train_graphs_path != "":
            for file in os.listdir(train_graphs_path):
                id_,L = read_graph_file(train_graphs_path+"/"+file)
                graph_points.append((id_,L))
        else:
            print("Graphs not found")
            return -1

    else:
        #create folder with graph info and unique ID
        train_graphs_path = train_dir+f"/N{N}_{rd.randint(1e7,1e8)}"
        os.makedirs(train_graphs_path)
    
        ids = [rd.randint(1e7,1e8) for i in range(n_graphs)]
        ids = list(dict.fromkeys(ids))
        while(len(ids)) < n_graphs:
            r = rd.randint(1e7,1e8)
            if r not in ids:
                ids.append(r)
        
        for _id in ids:
            L = [(rd.random()*max_coord,rd.random()*max_coord) for i in range(n_nodes)]
            graph_points.append((_id,L))
            
            with open(f"{train_graphs_path}/{_id}.tsp", "w") as f:
                f.write(f"NAME : {_id}\nTYPE : TSP\nDIMENSION : {n_nodes}\nEDGE_WEIGHT_TYPE : EUC_2D\nNODE_COORD_SECTION\n")
                Ltxt = [f'{i+1} {a[0]} {a[1]}\n' for i,a in enumerate(L)]+["EOF"]
                f.writelines(Ltxt)
    
   
    for _id,L in graph_points:
        
        UB,UB_tour = ortoolstsp(L)
        UB = int(1+UB)
        G = nx.Graph()
        G.add_nodes_from([i for i in range(n_nodes)])
        for i in range(n_nodes):
            for j in range(i+1,n_nodes):
                G.add_edge(i,j, weight=np.linalg.norm(np.array(L[j])-np.array(L[i])))
        Gnx = G
            

        #get thetas
        LB_,true_t = thetas_C(f"{train_graphs_path}/{_id}.tsp",N,UB+1)
        if LB_ <= 1:
            continue
        onetree, LB = one_tree(G,theta=torch.tensor(true_t))
        t = true_t.copy() #normalized thetas
        print(f'\n_id {_id}.tsp')
        print(f'\n--LB_ {LB_} LB {LB} UB {UB}')
        if LB_ > UB:
            UB = LB_
        RIe = one_tree_costs_fast(G)
        
        
        print(f'LB {LB} UB {UB}')
        max_theta = max(max(true_t),abs(min(true_t)))
        if LB <= 1:
            continue
            for _t in range(len(t)):
                t[_t] = 0.0
        else:
            for _t in range(len(t)):
                t[_t] /= max_theta
                    
        #create dgl graph G   
        G = dgl.graph(([],[]))
        for i in range(N):
            G.add_edges(i, [j for j in range(i+1,N)])
        G = dgl.add_reverse_edges(G)
        
        #node feature
        n_feat = torch.ones(len(G.nodes()), 32)
        n_feat[0][0] = 0.0 #indicate node 1
        mean_dist = [0 for i in G.nodes()]
        min_dist = [10000 for i in G.nodes()]
        max_dist = 0
        for i,e1 in enumerate(G.edges()[0]):
            e2 = G.edges()[1][i]
            dist = np.linalg.norm( np.array(L[e1])-np.array(L[e2]) )
            if dist > max_dist:
                max_dist = dist
            if dist < min_dist[e1]:
                min_dist[e1] = dist
            if dist < min_dist[e2]:
                min_dist[e2] = dist
            mean_dist[e1] += dist
            mean_dist[e2] += dist
        for i in range(len(mean_dist)):
            mean_dist[i] = mean_dist[i]/((len(G.edges()[0])-1)*max_dist)
        for i in range(len(mean_dist)):
            min_dist[i] /= max_dist
        tree, tw = one_tree(Gnx)
        deg = [0 for i in Gnx.nodes()]
        for e in tree.edges():
            deg[e[0]] += 1
            deg[e[1]] += 1
        for i,_ in enumerate(n_feat):
            idx, rep = 2,9
            for j in range(rep):
                n_feat[i][idx] = mean_dist[i]
                n_feat[i][idx+1] = min_dist[i]
                n_feat[i][idx+2] = deg[i]/len(Gnx.nodes())
                idx += 3
            
        #edge feature
        e_feat = torch.ones(len(G.edges()[0]), 1)
        for i,e1 in enumerate(G.edges()[0]):
            e2 = G.edges()[1][i]
            e_feat[i][0] = np.linalg.norm( np.array(L[e1])-np.array(L[e2]) ) /max_dist

        G.edata['weight'] = e_feat
        G.ndata['x'] = n_feat
        
        if LB > 1:
            t = torch.tensor(t)
            training_set.append((G,Gnx,t,max_theta,(LB,UB,LB_,UB_tour),true_t,RIe))
        print((G,t))
    
    rd.shuffle(training_set)
    training = training_set[:int((1-validation_set_p)*len(training_set))]
    validation = training_set[int((1-validation_set_p)*len(training_set)):]
    with open(f'{train_dir}/TSPdataset_s{n_graphs}_n{n_nodes}_c{max_coord}.pkl', 'wb') as f:
        pickle.dump((training,validation), f)
    return training,validation
  

def read_graph_file(file):
    L = []
    id_ = file.split('.')[0].split("/")[-1]
    with open(file,'r') as f:
        lines = f.readlines()
        for l in lines:
            l = l.split()
            if l[0].isnumeric():
                L.append((float(l[-2]),float(l[-1])))
    return id_,L


class oneTreeLoss(Function):
    @staticmethod
    def forward(G,logits,max_theta):
        nnt, nnw = one_tree(G,theta=(logits*max_theta).tolist())
        return -nnw

    @staticmethod
    def backward(g,th_true,th_est,mt):
        return compute_f_gradient(g,th_true,th_est,mt),None
    
class degreeLoss(Function):
    @staticmethod
    def forward(G,logits,max_theta):
        nnt, nnw = one_tree(G,theta=(logits*max_theta).tolist())
        return -nnw
    @staticmethod
    def backward(g,th_true,th_est,mt):
        return compute_deg_gradient(g,th_true,th_est,mt),None
      
 
def compute_deg_gradient(g,th_true,th_est,mt):
    
    C = 1.0
    grad_f = np.zeros(len(th_true)) #vector size theta
    
    t, tw = one_tree(g,theta=th_est)
    deg = [0 for i in th_est]
    for e in t.edges():
        deg[e[0]] += 1
        deg[e[1]] += 1
    
    for i in g.nodes():
        grad_f[i] = -C*(2-deg[i])
    print(f"grad_f {grad_f}")
    for i in g.nodes():
        grad_f[i] -= 0.0001*th_est[i]
    
    grad_f = torch.from_numpy(grad_f)
    print(f"grad_f' {grad_f}")
    return grad_f, None

grad_clip = 10.0
def compute_f_gradient(G,th_true,th_est,mt):
    
    grad_f = np.zeros(len(th_true)) #vector size theta

    edata_true,edata_est = costs_derivative(G,th_true,th_est)
    R,I,Rhat,Ihat = {},{},{},{}
    for k in edata_true:
        R[k] = edata_true[k][0]
        I[k] = edata_true[k][1]
    for k in edata_est:
        Rhat[k] = edata_est[k][0]
        Ihat[k] = edata_est[k][1]
    for i in G.nodes():
        
        ep,epval,epp,eppval = None,M,None,M
        for j in G.nodes():
            if i != j:
                e = (i,j)
                
                if e not in Rhat:
                    e = e[::-1]
                if Rhat[e] > 0 and Rhat[e] < eppval:
                    eppval = Rhat[e]
                    epp = e
                    
                if e not in Ihat:
                    e = e[::-1]
                if Ihat[e] > 0 and Ihat[e] < epval:
                    epval = Ihat[e]
                    ep = e
                    

        if ep == None or epp == None:
            if ep == None:
                grad_f[i] = grad_clip
                print('No edge to insert')
            elif epp == None:
                grad_f[i] = -grad_clip
                print('No edge to replace, ERROR')
                sys.exit()
        else:
            grad_f[i] = 2*(I[ep]-R[epp])/(Ihat[ep]+Rhat[epp])
    print(f"grad_f {grad_f}")
    
    #degree
    t, tw = one_tree(g,theta=th_est)
    deg = [0 for i in th_true]
    for e in t.edges():
        deg[e[0]] += 1
        deg[e[1]] += 1
    
    for i in g.nodes():
        grad_f[i] -= -0.1*(2-deg[i])
    print(f"grad_f' {grad_f}")
    
    for i in G.nodes():
            if math.isnan(grad_f[i]):
                grad_f[i] = 0.0
            elif grad_f[i] > grad_clip:
                grad_f[i] = grad_clip
            elif grad_f[i] < -grad_clip:
                grad_f[i] = -grad_clip
    for i in range(len(th_true)):
        grad_f[i] -= 0.1*th_est[i] #prevent growth
        grad_f[i] -= 0.01 #remove sum thetas
        
    grad_f = torch.from_numpy(grad_f)
    print(f"grad_f'' {grad_f}")
    return grad_f, None
    
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


TS,TS_ = [],[]
k = 11
folder_path = f'./training_graphs/subgraphs_n{size_graphs}_{gtype}'
M = 1000000
ptoload = 1
print('loading dataset...',flush=True)
for file_name in os.listdir(folder_path): #load primary graphs
    if file_name.endswith('.pkl'):
        if '_' not in file_name:
            #print(f'loading {file_name}',flush=True)
            file_path = os.path.join(folder_path, file_name)
            numbers = re.findall(r'\d+', file_name)
            gid = int(numbers[0])
            with open(file_path, 'rb') as file:
                pklG = pickle.load(file) 
            G = nx.Graph()
            G1 = nx.Graph()
            n_nodes = pklG.num_nodes()
            G.add_nodes_from([i for i in range(n_nodes)])
            G1.add_nodes_from([i for i in range(n_nodes)])
            edges = [(i,j) for i in range(n_nodes) for j in range(i+1,n_nodes)]
            G.add_edges_from(edges)
            G1.add_edges_from(edges)
            for i in range(n_nodes):
                for j in range(i+1,n_nodes):
                    G[i][j]['weight'] = pklG.edata['weight'][pklG.edge_ids(i,j)][0].item()
                    G1[i][j]['weight'] = pklG.edata['weight'][pklG.edge_ids(i,j)][0].item()
            Gnx = G
            TS_.append((pklG,Gnx,G1,gid))

for file_name in os.listdir(folder_path): #load subgraphs
    if file_name.endswith('.pkl'):
        if '_' in file_name:
            file_path = os.path.join(folder_path, file_name) 
            numbers = re.findall(r'\d+', file_name)
            gid = int(numbers[0])
            kg = int(numbers[-1])
            if kg <= k:
                if rd.random() < ptoload:
                    #print(f'loading {file_name}',flush=True)
                    with open(file_path, 'rb') as file:
                        edges = pickle.load(file) 
                    for ts in TS_:
                        if ts[-1] == gid:
                            pklG = copy.deepcopy(ts[0])
                            G = copy.deepcopy(ts[1])
                            G1 = copy.deepcopy(ts[2])
                            
                            
                    E1 = [edges[4*i] for i in range(len(edges)//4)] 
                    E2 = [edges[4*i+1] for i in range(len(edges)//4)]
                    isForced = [edges[4*i+3] for i in range(len(edges)//4)]
                    for i in range(len(edges)//4):
                        for _id in [pklG.edge_ids(E1[i],E2[i]),pklG.edge_ids(E2[i],E1[i])]:
                            dist_ = pklG.edata['weight'][_id][0]
                            for i in range(5):
                                pklG.edata['weight'][_id][i*4] = dist_
                                pklG.edata['weight'][_id][1+i*4] = 1.0
                                pklG.edata['weight'][_id][2+i*4] = isForced[i]
                    for i,e1 in enumerate(pklG.edges()[0]):
                        e2 = pklG.edges()[1][i].item()
                        e1 = e1.item()
                        offs = 0
                        if not pklG.edata['weight'][i][1].item():#isPresent
                            G[e1][e2]['weight'] += M
                        elif pklG.edata['weight'][i][2].item():#isForced
                            G[e1][e2]['weight'] -= M
                    Gnx = G
                    TS.append((pklG,Gnx,G1,gid))

print(f'loaded {len(TS)} graphs')
rd.shuffle(TS)
validation_set_p = 0.1
training_set = TS[:int((1-validation_set_p)*len(TS))]
validation_set = TS[int((1-validation_set_p)*len(TS)):]

lr = 0.002
dfW,dfWval = [],[]
from onetree import model
gnnmodel = model(32,32,1)
opt = torch.optim.Adam(gnnmodel.parameters(),lr=lr,maximize=True,weight_decay=1e-5)

degloss = degreeLoss()

batch_s, itr = 5,200
for i in range(itr+1):
    torch.save(gnnmodel.state_dict(), f"./trained_models/ULmodeln{size_graphs}_{gtype}_k{k-1}.pt")
    dfw, dfwval = 0,0
    rd.shuffle(training_set)
    
        
    for gtrain in validation_set:

        graph,Gnx,G1,_ = gtrain
        
        feats = graph.ndata['x']
        e_weights = graph.edata['weight']
        

        logits = gnnmodel(graph, feats, e_weights)   
        th_DF = (logits).tolist()
        g = Gnx        

        print(f'--- {i+1}/{itr} validation')
        print(logits)

        if False:
            pass
        else:
            dft, dfw_ = one_tree(g,theta=th_DF,gnopen=G1)
            dfwval += dfw_
            print(f'1Tree DF thetas {dfw_}')
    
    batch_grad = []
    for gtrain in training_set:

        graph,Gnx,G1,_ = gtrain
        
        feats = graph.ndata['x']
        e_weights = graph.edata['weight']
        logits = gnnmodel(graph, feats, e_weights)
        
        print(f'--- {i+1}/{itr}')
        
        th_DF = (logits).tolist()
            
        g = Gnx
        opt.zero_grad()
        grad, _ = degloss.backward(g,th_DF,th_DF,1)
        
        b_grad = False
        if b_grad:
            batch_grad.append(grad[0])
            if len(batch_grad) == batch_s:
                grad = torch.zeros(len(batch_grad[0]))
                for grad_ in batch_grad:
                    grad += grad_
                logits.backward(grad)
                opt.step()
                batch_grad = []
        else:
            logits.backward(grad)
            opt.step()

        

        if False:
            pass
        else:
            dft, dfw_ = one_tree(g,theta=th_DF,gnopen=G1)
            dfw += dfw_
            print(f'1Tree DF thetas {dfw_}')

    dfW.append(dfw/len(training_set))
    dfWval.append(dfwval/len(validation_set))
    print(f"training bound dfW {dfW}")
    print(f"validation bound dfWval {dfWval}")


torch.save(gnnmodel.state_dict(), f"./trained_models/ULmodeln{size_graphs}_{gtype}_k{k-1}.pt")
sys.exit()
