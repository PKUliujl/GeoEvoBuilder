# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 13:51:21 2020

@author: Administrator
"""



import torch
from torch_geometric.nn import MessagePassing
from torch_geometric.utils import add_self_loops,degree

class GCNconv(MessagePassing):
    def __init__(self,in_channels, out_channels,edge_length=None):
        super(GCNconv, self).__init__(aggr='add')   ## 'add' aggregation
        self.in_channels = in_channels
        self.out_channels = out_channels
        self.lin1 = torch.nn.Linear(out_channels,out_channels)
        self.BN1 = torch.nn.BatchNorm1d(out_channels)
        self.lin2 = torch.nn.Linear(12,out_channels**2)
        self.BN2 = torch.nn.BatchNorm1d(out_channels)
        self.act = torch.nn.ELU()
        self.linCAT = torch.nn.Linear( in_channels+out_channels,out_channels)
        self.conv1 = torch.nn.Sequential(
                torch.nn.Conv1d(in_channels=in_channels+out_channels,out_channels=out_channels,kernel_size=3,padding=1),
                torch.nn.BatchNorm1d(out_channels),
                torch.nn.ELU(),
                torch.nn.Conv1d(in_channels=out_channels,out_channels=out_channels,kernel_size=3,padding=1),
                #torch.nn.BatchNorm1d(out_channels),
                #torch.nn.ELU(),
                )
        if edge_length:
            self.mlp = torch.nn.Sequential(
                torch.nn.Linear( edge_length + in_channels,out_channels),
                torch.nn.BatchNorm1d(out_channels),
                torch.nn.ELU(),
                torch.nn.Linear(out_channels,in_channels * out_channels),
                )
        else:
            self.mlp = torch.nn.Sequential(
                torch.nn.Linear(out_channels + in_channels,out_channels),
                torch.nn.BatchNorm1d(out_channels),
                torch.nn.ELU(),
                torch.nn.Linear(out_channels,in_channels * out_channels),
                torch.nn.BatchNorm1d(in_channels*out_channels),
                torch.nn.ELU(),
                )
        self.mlp2 = torch.nn.Sequential(
                torch.nn.Linear(in_channels * out_channels,out_channels),
                torch.nn.BatchNorm1d(out_channels),
                torch.nn.ELU(),
                )
        
    def forward(self,x,edge_index,edge_attr):
        
        # x has shape [N, in_channels]
        # edge_index has shape [2, E]
        
        # Step 1: Add self-loops to the adjacency matrix.
        #edge_index , _ =add_self_loops( edge_index,num_nodes=x.size(0) )

        x_j = x[edge_index[1]]
        x_i = x[edge_index[0]]
        edge_attr = torch.cat((x_j-x_i,edge_attr),1)
        #print(x_j-x_i)
        edge_attr = self.mlp(edge_attr)

        out = self.propagate(edge_index, x=x, edge_attr = edge_attr)
        out = self.BN2(out)
        out = self.act(out)
        #step 2: Linearly transform node feature matrix
        out = torch.cat((out,x),1)
        out = self.linCAT(out) 
        #print('forward:', edge_attr.size())
        #step 3-5: Start propagating messages.
        #return self.propagate(edge_index, x=x, edge_attr = edge_attr) + x
        edge_attr = self.mlp2(edge_attr)
        return out,edge_attr
 

    
    def message(self,x_j,edge_index,x,edge_attr):
        row, col = edge_index
        deg = degree(row, x.size(0),dtype=x_j.dtype )
        deg_inv_sqrt = deg.pow(-0.5)
        norm = deg_inv_sqrt[row] * deg_inv_sqrt[col]
        '''
        edge_attr = torch.cat((x_j-x_i,edge_attr),1)
        #print(x_j-x_i)
        edge_attr = self.mlp(edge_attr)
        '''
        weight = edge_attr.view(-1,self.in_channels,self.out_channels)
        return norm.view(-1,1) * (torch.matmul(x_j.unsqueeze(1), weight).squeeze(1))
    
    def update(self,aggr_out):
        # aggr_out has shape [N, out_channels]
        
        # step 5: Return new node embeddings.
        return aggr_out

