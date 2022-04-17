import collections
# import imp
from platform import node
import random
import copy
from re import I
# import re
# from subprocess import list2cmdline
from time import time
from tkinter.messagebox import NO
# import math
# from turtle import pos
import numpy as np
import networkx as nx
from matplotlib import pyplot as plt
# from scipy import rand
import warnings
import pandas as pd

from allocation import Node
from allocation import Net
from allocation import PopGene
from allocation import Initpre

nodelist = [['g8', 'gp2', 'g10', 'g11', 'g13'], ['g7', 'gp3', 'g10', 'g11'],
            ['g6', 'gp4', 'g11'], ['g14', 'gp5', 'g1', 'g18'], ['g5', 'gp6'],
            ['g5', 'gp7'], ['g17', 'gp8', 'g19', 'g20', 'g22'],
            ['g16', 'gp9', 'g19', 'g20'], ['g15', 'gp10','g20'], ['g0', 'gp11'],
            ['g0', 'gp12'], ['g0', 'gp13'], ['g2', 'gp14'], ['g2', 'gp15'],
            ['g2', 'gp16'], ['g25', 'gp17', 'g27', 'g28', 'g30'],
            ['g24', 'gp18', 'g27', 'g28'], ['g23', 'gp19', 'g28'],
            ['gp20', 'g10', 'g5', 'g19', 'g0', 'g2', 'g27', 'g11', 'g12', 'g20', 'g21', 'g28', 'g29'],
            ['gp21', 'g10', 'g5', 'g19', 'g0', 'g2', 'g27', 'g11', 'g12', 'g20', 'g21', 'g28', 'g29'],
            ['gp0', 'g4', 'g5', 'g1', 'g0', 'g2', 'g3'],
            ['gp1', 'g5', 'g0', 'g2', 'g3', 'g26'], ['g4', 'g6', 'g7', 'g8'],
            ['g11', 'g6'], ['g9', 'g6', 'g7', 'g8'], ['g12', 'g6', 'g7', 'g8'],
            ['g10', 'g7'], ['g13', 'g8'], ['g1', 'g15', 'g16', 'g17'],
            ['g20', 'g15'], ['g18', 'g15', 'g16', 'g17'],
            ['g21', 'g15', 'g16', 'g17'], ['g19', 'g16'], ['g22', 'g17'],
            ['g3', 'g23', 'g24', 'g25'], ['g28', 'g23'],
            ['g26', 'g23', 'g24', 'g25'], ['g29', 'g23', 'g24', 'g25'],
            ['g27', 'g24'], ['g30', 'g25']]
quadlist = [{'g3': 53, 'g4': 52, 'g12': 4, 'g13': 3, 'g17': 6, 'gp0': 1,
             'gp5': 1, 'gp6': 1, 'gp7': 1, 'gp12': 1, 'gp17': 1, 'gp18': 1},
            {'g2': 65, 'g6': 6, 'g7': 6, 'g9': 2, 'g10': 6, 'g11': 7, 'g18': 3,
             'g22': 3, 'g27': 6, 'g28': 7, 'g29': 4, 'g30': 3, 'gp2': 1, 'gp3': 1,
             'gp8': 1, 'gp15': 1, 'gp16': 1, 'gp21': 1},
            {'g0': 65, 'g1': 53, 'g26': 3, 'gp9': 1, 'gp10': 1, 'gp14': 1, 'gp19': 1, 'gp20': 1},
            {'g5': 64, 'g8': 6, 'g14': 1, 'g15': 6, 'g16': 6, 'g19': 6, 'g20': 7,
             'g21': 4, 'g23': 6, 'g24': 6, 'g25': 6, 'gp1': 1, 'gp4': 1, 'gp11': 1, 'gp13': 1}]
exp_list = [{'name':1}, {'name':2}]
gene_exp = {'name':'g0', 'weight':1, 'coord':[0.1, 0.2]}
expdict = {
    'g0': 65,
    'g1': 53,
    'g2': 65,
    'g3': 53,
    'g4': 52,
    'g5': 64,
    'g6': 6,
    'g7': 6,
    'g8': 6,
    'g9': 2,
    'g10': 6,
    'g11': 7,
    'g12': 4,
    'g13': 3,
    'g14': 1,
    'g15': 6,
    'g16': 6,
    'g17': 6,
    'g18': 3,
    'g19': 6,
    'g20': 7,
    'g21': 4,
    'g22': 3,
    'g23': 6,
    'g24': 6,
    'g25': 6,
    'g26': 3,
    'g27': 6,
    'g28': 7,
    'g29': 4,
    'g30': 3,
    'gp0': 1,
    'gp1': 1,
    'gp2': 1,
    'gp3': 1,
    'gp4': 1,
    'gp5': 1,
    'gp6': 1,
    'gp7': 1,
    'gp8': 1,
    'gp9': 1,
    'gp10': 1,
    'gp11': 1,
    'gp12': 1,
    'gp13': 1,
    'gp14': 1,
    'gp15': 1,
    'gp16': 1,
    'gp17': 1,
    'gp18': 1,
    'gp19': 1,
    'gp20': 1,
    'gp21': 1
}
net_exp = []
set_exp = []





"""
# 每个列表第一个元素为驱动节点 该节点连接该线网内所有的驱动节点
# 先把所有驱动节点放到一个列表
# 再建立一个驱动节点与其所有的负载节点的关系
# 然后检查一个象限内的驱动节点
# 若该驱动节点的某个负载节点不在该象限 则连接线值+1
"""

"""
修改为
不再考虑节点类型 只考虑线网整体
一个线网内的若干个节点在某一象限 则连线值+1
建立象限集合与线网集合 若可取到交集 则连线值+1
"""

# 节点的具名元组 '节点名 分配区域 资源权重 坐标'
Vertex = collections.namedtuple('Vertex', 'node fpganum weight coordinates')

def read_net(filepath):
    
    netdata = []
    with open(filepath, "r") as f:
        for line in f.readlines():
            # 去掉列表中每一个元素的换行符
            netdata.append(line.strip('\n'))

    # 返回驱动节点的索引
    snode = []
    for i in range(len(netdata)):

        if netdata[i].count(' ') == 2:
            snode.append(i)
    print(snode)

    nodegroup = []
    for i in range(len(snode)):
        
        # i + 1 <= len(snode) - 1
        if i <= len(snode) - 2:
            # 同一个线网的节点
            nodegroup.append(netdata[snode[i]:snode[i + 1]])
        else:
            # 最后一个线网
            nodegroup.append(netdata[snode[i]:])
    print(nodegroup)
    """
    nodegroup[]
    [['g8 s 1', 'gp2 l', 'g10 l', 'g11 l', 'g13 l']
     ['g7 s 1', 'gp3 l', 'g10 l', 'g11 l']
     ['g6 s 1', 'gp4 l', 'g11 l']]
    """
    return nodegroup




# # 测试语句 # #
# netpath = "data\design.net"
# read_net(netpath)
"""
输入：线网列表 节点分配的字典列表
输出：连线值列表 
"""
def linkCal(nodelist, quadlist):
    
    nodeset = []  # 线网集合的list
    for i in nodelist:
        j = set(i)
        nodeset.append(j)
    # print(nodeset)

    quadset = []  # 象限集合的list
    for i in quadlist:
        j = set(list(i.keys()))
        quadset.append(j)

    netlink = []
    for i in quadset:  # 循环象限
        
        link = 0
        for j in nodeset:  # 循环所有线网
            if len(i & j) != 0:  # 有交集 线连接数+1
                link += 1
        netlink.append(link)
    # print(netlink)
    return netlink
        
