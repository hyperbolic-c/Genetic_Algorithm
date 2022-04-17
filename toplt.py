from matplotlib import pyplot as plt
import numpy as np
import networkx as nx

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

pos_copy = {
    'g0': [-1., 0.23551616],
    'g1': [-0.76203491, -0.03619811],
    'g2': [-0.07089636, -0.27689165],
    'g3': [-0.17387499, 0.11831481],
    'g4': [0.48751594, -0.4377812],
    'g5': [0.96534748, 0.02099123],
    'g6': [0.55394284, 0.37604875]
}
# print(pos_copy)
poskeys = pos_copy.keys()
# print(poskeys)
for i in poskeys:
    pos_copy[i] = np.array(pos_copy[i])
# print(pos_copy)

G = nx.Graph()
G.add_nodes_from(list(poskeys))
nx.draw(G, pos_copy, with_labels=True)

plt.show()


"""
确认节点位于哪一象限
返回各个象限的节点-权重字典
"""
def posCal(posdict):
    Z1_dict = {}
    Z2_dict = {}
    Z3_dict = {}
    Z4_dict = {}
    for key in posdict.keys():

        val = posdict[key]
        x = val[0]
        y = val[1]

        if x>0 and y>0 :
            Z1_dict[key] = expdict[key]
            print("1")
        elif x>0 and y<0 :
            Z4_dict[key] = expdict[key]
            print("4")
        elif x<0 and y>0 :
            Z2_dict[key] = expdict[key]
            print("2")
        else :
            Z3_dict[key] = expdict[key]
            print("3")
    print(Z1_dict, Z2_dict, Z3_dict, Z4_dict)

    return Z1_dict, Z2_dict, Z3_dict, Z4_dict


"""
输入节点-权值字典
返回字典中权值总和
返回目标节点用作输出信息
"""
def weightCal(targdict):

    targkey = list(targdict.keys())
    weightsum = 0
    for i in targkey:
        val = targdict[i]
        weightsum = weightsum + val
    return targkey, weightsum
