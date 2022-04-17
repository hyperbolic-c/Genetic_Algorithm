import networkx as nx
from matplotlib import pyplot as plt

g = nx.Graph()
nodetext = [['g8 s 1', 'gp2 l', 'g10 l', 'g11 l', 'g13 l'],
            ['g7 s 1', 'gp3 l', 'g10 l', 'g11 l'],
            ['g6 s 1', 'gp4 l', 'g11 l']]

netnodes = []
for nodei in nodetext:  # 所有子图
    netnode = []  # 只含有节点的列表
    for nodej in nodei:  # 一个子图
        if nodej.find('s') > 0:
            netnode.append(nodej[:nodej.index(' ')])
        else:
            netnode.append(nodej[:nodej.index(' ')])
    print(netnode)

    g.add_node(netnode[0])
    for i in range(1, len(netnode)):
        g.add_node(netnode[i])
        g.add_edge(netnode[0], netnode[i])

    netnodes.append(netnode)

print(netnodes)
"""
[['g8', 'gp2', 'g10', 'g11', 'g13'],
 ['g7', 'gp3', 'g10', 'g11'],
 ['g6', 'gp4', 'g11']]
"""
nx.draw_networkx(g)
plt.show()
