"""
集成电路计算机辅助
实验二
2022.3.28 21:00
"""


import random
import copy
from time import time
import numpy as np
import networkx as nx
from matplotlib import pyplot as plt
import warnings
import datetime
import sys

warnings.filterwarnings('ignore')


def weightCal(targlist):
    """
    计算各象限的资源权重和
    输入节点-权值字典列表
    返回字典中权值总和列表
    """
    weightlist = []
    for targdict in targlist:  # 遍历列表
        targkey = list(targdict.keys())
        weightsum = 0
        for i in targkey:  # 遍历其中一个象限
            val = targdict[i]
            weightsum = weightsum + val
        weightlist.append(weightsum)
    return weightlist


def read_are(filepath):
    """
    读取are文件 返回字典
    键为节点名 值为节点的资源和
    {'g0': 65, 'g1': 53}
    """
    aredata = []
    aredict = {}
    aredata = np.genfromtxt(filepath, dtype=str)

    for i in aredata:
        c = np.array(i[1:].astype(np.int).tolist())
        are_not_zero = c.ravel()[np.flatnonzero(c)]  # 筛选非零
        aredict[i[0]] = sum(are_not_zero.tolist())
    return aredict


def initPop(N, nc):
    """
    产生N个染色体的初始群体，保存在pop
    将染色体编码为范围为(-1, 1)的浮点数 也就是53个节点的位置pos信息
    前nc个随机数为x 后nc个随机数为y
    编码为节点会方便一些
    引入坐标为计算布局 连线做准备
    """
    pop = -1 + 2 * np.random.random((N, nc * 2))
    return pop


def posCal(posdict, aredict):
    """
    确认节点位于哪一象限
    aredict为所有节点的资源权重字典
    输入节点-位置字典 返回各个象限的节点-权重字典
    """
    Z1_dict = {}
    Z2_dict = {}
    Z3_dict = {}
    Z4_dict = {}

    for key in posdict.keys():
        val = posdict[key]
        x = val[0]
        y = val[1]

        if x > 0 and y > 0:
            Z1_dict[key] = aredict[key]
        elif x > 0 and y < 0:
            Z4_dict[key] = aredict[key]
        elif x < 0 and y > 0:
            Z2_dict[key] = aredict[key]
        else:
            Z3_dict[key] = aredict[key]
    return [Z1_dict, Z2_dict, Z3_dict, Z4_dict]


def varCal(pop):
    """
    计算方差 均值
    """
    N, nc = np.shape(pop)
    lstvar = np.zeros([N, 1])
    lstmean = np.zeros([N, 1])
    lstweight = []

    for i in range(N):
        aredict = read_are(arepath)
        gene = pop[i, :]  # 取一个染色体
        pos = geneDict(gene, aredict)  # 将gene中的坐标绑定节点
        quadict = posCal(pos, aredict)  # 判定节点位于哪个区域
        weightlist = weightCal(quadict)  # 计算权重
        lstvar[i] = np.var(weightlist)  # 计算方差
        lstmean[i] = np.mean(weightlist) # 计算平均值
        lstweight.append(weightlist)
    return lstvar, lstmean, lstweight


def geneDict(gene, aredict):
    """
    将gene坐标与节点绑定
    前nc个gene为x 后nc个gene为y
    """
    pos = {}
    nc = int(len(gene) / 2)
    targkey = list(aredict.keys())
    i = 0
    for key in targkey:
        pos[key] = [gene[i], gene[i + nc]]
        i = i + 1
    return pos


def fitnessCal(pop, mode):
    """
    根据方差/连线和 计算各个染色体的适应值fitness
    方差lstvar = np.zeros([N, 1])
    连线和lstsum = np.zeros([N, 1])
    
    3个模式 mode == 1方差最小 mode == 2连线和最小 mode == 3满足条件使连线和最小
    传入pop
    再调用函数计算var sum
    """
    lstvar, lstmean, listweight = varCal(pop)
    lstsum = linkCal(pop)

    # 方差最小
    if(mode == 1):
        fitness = 1.0 / lstvar
    
    # 连线和最小
    elif(mode == 2):
        fitness = 1.0 / lstsum

    elif(mode == 3):
        # fitlink = 1.0 / lstsum
        """
        0.9M<Zi<1.1M
        abs(Zi)<0.2M
        """

        fitlink = np.zeros(lstvar.shape)
        for i in range(len(listweight)):  # 循环整个群体
            zmax = 0
            for z in listweight[i]:  # 循环一个gene
                zabs = abs(lstmean[i] - z)
                if zabs > zmax:
                    zmax = zabs
            if(zmax < lstmean[i] * 0.2):
                fitlink[i] = 1.0 / lstsum[i]
            else:
                fitlink[i] = (1.0 / lstsum[i]) / zmax
        fitness = fitlink
    else:
        print('error in fitnessCal')
    return fitness


def linkCal(pop):
    """
    计算染色体群连线和
    """
    N, nc = np.shape(pop)
    lstsum = np.zeros([N, 1])

    for i in range(N):
        aredict = read_are(arepath)
        gene = pop[i, :]  # 取一个染色体
        pos = geneDict(gene, aredict)  # 将gene中的坐标绑定节点
        quadict = posCal(pos, aredict)  # 判定节点位于哪个区域
        nodelist = read_net(netpath)  # 所有线网
        netlink = linkSum(nodelist, quadict)
        lstsum[i] = sum(netlink)
    return lstsum


def linkSum(nodelist, quadlist):
    """
    计算一个分布中四个象限的连线和
    输入：所有线网列表 节点分配的字典列表
    输出：连线值列表 
    """
    maxnum = len(nodelist)
    nodeset = []  # 线网集合的list
    for i in nodelist:
        j = set(i)
        nodeset.append(j)

    quadset = []  # 象限集合的list
    for i in quadlist:
        j = set(list(i.keys()))
        quadset.append(j)

    netlink = []  # 四个象限的连接数
    for i in quadset:  # 循环象限
        link = 0
        for j in nodeset:  # 循环所有线网
            # if (len(a & b) != 0) & (len(a & b) != len(b)):
            if (len(i & j) != 0) & (len(i & j) != len(j)):
                link += 1
        # 如果某一区域的连线和为0
        if link == 0:
            link = maxnum
        netlink.append(link)
    return netlink


def findBest(pop, fitness):
    """
    根据染色体群体pop已经对应的适应值fitness
    找到最高的适应值f，以及对应的染色体bst和其在pop中的编号/下标ind
    """
    f = np.max(fitness)
    ind = np.asarray(np.where(fitness == f)).flatten()
    bst = pop[ind, :]
    return [bst, f, ind]


def satified(fitness):
    return 0


def chooseNewP(pop, fitness):
    """
    根据染色体的适应值 按照一定的概率 生成新一代染色体群体newpop
    取得原群体的某些染色体
    """
    N, nc = np.shape(pop)
    fitness = np.cumsum(fitness)  # 转为一维数组
    selelst = np.zeros([N, 1])  # 原群体中的某些染色体的编号
    for i in range(N):
        rval = np.random.rand()
        selelst[i] = 0
        for j in range(N - 1, -1, -1):
            if rval > fitness[j]:
                selelst[i] = j
                break
    newpop = pop[list(selelst.flatten().astype(np.uint8)), :]
    return newpop


def crossPop(pop, pc, fitness, SelfAdj):
    """
    根据交叉频率pc 以及各染色体的适应值fitness
    通过交叉的方式生成新群体
    SelfAdj = 1时为自适应，否则取固定的交叉概率pc
    """
    N, nc = np.shape(pop)
    crosslst = list(range(N))  # 原群体中的某些染色体的编号
    np.random.shuffle(crosslst)

    fmax = np.max(fitness)
    fmean = np.mean(fitness)
    k1 = pc
    k2 = pc
    i = 0

    while i < int(N / 2):
        rval = np.random.rand()
        j = np.random.randint(int(N / 2), N)

        if SelfAdj == 1:
            fprime = np.max(fitness[i], fitness[j])
            if fprime > fmean:
                pc = k1 * (fmax - fprime) / (fmax - fmean)
            else:
                pc = k2

        if rval > pc:
            continue

        partner1 = copy.copy(pop[crosslst[i], :])  # 浅复制
        partner2 = copy.copy(pop[crosslst[j], :])

        if (partner1 == partner2).all():  # 迭代 全相等
            continue

        child1, child2 = genecross(partner1, partner2)
        pop[crosslst[i], :] = copy.copy(child1)
        pop[crosslst[j], :] = copy.copy(child2)
        i = i + 1


def genecross(partner1, partner2):
    """
    父染色体partner1,partner2
    通过交叉方式
    生成两个子染色体child1,child2
    """
    length = len(partner1)
    idx1 = 0
    idx2 = 0
    while idx1 == idx2:
        idx1 = random.randint(0, length - 1)
        idx2 = random.randint(0, length - 1)
    ind1 = min(idx1, idx2)
    ind2 = max(idx1, idx2)

    child1 = copy.copy(partner1)
    child2 = copy.copy(partner2)

    tem1 = copy.copy(child1[ind1:ind2])
    tem2 = copy.copy(child2[ind1:ind2])

    if (tem1 == tem2).all():
        return [child1, child2]
    if set(tem1) == set(tem2):
        child1[ind1:ind2] = tem2
        child2[ind1:ind2] = tem1
        return [child1, child2]

    temdff1 = set(tem1).difference(set(tem2))
    temdff2 = set(tem2).difference(set(tem1))

    for i in range(len(temdff1)):
        child1[np.where(child1 == list(temdff2)[i])] = list(temdff1)[i]
        child2[np.where(child2 == list(temdff1)[i])] = list(temdff2)[i]

    child1[ind1:ind2] = tem2
    child2[ind1:ind2] = tem1
    return [child1, child2]


def mutPop(pop, pw, fitness, SelfAdj):
    """
    根据变异概率pw 以及各染色体的适应值fitness
    通过编译的方式生成新群体
    SelfAdj = 1时为自适应，否则取固定的交叉概率pw
    """
    N, nc = np.shape(pop)

    fm = max(fitness)
    fa = np.mean(fitness)
    k3 = pw
    k4 = pw
    for i in range(N):
        rval = random.random()
        if SelfAdj == 1:  # 自适应，参考149页
            if fitness[i] > fa:
                pw = k3 * (fm - fitness[i]) / (fm - fa)
            else:
                pw = k4
        if rval > pw:
            continue
        pop[i, :] = np.asarray(mutWeightGene(pop[i, :]))


def mutWeightGene(gene):
    """
    资源变异
    找出资源权重和最大与最小的区域
    交换某些节点
    """
    aredict = read_are(arepath)  # 所有节点
    pos = geneDict(gene, aredict)  # {'节点':坐标}
    quadict = posCal(pos, aredict)
    weightlist = weightCal(quadict)

    # 先取得最大权重和的下标 再在象限列表中取得对应的字典
    quadmin = quadict[weightlist.index(min(weightlist))]
    quadmax = quadict[weightlist.index(max(weightlist))]

    # 权重最值的键 也就是节点
    keymax = max(quadmax, key=quadmax.get)
    keymin = min(quadmin, key=quadmin.get)

    # 将 节点-权重 字典更新为 节点-坐标 字典
    for i in quadmax.keys():
        quadmax[i] = pos[i]

    for i in quadmin.keys():
        quadmin[i] = pos[i]

    # 返回权重最值的节点的坐标
    coormax = quadmax[keymax]
    coormin = quadmin[keymin]

    # 分为x y来交换
    pos1x = list(gene).index(coormax[0])
    pos1y = list(gene).index(coormax[1])
    pos2x = list(gene).index(coormin[0])
    pos2y = list(gene).index(coormin[1])

    listgene = swaPosition(list(gene), pos1x, pos2x)
    listgene = swaPosition(list(gene), pos1y, pos2y)

    return listgene


def swaPosition(list, pos1, pos2):
    """
    对调列表两个元素
    """
    list[pos1], list[pos2] = list[pos2], list[pos1]
    return list


def mutLinkPop(pop, pw, fitness, SelfAdj):
    """
    连接数变异
    """
    N, nc = np.shape(pop)

    fm = max(fitness)
    fa = np.mean(fitness)
    k3 = pw
    k4 = pw
    for i in range(N):
        rval = random.random()
        if SelfAdj == 1:  # 自适应，参考149页
            if fitness[i] > fa:
                pw = k3 * (fm - fitness[i]) / (fm - fa)
            else:
                pw = k4
        if rval > pw:
            continue
        pop[i, :] = mutNetLink(pop[i, :])


def mutNetLink(gene):
    """
    对连接数最大的区域进行变异
    (x, y)以一定的概率->(x, -y) (-x, y) (-x, -y)
    """
    aredict = read_are(arepath)  # 所有节点
    pos = geneDict(gene, aredict)  # {'节点':坐标}
    quadict = posCal(pos, aredict)  # {'节点':权重}
    nodelist = read_net(netpath)
    netlink = linkSum(nodelist, quadict)  # 四个象限的连接数
    quadnum = netlink.index(max(netlink))  # 连接数最大的象限

    # 将节点-权重更新为节点-坐标
    quadcopy = copy.deepcopy(quadict)
    for dict in quadcopy:
        for i in dict.keys():
            dict[i] = pos[i]

    mutquad = quadcopy[quadnum]
    for i in mutquad.keys():
        rval = np.random.rand()
        if rval < 0.25:
            mutquad[i] = [-mutquad[i][0], mutquad[i][1]]
        elif 0.25 < rval < 0.5:
            mutquad[i] = [-mutquad[i][0], -mutquad[i][1]]
        elif 0.5 < rval < 0.75:
            mutquad[i] = [mutquad[i][0], -mutquad[i][1]]
        else:
            mutquad[i] = [mutquad[i][0], mutquad[i][1]]
    # 将象限内节点的坐标信息传递给aredict
    mutaredict = copy.deepcopy(aredict)
    for dict in quadcopy:
        for i in dict.keys():
            mutaredict[i] = dict[i]

    # 更新gene中的坐标
    nc = int(len(gene) / 2)
    i = 0
    for key in mutaredict.keys():
        gene[i] = mutaredict[key][0]
        gene[i + nc] = mutaredict[key][1]
        i = i + 1

    return gene


def read_net(filepath):
    """
    读取net文件
    返回一个列表
    列表元素为同一个线网的节点列表
    """
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

    nodegroup = []
    for i in range(len(snode)):

        # i + 1 <= len(snode) - 1
        if i <= len(snode) - 2:
            # 同一个线网的节点
            nodegroup.append(netdata[snode[i]:snode[i + 1]])
        else:
            # 最后一个线网
            nodegroup.append(netdata[snode[i]:])
    """
    nodegroup[]
    [['g8 s 1', 'gp2 l', 'g10 l', 'g11 l', 'g13 l']
     ['g7 s 1', 'gp3 l', 'g10 l', 'g11 l']
     ['g6 s 1', 'gp4 l', 'g11 l']]
    """
    netnodes = []
    for nodei in nodegroup:  # 所有子图
        netnode = []  # 只含有节点的列表
        for nodej in nodei:  # 一个子图
            if nodej.find('s') > 0:
                # 找到s 第一个空格前的字符就是驱动节点
                netnode.append(nodej[:nodej.index(' ')])
            else:
                # 第一个空格前的字符就是负载节点
                netnode.append(nodej[:nodej.index(' ')])
        netnodes.append(netnode)
    """
    netnodes
    [['g8', 'gp2', 'g10', 'g11', 'g13'],
    ['g7', 'gp3', 'g10', 'g11'],
    ['g6', 'gp4', 'g11']]
    """
    return netnodes


def draw_net(netnodes, pos):
    """
    再传入一个pos 匹配线网后 进行绘图
    """
    g = nx.Graph()
    
    # 提取节点名
    for netnode in netnodes:
        # g.add_node(netnode[0])  # 添加驱动节点
        # for i in range(1, len(netnode)):
        #     g.add_node(netnode[i])  # 添加负载节点
        #     g.add_edge(netnode[0], netnode[i])  # 将该线网中的负载节点与驱动节点相连
        g.add_nodes_from(netnode)

    nx.draw(g, pos=pos, with_labels=True)
    plt.show()


def drawGraph(gene):
    """
    绘图 参考networkx的绘图
    """
    aredict = read_are(arepath)
    pos = geneDict(gene, aredict)
    for i in pos.keys():
        pos[i] = np.array(pos[i])
    netnodes = read_net(netpath)
    draw_net(netnodes, pos)


def writeResTxt(quadlist):
    """
    将结果写进txt文件
    """
    res.write('\n')
    now_time = datetime.datetime.now()
    res.write('- - -This result is at ' + str(now_time) + '- - -')
    res.write('\n')
    nosFPGA = len(quadlist)
    for fpga in range(nosFPGA):
        res.write("F"+str(fpga))
        dict = quadlist[fpga]
        res.write(" " + str(list(dict.keys())))
        res.write('\n')
"""
1 计算四块板各自的资源总和 Z1 Z2 Z3 Z4 使得它们的方差最小 M为平均数
  s^2 = [(M-Z1)^2 + (M-Z2)^2 + (M-Z3)^2 + (M-Z4)^2]/4
  arr = [Z1, Z2, Z3, Z4]
  M = np.mean(arr)
  s^2 = np.var(arr)

2 计算每块板与其他各个板的连接线 L1 L2 L3 L4 使得它们的和最小

3 计算平均值M 在所有Zi满足 0.9M<Zi<1.1M 的情况下 使得L1 L2 L3 L4的和最小
节点类似城市
计算距离类似计算资源总值
"""
# 超参数定义
arepath = "data\design.are"
netpath = "data\design.net"
showPlot = 1  # showPlot=1，显示最后的最佳分布

N = 100  # 染色体群体规模

MAXITER = 10**3  # 最大循环次数
SelfAdj = 0  #SelfAdj = 1时为自适应
randAccordDist = 0  # 在mutate的时候根据位置值分配对应的概率，此概率决定哪个城市参与mutate

pc = 0.85  # 交叉概率
pw = 0.15  # 变异概率
mode = 3   # 计算模式
# # # # # 遗传算法 # # # # #

# 步骤1 产生N个染色体的初始群体,保存在pop里面
aredict = read_are(arepath)
nc = len(aredict)  # 节点的数量
pop = initPop(N, nc)  # pop 遗传算法中的种群
nodelist = read_net(netpath)

iter = 0
tstart = time()

# # # 主循环 # # #
while iter < MAXITER:
    iter = iter + 1

    # 步骤2
    lstvar, _, _ = varCal(pop)  # 计算方差
    lstsum = linkCal(pop)  # 计算连线和

    # 分模式计算pop的fitness 
    fitness = fitnessCal(pop, mode)
    b, f, ind = findBest(pop, fitness)  # 当前种群中 适应度最高的个体

    # 步骤3 如果满足某些标准 算法停止
    if satified(fitness):
        break

    chooseNewP(pop, fitness)  # 否则 选取出一个新的群体
    crossPop(pop, pc, fitness, SelfAdj)  # 步骤4 交叉产生新染色体，得到新群体
    
    # 步骤5 基因变异
    # mutPop(pop, pw, fitness, SelfAdj)  # 权重变异 mode=1
    mutLinkPop(pop, pw, fitness, SelfAdj)  # 连接数变异 mode=2
    pop[0, :] = b[0, :]  # 保留上一代中适应值最高的染色体

    posdict = geneDict(gene=b[0, :], aredict=aredict)
    quadlist = posCal(posdict, aredict)  # 计算节点的象限分布
    weightlist = weightCal(quadlist)  # 计算各象限的权重和
    linklist = linkSum(nodelist, quadlist)  # 计算连线和

    if np.mod(iter, MAXITER / 10) == 1:
        titer = time()
        print("iter = %d after running %d seconds with variance = %f" %
              (iter, int(titer - tstart), min(lstvar)))

        print("The sum of the weights of each quadrant is ")
        print(weightlist)

        print('The link num of each quadrant is')
        print(linklist)

# # # # # 迭代结束 # # # # #
# 输出节点数量 种群大小 迭代次数 交叉概率 变异概率
print('- - - - -')
print('nc=%d,N=%d, MAXITER=%d, pc=%f, pw=%f ' % (nc, N, MAXITER, pc, pw))

# 输出最优数据
print('适应值为 %f, 最优的方差为 %f' % (f, min(lstvar)))

#计算运行时间
tend = time()
print('The program need :  %f seconds.' % (tend - tstart))
print('- - - - -')
# 最优的分布
print('最优的分布为')
posdict = geneDict(gene=b[0, :], aredict=aredict)
quadlist = posCal(posdict, aredict)
print('F0:' + str(list(quadlist[0].keys())))
print('F1:' + str(list(quadlist[1].keys())))
print('F2:' + str(list(quadlist[2].keys())))
print('F3:' + str(list(quadlist[3].keys())))

weightlist = weightCal(quadlist)
print("The sum of the weights of each quadrant is ")
print(weightlist)
nodelist = read_net(netpath)
linklist = linkSum(nodelist, quadlist)
print('The link num of each quadrant is')
print(linklist)
print('The sum of link %d' % (sum(linklist)))

# 写进结果文件
filename = "result.res"
res = open(filename,'a')
res.write('The result to algorithm %d' % (mode))
writeResTxt(quadlist)
res.write('variance: %f, linksum: %d' % (min(lstvar), sum(linklist)))
res.write('\n')
# 显示节点分布
if showPlot == 1:
    drawGraph(b[0, :])

