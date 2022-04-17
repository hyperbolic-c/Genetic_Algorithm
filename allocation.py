from operator import itemgetter
from os import popen
from platform import node
import random
import copy
from time import time
import numpy as np
import networkx as nx
from matplotlib import pyplot as plt
import warnings
import datetime


# from allocation import Node
class Node(object):
    """
    输入字典{'name':g0, 'weight':1, 'coord':[0.1, 0.2]}
    实例化之后[{'name':g0, 'weight':1, 'coord':[0.1, 0.2]},
              {'name':g1, 'weight':2, 'coord':[0.2, 0.2]]
              type = <class 'allocation.Node'>
    exp_node = Node(**exp_dict)
    """

    def __init__(self, **dict):
        self.__dict__.update(dict)
    
    def todone(self):
        print('done')
    
    def quad(self):
        x = self.coord[0]
        y = self.coord[1]
        if x > 0 and y > 0:
            self.rant = 1
            # self.
        elif x > 0 and y < 0:
            self.rant = 4  
        elif x < 0 and y > 0:
            self.rant = 2
        else:
            self.rant = 3

    # def node_quad(self):
        

class Net(object):
    def __init__(self, nodelist=None):
        self.nodelist = nodelist
        self.quadment = set()

    def todone(self):
        print('done')

    def calMent(self, node_list):
        b = PopGene.calQuad(node_list)

        set_list = []
        for i in b:
            n_list = []
            for exp in i:
                n_list.append(exp.name)
            set_list.append(set(n_list))
        j = 0
        for s in set_list:
            if 0 < len(set(self.nodelist) & s) < len(set(self.nodelist)):
                self.quadment.add(j)
            j = j + 1

        




class Initpre:
    """
    将关于节点的一些信息集成在一个字典
    """
    @classmethod
    def initPop(cls, N, nc):
        """
        产生N个染色体的初始群体，保存在pop
        将染色体编码为范围为(-1, 1)的浮点数 也就是53个节点的位置pos信息
        前nc个随机数为x 后nc个随机数为y
        编码为节点会方便一些
        引入坐标为计算布局 连线做准备
        """
        pop = -1 + 2 * np.random.random((N, nc * 2))
        return pop


    @classmethod
    def read_are(cls, filepath):
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

    
    @classmethod
    def read_net(cls, filepath):
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
    # @classmethod
    # @classmethod
    # @classmethod



class PopGene:
    
    @classmethod
    def fitnessCal(cls, pop, net_exp, mode):
        """
        根据方差/连线和 计算各个染色体的适应值fitness
        方差lstvar = np.zeros([N, 1])
        连线和lstsum = np.zeros([N, 1])
        
        3个模式 mode == 1方差最小 mode == 2连线和最小 mode == 3满足条件使连线和最小
        传入pop
        再调用函数计算var sum
        """
        lstvar, lstmean, listweight = cls().calVar(pop)
        lstsum = cls().calLink(pop, net_exp)

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


    @classmethod
    def calLink(cls, pop, net_exp):
        N, _ = np.shape(pop)
        lstsum = np.zeros([N ,1])
        # pop_exp = cls().popExp(pop)
        allgene = cls().geneCoord(pop)
        for i in range(N):
            # gene_exp = pop_exp[i]
            # gene = pop[i, :]
            g_exp = allgene[i]
            gene_exp = []
            for dict in g_exp:
                n_exp = Node(**dict)
                n_exp.quad()
                gene_exp.append(n_exp)
            # print(gene_exp[0].name)
            L = [0, 0, 0, 0]
            # quadlist = cls().calQuad(gene_exp)
            for net in net_exp:
                net.calMent(gene_exp)
                # print(net.quadment)
                if len(net.quadment) == 1:
                    continue
                for fpga in net.quadment:
                    L[fpga] += 1
            # print(L)
            lstsum[i] = sum(L)
        return lstsum
    
    # 计算方差
    @classmethod
    def calVar(cls, pop):
        N, _ = np.shape(pop)
        lstvar = np.zeros([N, 1])
        lstmean = np.zeros([N, 1])
        lstweight = []
        pop_exp = cls().popExp(pop)
        for i in range(N):
            gene_exp = pop_exp[i]
            quadlist = cls().calQuad(gene_exp)
            wlist = cls().calWeight(quadlist)
            lstvar[i] = np.var(wlist)
            lstmean[i] = np.mean(wlist)
            lstweight.append(wlist)
        return lstvar, lstmean, lstweight


    # 计算权重
    @classmethod
    def calWeight(cls, quadlist):
        wlist = []
        for nodelist in quadlist:
            wsum = 0
            for node_exp in nodelist:
                w = node_exp.weight
                wsum = w + wsum
            wlist.append(wsum)
        return wlist

    # 节点象限
    @classmethod
    def calQuad(cls, gene_exp):
        Z1, Z2, Z3, Z4 = [], [], [], []
        for node_exp in gene_exp:
            if node_exp.quad() == 1:
                Z1.append(node_exp)
            if node_exp.quad() == 2:
                Z2.append(node_exp)
            if node_exp.quad() == 3:
                Z3.append(node_exp)
            if node_exp.quad() == 4:
                Z4.append(node_exp)
                # print(node_exp)
        # print([Z1, Z2, Z3, Z4])
        return [Z1, Z2, Z3, Z4]

    # pop实例化
    @classmethod
    def popExp(cls, pop):
        allgene = cls().geneCoord(pop)
        pop_exp = []

        for allnode in allgene:
            gene_exp = []
            for nodedict in allnode:
                node_exp = Node(**nodedict)
                node_exp.quad()
                gene_exp.append(node_exp)
            pop_exp.append(gene_exp)
        # print(pop_exp)
        return pop_exp
    
    @classmethod
    def geneCoord(cls, pop):
        """
        将节点的所有信息汇总到一个字典
        返回pop中所有gene对应的节点信息
        allgene = [[{'name': 'g0', 'weight': 65, 'coord':[0.1, 0.2]}, {...], [{...]]
        """
        N, nc = np.shape(pop)
        n = int(nc / 2)
        aredict = Initpre.read_are(arepath)
        allgene = []
        for i in range(N):
            allnode = []
            gene = pop[i, :]
            # allnode = cls().geneDict(gene, aredict)
            j = 0
            for key, item in aredict.items():
                pos = {}
                pos['name'] = key
                pos['weight'] = item
                pos['coord'] = [gene[j], gene[j + n]]
                pos['rant'] = 0
                allnode.append(pos)
                j += 1
            allgene.append(allnode)
        # print(allgene)
        return allgene

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
    pop_exp = PopGene.popExp(pop)
    for i in range(N):
        rval = random.random()
        if SelfAdj == 1:  # 自适应，参考149页
            if fitness[i] > fa:
                pw = k3 * (fm - fitness[i]) / (fm - fa)
            else:
                pw = k4
        if rval > pw:
            continue
        # pop[i, :] = np.asarray(mutGene(pop[i, :]))
        pop[i, :] = mutGene(pop_exp[i])


def mutGene(gene_exp):
    nc = len(gene_exp)
    gene = np.zeros(nc*2)
    i = 0
    for node in gene_exp:
        rval = np.random.rand()
        if rval < 0.25:
            node.coord = [node.coord[0], node.coord[1]]
        elif 0.25 < rval < 0.5:
            node.coord = [-node.coord[0], node.coord[1]]
        elif 0.5 < rval < 0.75:
            node.coord = [node.coord[0], -node.coord[1]]
        else:
            node.coord = [-node.coord[0], -node.coord[1]]
        gene[i] = node.coord[0]
        gene[i + nc] = node.coord[1]
        i = i + 1
    return gene

if __name__ == '__main__':

    arepath = "data\design.are"
    netpath = "data\design.net"
    showPlot = 1  # showPlot=1，显示最后的最佳分布
    N = 100  # 染色体群体规模

    MAXITER = 10**3  # 最大循环次数
    SelfAdj = 0  #SelfAdj = 1时为自适应
    randAccordDist = 0  # 在mutate的时候根据位置值分配对应的概率，此概率决定哪个城市参与mutate

    pc = 0.85  # 交叉概率
    pw = 0.15  # 变异概率

    # 步骤1 产生N个染色体的初始群体,保存在pop里面
    aredict = Initpre.read_are(arepath)
    netnodes = Initpre.read_net(netpath)
    nc = len(aredict)  # 节点的数量
    pop = Initpre.initPop(N, nc)
    # allgene = PopGene.geneCoord(pop)

    # Node类实例化
    # pop_exp = PopGene.popExp(pop)

    # Net类实例化
    net_exp = []
    for i in netnodes:
        n_exp = Net(i)
        net_exp.append(n_exp)

    iter = 0
    tstart = time()
    while iter < MAXITER:
        iter = iter + 1

        # Node类实例化
        pop_exp = PopGene.popExp(pop)

        # Net类实例化
        net_exp = []
        for i in netnodes:
            n_exp = Net(i)
            net_exp.append(n_exp)
        
        # 计算方差 计算连线和
        lstvar, _, _ = PopGene.calVar(pop)
        lstsum = PopGene.calLink(pop, net_exp)
        # 分模式计算pop的fitness 
        fitness = PopGene.fitnessCal(pop, net_exp, 2)
        # print(pop)
        # for i in net_exp:
        #     print(i.quadment)
        b, f, ind = findBest(pop, fitness)  # 当前种群中 适应度最高的个体

        # 步骤3 如果满足某些标准 算法停止
        if satified(fitness):
            break

        chooseNewP(pop, fitness)  # 否则 选取出一个新的群体
        crossPop(pop, pc, fitness, SelfAdj)  # 步骤4 交叉产生新染色体，得到新群体
        mutPop(pop, pw, fitness, SelfAdj)  # 权重变异 mode=1
        # mutLinkPop(pop, pw, fitness, SelfAdj)  # 连接数变异 mode=2
        pop[0, :] = b[0, :]  # 保留上一代中适应值最高的染色体

        b_exp = pop_exp[ind[0]]
        quadlist = PopGene.calQuad(b_exp)
        weightlist = PopGene.calWeight(quadlist)

        if np.mod(iter, MAXITER / 10) == 1:
            titer = time()
            print("iter = %d after running %d seconds with variance = %f" %
                (iter, int(titer - tstart), min(lstvar)))

            print("The sum of the weights of each quadrant is ")
            print(weightlist)

            print('The link num of each quadrant is')
            print(min(lstsum))

            # print(fitness)
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