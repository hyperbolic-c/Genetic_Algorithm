
import numpy as np
import math
import re
import matplotlib.pyplot as plt
import matplotlib.patches as pth
import shapely.geometry.polygon as shplyg
from time import time
import random



class area():
    def __init__(self, ist, ru):
        self._pntList = ist
        self.SD = ru[0]
        self.GATE = ru[1]
        self.SD_GATE = ru[2]
        self.SD_ITO = ru[3]
        self.GATE_ITO = ru[4]


class polygon():
	def __init__(self, lst, id, port):
		n = len(lst)
		self._id = id                        # 多边形的 id 号
		self._nospnt = n                     # 多边形所包含的顶点个数
		self._color = np.random.rand(3)      # 多边形所填充的颜色（方便查看）

		self._angle = 0                      # 初始化，多边形的角度方向[0,1,2,3]分别表示[0度，90度，180度，270度]
		self._mirror = 0                     # 初始化，多边形的镜像[0,1]分别表示 [没有镜像、以center为中心上下镜像]
		self._shift = (0,0)                  # 初始化，多边形的位置偏置

		ar = np.asarray(lst)
		self._ld = np.min(lst,0)             # 把多边形围起来的矩阵的左下角坐标
		self._ru = np.max(lst,0)             # 把多边形围起来的矩阵的右上角坐标
		self._center = np.mean([self._ld,self._ru],0)		  # 多边形的中心点
		self._pntList = lst                  # 多边形的顶点坐标列表（根据angle/mirror/shift参数更新后的坐标）
		# self._pntTuple 为原始顶点坐标，不可变更
		pntLst = []
		for i in range(n):
			pntLst.append(lst[i])
		self._pntTuple = tuple(pntLst)
        
		self._port = port  # port列表中第1、3个是SD,第2个是GATE

	def getCorner(self):
		return [self._ld, self._ru]
		#return [self._ld[0],self._ld[1],self._ru[0],self._ru[1]]
	def getColor(self):
		return self._color
	def getID(self):
		return self._id
	def getCenter(self):
		return self._center
	def getCurrentPointsLst(self):
		return self._pntList

	# 根据给出来的参数，从***初始状态***下应用上述参数,更新多边形数据
	def updateParameters(self,x,y,mirror,angle):
		pntLst = list(self._pntTuple)  #从***初始状态***开始
		# 把多边形整体偏移，偏移量=[x,y]
		if x!=0 or y!=0:
			for i in range(self._nospnt):
				pntLst[i] = [pntLst[i][0]+x, pntLst[i][1]+y]
		self._shift = [x,y] 			 	# 更新
		self._center = np.mean([np.min(pntLst,0),np.max(pntLst,0)],0)	# 更新
		# 把多边形沿着 self.center 横轴做镜像
		self._mirror=mirror					# 更新
		if self._mirror!=0:
			doubleCenterY = 2*self._center[1]
			for i in range(self._nospnt):
				pntLst[i] = [pntLst[i][0],doubleCenterY-pntLst[i][1]]
		# 把多边形围绕着 self.center 旋转，旋转角度为=angle
		self._angle = angle *90				# 更新
		if self._angle!=0:
			for i in range(self._nospnt):
				pntLst[i] = self._rotatePoint(pntLst[i], self._center, self._angle)

		self._ld = np.min(pntLst,0)
		self._ru = np.max(pntLst,0)
		self._pntList = pntLst
	def _rotatePoint(self, mvpoint, pivotpoint, angle):
		# 顶点 mvpoint 围绕顶点 pivotpoint 旋转 angle 角度
		# 参考了 https://community.esri.com/t5/arcgis-pro-questions/rotate-polygon-using-attribute/td-p/1055568
		angle_rad = -math.radians(angle)
		ox, oy = pivotpoint
		px, py = mvpoint
		qx = ox + math.cos(angle_rad) * (px - ox) - math.sin(angle_rad) * (py - oy)
		qy = oy + math.sin(angle_rad) * (px - ox) + math.cos(angle_rad) * (py - oy)
		mvpoint = [round(qx,1), round(qy,1)]
		return mvpoint

# ------------------------------------------------------
# ----------------- polygon 类定义结束------------------
# ------------------------------------------------------



def implement(caseNum):
	# 根据案例号 caseNum 读取对应文件 filename ，并获得多边形列表 polyLst
	filename = "contest_case/Ports_area_etc_input_"+str(caseNum)+".txt"
	_,polyLst = readFile(filename)

	# 设置图像
	fig = plt.figure()
	# 设置子图121画布 ax
	ax = fig.add_subplot(121,aspect='equal')
	# 把多边形列表 polyLst 中所有多边形显示在 ax 中
	[ox, oy, dx, dy] = plotPolys(ax,polyLst)
	# 计算边框面积大小
	orgsize = (dx-ox)*(dy-oy)
	orgratio = (dx-ox)/(dy-oy)
	# 打印初始的box大小
	print("original box-size = ", orgsize, " , ratio = {:.2f}".format(orgratio))



	"""
	主要修改 changePolygons 函数和相关内容
	"""
	# 修改多边形列表 polyLst 中多边形的位置/镜像/角度
	changePolygons(polyLst)


	# 检查 polyLst
	checkPolygonsOverlap(polyLst)


	# 设置子图122画布 ax
	ax = fig.add_subplot(122,aspect='equal')
	# 把多边形列表 polyLst 中所有多边形显示在 ax 中
	[ox, oy, dx, dy] = plotPolys(ax,polyLst)
	# 计算边框面积大小
	newsize = (dx-ox)*(dy-oy)
	newratio = (dx-ox)/(dy-oy)
	# 打印更新后box大小
	print("new box-size = ", newsize, " , ratio = {:.2f}".format(newratio))
	print("----- new-size / org-size : ", int(newsize/orgsize*100), "%")

	# 显示plot
	plt.show()
	return [polyLst,newsize,newratio]

# 检查多边形之间的重叠
def checkPolygonsOverlap(polyLst):
	# 采用 shapely 库中polygon类，把当前所有多边形封装在列表 shPlygonList[]
	shPlygonList = []
	for i in range(len(polyLst)):
		shPlygonList.append(shplyg.Polygon(polyLst[i].getCurrentPointsLst()))
	# 两两多边形进行检查是否存在重叠
	for i in range(len(shPlygonList)):
		for j in range(i+1,len(shPlygonList)):
			if isOverlap(shPlygonList[i],shPlygonList[j]):
				print('overlap: %d , %d'%(i+1,j+1))

# 判断i，j两个多边形是否重叠
def isOverlap(shplgi,shplgj):
	itstn=shplgi.intersection(shplgj)
	if itstn.area!=0:
		return True
	else:
		return False

"""
函数 changePolygons 及相关代码是同学最需要修改的代码，
考核的核心就在这部分，即，如何变更多边形的状况，
使得外边框面积最小
"""
# 改变所有多边形的状况
def changePolygons(polyLst):
	prosspolys = polyLst
	autoCombine(prosspolys)

# # 读取文件filename里面的多边形，并把所有多边形放进列表polyLst并返回
# def readFile(filename):
# 	polyLst = []	#polyLst包含所有的多边形顶点坐标
# 	countPolyNos = 0
# 	fo = open(filename,"r")
# 	for line in fo:
# 		if line!="Polygon:\n":
# 			lst = []	#lst是某一个单独多边形顶点坐标
# 			while len(line)!=0:
# 				# result = re.match('\((\d+),(\d+)\)(.*)',line)
# 				result = re.search('\((\d+),(\d+)\)(.*)', line)
# 				lst.append([int(result.group(1)),int(result.group(2))])
# 				line = result.group(3)
# 			countPolyNos+=1   #计算多边形的个数，并作为当前多边形的ID号
# 			polyLst.append(polygon(lst,countPolyNos))
# 			pass
# 	fo.close()
# 	return  polyLst

#读取文件
def readFile(filename):
	#filename = "contest_case/Ports_area_etc_input_1.txt"
	fo = open(filename, "r")
	Area = []
	ist = []
	rul = []

	polyLst=[]
	id = ""
	lst = []
	port = []
	cnt=0#记录一个module的port的个数
	for line in fo:

		if line[0] == "A":
			result = re.findall(r'-?\d+\.?\d*e?-?\d*?', line)
			pntNum = int(len(result) / 2)  # area点的个数
			for i in range(0, pntNum * 2 - 1, 2):
				ist.append((float(result[i]), float(result[i + 1])))

		elif line[0] == "R":
			result = re.findall(r'-?\d+\.?\d*e?-?\d*?', line)
			i = 0
			while (i < 7):
				if i < 4:
					rul.append((float(result[i]), float(result[i + 1])))
					i = i + 2
				else:
					rul.append((float(result[i])))
					i = i + 1
		elif line[0] == "M":
			result = re.findall(r'M-?\d+', line)  # result为list
			id = ' '.join(result)  # 转化为字符串

		elif line[0] == "B":
			result = re.findall(r'-?\d+\.?\d*e?-?\d*?', line)
			pntNum = int(len(result) / 2)  # module点的个数
			for i in range(0, pntNum * 2 - 1, 2):
				lst.append((float(result[i]), float(result[i + 1])))
				print("读取的Boundary数据是",(float(result[i]), float(result[i + 1])))
		elif line[0] == "P":
			po=[]
			result = re.findall(r'-?\d+\.?\d*e?-?\d*?', line)
			pntNum = int(len(result) / 2)  # module点的个数
			for i in range(0, pntNum * 2 - 1, 2):
				po.append((float(result[i]), float(result[i + 1])))
			port.append(po)
			cnt=cnt+1
			if(cnt==3):
				polyLst.append(polygon(lst, id, port))
				cnt=0
				id=""
				lst=[]
				port=[]

	Area.append(area(ist, rul))
	return Area,polyLst


#根据给出的polygon类对象onePolygon，在画布ax显示出来
def plotonepoly(ax, onePolygon):
	# 把多边形的顶点列表围成一个patch，并放在画布ax中
	tupleLst = onePolygon.getCurrentPointsLst()
	patchesPoly = pth.Polygon(tupleLst, alpha=0.7)
	patchesPoly.set_color(onePolygon.getColor())
	ax.add_patch(patchesPoly)
	# 画多边形的中心点
	'''
	x,y = onePolygon.getCenter()
	patchesPolyOrg = pth.Circle((x,y),0.1)
	patchesPolyOrg.set_color('black')
	ax.add_patch(patchesPolyOrg)
	'''
	# 把id号显示在中心位置
	'''
	plt.text(x,y,onePolygon.getID(),fontsize=14)
	'''

# 把polyLst列表中的所有polygon类对象在画布ax中显示出来
# （在此代码中，一定要先执行 plotonepoly 函数，再执行 minBox 函数，
#  因为 minBox 中需要使用更新后的 ld 和 ru 属性，而两个属性在
# getCurrentPointsLst 函数中被更新。plotonepoly 函数则调用了
# getCurrentPointsLst 函数。否则，外框的显示不正确）
def plotPolys(ax,polyLst):
	# 找到所有polyLst的最小矩形外框
	[LDx, LDy, RUx, RUy] = minBox(polyLst)
	# 用Rectangle类对象rec表达minBox结果
	p,w,h = (LDx, LDy),RUx-LDx,RUy-LDy
	rec = pth.Rectangle(p,w,h,alpha=1,facecolor='white',edgecolor='red',linewidth=2)
	# 把矩形rec放进画布ax
	ax.add_patch(rec)

	# 把所有多边形逐个放进画布ax
	for i in range(len(polyLst)):
		plotonepoly(ax, polyLst[i])

	plt.grid()

	# 设置画布ax大小
	ax.set_xbound(LDx-10,RUx+10)
	ax.set_ybound(LDy-10,RUy+10)
	return [LDx, LDy, RUx, RUy]

# 找到一个最小的矩形，可以容纳 polyLst 里面所有多边形
# （此函数调用了 ld / ru 两个属性，使用前需确保属性已经被更新）
def minBox(polyLst):
	LDx, LDy = 10**6,10**6
	RUx, RUy = -10**6,-10**6

	for i in range(len(polyLst)):
		ld, ru = polyLst[i].getCorner()
		ldx, ldy = ld
		rux, ruy = ru

		if ldx<LDx:
			LDx = ldx
		if ldy<LDy:
			LDy = ldy
		if rux>RUx:
			RUx = rux
		if ruy>RUy:
			RUy = ruy
	return [LDx, LDy, RUx, RUy]



'''———————————add———————————'''
# 检测MOVE多边形与KEEP是否重叠
def checkPolysOverlap(lcorner,polyLst,MOVE,need):
	global check
	if check == FUSE:
		overlap = checkPolysOverFUSE(lcorner,polyLst[MOVE]._pntList)
		if overlap == INVALID: # FUSE检测失效，启动AROUND检测
			check = AROUND
			print('------FUSE check pattern is invaild!------')
			print('------------enter AROUND pattern!-----------')
	if check == AROUND:
		overlap = checkPolysOverAROUND(polyLst,MOVE,need)
	return overlap

# 检测MOVE多边形与KEEP是否重叠（FUSE）
def checkPolysOverFUSE(pntlst1,pntlst2):
	check1 = shplyg.Polygon(pntlst1) # 类型转换
	if not check1.is_valid:
		return INVALID
	check2 = shplyg.Polygon(pntlst2)
	return isOverlap(check1,check2)

# 检测MOVE多边形与KEEP是否重叠（AROUND）
def checkPolysOverAROUND(polyLst,MOVE,need):
	move = shplyg.Polygon(polyLst[MOVE].getCurrentPointsLst())
	for i in range(len(need)):
		checkpoly = shplyg.Polygon(polyLst[need[i]].getCurrentPointsLst()) # 类型转换
		if isOverlap(move,checkpoly):
			return True
	return False

# 检测多边形是否超出边框
def checkOverframe(lframe,Corner):
	if Corner[RU][Y]>lframe[RU][Y] or Corner[LD][Y]<lframe[LD][Y] or Corner[LD][X]<lframe[LD][X] or Corner[RU][X]>lframe[RU][X]:
		return True
	else:
		return False

# 判断拼接向量是否同向（处理特殊情况）
def checkVectorSameDirect(kpindex,mpindex,lcorner,lmovepnt,d1d2,kpflag):
	'''
	# 如果根本没有拼接到
	if lcorner[kpindex] != lmovepnt[mpindex]:
		return False
	'''
	# 第一个角点
	IKp1 = kpindex+1
	IMp1 = mpindex-d1d2
	if IKp1 > len(lcorner)-1:
		IKp1=0
	if IMp1 > len(lmovepnt)-1:
		IMp1=0
	elif IMp1 < 0:
		IMp1=len(lmovepnt)-1
	# 第二个角点
	IKp2 = kpindex-1
	IMp2 = mpindex+d1d2
	if IKp2 < 0:
		IKp2=len(lcorner)-1
	if IMp2 > len(lmovepnt)-1:
		IMp2=0
	elif IMp2 < 0:
		IMp2=len(lmovepnt)-1
	# 拼接向量
	dotmul1 = (lcorner[IKp1][X]-lcorner[kpindex][X])*(lmovepnt[IMp1][X]-lcorner[kpindex][X]) + (lcorner[IKp1][Y]-lcorner[kpindex][Y])*(lmovepnt[IMp1][Y]-lcorner[kpindex][Y]) # vK1·vM1
	dotmul2 = (lcorner[IKp2][X]-lcorner[kpindex][X])*(lmovepnt[IMp2][X]-lcorner[kpindex][X]) + (lcorner[IKp2][Y]-lcorner[kpindex][Y])*(lmovepnt[IMp2][Y]-lcorner[kpindex][Y]) # vK2·vM2

	if kpflag == COMMON:
		# 是否都同向
		if dotmul1 > 0 and dotmul2 > 0:
			return True
		else:
			return False
	else:
		# 一个同向即可
		if dotmul1 > 0 or dotmul2 > 0:
			return True
		else:
			return False

# 检测长宽比是否满足[0.9~1.1]
def checkRadioStandard(lminBoxcorner):
	if lminBoxcorner[RU][Y]-lminBoxcorner[LD][Y] > 1.1*(lminBoxcorner[RU][X]-lminBoxcorner[LD][X]) or lminBoxcorner[RU][Y]-lminBoxcorner[LD][Y] < 0.9*(lminBoxcorner[RU][X]-lminBoxcorner[LD][X]):
		print('ratio does not meet the requirements!')
	else:
		print('ratio meets the requirements~')

# 计算多边形的面积（return面积列表）
def computeArea(polyLst):
	larea = []
	for i in range(len(polyLst)):
		shPlygon = shplyg.Polygon(polyLst[i].getCurrentPointsLst())
		larea.append(shPlygon.area)
	return larea

# 计算最小外接矩形的面积
def computeMinBoxArea(Corner):
	return (Corner[RU][X]-Corner[LD][X]) * (Corner[RU][Y]-Corner[LD][Y])

# 计算多边形拼接的顺序及面积
def computeWait(polyLst):
	larea = computeArea(polyLst)
	if order == SEQUENT: 
		lwait = sorted(range(len(larea)),key=lambda i:larea[i],reverse=True)
	else: # order == RANDOM
		lwait = list(range(len(polyLst)))
		random.shuffle(lwait)
	return lwait,larea

# 计算边框
def computeFrame(larea):
	totalArea = 0
	for i in range(len(larea)):
		totalArea += larea[i]
	width = height = math.ceil(mult * np.sqrt(totalArea))
	return [[0+X0,0+Y0],[width+X0,height+Y0]],width

# 计算所有多边形最小外接矩形的对角线长度
def computeDiagonal(polyLst):
	ldiagdist = []
	for i in range(len(polyLst)):
		ld, ru = polyLst[i].getCorner()
		ldiagdist.append(np.linalg.norm(ld-ru))
	return ldiagdist

# 计算old状态
def computeOldState(polyLst,MOVE):
	cen_old = polyLst[MOVE]._center
	mir_old = polyLst[MOVE]._mirror
	ang_old = polyLst[MOVE]._angle
	dx0 = cen_old[X] - polyLst[MOVE]._center[X]
	dy0 = cen_old[Y] - polyLst[MOVE]._center[Y]
	return dx0,dy0,mir_old,ang_old

# 计算需要检测重叠的多边形
def computeNeedToCheck(ldiagdist,MOVE,combined,lcombcen,keypnt):
	lneed = []
	r1 = ldiagdist[MOVE]
	for i in range(len(combined)):
		dist = np.linalg.norm(lcombcen[i]-keypnt)
		r2 = ldiagdist[combined[i]]/2
		if dist < r1+r2:
			lneed.append(combined[i])
	return lneed

# 计算多边形点集物理下标递增的方向
def computeDirect(pntlst):
	ring = shplyg.LinearRing(pntlst)
	if ring.is_ccw:
		return 1 # 顺时针
	else:
		return -1

# 计算多边形朝向（空缺角点）
def computeOrient(polygon,pattern):
	nospnt = polygon._nospnt
	if nospnt == 4:
		if pattern == FULL:
			return 0,0,0,0
		return 0
	
	Corner,pntList = polygon.getCorner(),polygon._pntList
	nlu = nru = nld = nrd = 1
	for i in range(len(pntList)):
		if pntList[i][Y] == Corner[RU][Y] and pntList[i][X] == Corner[LD][X]:
			nlu-=1
		elif pntList[i][Y] == Corner[RU][Y] and pntList[i][X] == Corner[RU][X]:
			nru-=1
		elif pntList[i][Y] == Corner[LD][Y] and pntList[i][X] == Corner[LD][X]:
			nld-=1
		elif pntList[i][Y] == Corner[LD][Y] and pntList[i][X] == Corner[RU][X]:
			nrd-=1
	
	if pattern == FULL:
		return nlu,nru,nld,nrd
	# pattern == PART
	n1,n2 = (nlu+nru),(nld+nrd)
	if n1==n2==0:
		return 0
	return (n1-n2)/(n1+n2)

# 计算KEEP和MOVE连续重合的点数，或者连接点
def computeOverandConnect(kpindex,mpindex,lcorner,lmovepnt,d1d2,pattern,kpflag):
	Kct1=Kct2=kpindex
	Mct1=Mct2=mpindex
	count = 1

	# 判断拼接向量是否同向
	if checkVectorSameDirect(kpindex,mpindex,lcorner,lmovepnt,d1d2,kpflag) == False:
		return False

	# 找到KEEP、MOVE的第1个连接点
	while True:
		Kct1+=1
		Mct1-=d1d2
		if Kct1 > len(lcorner)-1:
			Kct1=0
		if Mct1 > len(lmovepnt)-1:
			Mct1=0
		elif Mct1 < 0:
			Mct1=len(lmovepnt)-1
		if lcorner[Kct1] != lmovepnt[Mct1]:
			break
		count+=1

	# 找到KEEP、MOVE的第2个连接点
	while True:
		Kct2-=1
		Mct2+=d1d2
		if Kct2 < 0:
			Kct2=len(lcorner)-1
		if Mct2 > len(lmovepnt)-1:
			Mct2=0
		elif Mct2 < 0:
			Mct2=len(lmovepnt)-1
		if lcorner[Kct2] != lmovepnt[Mct2]:
			break
		count+=1
	
	if pattern == OVERLAP:
		return count
	else:
		return [Kct1,Kct2,Mct1,Mct2]

# 计算coc
def computeCOC(kpindex,mpindex,lcorner,lmovepnt,kpflag):
	d1d2=computeDirect(lcorner)*computeDirect(lmovepnt)
	return computeOverandConnect(kpindex,mpindex,lcorner,lmovepnt,d1d2,OVERLAP,kpflag)

# 计算适应度
def computeFitness(kpindex,mpindex,lcorner,polyLst,MOVE,lminBoxcorner,nc,sm1,kpflag):
	# ncoc/nc
	ncoc = computeCOC(kpindex,mpindex,lcorner,polyLst[MOVE]._pntList,kpflag)
	if ncoc == False: # 不合格状态
		return False
	# sm1/sm3
	tmp = updateMinBox(lminBoxcorner,polyLst[MOVE]._pntList)
	sm3 = computeMinBoxArea(tmp)
	# (n1-n2)/(n1+n2)
	dn = computeOrient(polyLst[MOVE],PART)
	# 总适应度
	fitness = A*(ncoc/nc) + B*(sm1/sm3) + C*(dn)
	return fitness

# 初始化KEEP的位姿
def initialTranAndPose(X0,Y0,polyLst,KEEP):
	if polyLst[KEEP]._nospnt == 4: # 矩形
		# 平移
		dx = X0 - polyLst[KEEP]._ld[X]
		dy = Y0 - polyLst[KEEP]._ld[Y]
		polyLst[KEEP].updateParameters(dx,dy,0,0)
	else: # 非矩形
		nlu,nru,nld,nrd = computeOrient(polyLst[KEEP],FULL)
		if nru == True:
			dx = X0 - polyLst[KEEP]._ld[X]
			dy = Y0 - polyLst[KEEP]._ld[Y]
			polyLst[KEEP].updateParameters(dx,dy,0,0)
		elif nlu == True: # 顺时针旋转90°
			polyLst[KEEP].updateParameters(0,0,0,1)
			dx = X0 - polyLst[KEEP]._ld[X]
			dy = Y0 - polyLst[KEEP]._ld[Y]
			polyLst[KEEP].updateParameters(dx,dy,0,1)
		elif nld == True: # 顺时针旋转180°
			polyLst[KEEP].updateParameters(0,0,0,2)
			dx = X0 - polyLst[KEEP]._ld[X]
			dy = Y0 - polyLst[KEEP]._ld[Y]
			polyLst[KEEP].updateParameters(dx,dy,0,2)
		elif nrd == True: # 镜像
			polyLst[KEEP].updateParameters(0,0,1,0)
			dx = X0 - polyLst[KEEP]._ld[X]
			dy = Y0 - polyLst[KEEP]._ld[Y]
			polyLst[KEEP].updateParameters(dx,dy,1,0)

# 更新最小外接矩形
def updateMinBox(lminBoxcorner,newpntlst):
	ldx,ldy = lminBoxcorner[LD][X],lminBoxcorner[LD][Y]
	rux,ruy = lminBoxcorner[RU][X],lminBoxcorner[RU][Y]
	for i in range(len(newpntlst)):
		lx,ry = newpntlst[i][X],newpntlst[i][Y]
		if lx < ldx:
			ldx = lx
		if ry < ldy:
			ldy = ry
		if lx > rux:
			rux = lx
		if ry > ruy:
			ruy = ry
	return [[ldx,ldy],[rux,ruy]]

# 更新角点
def updateCorners(kpindex,mpindex,lcorner,lmovepnt,kpflag):
	# 2个多边形方向乘积
	d1d2=computeDirect(lcorner)*computeDirect(lmovepnt)

	# 计算连接点
	tmp = computeOverandConnect(kpindex,mpindex,lcorner,lmovepnt,d1d2,CONNECT,kpflag)
	if tmp == False:
		return False
	Kct1,Kct2,Mct1,Mct2 = tmp[0],tmp[1],tmp[2],tmp[3]

	# 更新角点
	newcorner=[]

	# 插入KEEP的角点
	while Kct1 != Kct2:
		newcorner.append(lcorner[Kct1])
		Kct1+=1
		if Kct1>len(lcorner)-1:
			Kct1 = 0
	newcorner.append(lcorner[Kct2])

	# 插入MOVE的角点
	while Mct2 != Mct1:
		newcorner.append(lmovepnt[Mct2])
		Mct2+=d1d2
		if Mct2 > len(lmovepnt)-1:
			Mct2=0
		elif Mct2 < 0:
			Mct2=len(lmovepnt)-1
	newcorner.append(lmovepnt[Mct1])

	return newcorner

# 更新关键点及拼接模式
def updateKeypoint(lminBoxcorner,newcorners):
	if len(newcorners) == 4:
		keypoints = newcorners
		kpindex = list(range(len(newcorners)))
		kpflag = SPECIAL
		print('updated keypoints do not exist kpflag == SPECIAL')
	else:
		keypoints,kpindex = [],[]
		kpflag = COMMON
		for i in range(len(newcorners)):
			if lminBoxcorner[LD][X] < newcorners[i][X] < lminBoxcorner[RU][X] and lminBoxcorner[LD][Y] < newcorners[i][Y] < lminBoxcorner[RU][Y]:
				keypoints.append(newcorners[i])
				kpindex.append(i)
		print('updated keypoints exist, kpflag == COMMON~')
	return keypoints,kpindex,kpflag

# 打印MEET选取或尾插入的信息
def printInform(times,nsuccess,polyLst,moveid,pattern):

	print('——————————————————————————————————————————')
	if pattern == MEET:
#		print('[ %d) 第 %d/%d 个多边形插入，序号为%d ]' % (times,nsuccess+1,len(polyLst),moveid))
		print('[ %d) 第 %d/%d 个多边形插入，序号为%s ]' % (times,nsuccess+1,len(polyLst),moveid))
	else: # pattern == TAIL
		print('[ %d) 第 %d/%d 个多边形插入，尾处理... ]' % (times,nsuccess+1,len(polyLst)))
	if check == FUSE:
		print('检测模式为：FUSE')
	else:
		print('检测模式为：AROUND')

# 生成第二类关键点
def createKeypoint(lcorner):
	kpflag = SPECIAL
	print('created keypoints exist, kpflag == SPECIAL!')
	kpindex = list(range(len(lcorner)))
	return lcorner,kpindex,kpflag

# 多边形随角点移至关键点处
def changePolygonTranAndPose(polyLst,MOVE,lkeypnt,combkeypnt,combcorner,mirror,angle,pattern):
	polyLst[MOVE].updateParameters(0,0,mirror,angle)
	dx = lkeypnt[combkeypnt][X] - polyLst[MOVE]._pntList[combcorner][X]
	dy = lkeypnt[combkeypnt][Y] - polyLst[MOVE]._pntList[combcorner][Y]
	polyLst[MOVE].updateParameters(dx,dy,mirror,angle)
	if pattern == FORCE:
		return dx,dy

# MEET选取一个多边形
def meetSelect(polyLst,lcorner,lkeypnt,lkpindex,lwait,ldiagdist,lcombined,lcombcen,MOVE,lminBoxcorner,kpflag,lframe):
	combkeypnt = combcorner = find = maxF = 0
	nc = polyLst[MOVE]._nospnt
	sm1 = computeMinBoxArea(lminBoxcorner)
	# 达到阈值即可
	for combkeypnt in range(len(lkeypnt)):
		if len(lwait) < thfail: # 尾处理
			break
		need = computeNeedToCheck(ldiagdist,MOVE,lcombined,lcombcen,lkeypnt[combkeypnt])
		for combcorner in range(polyLst[MOVE]._nospnt):
			for mirror in range(2):
				for angle in range(4):
					# 多边形随角点移至关键点处
					changePolygonTranAndPose(polyLst,MOVE,lkeypnt,combkeypnt,combcorner,mirror,angle,MEET)
					# 判断是否重叠及出界
					overlap = checkPolysOverlap(lcorner,polyLst,MOVE,need)
					overframe = checkOverframe(lframe,polyLst[MOVE].getCorner())
					if overlap or overframe:
						continue
					# 计算适应度
					fitness = computeFitness(lkpindex[combkeypnt],combcorner,lcorner,polyLst,MOVE,lminBoxcorner,nc,sm1,kpflag)
					if fitness == False: # 不合格状态
						continue
					if fitness > maxF:
						maxF = fitness
					if fitness > thF: # 适应度是否达标
						find = 1
						print('MOVE = %d find~' % polyLst[MOVE]._id)
						break
				if find == True:
					break
			if find == True:
				break
		if find == True:
			break
	return find,maxF,combkeypnt,combcorner

# FORCE选取一个多边形
def forceSelect(lwait,polyLst,lkeypnt,lkpindex,ldiagdist,lcombined,lcombcen,lcorner,lframe,lminBoxcorner,kpflag):
	maxfit = maxid = maxdx = maxdy = maxmir = maxang = maxcombkp = maxcombcorn = 0
	sm1 = computeMinBoxArea(lminBoxcorner)
	print('----------选取适应度最高的多边形----------')
	for i in range(min(thfail,len(lwait))):
		TMOVE = lwait[i]
		nc = polyLst[TMOVE]._nospnt
		print('polygon currently detected is',TMOVE+1)
		# 保存旧状态
		dx0,dy0,mir_old,ang_old = computeOldState(polyLst,TMOVE)
		for combkeypnt in range(len(lkeypnt)):
			need = computeNeedToCheck(ldiagdist,TMOVE,lcombined,lcombcen,lkeypnt[combkeypnt])
			for combcorner in range(polyLst[TMOVE]._nospnt):
				for mirror in range(2):
					for angle in range(4):
						# 多边形随角点移至关键点处
						dx,dy = changePolygonTranAndPose(polyLst,TMOVE,lkeypnt,combkeypnt,combcorner,mirror,angle,FORCE)
						# 判断是否重叠及出界
						overlap = checkPolysOverlap(lcorner,polyLst,TMOVE,need)
						overframe = checkOverframe(lframe,polyLst[TMOVE].getCorner())
						if overlap or overframe:
							continue
						# 计算适应度
						fitness = computeFitness(lkpindex[combkeypnt],combcorner,lcorner,polyLst,TMOVE,lminBoxcorner,nc,sm1,kpflag)
						if fitness == False:
							continue
						if fitness > maxfit:
							maxfit,maxid,maxdx,maxdy,maxmir,maxang,maxcombkp,maxcombcorn = fitness,i,dx,dy,mirror,angle,combkeypnt,combcorner		
		# 恢复旧状态
		polyLst[TMOVE].updateParameters(dx0,dy0,mir_old,ang_old)
	return maxfit,maxid,maxdx,maxdy,maxmir,maxang,maxcombkp,maxcombcorn

# 将多边形移至等待列表的随机位置
def shuffleWait(thfail,lwait):
	nmove = min(thfail,len(lwait))
	moveid = [i for i in range(len(lwait))]
	random.shuffle(moveid)
	for i in range(nmove):
		lwait.insert(moveid[i],lwait.pop(0))

# 扩展边框
def extendFrame(lframe,length):
	for i in range(2):
		lframe[RU][i] += math.ceil(exstep*length)
	for i in range(2):
		lframe[LD][i] -= math.ceil(exstep*length)

# 修复边框
def restoreFrame(lframe,lastframe,lminBoxcorner):
	for i in range(2): 
		lframe[RU][i] = max(lastframe[RU][i],lminBoxcorner[RU][i])
	for i in range(2):
		lframe[LD][i] = min(lastframe[LD][i],lminBoxcorner[LD][i])

# 自动拼接
def autoCombine(polyLst):
	'''——————initial——————'''
	start_time = time()
	# KEEP多边形
	lwait,larea = computeWait(polyLst)
	KEEP = lwait.pop(0)
	lcombined = [KEEP]
	initialTranAndPose(X0,Y0,polyLst,KEEP)
	lminBoxcorner = [polyLst[KEEP]._ld,polyLst[KEEP]._ru]
	lcorner = polyLst[KEEP]._pntList
	lkeypnt,lkpindex,kpflag = updateKeypoint(lminBoxcorner,lcorner)
	lcombcen = [polyLst[KEEP]._center]
	# 边框
	lframe,length = computeFrame(larea)
	lastframe = []
	# MOVE多边形
	times,nfail,nsuccess,forceflag,danger = 0,0,1,OFF,False
	ldiagdist = computeDiagonal(polyLst)

	'''——————依次拼接——————'''
	while nsuccess < len(polyLst):
		if len(lwait) <= thfail: # 尾处理
			forceflag = ON

		'''——————MEET拼接——————'''
		if forceflag == OFF: # MEET选取
			# 选取MOVE
			times+=1
			MOVE = lwait.pop(0)
			printInform(times,nsuccess,polyLst,polyLst[MOVE]._id,MEET)
			# 保存旧状态
			dx0,dy0,mir_old,ang_old = computeOldState(polyLst,MOVE)
			# MEET选取
			find,maxF,combkeypnt,combcorner = meetSelect(polyLst,lcorner,lkeypnt,lkpindex,lwait,ldiagdist,lcombined,lcombcen,MOVE,lminBoxcorner,kpflag,lframe)
			print('max fitness：%.3f' % maxF)
			
			# 如果找到了，为下一次MEET选取做准备（更新最小外接矩形、融合角点、融合关键点）
			if find == True:
				lcombined.append(MOVE)
				lcombcen.append(polyLst[MOVE]._center)
				nfail = 0
				# forceflag = ON
				nsuccess+=1
				lmovepnt = polyLst[MOVE]._pntList
				lminBoxcorner = updateMinBox(lminBoxcorner,lmovepnt)
				if danger == True:
					restoreFrame(lframe,lastframe,lminBoxcorner)
					danger = False
				lcorner = updateCorners(lkpindex[combkeypnt],combcorner,lcorner,lmovepnt,kpflag)
				lkeypnt,lkpindex,kpflag = updateKeypoint(lminBoxcorner,lcorner)
			
			# 如果没找到
			else:
				nfail+=1
				# 恢复到old状态
				polyLst[MOVE].updateParameters(dx0,dy0,mir_old,ang_old)
				if len(lwait) >= thfail: # 尾处理
#					print(type(polyLst[MOVE]._id))
					print('(≧ ︿ ≦) MOVE = %s Not found!' % polyLst[MOVE]._id)
#					print('(≧ ︿ ≦) MOVE = %d Not found!' % polyLst[MOVE]._id)
				# 将拼接失败的多边形插入到等待队列，开启FORCE选取
				lwait.insert(min(thfail-1,len(lwait)),MOVE)
				if nfail >= thfail:
					forceflag = ON

		'''——————FORCE拼接——————'''
		# 如果连续thfail次，或需要尾处理
		# 连续thfail次：情况一： 处于无效关键点状态；情况二：处于阈值过高状态
		if forceflag == ON:
			if len(lwait) <= thfail:
				printInform(times,nsuccess,polyLst,-1,-1)
			# 生成足够多的关键点
			if kpflag == COMMON:
				lkeypnt,lkpindex,kpflag = createKeypoint(lcorner)
			# FORCE选取
			maxfit,maxid,maxdx,maxdy,maxmir,maxang,maxcombkp,maxcombcorn = forceSelect(lwait,polyLst,lkeypnt,lkpindex,ldiagdist,lcombined,lcombcen,lcorner,lframe,lminBoxcorner,kpflag)

			# 如果找到了，为下一次MEET选取做准备
			if maxfit > 0:
				nfail = 0
				forceflag = OFF
				nsuccess+=1
				MOVE = lwait.pop(maxid)
				print('(ง •_•)ง选取适应度最高',round(maxfit,3),'的多边形 %d 插入' % (MOVE+1))
				lcombined.append(MOVE)
				polyLst[MOVE].updateParameters(maxdx,maxdy,maxmir,maxang)
				lcombcen.append(polyLst[MOVE]._center)
				lmovepnt = polyLst[MOVE]._pntList
				lminBoxcorner = updateMinBox(lminBoxcorner,lmovepnt)
				if danger == True:
					restoreFrame(lframe,lastframe,lminBoxcorner)
					danger = False
				lcorner = updateCorners(lkpindex[maxcombkp],maxcombcorn,lcorner,lmovepnt,kpflag)
				if len(lwait) <= thfail:
					lkeypnt,lkpindex,kpflag = createKeypoint(lcorner)
				else:
					lkeypnt,lkpindex,kpflag = updateKeypoint(lminBoxcorner,lcorner)

			# 仍然找不到
			else:
				print('----------------仍然失败！----------------')
				if not danger and len(lwait)>thfail: # 移至随机位置，开启危险预警，继续FORCE选取
					print('---------------移至随机位置---------------')
					shuffleWait(thfail,lwait)
					danger = True
					lastframe = [[lframe[LD][X],lframe[LD][Y]],[lframe[RU][X],lframe[RU][Y]]]
				else: # 危险状态下扩展边框，继续FORCE选取
					print('----------(ﾟДﾟ*)ﾉ extend frame!----------')
					extendFrame(lframe,length)

	# 结束
	print('\n——————————————————————————————————————————————————')
	print('总共需要拼接的多边形：%d个，成功拼接的多边形：%d个' % (len(polyLst),nsuccess))
	checkRadioStandard(lminBoxcorner)
	end_time = time()
	print('The program need %.3f seconds' % (end_time-start_time))
	return 0

'''———————字符常量———————'''
RANDOM  = MEET  = FUSE   = OVERLAP = PART = COMMON  = ON  = LD = X = 0
SEQUENT = FORCE = AROUND = CONNECT = FULL = SPECIAL = OFF = RU = Y = 1
INVALID = -1

'''——————超参数设置——————'''
check  = FUSE	 # 检测重叠模式（FUSE AROUND）
order  = SEQUENT # 多边形拼接的顺序（SEQUENT RANDOM）
A,B,C  = 1.3,0.45,0.1 # 适应度系数
thF    = 0.72	 # 适应度阈值
X0,Y0  = 0,0	 # 原点
mult   = 1.042   # 边框大小（SEQUENT:1.042 RANDOM:1.105）
exstep = 0.019   # 边框扩展步长
thfail = 1		 # 连续失败次数阈值（SEQUENT:1 RANDOM:3）

#———————————main———————————#
if __name__ == '__main__':
	iter = 1 # 迭代多少次
	caseNum = 1 # 第几个文件
	polyLst,a,b = implement(caseNum)

	# 保存 polyLst 到本地文件
	npyfile='savePolygonsInCase'+str(caseNum)+'InIter'+str(iter)+'.npy'
	np.save(npyfile, polyLst)
	print("data has been saved.")
	# polyLst=np.load(npyfile)      # 读取文件


