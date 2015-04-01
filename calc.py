#coding=utf-8
import xml.dom.minidom
import pylab as plb
import numpy as np
SEQ="MKAILVVLLYTFATANADTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKLRGVAPLHLGKCNIAGWILGNPECESLSTASSWSYIVETPSSDNGTCYPGDFIDYEELREQLSSVSSFERFEIFPKTSSWPNHDSNKGVTAACPHAGAKSFYKNLIWLVKKGNSYPKLSKSYINDKGKEVLVLWGIHHPSTSADQQSLYQNADTYVFVGSSRYSKKFKPEIAIRPKVRDQEGRMNYYWTLVEPGDKITFEATGNLVVPRYAFAMERNAGSGIIISDTPVHDCNTTCQTPKGAINTSLPFQNIHPITIGKCPKYVKSTKLRLATGLRNIPSIQSRGLFGAIAGFIEGGWTGMVDGWYGYHHQNEQGSGYAADLKSTQNAIDEITNKVNSVIEKMNTQFTAVGKEFNHLEKRIENLNKKVDDGFLDIWTYNAELLVLLENERTLDYHDSNVKNLYEKVRSQLKNNAKEIGNGCFEFYHKCDNTCMESVKNGTYDYPKYSEEAKLNREEIDGVKLESTRIYQILAIYSTVASSLVLVVSLGAISFWMCSNGSLQCRICI"
# node为DOM的叶节点
def getValue(leafNode):
	value=leafNode.childNodes[0].nodeValue
	return str(value)
# 获取info.xml中的scale，以词典的形式返回
def getScale(filename="info.xml",propertyID="Hydrophlicity"):#,subpropertyID="0"):
	scaleDict={}
	doc=xml.dom.minidom.parse(filename) 		#加载xml文件
	root=doc.documentElement			#获取xml文档对象	
	nodes=root.getElementsByTagName("property")	#获取property节点集合
	for node in nodes:
		if str(node.getAttribute("id"))==propertyID:#获得每个节点的id,并转换为string
			scaleNode=node.getElementsByTagName("scale")[0] #获得scale节点
			lis=scaleNode.getElementsByTagName("li")	#获取property节点集合
			for li in lis:
				aa=getValue(li.getElementsByTagName("aa")[0])
				value=getValue(li.getElementsByTagName("value")[0])
				scaleDict[aa]=value
	return scaleDict
# Check if the amino acid sequnce include some worng symbols
def checkSeq(seq):
	stdSet=set("ACDEFGHIKLMNPQRSTVWY")
	seq=seq.upper()
	checkSet=set(seq)
	if not(checkSet.issubset(stdSet)):
		raise Exception("Wrong amino acid sequnce!!!")
# Calculate the profile of a given window using "Hydrophlicity","Antigenicity","Beta_turn"
def windowCalc(window,scaleDict=getScale()):
	checkSeq(window)
	winSum=0.0;
	window=window.upper()
	for index,aa in enumerate(window):
		winSum+=float(scaleDict[aa])
	return winSum/len(window)
# Calculate the profile of a given window using "Accessibility"
def windowCalcAcc(window,scaleDict=getScale(propertyID="Accessibility")):
	checkSeq(window)
	window=window.upper()
	product=1.0
	for index,aa in enumerate(window):
		product*=float(scaleDict[aa])
	product*=pow(0.37,-6)
	return product
# Split amino acid sequnce and return peptide list
def seqSplit(seq,windowSize):
	windowList=[]
	for i in range(len(seq)-windowSize+1):
		windowList.append(seq[i:windowSize+i])
	return windowList
# Get profile matrix of a window list using "Hydrophlicity","Antigenicity","Beta_turn" 
def profileMat(windowList,*scaleDictTuple):
	row=len(windowList)
	column=len(scaleDictTuple)
	proMat=np.zeros((row,column))	# 先定义一个每个元素都为0.的多为数组
	for index in range(column):
		windowScaleList=[]
		for window in windowList:
			windowScaleList.append(windowCalc(window,scaleDictTuple[index]))
		proMat[:,index]=windowScaleList
	return proMat
# Get the profile of a give sequence using "Accessibility"
def profileAcc(sequence, windowSize):
	windowList=seqSplit(sequence,windowSize)
	windowProfileList=[]
	for window in windowList:
		windowProfileList.append(windowCalcAcc(window))
	average=sum(windowProfileList)/len(windowProfileList)
	return windowList,[windowPro/average for windowPro in windowProfileList]
# Get the profile of a give sequence using "Flexibility"
def profileFlex(sequence, windowSize):    
        AA = ["K","S","G","P","D","E","Q","T","N","R","A","L","H","V","Y","I","F","C","W","M"]
        BNORM0 = [1.093,1.169,1.142,1.055,1.033,1.094,1.165,1.073,1.117,1.038,1.041,0.967,0.982,0.982,0.961,1.002,0.930,0.960,0.925,0.947]
        BNORM1= [1.082,1.048,1.042,1.085,1.089,1.036,1.028,1.051,1.006,1.028,0.946,0.961,0.952,0.927,0.930,0.892,0.912,0.878,0.917,0.862]
        BNORM2 = [1.057,0.923,0.923,0.932,0.932,0.933,0.885,0.934,0.930,0.901,0.892,0.921,0.894,0.913,0.837,0.872,0.914,0.925,0.803,0.804]
        WT = [0.25,0.50,0.75,1.00,0.75,0.50,0.25]
        l= len(sequence)
	RR=[0 for x in range(l)]
	for i in range(l):
		if AA.index(sequence[i])>9:
			RR[i]=1
	NAYB=[RR[x-1]+RR[x+1] for x in range(1,l-1)]
        NAYB.insert(0, 0)
        NAYB.append(0)
	windowList=[]
        scaleFlex= []
        for i in range(l-windowSize+1):
            sum = 0
            peptide = ""
            for j in range(windowSize):
                res = sequence[i+j:i+j+1]
                peptide += res
                index = AA.index(res)
                if NAYB[i+j] == 0:
                    sum += BNORM0[index] * WT[j] / 4.0
                if NAYB[i+j] == 1:
                    sum += BNORM1[index] * WT[j] / 4.0
                if NAYB[i+j] == 2:
                    sum += BNORM2[index] * WT[j] / 4.0
            windowList.append(peptide)
	    scaleFlex.append(sum)
        return  windowList,scaleFlex
#
def drowProfile(winScaleList,windowSize=1):
	x=range(len(winScaleList))
	x=[i+windowSize/2 for i in x]
	fig=plb.figure(1,facecolor='white')
	plb.plot(x,winScaleList)
	plb.show()
def scatterMat(mat,labels):
	l=len(mat[0])
	plb.figure(figsize=(10,10),dpi=80,facecolor='white')
	plb.suptitle("Scatter Matrix",fontsize=40)		# 设置图标标题
	for i in range(l):
		for j in range(l):
			plb.subplot(l,l,(l-i-1)*l+j+1)
			plb.scatter(mat[:,j],mat[:,i])
			if j==0:
				plb.ylabel(labels[i],fontsize=20)
			if i==0:
				plb.xlabel(labels[j],fontsize=20)
	plb.savefig("cheng.tif",dpi=200)  
	plb.show()
def drowScatter(mat,labels):
	x=range(len(mat))
	l=len(mat[0])
	plb.figure(figsize=(10,10),dpi=80,facecolor='white')
	plb.suptitle("Scatter",fontsize=40)		# 设置图标标题
	for i in range(l):
		plb.subplot(l,1,i+1)
		plb.scatter(x,mat[:,i])
		plb.ylabel(labels[i],fontsize=20)
	plb.savefig("cheng2.tif",dpi=200)  
	plb.show()
if __name__=="__main__":
	scale1=getScale(propertyID="Hydrophlicity")
	scale2=getScale(propertyID="Antigenicity")
	scale3=getScale(propertyID="Beta_turn")
	windowList=seqSplit(SEQ,7)
	matPro=profileMat(windowList,scale1,scale2,scale3)
	wl,accPro=profileAcc(SEQ,7)		
	wl,flexPro=profileFlex(SEQ,7)
	profileMat=np.zeros((len(matPro),5))	# 先定义一个每个元素都为0.的多为数组
	profileMat[:,0:3]=matPro
	profileMat[:,3]=accPro
	profileMat[:,4]=flexPro
	labels=["Hydrophlicity","Antigenicity","Beta_turn","Accessibility","Flexibility"]
	scatterMat(profileMat,labels)
	drowScatter(profileMat,labels)
