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
# Calcutate the average scale of the "window"s
def windowCalc(window,scaleDict=getScale()):
	checkSeq(window)
	winSum=0.0;
	window=window.upper()
	for index,aa in enumerate(window):
		winSum+=float(scaleDict[aa])
	return winSum/len(window)
# Check if the amino acid sequnce include some worng symbols
def checkSeq(seq):
	stdSet=set("ACDEFGHIKLMNPQRSTVWY")
	seq=seq.upper()
	checkSet=set(seq)
	if not(checkSet.issubset(stdSet)):
		raise Exception("Wrong amino acid sequnce!!!")
# Split amino acid sequnce and return peptide list
def seqSplit(seq,windowSize):
	windowList=[]
	for i in range(len(seq)-windowSize+1):
		windowList.append(seq[i:windowSize+i])
	return windowList
# Get the scale matrix 
def scaleMat(windowList,*scaleDictTuple):
	row=len(windowList)
	column=len(scaleDictTuple)
	scaleMat=np.zeros((row,column))	# 先定义一个每个元素都为0.的多为数组
	for index in range(column):
		windowScaleList=[]
		for window in windowList:
			windowScaleList.append(windowCalc(window,scaleDictTuple[index]))
		scaleMat[:,index]=windowScaleList
	return scaleMat
#
def drowProfile(winScaleList,windowSize=1):
	x=range(len(winScaleList))
	x=[i+windowSize/2 for i in x]
	fig=plb.figure(1,facecolor='white')
	plb.plot(x,winScaleList)
	plb.show()
def drowScatter(mat):
	l=len(mat[0])
	plb.figure(figsize=(10,10),dpi=80,facecolor='white')
	for i in range(l):
		for j in range(l):
			plb.subplot(l,l,(l-i-1)*l+j+1)
			plb.scatter(mat[:,j],mat[:,i])
	plb.savefig("cheng.tif",dpi=200)  
	plb.show()
		

if __name__=="__main__":
	scale1=getScale(propertyID="Hydrophlicity")
	scale2=getScale(propertyID="Antigenicity")
	scale3=getScale(propertyID="Accessibility")
	scale4=getScale(propertyID="Beta_turn")
	windowList=seqSplit(SEQ,7)
	mat=scaleMat(windowList,scale1,scale2,scale3,scale4)
	drowScatter(mat)
