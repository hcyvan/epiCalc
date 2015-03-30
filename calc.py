#coding=utf-8
import xml.dom.minidom
# node为DOM的叶节点
def getValue(leafNode):
	value=leafNode.childNodes[0].nodeValue
	return str(value)
# 获取info.xml中的scale，以词典的形式返回
def getScale(filename="info.xml",propertyID="Hydrophlicity",subpropertyID="0"):
	scaleDict={}
	doc=xml.dom.minidom.parse(filename) 		#加载xml文件
	root=doc.documentElement			#获取xml文档对象	
	nodes=root.getElementsByTagName("property")	#获取property节点集合
	for node in nodes:
		if str(node.getAttribute("id"))==propertyID:#获得每个节点的id,并转换为string
			subNodes=node.getElementsByTagName("subproperty") #获取subproperty节点集合
			for subNode in subNodes:
				if str(subNode.getAttribute("id"))==subpropertyID:
					scaleNode=subNode.getElementsByTagName("scale")[0] #获得scale节点
					lis=scaleNode.getElementsByTagName("li")#获取subproperty节点集合
					for li in lis:
						aa=getValue(li.getElementsByTagName("aa")[0])
						value=getValue(li.getElementsByTagName("value")[0])
						scaleDict[aa]=value
	return scaleDict

if __name__=="__main__":
	scale=getScale()
	for aa in scale.keys():
		print aa+" "+scale[aa]
