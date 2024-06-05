from cifFileParser import CIFParser
import os
import numpy as np
from graph import Graph
from render import Render,ExportUnit
from heatmap import InteractivePlot

def magnitude(vector):
  mag=0
  for i in vector:
    mag+=(i**2)
  
  return mag**0.5

def getAngle(left,center,right):
  v1=[]
  v2=[]
  for i in range(3):
    v1.append(left[i]-center[i])
    v2.append(right[i]-center[i])

  v1=np.array(v1)
  v2=np.array(v2)
  dot=np.dot(v1,v2)
  cosTheta=dot/(magnitude(v1) * magnitude(v2))
  angle=np.arccos(cosTheta)
  return round(np.degrees(angle),6)

######################################################SNSManipulation

def IdentifySNSAngles(parser,lowerLimit,upperLimit,invalidFiles,distanceValues,ExportData):
      sulphurs=[]
      nitrogens=[]

      nitrogens=parser.getElementAtoms("N")
      sulphurs=parser.getElementAtoms("S")
      
      
      for n in nitrogens:
        for s in sulphurs:
          distance=n.getDistance(s)

          if(distance>=lowerLimit and distance<=upperLimit):
            distanceValues[(n,s)]=distance
      
      if(len(distanceValues)==0):
        invalidFiles.append([parser.fileName,"No S-N bonds found"])
        
      # print(f"S-N-S bonds for {file}")
      occurences={}
      for (i,j) in distanceValues:
        if(i in occurences.keys()):
          occurences[i]+=1
        else:
          occurences[i]=1
        
        if(j in occurences.keys()):
          occurences[j]+=1
        else:
          occurences[j]=1
        
      SNSBonds=[]
      for key in occurences.keys():#Cycles over each nitrogen atom in the occurence list
        left,right,center=None,None,None
        if(occurences[key]>=2 and key.symbol=="N"):
            center=key
            bonds=list(distanceValues.keys())
            for n,s in bonds:
              if(n==center):
                if(left is None):
                  left=s
                elif(right is None):
                  right=s

            angle = getAngle(left.positionVector,center.positionVector,right.positionVector)
            g=Graph([left,center,right],angle)
            #Separate the Export stuff from the SNS Bonds Identification
            leftd=distanceValues[(g.center,g.left)]
            rightd=distanceValues[(g.center,g.right)]
            # Extremes(maxProps,minProps,parser.fileName,g.center.symbol,g.left.symbol,leftd)
            # Extremes(maxProps,minProps,parser.fileName,g.center.symbol,g.left.symbol,rightd)
            ExportData.append(ExportUnit(parser.fileName,g.bondAngle,[center,left],leftd,[center,right],rightd))
            SNSBonds.append(g)
      return SNSBonds


def SNSManipulation():
  folder="Tej2"
  allFileNames=os.listdir(folder)
  renderModule=Render()
  lowerLimit=1.5
  upperLimit=1.7
  AnglePlotValues=[]
  distanceDiff=[]
  distanceplotValues=[]
  invalidFiles=[]
  ExportData=[]
  maxProps=[0,"",""]
  minProps=[100,"",""]
  progress=0
  
  for file in allFileNames:
    print(f"Progress: {progress}/{len(allFileNames)}")
    try:
        SNSBonds=[]
        parser=CIFParser(f"{folder}\{file}")
        if(not parser.validFile):
          invalidFiles.append([file,"Invalid Input"])
          continue
        distanceValues={}
        SNSBonds=IdentifySNSAngles(parser,lowerLimit,upperLimit,invalidFiles,distanceValues,ExportData)
        
        
        for bond in SNSBonds:
          AnglePlotValues.append(bond.bondAngle)
          for j in bond.structure: #For each atom in the bond
            temp=[]
            if(j.symbol=="N"):
              for k in bond.structure[j]:
                distanceplotValues.append(k[1])
                temp.append(k[1])
              distanceDiff.append((temp[0]+temp[1])/2)

    except Exception as e:
      print(distanceValues)
      invalidFiles.append(file)
      raise e
    progress+=1
  

  # renderModule.plotHistogram(distanceplotValues,"S-N Distance","Distance (10^-10 m)","Frequency")
  # renderModule.plotHistogram(AnglePlotValues,"S-N-S Angle","Angle (deg)","Frequency")
  # renderModule.scatterPlot(AnglePlotValues,distanceDiff,"S-N-S Angle vs S-N-S Distance Avg ","Angle (deg)","Difference (10^-10 m)")
  x,y,names=[],[],[]
  for unit in ExportData:
      x.append(unit.angle)
      y.append((unit.distance1+unit.distance2)/2)
      names.append([unit.file,unit.bond1,unit.bond2])

  # x1,names1=[],[]
  # for unit in ExportData:
  #   x1.append(unit.angle)
  #   names1.append(unit.file)
  
  InteractivePlot.plotInteractivePlot(x,y,names,"Angle (deg)","Average Distance(10^-10 m)",f"S-N-S Angle vs S-N-S Distance Avg for {folder}")
  #InteractivePlot.InteractiveHistogram(x1,names1,"Angle (deg)","S-N-S Angle Frequency")
  #ExportUnit.Export(ExportData,maxProps,minProps)
  print(f"Invalid Files:")
  CIFParser.printNicely(invalidFiles)
  print(len(invalidFiles))

def Extremes(maxProps, minProps, file, n, s, distance):
    if(distance>maxProps[0]):
      maxProps[0]=distance
      maxProps[1]=file
      maxProps[2]=f"{n} -- {s}"
    if(distance<minProps[0]):
      minProps[0]=distance
      minProps[1]=file
      minProps[2]=f"{n} -- {s}"
######################################################End of SNSManipulation/Start of Torison Angle
def TorisonAngle():
  folder="testcif"
  allFileNames=os.listdir(folder)
  renderModule=Render()
  AnglePlotValues=[]
  invalidFiles=[]
  progress=0
  ExportDataTorsion=[]
  for file in allFileNames:
    print(f"Progress: {progress}/{len(allFileNames)}")
    try:
        parser=CIFParser(f"{folder}\{file}")
        if(not parser.validFile):
          invalidFiles.append([file,"Invalid Input"])
          progress+=1
          continue
        distanceValues={}
        SNSBonds=IdentifySNSAngles(parser,1.5,1.7,invalidFiles,distanceValues,ExportDataTorsion)
        for bond in SNSBonds:
          carbons = parser.getElementAtoms("C")
          cDistLeft,cDistRight=100,100
          cLeft,cRight=None,None
          for c in carbons:
            distL=bond.left.getDistance(c)
            distR=bond.right.getDistance(c)
            if(distL<cDistLeft):
              cDistLeft=distL
              cLeft=c
            if(distR<cDistRight):
              cDistRight=distR
              cRight=c
          if(cLeft is not None and cRight is not None):
            bond.addBond(bond.left,cLeft,cDistLeft)
            bond.addBond(bond.right,cRight,cDistRight)
            print(bond)
          print("=====================================")
    except Exception as e:
      invalidFiles.append(file)
      print(file)
      raise e
    progress+=1
        
        
######################################################End of Torison Angle/Start of Main
def main():
  TorisonAngle()
main()


