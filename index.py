from cifFileParser import CIFParser
import os
import numpy as np
from graph import Graph
from render import Render,ExportUnit

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

def main():
  folder="Tej"
  allFileNames=os.listdir(folder)
  renderModule=Render()
  lowerLimit=1.5
  upperLimit=1.7
  AnglePlotValues=[]
  distanceplotValues=[]
  invalidFiles=[]
  ExportData=[]
  maxProps=[0,"",""]
  minProps=[100,"",""]
  for file in allFileNames:
    sulphurs=[]
    nitrogens=[]

    parser=CIFParser(f"{folder}\{file}")
    if(not parser.validFile):
      invalidFiles.append(file)
      continue
    nitrogens=parser.getElementAtoms("N")
    sulphurs=parser.getElementAtoms("S")
    distanceValues={}
    for n in nitrogens:
      for s in sulphurs:
        distance=n.getDistance(s)
        if(distance>=lowerLimit and distance<=upperLimit):
          if(distance>maxProps[0]):
            maxProps[0]=distance
            maxProps[1]=file
            maxProps[2]=f"{n} -- {s}"
          if(distance<minProps[0]):
            minProps[0]=distance
            minProps[1]=file
            minProps[2]=f"{n} -- {s}"
          distanceValues[(n,s)]=distance

    print(f"S-N-S bonds for {file}")
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
    for key in occurences.keys():#Cycles over each atom in the occurence list
      left,right,center=None,None,None
      if(occurences[key]==2):
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
          SNSBonds.append(g)
          ExportData.append(ExportUnit(file,angle,[center,left],distanceValues[(center,left)],[center,right],distanceValues[(center,right)]))

    parser.printNicely(SNSBonds)
    for bond in SNSBonds:
      AnglePlotValues.append(bond.bondAngle)
      for j in bond.structure:
        if(j.symbol=="N"):
          for k in bond.structure[j]:
            distanceplotValues.append(k[1])
    print("="*20)
  

  # renderModule.plotHistogram(distanceplotValues,"S-N Distance","Distance (10^-10 m)","Frequency")
  # renderModule.plotHistogram(AnglePlotValues,"S-N-S Angle","Angle (deg)","Frequency")
  ExportUnit.Export(ExportData,maxProps,minProps)
  print(f"Invalid Files: {invalidFiles}")



main()


