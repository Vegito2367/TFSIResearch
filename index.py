from cifFileParser import CIFParser
import os
import numpy as np

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

  print("=======")
  v1=np.array(v1)
  v2=np.array(v2)
  dot=np.dot(v1,v2)
  cosTheta=dot/(magnitude(v1) * magnitude(v2))
  angle=np.arccos(cosTheta)
  return round(np.degrees(angle),6)

def main():
  folder="testcif"
  allFileNames=os.listdir(folder)

  lowerLimit=1.5
  upperLimit=1.7

  for file in allFileNames:
    sulphurs=[]
    nitrogens=[]

    parser=CIFParser(f"{folder}\{file}")
    nitrogens=parser.getElementAtoms("N")
    sulphurs=parser.getElementAtoms("S")
    distanceValues={}
    for n in nitrogens:
      for s in sulphurs:
        distance=n.getDistance(s)
        if(distance>=lowerLimit and distance<=upperLimit):
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
    left,center,right=None,None,None
    for key in occurences.keys():
      if(occurences[key]==2):
          center=key
      else:
        if(left is None):
          left=key
        else:
          right=key
    
    print(occurences)
    print("=============Test==========")
    print(getAngle([0,1,0],[0,0,0],[1,1,0]))
    print("===============")
    print("left", left.positionVector)
    print("center", center.positionVector)
    print("right", right.positionVector)
    print(getAngle(left.positionVector,center.positionVector,right.positionVector))

    print("="*20)



main()


