from cifFileParser import CIFParser
import os


folder="cifstructures"
allFileNames=os.listdir(folder)
nitrogens=[]
firstfile=allFileNames[0]
# for file in allFileNames:
#   parser=CIFParser(f"{folder}\{file}")
#   atomsN=parser.getElementAtoms("N")
#   allNitrogens.append(atomsN)
sulphurs=[]
distanceLimit=0.3
parser=CIFParser(f"{folder}\{firstfile}")
nitrogens=parser.getElementAtoms("N")
sulphurs=parser.getElementAtoms("S")
distanceValues={}
for n in nitrogens:
  for s in sulphurs:
    distance=n.getDistance(s)
    if(distance<distanceLimit):
      distanceValues[(n,s)]=distance


for j in distanceValues:
  print(distanceValues[j],j)
