from cifFileParser import CIFParser
import os

def main():
  folder="cifstructures"
  allFileNames=os.listdir(folder)
  # firstfile=allFileNames[0]
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
    for j in distanceValues:
      print(distanceValues[j],j)
    print("="*20)



main()


