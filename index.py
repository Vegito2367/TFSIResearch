from cifFileParser import CIFParser
import os
def printNicely(myList):
    for j in myList:
      print(j)

folder="cifstructures"
allFileNames=os.listdir(folder)
allNitrogens=[]
for file in allFileNames:
  parser=CIFParser(f"{folder}\{file}")
  atomsN=parser.getElementAtoms("N")
  allNitrogens.append(atomsN)
  atomsN=[]

printNicely(allNitrogens)
