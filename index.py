from cifFileParser import CIFParser
import os


folder="cifstructures"
allFileNames=os.listdir(folder)
allNitrogens=[]
for file in allFileNames:
  parser=CIFParser(f"{folder}\{file}")
  atomsN=parser.getElementAtoms("N")
  allNitrogens.append(atomsN)
  atomsN=[]

CIFParser.printNicely(allNitrogens)
