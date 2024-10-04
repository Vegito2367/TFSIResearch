from cifFileParser import CIFParser
import os
import numpy as np
from graph import Graph
from render import Render,ExportUnit
from heatmap import InteractivePlot
import pandas as pd
from collections import Counter

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

def getNormalVector(x,y,z):
  pass
def getTorsionAngle(A,B,C,D):
  #BC is the common edge
  #ABC is first plane and BCD is the second plane
  normABC=np.cross(A.positionVector-B.positionVector,C.positionVector-B.positionVector)
  normBCD=np.cross(B.positionVector-C.positionVector,D.positionVector-C.positionVector)
  return np.degrees(np.arccos(np.dot(normABC,normBCD)/(magnitude(normABC)*magnitude(normBCD))))


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


# def SNSManipulation():
#   folder="TFSI_NoDisorder"
#   allFileNames=os.listdir(folder)
#   renderModule=Render()
#   lowerLimit=1.5
#   upperLimit=1.7
#   AnglePlotValues=[]
#   distanceDiff=[]
#   distanceplotAvg=[]
#   checkArray=[]
#   invalidFiles=[]
#   ExportData=[]
#   maxProps=[0,"",""]
#   minProps=[100,"",""]
#   progress=0
#   isCobaltPresent=[]
#   atomToLookFor=["Fe","Co","Au","Mn"]
#   fudgeFactor=0.5
#   for file in allFileNames:
#     print(f"Progress: {progress}/{len(allFileNames)}")
#     try:
#         SNSBonds=[]
#         parser=CIFParser(f"{folder}\{file}")
#         if(not parser.validFile):
#           invalidFiles.append([file,"Invalid Input"])
#           continue
#         distanceValues={}
#         SNSBonds=IdentifySNSAngles(parser,lowerLimit,upperLimit,invalidFiles,distanceValues,ExportData)
        
        
#         for bond in SNSBonds:
          
#           AnglePlotValues.append(bond.bondAngle)
          
#           temp=[]
#           for j,dist in bond.structure[bond.center]:
#             temp.append(dist)
          
#           distanceDiff.append(abs(temp[0]-temp[1]))
#           distanceplotAvg.append(sum(temp)/2)
#           surroundingAtoms = parser.getAtomsInARadius(bond.center,10)
#           cobalt=False
#           for atom,distance in surroundingAtoms:
#             if(atom.symbol in atomToLookFor):
#               if(distance<=(bond.center.covalentRadius+ atom.covalentRadius+fudgeFactor) and not cobalt):
#                 isCobaltPresent.append(1)
#                 cobalt=True
          
#           if(not cobalt):
#             isCobaltPresent.append(0)
          

#     except Exception as e:
#       print(distanceValues)
#       invalidFiles.append(file)
#       raise e
#     progress+=1
  

  
#   x,y,names=[],[],[]
#   for unit in ExportData:
#       x.append(unit.singleProperty)
#       y.append((unit.doubleProperty1+unit.doubleProperty1)/2)
#       names.append([unit.file,unit.atoms1,unit.atoms2])


  
#   # InteractivePlot.plotInteractivePlotColorArray(AnglePlotValues,distanceplotAvg,names,isCobaltPresent,"S-N-S Angle (deg)", "S-N Angle Average (Angstrom)",f"Angle vs Distance Avg for {folder} with {atomToLookFor} presence", False)
#   InteractivePlot.plotInteractivePlotColorArray(AnglePlotValues,distanceplotAvg,names,isCobaltPresent,"S-N-S Angle (deg)", "S-N Angle Average (Angstrom)",f"Angle vs Distance Avg for {folder} with {atomToLookFor} presence", False)
#   #InteractivePlot.InteractiveHistogram(x1,names1,"Angle (deg)","S-N-S Angle Frequency")
#   #ExportUnit.Export(ExportData,maxProps,minProps)
#   print(f"Invalid Files:")
#   CIFParser.printNicely(invalidFiles)
#   print(len(invalidFiles))

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
  folder="LatestTFSI_NoDisorder"
  allFileNames=os.listdir(folder)
  renderModule=Render()
  AnglePlotValues=[]
  invalidFiles=[]
  progress=0
  ExportDataTorsion=[]
  TorisonAngleAverages=[]
  names=[]
  bondlengthDeltas=[]
  bondlengthAverage=[]
  torsionDeltas=[]
  fudgeFactorMetalsBound={}
  testParse=CIFParser(f"{folder}\{allFileNames[0]}")
  atomToLookFor=list(testParse.covalentRadii.keys())
  atomDictionaryList={} #Stores a list of atoms that are present in the compound for each atom
  for atom in atomToLookFor:
    atomDictionaryList[atom]=[]
  fudgeFactor = [round(i * 0.1,1) for i in range(31)]
  currentFudgeFactor=1
  for factor in fudgeFactor:
    fudgeFactorMetalsBound[factor]=0
    

  totalSNS=0
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
        totalSNS+=len(SNSBonds)
        for bond in SNSBonds:
          temp=[]
          for s,dist in bond.structure[bond.center]:
            temp.append(dist)
          
          
          bondlengthAverage.append((temp[0]+temp[1])/2)
          bondlengthDeltas.append(abs(temp[0]-temp[1]))
          for atomSymbol in atomToLookFor:
            surroundingAtoms = parser.getAtomsInARadius(bond.center,3)
            cobalt=False
            for atom,distance in surroundingAtoms:
              if(atom.symbol==atomSymbol):#Checks if the given atom is one of the metals we are looking for
                if(not cobalt):

                  for factor in fudgeFactor:
                    if(distance<=(bond.center.covalentRadius + atom.covalentRadius+factor)):
                      fudgeFactorMetalsBound[factor]+=1

                  if(distance<=(bond.center.covalentRadius + atom.covalentRadius+currentFudgeFactor)):
                    atomDictionaryList[atomSymbol].append("Metal present and N-bound to TFSI")  #The atom symbol is counted as a metal connected to the TFSI
                    cobalt=True
            
            if(not cobalt):
              if(parser.containsAtom(atomSymbol)):
                atomDictionaryList[atomSymbol].append("Metal present but not N–bound to TFSI")  #The atom symbol is present in the compound but not connected to the TFSI
              else:
                atomDictionaryList[atomSymbol].append("No metal present in the structure")

          AnglePlotValues.append(bond.bondAngle)
          # carbons = parser.getElementAtoms("C")
          # cDistLeft,cDistRight=100,100
          # cLeft,cRight=None,None
          # for c in carbons:
          #   distL=bond.left.getDistance(c)
          #   distR=bond.right.getDistance(c)
          #   if(distL<cDistLeft):
          #     cDistLeft=distL
          #     cLeft=c
          #   if(distR<cDistRight):
          #     cDistRight=distR
          #     cRight=c
          # if(cLeft is not None and cRight is not None):
          #   bond.addBond(bond.left,cLeft,cDistLeft)
          #   bond.addBond(bond.right,cRight,cDistRight)
            
          
          # #Right Torsion Angle
          # A1,B1,C1,D1 = bond.left,bond.center,bond.right,cRight
          # TorisonAngleRight=getTorsionAngle(A1,B1,C1,D1)
          # TorisonAngleLeft=getTorsionAngle(C1,B1,A1,cLeft)
          
          # TorisonAngleAverages.append((TorisonAngleRight+TorisonAngleLeft)/2)
          # torsionDeltas.append(abs(TorisonAngleRight-TorisonAngleLeft))
          names.append(f"{file} {bond.left} {bond.center} {bond.right}")
          # ExportDataTorsion.append(ExportUnit(file,bond.bondAngle,[A1,B1,C1,D1],TorisonAngleRight,[C1,B1,A1,cLeft],TorisonAngleLeft))
  
    except Exception as e:
      invalidFiles.append(file)
      print(file)
      raise e
    
    progress+=1
  
  GenerateAllGraphs(folder, AnglePlotValues, names, bondlengthAverage, atomDictionaryList)
  # GenerateAllGraphs(folder, AnglePlotValues, names, bondlengthAverage, {"Au":atomDictionaryList["Au"],"Ag":atomDictionaryList["Ag"]})
  totalLen=[]
  metalBoundCount=0
  metalPresent=0
  MetalPresenceCount=0
  for atomSymbol in atomDictionaryList:
    metalCount=atomDictionaryList[atomSymbol].count("Metal present and N-bound to TFSI")
    metalPresent=atomDictionaryList[atomSymbol].count("Metal present but not N–bound to TFSI")
    totalLen.append([metalCount,atomSymbol])
    metalBoundCount+=metalCount
    MetalPresenceCount+=(metalCount+metalPresent)
    print(f"Percentage presence of {atomSymbol} {(metalCount/len(atomDictionaryList[atomSymbol])*100):0.6f}")
    temp.append((atomDictionaryList[atomSymbol].count(atomSymbol)/len(atomDictionaryList[atomSymbol])*100))
  
  #Each list in atomDictionaryList has length = total number of compounds
  histX=[]
  histY=[]
  for num,symbol in totalLen:
    histX.append(symbol)
    histY.append(num)
  print(totalLen)
  print("Metal Presence: ",MetalPresenceCount)
  print("Metal Bound: ",metalBoundCount)
  print("Total SNS Bonds: ",totalSNS)
  print("Number of Avgs",len(bondlengthAverage))
  val = (metalBoundCount/totalSNS)*100
  print(f"The number of metals for fudge Factor 1 = {fudgeFactorMetalsBound[currentFudgeFactor]}")
  renderModule.barGraphFrequencies(histX,histY,totalSNS,"Metal Presence in Compound","Element","Frequency","Percentage")
  # renderModule.barGraphFrequencies(["Metal bound to\nnitrogen in structure","Metal not bound\nbut present in structure"],[metalBoundCount,MetalPresenceCount],totalSNS,"","","Frequency","Percentage of all\nmetal-containing structures")
  #renderModule.barGraph(["Metal-containing","Metal-free"],[MetalPresenceCount,totalSNS-MetalPresenceCount],"","","Frequency")
  # renderModule.plotHistogram(AnglePlotValues,"SNS Angle Spread","Angle (deg)","Frequency")
  # renderModule.plotHistogram(bondlengthAverage,"SN Distance Spread","SN Distance (Ang)","Frequency")
  # renderModule.plotLine(list(fudgeFactorMetalsBound.keys()),list(fudgeFactorMetalsBound.values()),"Sensitivity Analysis the Fudge Factor","Fudge Factor","Frequency of N bound metals")



def GenerateAllGraphs(folder, AnglePlotValues, names, bondlengthAverage, atomDictionaryList):
    for atom in atomDictionaryList:
      InteractivePlot.plotInteractivePlot(AnglePlotValues, bondlengthAverage,names,atomDictionaryList[atom],"SNS Angle (°)","Bond Length Avg (Å)",f"S—N—S angle vs. average S—N bond length for structures containing {atom}", True, f"WithCaption_{folder}_{atom}_SNS_Angle_BondLength")        
        
######################################################End of Torison Angle
def getDataFromExcelFiles():
  # Define variable to load the dataframe
  dataframe = pd.read_excel("MetalOxygenSNSAngles.xlsx")
  SNSAngles=list(dataframe.get("Angle"))
  RefCodes=list(dataframe.get("Refcode"))
  
  bondLenDataFrame=pd.read_excel("MetalOxygenSNSBonds.xlsx")
  SNSAvgExcel=list(bondLenDataFrame.get("AVG"))
  RefCodesBond=list(bondLenDataFrame.get("Refcode"))
  SNSAvgFinal=[]
  FinalRefcodes=[]

  i=1
  prevRef=RefCodesBond[0]
  prevAvg=SNSAvgExcel[0]
  FinalRefcodes.append(RefCodesBond[0])
  SNSAvgFinal.append(SNSAvgExcel[0])
  while(i<len(SNSAvgExcel)):
    if(RefCodesBond[i]==prevRef and prevAvg==SNSAvgExcel[i]):
      #print(RefCodesBond[i],SNSAvgExcel[i])
      i+=1
      continue
    else:
      SNSAvgFinal.append(SNSAvgExcel[i])
      FinalRefcodes.append(RefCodesBond[i])
      prevRef=RefCodesBond[i]
      prevAvg=SNSAvgExcel[i]
      i+=1
  
  print(len(SNSAvgFinal),len(SNSAngles))

  counterAvg=Counter(FinalRefcodes)
  counterAngle=Counter(RefCodes)
  for key in counterAvg:
    if not (counterAvg[key]==counterAngle[key]):
      print(key,counterAvg[key],counterAngle[key])

    
  
  

  

def main():
  TorisonAngle()
  
main()


