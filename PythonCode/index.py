from cifFileParser import CIFParser
import os
import numpy as np
from graph import Graph
from render import Render,ExportUnit
from interactiveGraphs import InteractivePlot
import pandas as pd
from collections import Counter


'''
Vector math utilities
'''
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

def getTorsionAngle(A,B,C,D):
  #BC is the common edge
  #ABC is first plane and BCD is the second plane
  normABC=np.cross(A.positionVector-B.positionVector,C.positionVector-B.positionVector)
  normBCD=np.cross(B.positionVector-C.positionVector,D.positionVector-C.positionVector)
  return np.degrees(np.arccos(np.dot(normABC,normBCD)/(magnitude(normABC)*magnitude(normBCD))))


######################################################SNSManipulation

def IdentifySNSBonds(parser,lowerLimit,upperLimit,invalidFiles):
      sulphurs=[]
      nitrogens=[]
      distanceValues={}
      occurences={}
      nitrogens=parser.getElementAtoms("N")
      sulphurs=parser.getElementAtoms("S")
      
      
      for n in nitrogens:
        for s in sulphurs:
          distance=n.getDistance(s)

          if(distance<=(n.covalentRadius+s.covalentRadius+1)):
            distanceValues[(n,s)]=distance
            if(n in occurences):
              occurences[n]+=1
            else:
              occurences[n]=1
      
      if(len(distanceValues)==0):
        invalidFiles.append([parser.fileName,"No S-N bonds found"])

        
      SNSBonds=[]
      for key in occurences:#Cycles over each nitrogen atom in the occurence list
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
            #Values used to create ExportUnit Objects that were used to create the excel sheet during the early days of the project
            # leftd=distanceValues[(g.center,g.left)]
            # rightd=distanceValues[(g.center,g.right)]
            
            #ExportData.append(ExportUnit(parser.fileName,g.bondAngle,[center,left],leftd,[center,right],rightd))
            SNSBonds.append(g)
      return SNSBonds



######################################################  End of SNSManipulation/Start of Metal Binding Analysis
def MetalBinding():
  folder="TFSI_NoDisorder" # Folder name containing all the CIF files
  datasetPath=os.listdir(folder)
  renderModule=Render()
  AnglePlotValues=[]
  invalidFiles=[]
  progress=0
  totalSNS=0
  

  names=[]
  bondlengthDeltas=[]
  bondlengthAverage=[]
  
  fudgeFactorMetalsBound={}
  testParse=CIFParser(os.path.join(folder,datasetPath[0]))
  atomToLookFor=list(testParse.covalentRadii.keys()) # Retrieves the list of colavent radii for each atom
  atomOccurences={} #Stores a list of atoms that are present in the compound for each atom
  fudgeFactor = [round(i * 0.1,1) for i in range(31)] #Cycles through various fudge factors from 0 till 3 angstrom for sensitivity analysis
  currentFudgeFactor=1 #Standard fudge factor used for the analysis in angstrom (10^-10m)

  for atom in atomToLookFor:
    atomOccurences[atom]=[]
  
  for factor in fudgeFactor:
    fudgeFactorMetalsBound[factor]=0
  
  
  #Variables for Torsion Angle block of code
  #ExportDataTorsion=[]
  TorisonAngleAverages=[]
  torsionDeltas=[]
    

  
  for file in datasetPath:
    
    print(f"Progress: {progress}/{len(datasetPath)}")
    try:
        parser=CIFParser(f"{folder}/{file}")
        if(not parser.validFile):
          invalidFiles.append([file,"Invalid Input"])
          progress+=1
          continue

        SNSBonds=IdentifySNSBonds(parser,1.5,1.7,invalidFiles)
        totalSNS+=len(SNSBonds)
        for bond in SNSBonds:
          temp=[]
          for s,dist in bond.structure[bond.center]:
            temp.append(dist)
          
          
          bondlengthAverage.append((temp[0]+temp[1])/2)
          bondlengthDeltas.append(abs(temp[0]-temp[1]))
          AnglePlotValues.append(bond.bondAngle)

          for atomSymbol in atomToLookFor:
            surroundingAtoms = parser.getAtomsInARadius(bond.center,3)
            metalBound=False
            for atom,distance in surroundingAtoms:
              if(atom.symbol==atomSymbol):#Checks if the given atom is one of the metals we are looking for
                if(not metalBound):

                  for factor in fudgeFactor:
                    if(distance<=(bond.center.covalentRadius + atom.covalentRadius+factor)):
                      fudgeFactorMetalsBound[factor]+=1

                  if(distance<=(bond.center.covalentRadius + atom.covalentRadius+currentFudgeFactor)):
                    atomOccurences[atomSymbol].append("Metal present and N-bound to TFSI")  #The atom symbol is counted as a metal connected to the TFSI
                    metalBound=True
            
            if(not metalBound):
              if(parser.containsAtom(atomSymbol)):
                atomOccurences[atomSymbol].append("Metal present but not N–bound to TFSI")  #The atom symbol is present in the compound but not connected to the TFSI
              else:
                atomOccurences[atomSymbol].append("No metal present in the structure")

          
          '''
          Below code calculates torsion angles and adds them to graph data
          '''
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
          #ExportDataTorsion.append(ExportUnit(file,bond.bondAngle,[A1,B1,C1,D1],TorisonAngleRight,[C1,B1,A1,cLeft],TorisonAngleLeft))
          #names.append(f"{file} {bond.left} {bond.center} {bond.right} {TorsionAngleAverages[-1]}") - Alternative names append
          
          names.append(f"{file} {bond.left} {bond.center} {bond.right}")
          
          
    except Exception as e:
      invalidFiles.append(file)
      print(file)
      raise e
    
    progress+=1
  
  #GenerateAllGraphs(folder, AnglePlotValues, names, bondlengthAverage, atomOccurences)
  #GenerateTorsionAngles(folder, TorisonAngleAverages, names, bondlengthAverage, atomOccurences)
  #GenerateAllGraphs(folder, AnglePlotValues, names, bondlengthAverage, {"Au":atomOccurences["Au"]})
  totalLen=[]
  metalBoundCount=0
  metalPresent=0
  MetalPresenceCount=0
  for atomSymbol in atomOccurences:
    metalCount=atomOccurences[atomSymbol].count("Metal present and N-bound to TFSI")
    metalPresent=atomOccurences[atomSymbol].count("Metal present but not N–bound to TFSI")
    totalLen.append([metalCount,atomSymbol])
    metalBoundCount+=metalCount
    MetalPresenceCount+=(metalCount+metalPresent)
    print(f"Percentage presence of {atomSymbol} {(metalCount/len(atomOccurences[atomSymbol])*100):0.6f}")
    temp.append((atomOccurences[atomSymbol].count(atomSymbol)/len(atomOccurences[atomSymbol])*100))
  
  #Each list in atomOccurences has length = total number of compounds
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
  '''
  Below code uses Render() class to create matplotlib graphs
  '''

  # renderModule.barGraphFrequencies(histX,histY,totalSNS,"Metal Presence in Compound","Element","Frequency","Percentage",20)
  # renderModule.barGraphFrequencies(["Metal bound to\nnitrogen in structure","Metal not bound\nbut present in structure"],[metalBoundCount,MetalPresenceCount],totalSNS,"","","Frequency","Percentage of all\nmetal-containing structures",40)
  # renderModule.barGraph(["Metal-containing","Metal-free"],[MetalPresenceCount,totalSNS-MetalPresenceCount],"","","Frequency")
  # renderModule.plotHistogram(AnglePlotValues,"SNS Angle Spread","Angle (deg)","Frequency")
  # renderModule.plotHistogram(bondlengthAverage,"SN Distance Spread","SN Distance (Ang)","Frequency")
  # renderModule.plotLine(list(fudgeFactorMetalsBound.keys()),list(fudgeFactorMetalsBound.values()),"Sensitivity Analysis the Fudge Factor","Fudge Factor","Frequency of N bound metals")
  ExportUnit.exportPairValues(AnglePlotValues,bondlengthAverage,atomOccurences)

def GenerateAllGraphs(folder, AnglePlotValues, names, bondlengthAverage, atomOccurences):
    for atom in atomOccurences:
      InteractivePlot.plotInteractivePlot(AnglePlotValues, bondlengthAverage,names,atomOccurences[atom],"SNS Angle (°)","Bond Length Avg (Å)",f"S—N—S angle vs. average S—N bond length for structures containing {atom}", False, f"new-{atom}-graph")

def GenerateTorsionAngles(folder, TorsionAngles, names, bondlengthAverage, atomOccurences):
    for atom in atomOccurences:
      InteractivePlot.plotInteractivePlot(TorsionAngles, bondlengthAverage,names,atomOccurences[atom],"Torsion Angle (°)","Bond Length Avg (Å)",f"Torsion angle average vs. average S—N bond length for structures containing {atom}", False, f"new-{atom}-graph")         
        
def GenerateExcelSheets(folder, AnglePlotValues, names, bondlengthAverage, atomOccurences):
    for atom in atomOccurences:
      InteractivePlot.plotInteractivePlot(AnglePlotValues, bondlengthAverage,names,atomOccurences[atom],"SNS Angle (°)","Bond Length Avg (Å)",f"S—N—S angle vs. average S—N bond length for structures containing {atom}", False, f"new-{atom}-graph")

######################################################End of Torison Angle
'''
Below code extracts data for several different moeities of TFSI from an excel sheet input
'''
# def getDataFromExcelFiles():
#   # Define variable to load the dataframe
#   dataframe = pd.read_excel("MetalOxygenSNSAngles.xlsx")
#   SNSAngles=list(dataframe.get("Angle"))
#   RefCodes=list(dataframe.get("Refcode"))
  
#   bondLenDataFrame=pd.read_excel("MetalOxygenSNSBonds.xlsx")
#   SNSAvgExcel=list(bondLenDataFrame.get("AVG"))
#   RefCodesBond=list(bondLenDataFrame.get("Refcode"))
#   SNSAvgFinal=[]
#   FinalRefcodes=[]

#   i=1
#   prevRef=RefCodesBond[0]
#   prevAvg=SNSAvgExcel[0]
#   FinalRefcodes.append(RefCodesBond[0])
#   SNSAvgFinal.append(SNSAvgExcel[0])
#   while(i<len(SNSAvgExcel)):
#     if(RefCodesBond[i]==prevRef and prevAvg==SNSAvgExcel[i]):
#       #print(RefCodesBond[i],SNSAvgExcel[i])
#       i+=1
#       continue
#     else:
#       SNSAvgFinal.append(SNSAvgExcel[i])
#       FinalRefcodes.append(RefCodesBond[i])
#       prevRef=RefCodesBond[i]
#       prevAvg=SNSAvgExcel[i]
#       i+=1
  
#   print(len(SNSAvgFinal),len(SNSAngles))

#   counterAvg=Counter(FinalRefcodes)
#   counterAngle=Counter(RefCodes)
#   for key in counterAvg:
#     if not (counterAvg[key]==counterAngle[key]):
#       print(key,counterAvg[key],counterAngle[key])


def main():
  MetalBinding()
  
main()


