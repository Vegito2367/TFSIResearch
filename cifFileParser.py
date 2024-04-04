import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from atom import Atom
import os

class CIFParser:
  @staticmethod
  def printNicely(myList):
    for j in myList:
      print(j)

  def __init__(self,filePath):

    myfile=open(filePath,'r')
    
    alllines=myfile.read().split("\n")
    startIndex=alllines.index("_atom_site_fract_z")+1
    self.cellvalues={}
    for i in range(startIndex): #Retrieveing cell properties
      if("_cell_length_a" in alllines[i]):
        val=(alllines[i].split(" ")[1])
        if("(" in val):
          val=val[:val.index("(")]
        self.cellvalues["cell_length_a"]=(float(val))

      if("_cell_length_b" in alllines[i]):
        val=(alllines[i].split(" ")[1])
        if("(" in val):
          val=val[:val.index("(")]
        self.cellvalues["cell_length_b"]=(float(val))

      if("_cell_length_c" in alllines[i]):
        val=(alllines[i].split(" ")[1])
        if("(" in val):
          val=val[:val.index("(")]
        self.cellvalues["cell_length_c"]=(float(val))

      if("_cell_angle_alpha" in alllines[i]):
        val=(alllines[i].split(" ")[1])
        if("(" in val):
          val=val[:val.index("(")]
        self.cellvalues["cell_angle_alpha"]=(float(val))

      if("_cell_angle_beta" in alllines[i]):
        val=(alllines[i].split(" ")[1])
        if("(" in val):
          val=val[:val.index("(")]
        self.cellvalues["cell_angle_beta"]=(float(val))

      if("_cell_angle_gamma" in alllines[i]):
        val=(alllines[i].split(" ")[1])
        if("(" in val):
          val=val[:val.index("(")]
        self.cellvalues["cell_angle_gamma"]=(float(val))
    
    textofInterest=alllines[startIndex:alllines.index("#END")]

    
    astarnum=(np.cos(np.radians(self.cellvalues["cell_angle_beta"])) * np.cos(np.radians(self.cellvalues["cell_angle_gamma"]))) - np.cos(np.radians(self.cellvalues["cell_angle_alpha"]))
    astardenom=np.cos(np.radians(self.cellvalues["cell_angle_beta"])) * np.sin(np.radians(self.cellvalues["cell_angle_gamma"]))
    self.cellvalues["cell_astar"]=np.arccos(astarnum/astardenom)
    ConversionMatrix=[
      [self.cellvalues["cell_length_a"],self.cellvalues["cell_length_b"]*np.cos(np.radians(self.cellvalues["cell_angle_gamma"])),self.cellvalues["cell_length_c"]*np.cos(np.radians(self.cellvalues["cell_angle_beta"]) )],
                      
      [0,self.cellvalues["cell_length_b"]*np.sin(np.radians(self.cellvalues["cell_angle_gamma"])),-1 * self.cellvalues["cell_length_c"] * np.sin(np.radians(self.cellvalues["cell_angle_beta"])) * np.cos(np.radians(self.cellvalues["cell_astar"]))],

      [0,0,self.cellvalues["cell_length_c"] * np.sin(np.radians(self.cellvalues["cell_angle_beta"])) * np.sin(np.radians(self.cellvalues["cell_astar"]))]
      
      ]
    self.ConversionMatrix=np.array(ConversionMatrix)
    self.Atoms=[]
    for j in textofInterest:
      self.Atoms.append(Atom(j.split(" "),self.ConversionMatrix))

  def getElementAtoms(self,symbol):
    output=[]
    for j in range(len(self.Atoms)):
      if(self.Atoms[j].symbol == symbol):
        output.append(self.Atoms[j])

    return output
