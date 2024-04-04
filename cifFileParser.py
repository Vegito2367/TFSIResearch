import pandas as pd
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
    
    for i in self.cellvalues:
      print(f"{i} : {self.cellvalues[i]}")
    textofInterest=alllines[startIndex:alllines.index("#END")]

    self.Atoms=[]
    for j in textofInterest:
      self.Atoms.append(Atom(j.split(" ")))

  def getElementAtoms(self,symbol):
    output=[]
    for j in range(len(self.Atoms)):
      if(self.Atoms[j].symbol == symbol):
        output.append(self.Atoms[j])

    return output
