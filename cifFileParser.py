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

    
    astarnum=self.alphaStarNumerator()
    astardenom=self.alphaStarDenominator()

    self.cellvalues["cell_astar"]=np.arccos(astarnum/astardenom)

    ConversionMatrix=[
      [self.a(),self.b() * self.cos(self.gamma()),self.c() * self.cos(self.beta())],
                      
      [0, self.b() * self.sin(self.gamma()), -1 * self.c() * self.sin(self.beta()) * self.cos(self.cellvalues["cell_astar"])],

      [0,0,self.c() * self.sin(self.beta()) * self.sin(self.cellvalues["cell_astar"])]
      
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
  
  def cos(self,x):
    return np.cos(x)
  
  def sin(self,x):
    return np.sin(x)
  
  
  def alphaStarNumerator(self):
    return (self.cos(self.beta()) * self.cos(self.gamma())) - self.cos(self.alpha())
  
  def alphaStarDenominator(self):
    return self.sin(self.beta()) * self.sin(self.gamma())
  

  def b(self):
    return self.cellvalues["cell_length_b"]
  
  def a(self):
    return self.cellvalues["cell_length_a"]
  
  def c(self):
    return self.cellvalues["cell_length_c"]
  
  def gamma(self):
    return np.radians(self.cellvalues["cell_angle_gamma"])
  
  def beta(self):
    return np.radians(self.cellvalues["cell_angle_beta"])
  
  def alpha(self):
    return np.radians(self.cellvalues["cell_angle_alpha"])
  


  
