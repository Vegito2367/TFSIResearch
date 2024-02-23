import pandas as pd
import matplotlib.pyplot as plt
from atom import Atom
import os

class CIFParser:
  def printNicely(myList):
    for j in myList:
      print(j)

  def __init__(self,filePath):

    myfile=open(filePath,'r')

    alllines=myfile.read().split("\n")
    textofInterest=alllines[alllines.index("_atom_site_fract_z")+1:alllines.index("#END")]

    self.Atoms=[]
    for j in textofInterest:
      self.Atoms.append(Atom(j.split(" ")))

    self.nitrogens=next((i for i,atom in enumerate(self.Atoms) if atom.symbol == "N"),None)
    self.sulphurs=next((i for i,atom in enumerate(self.Atoms) if atom.symbol == "S"),None)
    print(next((i for i,atom in enumerate(self.Atoms) if atom.symbol == "H"),None))

  def getElementAtoms(self,symbol):
    output=[]
    for j in range(len(self.Atoms)):
      if(self.Atoms[j].symbol == symbol):
        output.append(j)

    return output


parser=CIFParser("testcif\AKONEU.cif")

print(parser.getElementAtoms("N"))