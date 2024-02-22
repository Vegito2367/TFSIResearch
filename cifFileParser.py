import pandas as pd
import matplotlib.pyplot as plt
from atom import Atom
import os

def printNicely(myList):
  for j in myList:
    print(j)

curDir=r"C:\Users\ASUS\Desktop\ProfJamesResources\testcif"
filename="AKONEU.cif"
allfiles=os.listdir(curDir)


myfile=open(f"{curDir}\{allfiles[0]}",'r')

alllines=myfile.read().split("\n")
textofInterest=alllines[alllines.index("_atom_site_fract_z")+1:alllines.index("#END")]

printNicely(textofInterest)