import pandas as pd
import matplotlib.pyplot as plt
from atom import Atom

textFile=open("something.txt")

alllines=textFile.read().split("\n")
eachLine=[]
i=0
while(i<len(alllines)):
  if("=" in alllines[i]):
    eachLine.append(f"{alllines[i]}{alllines[i+1]}")
    i+=1
  else:
    if(alllines[i].startswith("AFIX")):
      eachLine[-1]=f"{eachLine[-1]} {alllines[i]}"
    else:
      eachLine.append(alllines[i])
  i+=1


eachelement=[]

for j in eachLine:
  fullList=j.split(" ")
  temp=[]
  for k in fullList:
    if not (k == "" or k=="="):
      temp.append(k)
  
  eachelement.append(temp)


MyAtoms=[]

for j in eachelement:
  
  print(j)
  if(len(j)>0):
    MyAtoms.append(Atom(j))


for i in MyAtoms:
  print(i)







  

# ax.scatter()

# plt.show()