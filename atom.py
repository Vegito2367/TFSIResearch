class Atom:
  def __init__(self, inputList):



    elemname=inputList.pop(0)

    ind1,ind2=0,0
    for i in range(len(elemname)):
      if(str.isdecimal(elemname[i])):
        ind1=i
        break
    for i in range(len(elemname)-1,-1,-1):
      if(str.isdecimal(elemname[i])):
        ind2=i
        break
    
    if(ind1==0 and ind2==0):
      self.atomLetter=""
      self.atomNum=1
      self.element=elemname
    else:
      self.atomLetter=elemname[ind2+1:]
      self.element=elemname[:ind1]
      self.atomNum=elemname[ind1:ind2+1]

    # if(len(elemname)>1 and str.isalpha(elemname[-1]) and str.isdigit(elemname[-2])):
    #   self.atomLetter=elemname[-1]
    #   ind1,ind2=0,0
    #   for i in range(len(elemname)):
    #     if(str.isdecimal(elemname[i])):
    #       ind=i
    #       break
    #   for i in range(len(elemname)-1,-1,-1):
    #     if(str.isdecimal(elemname[i])):
    #       ind=i
    #       break
    #   self.element=elemname[0:ind]
    #   self.atomNum=elemname[ind:-1]
    # elif (str.isdigit(elemname[-1])):
    #   ind=int(0)
    #   for i in range(len(elemname)):
    #     if(str.isdecimal(elemname[i])):
    #       ind=i
    #       break
    #   self.atomLetter=""
    #   self.atomNum=elemname[ind:]
    #   self.element=elemname[0:ind]
    # else:
    #   self.element=elemname
    #   self.atomNum=1
    #   self.atomLetter=""


    self.positionVector=[inputList.pop(1),inputList.pop(1),inputList.pop(1)]
    self.remainingNumbers=inputList
  
  def __str__(self):
    atomletter="No letter" if self.atomLetter=="" else self.atomLetter
    return f"Name: {self.element}, position : {self.positionVector}, elemNum: {self.atomNum}, atom letter: {atomletter}\n"
  


