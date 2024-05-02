class Graph:
  def __init__(self,atomList):
    self.structure={}
    for atom in atomList:
      self.structure[atom]=[]
      for j in atomList:
        if(atom!=j):
          self.structure[atom].append((j,5))
    
  def __str__(self):
    output=""
    for key in self.structure.keys():
      output+=f"{key} : {self.structure[key]}\n"
    return output