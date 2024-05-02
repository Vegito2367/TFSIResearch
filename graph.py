class Graph:
  def __init__(self,atomList,angle):
    self.bondAngle=angle
    self.structure={}
    for atom in atomList:
      self.structure[atom]=[]
      for j in atomList:
        if(atom!=j):
          self.structure[atom].append((j,atom.getDistance(j)))
    
  def __str__(self):
    output=""
    for key in self.structure.keys():
      output+=f"{key}|{key.symbol} : {self.structure[key]}\n"
    output+=f"Angle: {self.bondAngle}"
    return output