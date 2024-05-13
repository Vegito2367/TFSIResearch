import matplotlib.pyplot as plt
import xlsxwriter as xl

class ExportUnit:
  def __init__(self,file,angle,bond1,distance1,bond2,distance2):
    self.file=file
    self.angle=angle
    self.bond1=f"{bond1[0]} -- {bond1[1]}"
    self.bond2=f"{bond2[0]} -- {bond2[1]}"
    self.distance1=distance1
    self.distance2=distance2
  
  @staticmethod
  def Export(dataArray,maxProps,minProps):
    dataArray.sort(key=lambda x: x.angle) # Sort dataArray by angle
    exportGrid=xl.Workbook("exportfullData.xlsx")
    exportSheet = exportGrid.add_worksheet("Angle_Sorted_Sheet")  # Define exportSheet variable

    row=1
    exportSheet.write(0,0,"File")
    exportSheet.write(0,1,"Angle")
    exportSheet.write(0,2,"S-N Distance1")
    exportSheet.write(0,3,"S-N Distance2")

    for data in dataArray:
      exportSheet.write(row,0,data.file)
      exportSheet.write(row,1,data.angle)
      exportSheet.write(row,2,data.bond1)
      exportSheet.write(row+1,2,data.distance1)
      exportSheet.write(row,3,data.bond2)
      exportSheet.write(row+1,3,data.distance2)
      row+=3
    exportSheet.write(row+1,0,"Max Distance")
    exportSheet.write(row+1,1,maxProps[0])
    exportSheet.write(row+2,0,"File")
    exportSheet.write(row+2,1,maxProps[1])
    exportSheet.write(row+3,0,"Atoms")
    exportSheet.write(row+3,1,maxProps[2])
    exportSheet.write(row+4,0,"Min Distance")
    exportSheet.write(row+4,1,minProps[0])
    exportSheet.write(row+5,0,"File")
    exportSheet.write(row+5,1,minProps[1])
    exportSheet.write(row+6,0,"Atoms")
    exportSheet.write(row+6,1,minProps[2])

    exportGrid.close()  # Close exportGrid workbook

class Render:
  def __init__(self):
    pass

  def plotHistogram(self, data, title, xLabel, ylabel):
    plt.hist(data, bins=10, alpha=0.5, color='b', edgecolor='black')
    plt.title(title)
    plt.xlabel(xLabel)
    plt.ylabel(ylabel)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.show()

  def scatterPlot(self, xData, yData, title, xLabel, yLabel):
    plt.scatter(xData, yData)
    plt.title(title)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.show()