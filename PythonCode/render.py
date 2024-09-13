import matplotlib.pyplot as plt
import xlsxwriter as xl

class ExportUnit:
  def __init__(self,file,singleProperty,atoms1,doubleProperty1,atoms2,doubleProperty2):
    self.file=file
    self.singleProperty=singleProperty
    self.atoms1=""
    for j in atoms1:
      self.atoms1+=f" {j} --"
    self.atoms2=""
    for j in atoms2:
      self.atoms2+=f" {j} --"
    self.doubleProperty1=doubleProperty1
    self.doubleProperty2 = doubleProperty2
  
  @staticmethod
  def ExportSingeProp_DoubleProp(dataArray, singlePropName,doublePropName, title):
    dataArray.sort(key=lambda x: x.singleProperty) # Sort dataArray by angle
    exportGrid=xl.Workbook("exportfullData.xlsx")
    exportSheet = exportGrid.add_worksheet(title)  # Define exportSheet variable

    row=1
    exportSheet.write(0,0,"File")
    exportSheet.write(0,1,singlePropName)
    exportSheet.write(0,2,f"{doublePropName}1")
    exportSheet.write(0,3,f"{doublePropName}2")

    for data in dataArray:
      exportSheet.write(row,0,data.file)
      exportSheet.write(row,1,data.singleProperty)
      exportSheet.write(row,2,data.atoms1)
      exportSheet.write(row+1,2,data.doubleProperty1)
      exportSheet.write(row,3,data.atoms2)
      exportSheet.write(row+1,3,data.doubleProperty2)
      row+=3
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

  def barGraph(self, xData,yData, title, xLabel, yLabel):
    plt.bar(xData, yData)
    plt.title(title)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.show()
  
  def barGraphFrequencies(self,xData,yData, totalCount,title,xLabel, yLabel):
    fig,ax1=plt.subplots()
    ax1.bar(xData,yData,color='blue',alpha=0.6)
    ax1.set_xlabel(xLabel)
    ax1.set_ylabel(yLabel,color='blue')
    ax2=ax1.twinx()
    ax2.set_ylabel('Percentage',color='red')
    percentages=[round(100*yElem/totalCount,3) for yElem in yData]
    ax2.plot(xData,percentages,'r--',marker='o')
    ax2.set_ylim(0,100)
    plt.title(title)
    plt.show()

  def barGraphHistogramPercentage(self,yData, totalCount,title,xLabel, yLabel):
    fig,ax1=plt.subplots()
    counts,bins,patches=ax1.hist(yData,bins=5,edgecolor="black",color='blue',alpha=0.6)
    ax1.set_xlabel(xLabel)
    ax1.set_ylabel(yLabel, color='blue')
    ax1.set_title(title)
    ax2 = ax1.twinx()
    ax2.set_ylabel('Relative Percentage', color='red')
    percentages=[round(100*count/totalCount,3) for count in counts]

    ax2.plot(bins[:-1], percentages, 'r--', marker='o')
    ax2.set_ylim(0,100)
    plt.title(title)
    plt.show()
