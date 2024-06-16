import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
import pandas as pd
class InteractivePlot:
  def __init__(self):
    pass
  @staticmethod
  def plotInteractivePlot(x,y,pointData, xlabel,ylabel,title):
    
  ########################Angle vs Average Distance
    
  #####################

    df = pd.DataFrame({
    xlabel: x,
    ylabel: y,
    "Data": pointData
    })


    fig = px.scatter(df, x=xlabel, y=ylabel, color=xlabel, hover_data=['Data'])
    fig.update_layout(title=title,
                  xaxis_title=xlabel,
                  yaxis_title=ylabel)
    fig.show()
    
    file = open(f"{title}.html","w")
    fig.write_html(f"{title}.html")
    file.close()


  def plotInteractivePlotColorArray(x,y,pointData,colorArray, xlabel,ylabel,title, makeFile):
    
  ########################Angle vs Average Distance
    
  #####################

    df = pd.DataFrame({
    xlabel: x,
    ylabel: y,
    "Data": pointData,
    "colorArray": colorArray
    })


    fig = px.scatter(df, x=xlabel, y=ylabel,color="colorArray", color_discrete_map={
      0:"red",
      1:"blue"
    }, symbol="colorArray", hover_data=['Data'])
    fig.update_layout(title=title,
                  xaxis_title=xlabel,
                  yaxis_title=ylabel)
    fig.show()
    
    if(makeFile):
      file = open(f"{title}.html","w")
      fig.write_html(f"{title}.html")
      file.close()


  def update_point(trace, points, selector,scatter):
      print("Hello")

  def InteractiveHistogram(x,pointData, xlabel,title):
    df = pd.DataFrame({
    'Angle': x,
    'Compound': pointData
    })

    fig=px.histogram(df,x='Angle',color='Angle',hover_data=['Compound'])
    #fig = px.scatter(df, x='Angle', y='Distance', color='Angle', hover_data=['Compound'])
    fig.update_layout(title=title,
                  xaxis_title=xlabel,
                  yaxis_title="Frequency")
    fig.write_html("exportfile.html")
    fig.show()
