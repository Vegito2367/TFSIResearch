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
    'Angle': x,
    'Distance': y,
    'Compound': pointData
    })


    fig = px.scatter(df, x='Angle', y='Distance', color='Angle', hover_data=['Compound'])
    fig.update_layout(title=title,
                  xaxis_title=xlabel,
                  yaxis_title=ylabel)
    fig.show()


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
    fig.show()
