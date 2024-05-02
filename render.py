import matplotlib.pyplot as plt

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