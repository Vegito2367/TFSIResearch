import matplotlib.pyplot as plt
import numpy as np

# Sample data
data = np.random.normal(0, 1, 1000)

# Create a histogram
fig, ax1 = plt.subplots()

# Plot the frequency histogram on the left Y-axis
counts, bins, patches = ax1.hist(data, bins=30, color='blue', alpha=0.6)
ax1.set_xlabel('Data')
ax1.set_ylabel('Frequency', color='blue')

# Create the second Y-axis for percentages
ax2 = ax1.twinx()
ax2.set_ylabel('Percentage', color='red')

# Convert counts to percentage
percentages = 100 * counts / counts.sum()
ax2.plot(bins[:-1], percentages, 'r--', marker='o')

# Adjust the ticks for the percentage axis
ax2.set_ylim(0, 100)  # Adjust this based on your data

plt.title('Histogram with Frequency and Percentage')
plt.show()