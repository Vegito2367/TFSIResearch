# import matplotlib.pyplot as plt
# import numpy as np

# # Sample data
# data = np.random.normal(size=1000)

# # Create a histogram
# fig, ax = plt.subplots()
# ax.hist(data, bins=30, color='blue', alpha=0.7)

# # Caption text with line breaks
# caption_text = (
#     "Interactive plot showing the distribution of structures as functions of S–N–S angle"
#     "and average S–N bond length. Moving the cursor over a specific point reveals the"
#     "properties of the given structure including the values of the S–N–S angle and the"
#     "S–N bond length as well as the CSD refcode for the corresponding entry in the"
#     "Cambridge Structural Database. The clustering of structural hits for {metal}"
#     "indicates that coordination of {metal} to the N atom of TFSI results in a significant"
#     "structural distortion of the geometric properties of the TFSI core structure."
#     "This is reasonable because..."
# )

# # Add caption below the plot using fig.text()
# fig.text(
#     0.1, -0.01, caption_text,wrap=True, horizontalalignment='left'
# )

# # Adjust the layout to make space for the caption
# plt.subplots_adjust(bottom=0.3)  # Increase bottom margin

# # Display the plot
# plt.show()
arr=[45, 32, 38, 47, 51, 33, 42, 28, 35, 44, 39, 29, 30, 34, 35, 40, 30, 48, 30, 40]
arr.sort()
print(arr)
print(sum(arr))
print(len(arr))