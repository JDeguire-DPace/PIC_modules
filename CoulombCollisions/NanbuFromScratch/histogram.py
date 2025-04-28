import numpy as np
import matplotlib.pyplot as plt

# Load the samples from the file
samples = np.loadtxt("samples.txt")

# Plot the histogram
plt.figure(figsize=(8, 5))
plt.hist(samples, bins=500, alpha=0.7, color='red')

# Labels and title
plt.xlabel("x")
plt.ylabel("Probability Density")
plt.title("Histogram of Samples from f(x)")

plt.grid(True)
plt.tight_layout()
plt.show()