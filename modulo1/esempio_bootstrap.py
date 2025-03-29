import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 10, 100)
y = np.sin(x)

plt.plot(x, y)
plt.xticks(np.arange(0, 10.5, 0.5))
plt.yticks(np.arange(-1, 1.1, 0.1))
plt.minorticks_on()
#plt.grid(which='both', axis='both', linestyle='--', alpha=0.5)

plt.show()

