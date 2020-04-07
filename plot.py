import json
import numpy as np
import matplotlib.pyplot as plt

f = open('fft.json')
data = np.array(json.loads(f.read())["data"]).T

plt.plot(*data)
plt.show()