import json
import numpy as np
import matplotlib.pyplot as plt

f = open('fft.json')
data = json.loads(f.read())
f = np.array(data["data"]).T
t = np.array(data["samples"]).T

plt.subplot(1,2,1)
plt.plot(*f)
plt.subplot(1,2,2)
plt.plot(*t)
plt.show()