import matplotlib.pyplot as plt
import numpy as np

title = r"Grafico del potenziale per y = 200"
folder = "./cond5"
file = "_final"

fig, ax = plt.subplots(figsize=(12*1.5, 5.5*1.5))
ax.set_box_aspect(1)
plt.suptitle(title)

# Potential
x, y, z = np.loadtxt("%s/potential%s.dat" % (folder, file), unpack=True)
x = np.arange(0, max(x)+1, 1.0)
y = np.arange(0, max(y)+1, 1.0)

x, y = np.meshgrid(x, y)
z = np.reshape(z, (len(x), len(y)))

pots = []
for i in range(len(x)):
    pots += [z[i, 200]]

xs = np.arange(0, len(x), 1.0)

plt.xlabel("x")
plt.ylabel("V")
plt.plot(xs, pots)

a = 1*pots[200]
plt.plot(xs, -a*np.log(abs(204-xs)))
plt.tight_layout()
plt.show()