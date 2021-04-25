import matplotlib.pyplot as plt
import numpy as np

#title = r"Conduttore a destra $Q_{TOT}=-27$, Coduttore al centro $Q_{TOT}=0$, Conduttore a sinistra $Q_{TOT}=27$"
#title = r"Conduttore a destra $Q_{TOT}=-27$, Conduttore a sinistra $Q_{TOT}=27$"
title = r"Conduttore esterno $Q_{TOT}=-27$, conduttore interno $Q_{TOT}=27$"
#title = r"Conduttore al centro $Q_{TOT}=30$, altri neutri"

folder = "./cond5"
file = "_final"

fig, ax = plt.subplots(1, 2, figsize=(12*1.5, 5.5*1.5))
ax1 = ax[0]
ax2 = ax[1]
ax1.set_box_aspect(1)
ax2.set_box_aspect(1)
plt.suptitle(title)

# Potential
x, y, z = np.loadtxt("%s/potential%s.dat" % (folder, file), unpack=True)
x = np.arange(0, max(x)+1, 1.0)
y = np.arange(0, max(y)+1, 1.0)

x, y = np.meshgrid(x, y)
z = np.reshape(z, (len(x), len(y)))

c = ax1.pcolormesh(y, x, z, cmap='Greys_r')
ax1.set_title("Potenziale")
fig.colorbar(c, ax=ax1)

# Charges
x, y, z = np.loadtxt("%s/density%s.dat" % (folder, file), unpack=True)
x = np.arange(0, max(x)+1, 1.0)
y = np.arange(0, max(y)+1, 1.0)
max_val = max(abs(z))

x, y = np.meshgrid(x, y)
z = np.reshape(z, (len(x), len(y)))

c = ax2.pcolormesh(y, x, z, cmap='RdBu_r', vmin=-max_val, vmax=max_val, shading='auto')
ax2.set_title("Distribuzione cariche")
# set the limits of the plot to the limits of the data
ax2.axis([x.min(), x.max(), y.min(), y.max()])
fig.colorbar(c, ax=ax2)
plt.tight_layout()
plt.savefig("grafici/conduttori1.png")
plt.show()