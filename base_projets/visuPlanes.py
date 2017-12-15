import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

normals = []
centers_x = []
centers_y = []
centers_z = []
dx = []
dy = []
dz = []

with open("planes.data") as f:
    content = f.readlines()

for x in content:
    tmp = x.replace("\n", "").split()
    x = [float(i) for i in tmp]

    centers_x.append(x[0])
    centers_y.append(x[1])
    centers_z.append(x[2])
    dx.append(x[3])
    dy.append(x[4])
    dz.append(x[5])

fig = plt.figure()
ax = Axes3D(fig)
ax.quiver(centers_x, centers_y, centers_z, dx, dy, dz, length=0.2, pivot="tail")
ax.scatter(centers_x, centers_y, centers_z)
plt.show()
