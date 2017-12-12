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

with open("example.txt") as f:
    content = f.readlines()

for x in content:
    tmp = x.replace("\n", "").split()
    x = [float(i) for i in tmp]
    # point  = np.array([x[0], x[1], x[2]])
    # normal = np.array([x[3], x[4], x[5]])
    # d = -point.dot(normal)
    # xx, yy = np.meshgrid(np.linspace(x[0]-0.1, x[0]+0.1, 5), np.linspace(x[1]-0.1, x[1]+0.1, 5))
    # z = (-normal[0] * xx - normal[1] * yy - d) * 1. /normal[2]
    # deltaz = max(1, np.max(z) - np.min(z))
    # xx, yy = np.meshgrid(np.linspace(x[0]-0.1/deltaz, x[0]+0.1/deltaz, 5), np.linspace(x[1]-0.1/deltaz, x[1]+0.1/deltaz, 5))
    # #for i in range()
    # plt3d.plot_surface(xx, yy, z)
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
