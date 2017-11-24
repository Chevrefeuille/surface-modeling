import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

normals = []
centers_x = []
centers_y = []
centers_z = []

with open("example.txt") as f:
    content = f.readlines()

for x in content:
    tmp = x.replace("\n", "").split()
    centers_x.append(float(tmp[0]))
    centers_y.append(float(tmp[1]))
    centers_z.append(float(tmp[2]))
    x = [float(i) for i in tmp]

#normals.append([0,0,0,1,1,1])
normals.append([1,1,1,0.1,2,0.1])


#print(normals)
soa = np.array(normals)
print(soa)
X, Y, Z, U, V, W = zip(*soa)
fig = plt.figure()
ax = Axes3D(fig)
ax.quiver(X, Y, Z, U, V, W)
ax.scatter(centers_x, centers_y, centers_z)
plt.show()
