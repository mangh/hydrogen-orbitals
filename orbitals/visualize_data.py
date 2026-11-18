import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import *

# Numbers of radial, polar, and azimuthal coordinates,
# maximum probability density (max_pd) and maximum real value (max_real)
radial_sz, azimuthal_sz, polar_sz, max_pd, max_real = np.genfromtxt('grid_data.txt', unpack=True)

radial_sz = int(radial_sz)
polar_sz = int(polar_sz)
azimuthal_sz = int(azimuthal_sz)

X, Z, Psi_sq = np.genfromtxt('data.txt', unpack=True)
X_r, Z_r, Psi_r = np.genfromtxt('data_re.txt', unpack=True)
X_im, Z_im, Psi_im = np.genfromtxt('data_im.txt', unpack=True)
X_r_phi_2, Y_r_phi_2, Psi_r_phi_2 = np.genfromtxt('data_re_phi_2.txt', unpack=True)
X_r_phi_2_im, Y_r_phi_2_im, Psi_r_phi_2_im = np.genfromtxt('data_im_phi_2.txt', unpack=True)

X_2d = np.zeros((radial_sz, azimuthal_sz))
Z_2d = np.zeros((radial_sz, azimuthal_sz))
Psi_sq_2d = np.zeros((radial_sz, azimuthal_sz))

X_2d_m = np.zeros((radial_sz, azimuthal_sz))
Z_2d_m = np.zeros((radial_sz, azimuthal_sz))
Psi_sq_2d_m = np.zeros((radial_sz, azimuthal_sz))

X_r_2d = np.zeros((radial_sz, azimuthal_sz))
Z_r_2d = np.zeros((radial_sz, azimuthal_sz))
Psi_r_2d = np.zeros((radial_sz, azimuthal_sz))

X_r_2d_m = np.zeros((radial_sz, azimuthal_sz))
Z_r_2d_m = np.zeros((radial_sz, azimuthal_sz))
Psi_r_2d_m = np.zeros((radial_sz, azimuthal_sz))

X_im_2d = np.zeros((radial_sz, azimuthal_sz))
Z_im_2d = np.zeros((radial_sz, azimuthal_sz))
Psi_im_2d = np.zeros((radial_sz, azimuthal_sz))

X_im_2d_m = np.zeros((radial_sz, azimuthal_sz))
Z_im_2d_m = np.zeros((radial_sz, azimuthal_sz))
Psi_im_2d_m = np.zeros((radial_sz, azimuthal_sz))

X_r_phi_2_2d = np.zeros((radial_sz, polar_sz))
Y_r_phi_2_2d = np.zeros((radial_sz, polar_sz))
Psi_r_phi_2_2d = np.zeros((radial_sz, polar_sz))

X_r_phi_2_2d_im = np.zeros((radial_sz, polar_sz))
Y_r_phi_2_2d_im = np.zeros((radial_sz, polar_sz))
Psi_r_phi_2_2d_im = np.zeros((radial_sz, polar_sz))

num = 0
for i in range(0, radial_sz):
    for j in range(0, azimuthal_sz):
        X_2d[i][j] = X[num]
        Z_2d[i][j] = Z[num]
        Psi_sq_2d[i][j] = Psi_sq[num]

        X_r_2d[i][j] = X_r[num]
        Z_r_2d[i][j] = Z_r[num]
        Psi_r_2d[i][j] = Psi_r[num]        

        X_im_2d[i][j] = X_im[num]
        Z_im_2d[i][j] = Z_im[num]
        Psi_im_2d[i][j] = Psi_im[num] 
        
        num = num + 1 

for i in range(0, radial_sz):
    for j in range(0, azimuthal_sz):
        X_2d_m[i][j] = X[num]
        Z_2d_m[i][j] = Z[num]
        Psi_sq_2d_m[i][j] = Psi_sq[num]

        X_r_2d_m[i][j] = X_r[num]
        Z_r_2d_m[i][j] = Z_r[num]
        Psi_r_2d_m[i][j] = Psi_r[num]

        X_im_2d_m[i][j] = X_im[num]
        Z_im_2d_m[i][j] = Z_im[num]
        Psi_im_2d_m[i][j] = Psi_im[num]        
        
        num = num + 1

X_t = np.zeros((2*radial_sz, azimuthal_sz))
Z_t = np.zeros((2*radial_sz, azimuthal_sz))
Psi_p = np.zeros((2*radial_sz, azimuthal_sz))

X_r_t = np.zeros((2*radial_sz, azimuthal_sz))
Z_r_t = np.zeros((2*radial_sz, azimuthal_sz))
Psi_r_t = np.zeros((2*radial_sz, azimuthal_sz))

X_im_t = np.zeros((2*radial_sz, azimuthal_sz))
Z_im_t = np.zeros((2*radial_sz, azimuthal_sz))
Psi_im_t = np.zeros((2*radial_sz, azimuthal_sz))

for i in range(0, radial_sz):
    for j in range(0, azimuthal_sz):
        X_t[i][j] = X_2d[i][j]
        Z_t[i][j] = Z_2d[i][j]
        Psi_p[i][j] = Psi_sq_2d[i][j]

        X_r_t[i][j] = X_r_2d[i][j]
        Z_r_t[i][j] = Z_r_2d[i][j]
        Psi_r_t[i][j] = Psi_r_2d[i][j]

        X_im_t[i][j] = X_im_2d[i][j]
        Z_im_t[i][j] = Z_im_2d[i][j]
        Psi_im_t[i][j] = Psi_im_2d[i][j]

for i in range(0, radial_sz):
    for j in range(0, azimuthal_sz):
        X_t[i+radial_sz][j] = X_2d_m[i][j]
        Z_t[i+radial_sz][j] = Z_2d_m[i][j]
        Psi_p[i+radial_sz][j] = Psi_sq_2d_m[i][j]

        X_r_t[i+radial_sz][j] = X_r_2d_m[i][j]
        Z_r_t[i+radial_sz][j] = Z_r_2d_m[i][j]
        Psi_r_t[i+radial_sz][j] = Psi_r_2d_m[i][j]

        X_im_t[i+radial_sz][j] = X_im_2d_m[i][j]
        Z_im_t[i+radial_sz][j] = Z_im_2d_m[i][j]
        Psi_im_t[i+radial_sz][j] = Psi_im_2d_m[i][j]        

num = 0
for i in range(0, radial_sz):
    for j in range(0, polar_sz):
        X_r_phi_2_2d[i][j] = X_r_phi_2[num]
        Y_r_phi_2_2d[i][j] = Y_r_phi_2[num]
        Psi_r_phi_2_2d[i][j] = Psi_r_phi_2[num]

        X_r_phi_2_2d_im[i][j] = X_r_phi_2_im[num]
        Y_r_phi_2_2d_im[i][j] = Y_r_phi_2_im[num]
        Psi_r_phi_2_2d_im[i][j] = Psi_r_phi_2_im[num]
        
        num = num + 1
        
        
cont = plt.contourf(X_t, Z_t, Psi_p, linspace(0, max_pd, 300), cmap = 'RdBu')
# cont = plt.contourf(X_r_t, Z_r_t, Psi_r_t, linspace(-max_real, max_real, 300), cmap = 'RdBu')
# cont = plt.contourf(X_im_t, Z_im_t, Psi_im_t, linspace(-max_real, max_real, 300), cmap = 'RdBu')
# cont = plt.contourf(X_r_phi_2_2d, Y_r_phi_2_2d, Psi_r_phi_2_2d, linspace(-max_real, max_real, 300), cmap = 'RdBu')
# cont = plt.contourf(X_r_phi_2_2d_im, Y_r_phi_2_2d_im, Psi_r_phi_2_2d_im, linspace(-max_real, max_real, 300), cmap = 'RdBu')

plt.colorbar()
plt.xlabel("X")
plt.ylabel("Z")

plt.show()
