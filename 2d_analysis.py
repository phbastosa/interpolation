import numpy as np
import matplotlib.pyplot as plt

original_Nx = 1001
original_Ny = 1001

rescaled_Nx = 51
rescaled_Ny = 51

original_N = original_Nx*original_Ny
rescaled_N = rescaled_Nx*rescaled_Ny

original_F = np.fromfile("outputs/original2D.bin", dtype = np.float32, count = 3*original_N).reshape([original_N, 3], order = "F")
rescaled_F = np.fromfile("outputs/rescaled2D.bin", dtype = np.float32, count = 3*rescaled_N).reshape([rescaled_N, 3], order = "F")

original_V = np.reshape(original_F[:,2], [original_Ny, original_Nx], order = "F")
rescaled_V = np.reshape(rescaled_F[:,2], [rescaled_Ny, rescaled_Nx], order = "F")

cubic_F =  np.fromfile("outputs/cubic2D.bin", dtype = np.float32, count = 3*original_N).reshape([original_N, 3], order = "F")
linear_F =  np.fromfile("outputs/linear2D.bin", dtype = np.float32, count = 3*original_N).reshape([original_N, 3], order = "F")

cubic_V = np.reshape(cubic_F[:,2], [original_Ny, original_Nx], order = "F")
linear_V = np.reshape(linear_F[:,2], [original_Ny, original_Nx], order = "F")

scale = 0.01*np.max(np.abs(original_V))

diff_cubic = np.zeros_like(original_V)
diff_linear = np.zeros_like(original_V)

diff_cubic[25:-25,25:-25] = original_V[25:-25,25:-25] - cubic_V[25:-25,25:-25]
diff_linear[25:-25,25:-25] = original_V[25:-25,25:-25] - linear_V[25:-25,25:-25]

plt.figure(figsize = (15,7))

plt.subplot(231)
plt.imshow(original_V, aspect = "auto", cmap = "jet", extent = [-10,10,10,-10])
plt.scatter(rescaled_F[:,0], rescaled_F[:,1], color = "black", s = 1)
plt.colorbar()

plt.subplot(232)
plt.imshow(linear_V, aspect = "auto", cmap = "jet", extent = [-10,10,10,-10])
plt.colorbar()

plt.subplot(233)
plt.imshow(cubic_V, aspect = "auto", cmap = "jet", extent = [-10,10,10,-10])
plt.colorbar()

plt.subplot(234)
plt.imshow(rescaled_V, aspect = "auto", cmap = "jet", extent = [-10,10,10,-10])
plt.colorbar()

plt.subplot(235)
plt.imshow(diff_linear, aspect = "auto", cmap = "jet", vmax = scale, vmin = -scale, extent = [-10,10,10,-10])
plt.colorbar()

plt.subplot(236)
plt.imshow(diff_cubic, aspect = "auto", cmap = "jet", vmax = scale, vmin = -scale, extent = [-10,10,10,-10])
plt.colorbar()

plt.tight_layout()
plt.show()
