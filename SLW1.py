import meep as mp
from meep import mpb
import save_class as sv

import numpy as np
#from numpy import linalg as LA
import matplotlib.pyplot as plt
#from IPython.display import Video

substrate = 1

r = 114          # hole radius in nm
h = 220          # thickness of Si slab in nm
alpha = 398      # array pitch in nm
num_periods = 12 # number of holes on either side of the waveguide

n = 3.48         # material index

# Slow Light modifications to the geometry:
r1 = 1.00        # radius multiplier for first row of holes
s1 = 0.00        # shift of first row of holes
x1 = 0.00        # horizontal shift of 1st row
w  = 1.00        # width of waveguide
     


#Define Geometry:
nSi = 3.48
Si = mp.Medium(index = nSi)
SiO2 =  mp.Medium(index = 1.45)
air = mp.air

Waveguide_Material = Si    # Material of the slab
Background_Material = air # Material of the holes 
Subs_Material = SiO2 # Material of the substrate

# ===================================================================================================

parameters = {"substrtate":substrate,
              "radius":r,
              "height":h,
              "alpha":alpha,
              "num_periods":num_periods,
              "r1":r1,
              "s1":s1,
              "x1":x1,
              "w":w,
              "Waveguide index":nSi,
              "Background index":1,
              "Subs_index":1.45}

name = 'W1_wg_'+str(100*np.random.randn())

for key, value in parameters.items() :
    print(key, " as key for ", value)
    name = name+str(round(value))

root = "/home/eric/Vincent/slow_waveguide"
#save = sv.Data_set(name,root)

# ===================================================================================================

#Normalize everything in terms of array pitch
r = r/alpha; h = h/alpha; r1 = r1*r
sc_z = 4*h


geometry = [mp.Block(size=(mp.inf,mp.inf,h),material = Waveguide_Material),
            mp.Cylinder(radius = r1, height = 1.1*h, center = (x1, s1+0.5*w*np.sqrt(3),0), material=Background_Material),
            mp.Cylinder(radius = r1, height = 1.1*h, center = (x1,-s1-0.5*w*np.sqrt(3),0), material=Background_Material)]
[geometry.append(mp.Cylinder(radius = r, height = 1.1*h, center = (0.5, 0.5*(w+ii)*np.sqrt(3),0), material=Background_Material)) for ii in range(1,num_periods,2)]
[geometry.append(mp.Cylinder(radius = r, height = 1.1*h, center = (0.5, -0.5*(w+ii)*np.sqrt(3),0), material=Background_Material)) for ii in range(1,num_periods,2)]
[geometry.append(mp.Cylinder(radius = r, height = 1.1*h, center = (0, 0.5*(w+ii)*np.sqrt(3),0), material=Background_Material)) for ii in range(2,num_periods,2)]
[geometry.append(mp.Cylinder(radius = r, height = 1.1*h, center = (0, -0.5*(w+ii)*np.sqrt(3),0), material=Background_Material)) for ii in range(2,num_periods,2)]
if substrate == 1:
    geometry += [mp.Block(material = Subs_Material, size = mp.Vector3(mp.inf,mp.inf,0.5*(sc_z-h)),center = mp.Vector3(z = 0.25*(sc_z+h)))]


geometry_lattice = mp.Lattice(size=mp.Vector3(1,(num_periods+w-1)*np.sqrt(3),4*h))
resolutionxy = 8
resolutionz = 30

resolution = mp.Vector3(resolutionxy,resolutionxy,resolutionz)
mesh_size = 8

#Set up and initialize simulation object:
ms = mpb.ModeSolver(geometry=geometry,
                    geometry_lattice=geometry_lattice,
                    resolution=resolution)
ms.init_params(mp.EVEN_Z,True)

#Plot index profile
md = mpb.MPBData(rectify=False, periods=3, resolution=42)
eps_data = ms.get_epsilon()
eps_plot = md.convert(eps_data)

#save.update_data_lib(eps_data,"eps_data")
 
plt.subplot(1,3,1)
plt.imshow(eps_data.T[np.int(eps_data.T[:,1,1].size/2),:,:], interpolation='spline36', cmap='binary')
plt.subplot(1,3,2)
plt.imshow(eps_data.T[:,np.int(eps_data.T[1,:,1].size/2),:], interpolation='spline36', cmap='binary')
plt.subplot(1,3,3)
plt.imshow(eps_data.T[:,:,np.int(eps_data.T[1,1,:].size/2)], interpolation='spline36', cmap='binary')
plt.show()

#plt.savefig(save.directory+"/"+save.name+"eps.png",transparent=False)

# ===================================================================================================
"""

k_points = [mp.Vector3(0.25),
            mp.Vector3(0.5)]
k_points = mp.interpolate(30, k_points)

num_bands = 30

ms = mpb.ModeSolver(num_bands=num_bands,
                    k_points=k_points,
                    geometry=geometry,
                    geometry_lattice=geometry_lattice,
                    resolution=resolution)
ms.run_zeven(mpb.fix_efield_phase)


Target_band = 26

freqs = ms.all_freqs

fig, ax = plt.subplots()

for k, omega in zip(k_points, ms.all_freqs):
    ax.scatter([k.x]*omega.size, omega, color='blue')

ax.plot([k.x for k in k_points],ms.all_freqs[:,Target_band], color='red')
ax.set_ylim([0.10, 0.35])
ax.set_xlim([0.25,0.5])
#plt.savefig(save.directory+"/"+save.name+"k_points.pdf",transparent=True)

#save.update_data_lib(k_points,"k_points")
#save.update_data_lib(freqs,"freqs")


# ================================================================================================

kx = [k.x for k in k_points]
omega = ms.all_freqs[:,Target_band]
wavelength = alpha/omega 
ng = -np.gradient(kx)/np.gradient(omega)

fig, ax = plt.subplots()
ax.plot(wavelength,ng)
ax.set_ylim([0.0, 200.0])
ax.set_xlim([1460,1510])

ax.set_xlabel('Wavelength (nm)');
ax.set_ylabel('Group Index');
plt.savefig(save.directory+"/"+save.name+"ng.pdf",transparent=True)

save.update_data_lib(wavelength,"wavelength")
save.update_data_lib(ng,"ng")

# ================================================================================================

field_profile = ms.get_efield(Target_band,bloch_phase=False)

efield = np.sum(np.abs(field_profile),axis=3)

plt.subplot(1,3,1)
plt.contour(eps_data.T[np.int(eps_data.T[:,0,0].size/2),:,:], cmap='binary')
plt.imshow(efield.T[np.int(efield.T[:,0,0].size/2)], interpolation='spline36', cmap='Reds')
plt.subplot(1,3,2)
plt.contour(eps_data.T[:,np.int(eps_data.T[0,:,0].size/2),:], cmap='binary')
plt.imshow(efield.T[:,np.int(efield.T[0,:,0].size/2),:], interpolation='spline36', cmap='Reds')
plt.subplot(1,3,3)
plt.contour(eps_data.T[:,:,np.int(eps_data.T[0,0,:].size/2)], cmap='binary')
plt.imshow(efield.T[:,:,np.int(efield.T[0,0,:].size/2)],interpolation='spline36', cmap='Reds')
#plt.savefig(save.directory+"/"+save.name+"efield.png",transparent=False)

#save.update_data_lib(efield,"efield")
"""
