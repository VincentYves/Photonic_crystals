import meep as mp
from meep import mpb
import save_class as sv

import numpy as np
#from numpy import linalg as LA
import matplotlib.pyplot as plt
#from IPython.display import Video

saving = 1
fig_show = 0

substrate = 1

r = 90           # hole radius in nm
h = 300          # thickness of Si slab in nm
alpha = 360      # array pitch in nm
num_periods = 12 # number of holes on either side of the waveguide

# Slow Light modifications to the geometry:
r1 = 1.00        # radius multiplier for first row of holes
s1 = 0.00        # shift of first row of holes
x1 = 0.00        # horizontal shift of 1st row
w  = 1.00        # width of waveguide
     
def find_band(freqs):
    edge = []
    for ii in range (np.shape(freqs)[1]):
        edge.append(freqs[-1,ii])

    edge = np.array(edge)
    gaps = np.diff(edge)
    loc = np.where(gaps>5*np.mean(gaps))
    Wband = int(loc[0][0])+1
    
    return Wband
# ===================================================================================================
def loop_SLW1(r,h,alpha,num_periods,r1,s1,x1,w,nwg=3.48,nsubs=1.45,nback=1,substrate=1):

    Waveguide_Material = mp.Medium(index = nwg)
    Background_Material =  mp.Medium(index = nback)
    Subs_Material = mp.Medium(index = nsubs)

    parameters = {"substrate":substrate,
                  "radius":r,
                  "height":h,
                  "alpha":alpha,
                  "num_periods":num_periods,
                  "r1":r1,
                  "s1":s1,
                  "x1":x1,
                  "w":w,
                  "Waveguide index":nwg,
                  "Background index":nback,
                  "Subs_index":nsubs}

    name = 'W1_wg_'+str(100*np.random.randn())
    

    
    for key, value in parameters.items() :
        print(key, " as key for ", value)
        name = name+str(round(value))
    

    root = "enter root"
    
    if saving == 1: 
        save = sv.Data_set(name,root)
        save.update_data_lib(parameters,"parameters")

    # ===========================================================================================

    #Normalize everything in terms of array pitch
    r = r/alpha; h = h/alpha; r1 = r1*r
    sc_z = 4*h

    geometry = [mp.Block(size=(mp.inf,mp.inf,h),material = Waveguide_Material),
                mp.Cylinder(radius = r1, height = 1.1*h, center = (x1, s1+0.5*w*np.sqrt(3),0), material=Background_Material),
                mp.Cylinder(radius = r1, height = 1.1*h, center = (x1,-s1-0.5*w*np.sqrt(3),0), material=Background_Material)]
    if substrate == 1:
        geometry += [mp.Block(material = Subs_Material, size = mp.Vector3(mp.inf,mp.inf,0.5*(sc_z-h)),center = mp.Vector3(z = 0.25*(sc_z+h)))]

    [geometry.append(mp.Cylinder(radius = r, height = 1.1*h, center = (0.5, 0.5*(w+ii)*np.sqrt(3),0), material=Background_Material)) for ii in range(1,num_periods,2)]
    [geometry.append(mp.Cylinder(radius = r, height = 1.1*h, center = (0.5, -0.5*(w+ii)*np.sqrt(3),0), material=Background_Material)) for ii in range(1,num_periods,2)]
    [geometry.append(mp.Cylinder(radius = r, height = 1.1*h, center = (0, 0.5*(w+ii)*np.sqrt(3),0), material=Background_Material)) for ii in range(2,num_periods,2)]
    [geometry.append(mp.Cylinder(radius = r, height = 1.1*h, center = (0, -0.5*(w+ii)*np.sqrt(3),0), material=Background_Material)) for ii in range(2,num_periods,2)]

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

    epsx = np.array(eps_data.T[np.int(eps_data.T[:,1,1].size/2),:,:])
    epsy = np.array(eps_data.T[:,np.int(eps_data.T[1,:,1].size/2),:])
    epsz = np.array(eps_data.T[:,:,np.int(eps_data.T[1,1,:].size/2)])

    eps_plot = [epsx,epsy,epsz]

    plt.subplot(1,3,1)
    plt.imshow(epsx, interpolation='spline36', cmap='binary')
    plt.subplot(1,3,2)
    plt.imshow(epsy, interpolation='spline36', cmap='binary')
    plt.subplot(1,3,3)
    plt.imshow(epsz, interpolation='spline36', cmap='binary')
    if fig_show == 1:
        plt.show()

    if saving == 1: 
        plt.savefig(save.directory+"/"+save.name+"eps.png",transparent=False)

    if saving == 1: 
        save.update_data_lib(eps_plot,"eps")

    # ===========================================================================================

    k_points = [mp.Vector3(0.25),
                mp.Vector3(0.5)]
    k_points = mp.interpolate(30, k_points) #mp.interpolate(10, k_points)

    num_bands =  30

    ms = mpb.ModeSolver(num_bands=num_bands,
                        k_points=k_points,
                        geometry=geometry,
                        geometry_lattice=geometry_lattice,
                        resolution=resolution)
    ms.run_zeven(mpb.fix_efield_phase)

    freqs = ms.all_freqs

    Target_band = find_band(freqs)

    fig, ax = plt.subplots()

    for k, omega in zip(k_points, ms.all_freqs):
        ax.scatter([k.x]*omega.size, omega, color='blue')

    ax.plot([k.x for k in k_points],ms.all_freqs[:,Target_band], color='red')
    ax.set_ylim([0.10, 0.35])
    ax.set_xlim([0.25,0.5])
    if fig_show == 1:
        plt.show()

    if saving == 1: 
        plt.savefig(save.directory+"/"+save.name+"k_points.pdf",transparent=True)

    freqs = np.array(freqs)
    k_points = np.array(k_points)

    if saving == 1: 
        save.update_data_lib(k_points,"k_points")
        save.update_data_lib(freqs,"freqs")


    # ===========================================================================================

    field_profile = ms.get_efield(Target_band+1,bloch_phase=False)

    efield = np.sum(np.abs(field_profile),axis=3)
    efieldx = np.array(efield.T[np.int(efield.T[:,0,0].size/2)])
    efieldy = np.array(efield.T[:,np.int(efield.T[0,:,0].size/2),:])
    efieldz = np.array(efield.T[:,:,np.int(efield.T[0,0,:].size/2)])

    efield_plot = [efieldx, efieldy, efieldz]

    plt.subplot(1,3,1)
    plt.contour(eps_data.T[np.int(eps_data.T[:,0,0].size/2),:,:], cmap='binary')
    plt.imshow(efieldx, interpolation='spline36', cmap='Reds')
    plt.subplot(1,3,2)
    plt.contour(eps_data.T[:,np.int(eps_data.T[0,:,0].size/2),:], cmap='binary')
    plt.imshow(efieldy, interpolation='spline36', cmap='Reds')
    plt.subplot(1,3,3)
    plt.contour(eps_data.T[:,:,np.int(eps_data.T[0,0,:].size/2)], cmap='binary')
    plt.imshow(efieldz,interpolation='spline36', cmap='Reds')
    if fig_show == 1:
        plt.show()

    if saving == 1: 
        plt.savefig(save.directory+"/"+save.name+"efield.png",transparent=False)
        save.update_data_lib(efield_plot,"efield")


r = 90 
alpha = 360
s1s = np.arange(0.02,0.16,0.02)
r1s = np.arange(0.80,1.0,0.02)

for r1 in r1s:
    for s1 in s1s:
        loop_SLW1(r,h,alpha,num_periods,r1,s1,x1,w,nwg=3.48,nsubs=1.45,nback=1,substrate=1)


"""
s1s = np.arange(0.02,0.22,0.02)
r1s = np.arange(0.9,1.11,0.02)
#r1=  0.1
s1 = 0.18 #s1s[4]#0#1#2#3#4
for r1 in r1s:
    loop_SLW1(r,h,alpha,num_periods,r1,s1,x1,w,nwg=3.48,nsubs=1.45,nback=1,substrate=1)

s1s = np.arange(0.02,0.22,0.02)
r1s = np.arange(0.9,1.11,0.02)
#r1=  0.1
s1 = 0.2 #s1s[4]#0#1#2#3#4
for r1 in r1s:
    loop_SLW1(r,h,alpha,num_periods,r1,s1,x1,w,nwg=3.48,nsubs=1.45,nback=1,substrate=1)

"""

"""
a = np.arange(330,440,10)

r = 82
for alpha in a:
    loop_SLW1(r,h,alpha,num_periods,r1,s1,x1,w,nwg=3.48,nsubs=1.45,nback=1,substrate=1)



"""
