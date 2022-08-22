import meep as mp 
import numpy as np
import matplotlib.pyplot as plt
from cavity_fct import *

nwg = 3.3 #3.48
nbg = 1.0
nsu = 1.55
Waveguide_Material = mp.Medium(index = nwg)
Background_Material =  mp.Medium(index = nbg)
Subs_Material = mp.Medium(index = nsu)



h = 250  # thickness of Si slab in nm
alpha = 190   # array pitch in nm
sz = 0#0.18*alpha
r1 = 55.18181818181818#0.29*alpha
r2 = r1#0.29*alpha
r3 = r1#0.29*alpha

def L34_find_modes(nwg,nbg,h,alpha,sz,r1,r2,r3,**kwargs):
    
    sz = sz/alpha ; r1 = r1/alpha ; r2 = r2/alpha ; r3 = r3/alpha ; h = h/alpha
    lcen = 650/alpha#1.3321
    fcen = 1/(lcen)
    df = 0.3#0.2*fcen
    name = "Sim_"+"h"+str(round(100*h))+"alpha"+str(round(alpha))+"sz"+str(round(sz*1e3))+"r1"+str(round(r1*1e3))+"r2"+str(round(r2*1e3))+"r3"+str(round(r3*1e3))
    data = {"h":h*alpha,
            "alpha":np.float(alpha),
            "nwg":np.float(nwg),
            "nbg":np.float(nbg),
            "sz":np.float(sz*alpha),
            "r1":np.float(r1*alpha),
            "r2":np.float(r2*alpha),
            "r3":np.float(r3*alpha),}
    L = L34_cavity(1,0,0,sz,r1,r2,r3,NL = 3)

    cell = mp.Vector3(20,15,5)

    geometry = [mp.Block(mp.Vector3(mp.inf,mp.inf,h),
                     center = mp.Vector3(),
                     material = Waveguide_Material)]

    for ii, cc in enumerate (L.maille.coord):
        r = L.R[ii]
        x,y = cc[1],cc[0]
        geometry.append(mp.Cylinder(radius = r, height = 1.1*h, center = (x,y,0), material=Background_Material)) 

        polar = mp.Ey

        sources = [mp.Source(src=mp.GaussianSource(fcen, fwidth= df),
                             component = polar,
                             center=mp.Vector3())]

    pml_layers = [mp.PML(1.0)]

    resolution = 10

    sim = mp.Simulation(cell_size = cell,
                    boundary_layers = pml_layers,
                    geometry = geometry,
                    sources = sources,
                    resolution = resolution)

    sim.symmetries.append(mp.Mirror(mp.Y, phase=-1))
    sim.symmetries.append(mp.Mirror(mp.Z, phase= 1))

    harminv = mp.Harminv(polar, mp.Vector3(), fcen, df)
    sim.run(mp.after_sources(harminv),
            until_after_sources=400)

    WL = []
    Q_factors = []
    number = []
    for ii,mode in enumerate (harminv.modes):
        number.append(int(ii))
        WL.append(alpha/mode.freq)
        Q_factors.append(mode.Q)
        print("========================================================")
        print("freq :",mode.freq/alpha)
        print("wl (um): ",alpha/mode.freq)
        print("Q factor: ",mode.Q)
        print("========================================================")
    
    data.update({"wl_modes":WL,
                 "Q":Q_factors})
 
    return name,data,h,alpha,sz,r1,r2,r3,number,WL,Q_factors

alphas = [190,200,210]#np.linspace(195,205,20)
R = np.linspace(55,58,400)
count = 0

"""
name,data,h,alpha,sz,r1,r2,r3,number,WL,Q = L34_find_modes(nwg,nbg,h,alpha,sz,r1,r2,r3)
print(WL)
print(Q)
print("-------------------------------------------------------------")
"""

for alpha in alphas:
    for r1 in R:
        r1 = r1
        r2 = r1
        r3 = r1
        print("-----------------------------------------------------------")
        print("COUNT", count)
        print("Launching sim :")
        print("alpha = ", alpha)
        print("r1 = ", r2)
        name,data,h,alpha,sz,r1,r2,r3,number,WL,Q = L34_find_modes(nwg,nbg,h,alpha,sz,r1,r2,r3)
        np.save(name+"_SIM_"+str(count)+".npy",data,allow_pickle = True)
        count += 1
        print("data saved ")
        print("closing sim :")
        print("alpha = ", alpha)
        print("r1 = ", r1)
        print("-----------------------------------------------------------")
"""
for r1 in R:
    r1 = r1
    r2 = r1
    r3 = r1
    print("-----------------------------------------------------------")
    print("Launching sim :")
    print("alpha = ", alpha)
    print("r1 = ", r2)
    name,data,h,alpha,sz,r1,r2,r3,number,WL,Q = L34_find_modes(nwg,nbg,h,alpha,sz,r1,r2,r3)
    np.save(name+"_SIM_"+str(count)+".npy",data,allow_pickle = True)
    count += 1
    print("data saved ")
    print("closing sim :")
    print("alpha = ", alpha)
    print("r1 = ", r1)
    print("-----------------------------------------------------------")
"""
