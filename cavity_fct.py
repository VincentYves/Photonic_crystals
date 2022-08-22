import numpy as np 
import matplotlib.pyplot as plt

def pair(N):
    m = int(N/2)
    N = m * 2
    return N


def ispair(N):
    m = pair(N)
    if m-N ==0:
        p = 1
    else:
        p = 0
    return p


# =============================================================================================
#
# =============================================================================================

class Maille_tri():
    
    def __init__(self,d,Nx,Nz,Ox = 0, Oz = 0,Inv_xs = False,**kwargs):
        
        self.d = d
        self.Nx = Nx
        self.Nz = Nz
        self.Ox = Ox
        self.Oz = Oz
        
        if ispair(self.Nx) == 1:
            self.Nx1 = int(self.Nx/2)
        else:
            self.Nx1 = int(self.Nx/2)+1
        if ispair(self.Nz) == 1:
            self.Nz1 = int(self.Nz/2)
        else:
            self.Nz1 = int(self.Nz/2)+1

            # Lattice 1
        self.dz = d*np.sqrt(3)

        self.xx1 = np.arange(-(self.Nx1/2)*self.dz+self.dz/2,(1+self.Nx1/2)*self.dz-self.dz/2,self.dz)
        self.zz1 = np.arange(-(self.Nz1/2)*self.d+self.d/2,(0.1+self.Nz1/2)*self.d-self.d/2,self.d)

        # Lattice 2
        self.xx2 = self.xx1[:self.Nx-self.Nx1] + self.dz/2
        self.zz2 = self.zz1[:self.Nz-self.Nz1] + self.d/2

        if Inv_xs == True:     
            self.xx1,self.xx2 = -self.xx1,-self.xx2

        # Size of the whole PhC
        self.length = abs(min(np.concatenate((self.zz1,self.zz2)))-max(np.concatenate((self.zz1,self.zz2))))
        self.width = abs(min(np.concatenate((self.xx1,self.xx2)))-max(np.concatenate((self.xx1,self.xx2))))

        self.center_z = (max((max(self.zz1),max(self.zz2)))+min((min(self.zz1),min(self.zz2))))/2
        self.center_x = (max((max(self.xx1),max(self.xx2)))+min((min(self.xx1),min(self.xx2))))/2
        
        # Translation to the origin
        self.xx1,self.xx2,self.zz1,self.zz2 = self.xx1-self.center_x+self.Ox,self.xx2-self.center_x+self.Ox,self.zz1-self.center_z+self.Oz,self.zz2-self.center_z+self.Oz
        
        # Creating tuples containing the coordinates of the lattice
        coord = []
        for x in self.xx1:
            for z in self.zz1:
                coord.append((x,z))

        for x in self.xx2:
            for z in self.zz2:
                coord.append((x,z))
                
        #self.coord = coord
        self.coord = sorted(coord, key=lambda t: (t[0], t[1]))
        
    def cast_coordinates(self):
        # Creating tuples containing the coordinates of the lattice
        coord = []
        for x in self.xx1:
            for z in self.zz1:
                coord.append((x,z))

        for x in self.xx2:
            for z in self.zz2:
                coord.append((x,z))
                
        self.coord = sorted(coord, key=lambda t: (t[0], t[1]))
             
    def plot_lattice(self,**kwargs):
        plt.figure()
        for cc in self.coord:
            plt.plot(cc[0],cc[1],'r o')

        
    def translate(self,dx,dz):   
        # Translation to the origin
        coord_temp = []
        for cc in self.coord:
            coord_temp.append((cc[0]+dx,cc[1]+dz))
        self.coord = coord_temp
        self.Ox += dx
        self.Oz += dz
          
    def rotate(self,theta):
        coord_temp = []
        for cc in self.coord:
            xtemp = cc[0]-self.Ox
            ztemp = cc[1]-self.Oz
            x = xtemp*np.cos(theta)+ztemp*np.sin(theta)+self.Ox
            z = ztemp*np.cos(theta)-xtemp*np.sin(theta)+self.Oz
            coord_temp.append((x,z))
        self.coord = coord_temp

    def copy(self,dx = 0, dz = 0, rot = 0):
        name = Maille_tri(self.d,self.Nx,self.Nz, Ox = self.Ox+dx, Oz = self.Oz+dz)
        name.rotate(rot)
        return name

# =============================================================================================
#
# =============================================================================================

class PhC_waveguide(Maille_tri):
    def __init__(self,L,w,d,Cx = 0,Cz = 0,Nx = 9,**kwargs):
        self.L = L 
        self.w = w
        self.d = d
        self.Cx = Cx
        self.Cz = Cz
        self.Nx = Nx
        self.Nz = 2*int(self.L/(d*np.sqrt(3)))+1
        
        self.maille1 = Maille_tri(d, self.Nx, self.Nz, Ox = self.Cx, Oz = self.Cz)
        self.maille2 = self.maille1.copy(dx = 0, dz = 0,rot = 0)# np.pi)
    
        self.maille1.translate((self.maille1.width+self.w)/2,0)
        self.maille2.translate(-(self.maille2.width+self.w)/2,0)
        
        self.width = self.w + self.maille1.width+self.maille2.width
        self.coord = sorted(self.maille1.coord+self.maille2.coord, key=lambda t: (t[0], t[1]))
        
        ind_midr = int(len(self.coord)/2)
        ind_midl = ind_midr - 1
        valr = self.coord[ind_midr][0]
        vall = self.coord[ind_midl][0]
        coord_mid_L = []
        coord_mid_R = []
        ind_sides = []
        
        for aa in range(len(self.coord)):         
            ii = len(self.coord)-aa-1
            cc = self.coord[ii]
            if round(10*cc[0]) == round(10*valr): 
                coord_mid_R.append(cc)
                ind_sides.append(ii)
            if round(10*cc[0]) == round(10*vall):
                coord_mid_L.append(cc)
                ind_sides.append(ii)
        for i in ind_sides:
            self.coord.pop(i)
            
        self.coord_mid_L = coord_mid_L
        self.coord_mid_R = coord_mid_R
        self.coord_mid = []
        for ii, cc in enumerate(self.coord_mid_L):
            xshift = cc[0]+self.w/2
            self.coord_mid.append((xshift,cc[1]))
           
    def plot_lattice(self,print_box = True,**kwargs):
        plt.figure()
        for cc in self.coord:
            plt.plot(cc[0],cc[1],'r o')
        for cc in self.coord_mid_L:
            plt.plot(cc[0],cc[1],'g o')
        for cc in self.coord_mid_R:
            plt.plot(cc[0],cc[1],'g o')
        for cc in self.coord_mid:
            plt.plot(cc[0],cc[1],'b o')
              

    def translate(self,dx,dz):   
        # Translation to the origin
        coord_temp = []
        for cc in self.coord:
            coord_temp.append((cc[0]+dx,cc[1]+dz))
        self.coord = coord_temp
        
        coord_temp = []
        for cc in self.coord_mid_L:
            coord_temp.append((cc[0]+dx,cc[1]+dz))
        self.coord_mid_L = coord_temp
        coord_temp = []
        for cc in self.coord_mid_R:
            coord_temp.append((cc[0]+dx,cc[1]+dz))
        self.coord_mid_R = coord_temp   
        coord_temp = []
        for cc in self.coord_mid:
            coord_temp.append((cc[0]+dx,cc[1]+dz))
        self.coord_mid = coord_temp    
        self.Cx += dx
        self.Cz += dz
        


    def central_row(self,dp):
        N = len(self.coord_mid)
        diff = abs(self.d - dp)
        x = -np.linspace(-N/2,N/2,N)

        for ii,cc in enumerate (self.coord_mid):
            shift = cc[1]-x[ii]*diff
            self.coord_mid[ii] = (cc[0],shift)
            
    def side_row(self,dp):                
        N = len(self.coord_mid_L)
        diff = abs(self.d - dp)
        x = -np.linspace(-N/2,N/2,N)

        for ii in range(N):
            ccl = self.coord_mid_L[ii]
            ccr = self.coord_mid_R[ii]
            shift = ccl[1]-x[ii]*diff
            self.coord_mid_L[ii] = (ccl[0],shift)
            self.coord_mid_R[ii] = (ccr[0],shift)


# =============================================================================================
#
# =============================================================================================


class L34_cavity(Maille_tri):
    
    def __init__(self,d,Cx,Cz,sz,r1,r2,r3,NL = 3,NX = 15, NZ = 31):
        
        self.d = d
        self.Cx = Cx 
        self.Cz = Cz
        self.sz = sz
        self.r1 = r1
        self.r2 = r2
        self.r3 = r3
        self.NX = NX
        self.NZ = NZ
        self.NL = NL
        self.maille = Maille_tri(d, self.NX,self.NZ, Ox = self.Cx, Oz = self.Cz)
        self.N = len(self.maille.coord)
        self.R = self.r1*np.ones(self.N)
        self.ind_rm = [int(self.N/2)]
        
        self.ind_n2 = [int(self.N/2)-int(self.NZ/2)-2,
                       int(self.N/2)+int(self.NZ/2)+2,
                       int(self.N/2)-int(self.NZ/2)+1,
                       int(self.N/2)+int(self.NZ/2)-1,
                       int(self.N/2)-2*int(self.NZ/2),
                       int(self.N/2)+2*int(self.NZ/2),
                       int(self.N/2)-2*int(self.NZ/2)-2,
                       int(self.N/2)+2*int(self.NZ/2)+2,]


        for kk in range (int((self.NL - 3)/2)+1):
            
            self.ind_rm.append(int(self.N/2)+1+kk)
            self.ind_rm.append(int(self.N/2)-1-kk)

            self.ind_n2.append(int(self.N/2)-int(self.NZ/2)-4-2*kk)
            self.ind_n2.append(int(self.N/2)+int(self.NZ/2)+4+2*kk)
            self.ind_n2.append(int(self.N/2)-int(self.NZ/2)+3+2*kk)
            self.ind_n2.append(int(self.N/2)+int(self.NZ/2)-3-2*kk)

            self.ind_n2.append(int(self.N/2)-2*int(self.NZ/2)-4-2*kk)
            self.ind_n2.append(int(self.N/2)+2*int(self.NZ/2)+4+2*kk)
            self.ind_n2.append(int(self.N/2)-2*int(self.NZ/2)+2+2*kk)
            self.ind_n2.append(int(self.N/2)+2*int(self.NZ/2)-2-2*kk)

        self.ind_n3 = [int(self.N/2)+(int((self.NL - 3)/2)+1)+1,int(self.N/2)-(int((self.NL - 3)/2)+1)-1]
        
        self.R[self.ind_rm] = 0
        self.R[self.ind_n2] = r2
        self.R[self.ind_n3] = r3

        cx1 = self.maille.coord[self.ind_n3[0]][0] 
        cz1 = self.maille.coord[self.ind_n3[0]][1] + self.sz
        
        cx2 = self.maille.coord[self.ind_n3[1]][0] 
        cz2 = self.maille.coord[self.ind_n3[1]][1] - self.sz
        
        self.maille.coord[self.ind_n3[0]] = (cx1,cz1)
        self.maille.coord[self.ind_n3[1]] = (cx2,cz2)   
        
    def plot_layout(self):
        coor_array = np.array(self.maille.coord)
        plt.figure()
        plt.scatter(coor_array[:,0],coor_array[:,1],(self.R),cmap = "viridis")
        plt.show()
