import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform

def uk(theta,phi): # unit vector in propagation direction
    THETA,PHI = np.meshgrid(theta,phi)
    ux = np.sin(np.ravel(THETA))*np.cos(np.ravel(PHI))
    uy = np.sin(np.ravel(THETA))*np.sin(np.ravel(PHI))
    uz = np.cos(np.ravel(THETA))
    return np.array([ux,uy,uz])

def dr(uk,epsilon): # solving anisotropic dispersion relation a_4 * n^4 + a_2 * n^2 + a_0 = 0
    a4 = uk[0]**2*epsilon[0] + uk[1]**2*epsilon[1] + uk[2]**2*epsilon[2]
    a2 = - (1-uk[0]**2)*epsilon[1]*epsilon[2] - (1-uk[1]**2)*epsilon[0]*epsilon[2] - (1-uk[2]**2)*epsilon[0]*epsilon[1]
    a0 = epsilon[0]*epsilon[1]*epsilon[2]
    p = a2/a4 
    q = a0/a4
    na = np.sqrt(- p/2 - np.real(np.sqrt(p**2/4-q+0j)))
    nb = np.sqrt(- p/2 + np.real(np.sqrt(p**2/4-q+0j)))
    return na, nb

def ns(theta,phi,epsilon): # produces coordinate vectors to plot normal surfaces 
    na,nb = dr(uk(theta,phi),epsilon)
    return na*uk(theta,phi), nb*uk(theta,phi)

def D(theta,phi,epsilon): # unit vectors of D field (to be finished)
    na,nb = dr(uk(theta,phi),epsilon)
    Dxa = epsilon[0]*uk(theta,phi)[0]*na/(epsilon[0]-na**2)
    Dya = epsilon[1]*uk(theta,phi)[1]*na/(epsilon[1]-na**2)
    Dza = epsilon[2]*uk(theta,phi)[2]*na/(epsilon[2]-na**2) 
    Dxb = epsilon[0]*uk(theta,phi)[0]*nb/(epsilon[0]-nb**2)
    Dyb = epsilon[1]*uk(theta,phi)[1]*nb/(epsilon[1]-nb**2)
    Dzb = epsilon[2]*uk(theta,phi)[2]*nb/(epsilon[2]-nb**2)
    # 
    # Dxa = 0
    # Dya = 0
    # Dza = 0
    # Dxb = 0
    # Dyb = 0
    # Dzb = 0
    # if na**2 == epsilon[0]:
    #     Dxa = 1
    # elif na**2 == epsilon[1]:
    #     Dya=1
    # elif na**2 == epsilon[2]:
    #     Dza=1
    # else:
    #     Dxa = epsilon[0]*uk(theta,phi)[0]*na/(epsilon[0]-na**2)
    #     Dya = epsilon[1]*uk(theta,phi)[1]*na/(epsilon[1]-na**2)
    #     Dza = epsilon[2]*uk(theta,phi)[2]*na/(epsilon[2]-na**2)        
    # if nb**2 == epsilon[0]:
    #     Dxb = 1
    # elif nb**2 == epsilon[1]:
    #     Dyb=1
    # elif nb**2 == epsilon[2]:
    #     Dzb=1
    # else:
    #     Dxb = epsilon[0]*uk(theta,phi)[0]*nb/(epsilon[0]-nb**2)
    #     Dyb = epsilon[1]*uk(theta,phi)[1]*nb/(epsilon[1]-nb**2)
    #     Dzb = epsilon[2]*uk(theta,phi)[2]*nb/(epsilon[2]-nb**2)
    return np.array([Dxa,Dya,Dza])/np.sqrt(Dxa**2+Dya**2+Dza**2),np.array([Dxb,Dyb,Dzb])/np.sqrt(Dxb**2+Dyb**2+Dzb**2)

class Arrow3D(FancyArrowPatch):

    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._xyz = (x, y, z)
        self._dxdydz = (dx, dy, dz)

    def draw(self, renderer):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)
        
    def do_3d_projection(self, renderer=None):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))

        return np.min(zs) 

def _arrow3D(ax, x, y, z, dx, dy, dz, *args, **kwargs):
    '''Add an 3d arrow to an `Axes3D` instance.'''

    arrow = Arrow3D(x, y, z, dx, dy, dz, *args, **kwargs)
    ax.add_artist(arrow)
    
def plot_ns(ax,theta0,phi0,epsilon):
    setattr(Axes3D, 'arrow3D', _arrow3D)
    
    if theta0 > np.pi/2:
        theta = np.linspace(np.pi/2, np.pi, 100)
    else:
        theta = np.linspace(0, np.pi/2, 100)
    phi = np.linspace((phi0//(np.pi/2))*np.pi/2, (phi0//(np.pi/2)+1)*np.pi/2, 100)

    nsxza,nsxzb = ns(theta,(np.sign(np.cos(phi[len(phi)//2]))-1)*np.pi/2,epsilon)
    nsyza,nsyzb = ns(theta,np.sign(np.sin(phi[len(phi)//2]))*np.pi/2,epsilon)
    nsxya,nsxyb = ns(np.pi/2,phi,epsilon)

    ax.plot(nsyza[0], nsyza[1], nsyza[2], label=r'$n_a$', color='b')
    ax.plot(nsyzb[0], nsyzb[1], nsyzb[2], label=r'$n_b$', color='r')
    ax.plot(nsxza[0], nsxza[1], nsxza[2], color='b')
    ax.plot(nsxzb[0], nsxzb[1], nsxzb[2], color='r')
    ax.plot(nsxya[0], nsxya[1], nsxya[2], color='b')
    ax.plot(nsxyb[0], nsxyb[1], nsxyb[2], color='r')
    ax.legend()

    ax.set_xlabel(r'$k_1c/\omega$')
    ax.set_ylabel(r'$k_2c/\omega$')
    ax.set_zlabel(r'$k_3c/\omega$')
    
    nsa,nsb = ns(theta0,phi0,epsilon)
    ax.arrow3D(0,0,0,nsa[0][0],nsa[1][0],nsa[2][0], mutation_scale=10, arrowstyle="-|>", color = 'b', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0)
    ax.arrow3D(0,0,0,nsb[0][0],nsb[1][0],nsb[2][0], mutation_scale=10, arrowstyle="-|>", color = 'r', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0)

    if np.abs(np.sin(phi0))>1.e-3 and np.abs(np.cos(phi0))>1.e-3 and np.abs(np.sin(theta0))>1.e-3:
        nsthetaa,nsthetab = ns(theta,phi0,epsilon)
        ax.plot(nsthetaa[0], nsthetaa[1], nsthetaa[2], color='b', linestyle = ':')
        ax.plot(nsthetab[0], nsthetab[1], nsthetab[2], color='r', linestyle = ':')

        nsphia,nsphib = ns(theta0,phi,epsilon)
        ax.plot(nsphia[0], nsphia[1], nsphia[2], color='b', linestyle = ':')
        ax.plot(nsphib[0], nsphib[1], nsphib[2], color='r', linestyle = ':')

    ax.arrow3D(0,0,0,1.1*ax.get_xlim()[np.argmax(np.abs(ax.get_xlim()))],0,0, mutation_scale=10, arrowstyle="-|>", color = 'k',  shrinkA=0,  shrinkB=0)
    ax.arrow3D(0,0,0,0,1.1*ax.get_ylim()[np.argmax(np.abs(ax.get_ylim()))],0, mutation_scale=10, arrowstyle="-|>", color = 'k',  shrinkA=0,  shrinkB=0)
    ax.arrow3D(0,0,0,0,0,1.1*ax.get_zlim()[np.argmax(np.abs(ax.get_zlim()))], mutation_scale=10, arrowstyle="-|>", color = 'k',  shrinkA=0,  shrinkB=0)
    
#    Da,Db = D(theta0,phi0,epsilon)
    
#    print(Da,Db)
#    print(np.sum(uk(theta0,phi0)*Da),np.sum(uk(theta0,phi0)*Db))
        
#    ax.arrow3D(0,0,0,Da[0][0],Da[1][0],Da[2][0], mutation_scale=10, arrowstyle="-|>", color = 'g', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0)
#    ax.arrow3D(0,0,0,Db[0][0],Db[1][0],Db[2][0], mutation_scale=10, arrowstyle="-|>", color = 'g', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0)