import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform

def uk(theta,phi): # unit vector in propagation direction (vectorial input possible)
    THETA,PHI = np.meshgrid(theta,phi)
    ux = np.sin(np.ravel(THETA))*np.cos(np.ravel(PHI))
    uy = np.sin(np.ravel(THETA))*np.sin(np.ravel(PHI))
    uz = np.cos(np.ravel(THETA))
    return np.array([ux,uy,uz])

def dr(uk,epsilon): # solving anisotropic dispersion relation a_4 * n^4 + a_2 * n^2 + a_0 = 0 (vectorial input possible)
    a4 = uk[0]**2*epsilon[0] + uk[1]**2*epsilon[1] + uk[2]**2*epsilon[2]
    a2 = - (1-uk[0]**2)*epsilon[1]*epsilon[2] - (1-uk[1]**2)*epsilon[0]*epsilon[2] - (1-uk[2]**2)*epsilon[0]*epsilon[1]
    a0 = epsilon[0]*epsilon[1]*epsilon[2]
    p = a2/a4 
    q = a0/a4
    na = np.sqrt(- p/2 - np.real(np.sqrt(p**2/4-q+0j)))
    nb = np.sqrt(- p/2 + np.real(np.sqrt(p**2/4-q+0j)))
    return na, nb

def ns(theta,phi,epsilon): # produces coordinate vectors to plot normal surfaces (vectorial input possible)
    na,nb = dr(uk(theta,phi),epsilon)
    return na*uk(theta,phi), nb*uk(theta,phi)

def D(theta,phi,epsilon): # unit vectors of D field (scalar input only)
    na,nb = dr(uk(theta,phi),epsilon)
    if na == nb: # propagation along optical axis
        Db = np.ravel(uk(theta+np.pi/2,phi))
        Da = np.cross(Db,np.ravel(uk(theta,phi)))
    else: 
        Dxa = np.array([0.])
        Dya = np.array([0.])
        Dza = np.array([0.])
        Dxb = np.array([0.])
        Dyb = np.array([0.])
        Dzb = np.array([0.])
        if np.abs(na**2 - epsilon[0]) < 1.e-6*epsilon[0]:
            Dxa = np.array([1.])
        elif np.abs(na**2 - epsilon[1]) < 1.e-6*epsilon[1]:
            Dya = np.array([1.])
        elif np.abs(na**2 - epsilon[2]) < 1.e-6*epsilon[2]:
            Dza = np.array([1.])
        else:
            Dxa = epsilon[0]*uk(theta,phi)[0]*na/(epsilon[0]-na**2)
            Dya = epsilon[1]*uk(theta,phi)[1]*na/(epsilon[1]-na**2)
            Dza = epsilon[2]*uk(theta,phi)[2]*na/(epsilon[2]-na**2)   
        if np.abs(nb**2 - epsilon[0]) < 1.e-6*epsilon[0]:
            Dxb = np.array([1.])
        elif np.abs(nb**2 - epsilon[1]) < 1.e-6*epsilon[1]:
            Dyb = np.array([1.])
        elif np.abs(nb**2 - epsilon[2]) < 1.e-6*epsilon[2]:
            Dzb = np.array([1.])
        else:
            Dxb = epsilon[0]*uk(theta,phi)[0]*nb/(epsilon[0]-nb**2)
            Dyb = epsilon[1]*uk(theta,phi)[1]*nb/(epsilon[1]-nb**2)
            Dzb = epsilon[2]*uk(theta,phi)[2]*nb/(epsilon[2]-nb**2)
        Da = np.array([Dxa[0],Dya[0],Dza[0]])/np.sqrt(Dxa**2+Dya**2+Dza**2)
        Db = np.array([Dxb[0],Dyb[0],Dzb[0]])/np.sqrt(Dxb**2+Dyb**2+Dzb**2)
        if epsilon[0] == epsilon[1]: # uniaxial
            if epsilon[2] > epsilon[0]: # check if positive uniaxial
                Da = np.cross(Db,np.ravel(uk(theta,phi)))
            else:
                Db = np.cross(np.ravel(uk(theta,phi)),Da)
        if epsilon[1] == epsilon[2]: # uniaxial
            if epsilon[0] > epsilon[1]: # check if positive uniaxial
                Da = np.cross(Db,np.ravel(uk(theta,phi)))
            else:
                Db = np.cross(np.ravel(uk(theta,phi)),Da) 
        if epsilon[0] == epsilon[2]: # uniaxial
            if epsilon[1] > epsilon[0]: # check if positive uniaxial
                Da = np.cross(Db,np.ravel(uk(theta,phi)))
            else:
                Db = np.cross(np.ravel(uk(theta,phi)),Da) 
    return Da,Db

def E(theta,phi,epsilon): # unit vectors of E field (scalar input only)
    Da,Db = D(theta,phi,epsilon)
    Ea = Da/epsilon
    Eb = Db/epsilon
    return Ea/np.sqrt(Ea[0]**2+Ea[1]**2+Ea[2]**2),Eb/np.sqrt(Eb[0]**2+Eb[1]**2+Eb[2]**2)

def H(theta,phi,epsilon): # unit vectors of H field (scalar input only)
    Da,Db = D(theta,phi,epsilon)
    k = np.ravel(uk(theta,phi))
    return np.cross(k,Da),np.cross(k,Db)

def S(theta,phi,epsilon): # unit Poynting vectors (scalar input only)
    Ea,Eb = E(theta,phi,epsilon)
    Ha,Hb = H(theta,phi,epsilon)
    return np.cross(Ea,Ha),np.cross(Eb,Hb)

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
    
def plot_ns(ax,theta0,phi0,epsilon,show_E,show_S):
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
    
    Da,Db = D(theta0,phi0,epsilon)
    
    Da=Da*0.2*np.sqrt(np.amax(epsilon))
    Db=Db*0.2*np.sqrt(np.amax(epsilon))
        
    ax.arrow3D(0,0,0,Da[0],Da[1],Da[2], mutation_scale=10, arrowstyle="-|>", color = 'b', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0)
    ax.text(Da[0],Da[1],Da[2],r'$\mathbf{D}^a$',color='b')
    ax.arrow3D(0,0,0,Db[0],Db[1],Db[2], mutation_scale=10, arrowstyle="-|>", color = 'r', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0)
    ax.text(Db[0],Db[1],Db[2],r'$\mathbf{D}^b$',color='r')
    
    if show_E == 'show':
        Ea,Eb = E(theta0,phi0,epsilon)
    
        Ea=Ea*0.2*np.sqrt(np.amax(epsilon))
        Eb=Eb*0.2*np.sqrt(np.amax(epsilon))
        
        ax.arrow3D(nsa[0][0],nsa[1][0],nsa[2][0],Ea[0],Ea[1],Ea[2], mutation_scale=10, arrowstyle="-|>", color = 'b', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0)
        ax.text(Ea[0]+nsa[0][0],Ea[1]+nsa[1][0],Ea[2]+nsa[2][0],r'$\mathbf{E}^a$',color='b')
        ax.arrow3D(nsb[0][0],nsb[1][0],nsb[2][0],Eb[0],Eb[1],Eb[2], mutation_scale=10, arrowstyle="-|>", color = 'r', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0)
        ax.text(Eb[0]+nsb[0][0],Eb[1]+nsb[1][0],Eb[2]+nsb[2][0],r'$\mathbf{E}^b$',color='r')

    if show_S == 'show':
        Sa,Sb = S(theta0,phi0,epsilon)
    
        Sa=Sa*0.2*np.sqrt(np.amax(epsilon))
        Sb=Sb*0.2*np.sqrt(np.amax(epsilon))
    
        ax.arrow3D(nsa[0][0],nsa[1][0],nsa[2][0],Sa[0],Sa[1],Sa[2], mutation_scale=10, arrowstyle="-|>", color = 'b', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0)
        ax.text(Sa[0]+nsa[0][0],Sa[1]+nsa[1][0],Sa[2]+nsa[2][0],r'$\mathbf{S}^a$',color='b')
        ax.arrow3D(nsb[0][0],nsb[1][0],nsb[2][0],Sb[0],Sb[1],Sb[2], mutation_scale=10, arrowstyle="-|>", color = 'r', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0)
        ax.text(Sb[0]+nsb[0][0],Sb[1]+nsb[1][0],Sb[2]+nsb[2][0],r'$\mathbf{S}^b$',color='r')
        
def plot_ns_uniaxial(ax,theta0,epsilon,show_E,show_S):
    
    theta = np.linspace(-np.pi, np.pi, 250)

    if epsilon[2] > epsilon[0]: # check if positive uniaxial
        nsyzor,nsyze = ns(theta,np.pi/2,epsilon)
        nsor,nse = ns(theta0,np.pi/2,epsilon)
        Dor,De = D(theta0,np.pi/2,epsilon)
        Eor,Ee = E(theta0,np.pi/2,epsilon)
        Sor,Se = S(theta0,np.pi/2,epsilon)
    else:
        nsyze,nsyzor = ns(theta,np.pi/2,epsilon)
        nse,nsor = ns(theta0,np.pi/2,epsilon)
        De,Dor = D(theta0,np.pi/2,epsilon)
        Ee,Eor = E(theta0,np.pi/2,epsilon)
        Se,Sor = S(theta0,np.pi/2,epsilon)
                
    ax.plot(nsyzor[1], nsyzor[2], label=r'$n_{\rm or}$', color='b')
    ax.plot(nsyze[1], nsyze[2], label=r'$n_{\rm e}$', color='r')
    ax.set_aspect('equal')
    ax.set_xlim([-1.5*np.sqrt(np.amax(epsilon)),1.5*np.sqrt(np.amax(epsilon))])
    ax.set_ylim([-1.5*np.sqrt(epsilon[0]),1.5*np.sqrt(epsilon[0])])
    ax.legend()

    ax.set_xlabel(r'$k_2c/\omega$')
    ax.set_ylabel(r'$k_3c/\omega$')
    
    ax.annotate("", xy=(0, 0), xytext=(1.1*np.sqrt(np.amax(epsilon)), 0), arrowprops=dict(arrowstyle="<-", shrinkA=0, shrinkB=0))
    ax.annotate(r"$k_2$", xy=(0, 0), xytext=(1.125*np.sqrt(np.amax(epsilon)), 0))
    ax.annotate("", xy=(0, 0), xytext=(0, np.sqrt(epsilon[0])+0.1*np.sqrt(np.amax(epsilon))), arrowprops=dict(arrowstyle="<-", shrinkA=0, shrinkB=0))
    ax.annotate(r"$k_3$", xy=(0, 0), xytext=(0, np.sqrt(epsilon[0])+0.125*np.sqrt(np.amax(epsilon))))   

    ax.annotate("", xy=(0, 0), xytext=(nsor[1][0],nsor[2][0]), arrowprops=dict(arrowstyle="<-", color = 'b', lw = 2, alpha = 0.8, shrinkA=0, shrinkB=0))    
    ax.annotate("", xy=(0, 0), xytext=(nse[1][0],nse[2][0]), arrowprops=dict(arrowstyle="<-", color = 'r', lw = 2, alpha = 0.8, shrinkA=0, shrinkB=0))    

    theta_circ = np.linspace(np.pi/2, theta0, 100)
    ax.plot(0.15*np.sqrt(epsilon[0])*np.sin(theta_circ), 0.15*np.sqrt(epsilon[0])*np.cos(theta_circ), color='k')
    ax.annotate(r"$\theta$", xy=(0, 0), xytext=(0.075*np.sqrt(epsilon[0])*np.sin(theta_circ[50]), 0.075*np.sqrt(epsilon[0])*np.cos(theta_circ[50])),horizontalalignment='center', verticalalignment='center')
    
    De=De*0.25*np.sqrt(np.amax(epsilon))
    ax.annotate("", xy=(0, 0), xytext=(De[1],De[2]), arrowprops=dict(arrowstyle="<-", color = 'r', lw = 2, alpha = 0.8, shrinkA=0, shrinkB=0))    
    ax.annotate(r'$\mathbf{D}^{\rm e}$', xy=(0, 0), xytext=(1.1*De[1],1.1*De[2]), color = 'r', horizontalalignment='center', verticalalignment='center')
        
    if show_E == 'show':
        Ee=Ee*0.25*np.sqrt(np.amax(epsilon))

        ax.annotate("", xy=(nse[1][0],nse[2][0]), xytext=(nse[1][0]+Ee[1],nse[2][0]+Ee[2]), arrowprops=dict(arrowstyle="<-", color = 'r', lw = 2, alpha = 0.8, shrinkA=0, shrinkB=0))    
        ax.annotate(r'$\mathbf{E}^{\rm e}$', xy=(nse[1][0],nse[2][0]), xytext=(nse[1][0]+1.1*Ee[1],nse[2][0]+1.1*Ee[2]), color = 'r', horizontalalignment='center', verticalalignment='center')

    if show_S == 'show':
        Se=Se*0.25*np.sqrt(np.amax(epsilon))

        ax.annotate("", xy=(nse[1][0],nse[2][0]), xytext=(nse[1][0]+Se[1],nse[2][0]+Se[2]), arrowprops=dict(arrowstyle="<-", color = 'r', lw = 2, alpha = 0.8, shrinkA=0, shrinkB=0))    
        ax.annotate(r'$\mathbf{S}^{\rm e}$', xy=(nse[1][0],nse[2][0]), xytext=(nse[1][0]+1.1*Se[1],nse[2][0]+1.1*Se[2]), color = 'r', horizontalalignment='center', verticalalignment='center') 