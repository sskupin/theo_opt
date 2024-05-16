import numpy as np
from numpy import linalg as LA
import scipy.optimize as spo
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

def ns(uk,epsilon): # produces coordinate vectors to plot normal surfaces (vectorial input possible)
    na,nb = dr(uk,epsilon)
    return na*uk, nb*uk

def Rinv(theta,phi,epsilon): # computes the change-of-basis matrix from lab to crystal frame, assuming x,y directions given by Da,Db
    Da, Db = D(theta,phi,epsilon)
    return np.transpose(np.array([Da,Db,uk(theta,phi)[:,0]]))    

def kz(kx,ky,epsilon,theta,phi): # computes kz for given kx, ky, and mean propagation direction given by theta, phi (scalar input only)
    M = Rinv(theta,phi,epsilon)
    def funca(x):
        k = np.sqrt(kx**2+ky**2+x**2)
        na, nb = dr(np.matmul(M,np.array([kx/k,ky/k,x/k])),epsilon)
        return na**2 - k**2
    def funcb(x):
        k = np.sqrt(kx**2+ky**2+x**2)
        na, nb = dr(np.matmul(M,np.array([kx/k,ky/k,x/k])),epsilon)
        return nb**2 - k**2    
    na, nb = dr(uk(theta,phi),epsilon)
    if (epsilon[0] == epsilon[1] and epsilon[0] == epsilon[2]):
        kza = np.sqrt(na**2-kx**2-ky**2)
        kzb = np.sqrt(nb**2-kx**2-ky**2)
    else:
        kza = spo.fsolve(funca, na)
        kzb = spo.fsolve(funcb, nb)
    return kza, kzb

def kz_taylor(epsilon,theta,phi): # computes Taylor coefficients up to second order for given kx, ky, and mean propagation direction given by theta, phi (scalar input only)
    kz0a, kz0b = kz(0,0,epsilon,theta,phi)
    na, nb = dr(uk(theta,phi),epsilon)
    deltakperp = 0.01*np.minimum(na,nb)
    kzpxa, kzpxb = kz(deltakperp,0,epsilon,theta,phi)
    kzmxa, kzmxb = kz(-deltakperp,0,epsilon,theta,phi)
    kzpya, kzpyb = kz(0,deltakperp,epsilon,theta,phi)
    kzmya, kzmyb = kz(0,-deltakperp,epsilon,theta,phi)
    deltaxa = (kzpxa - kzmxa)/(2*deltakperp)
    deltaxb = (kzpxb - kzmxb)/(2*deltakperp)
    deltaya = (kzpya - kzmya)/(2*deltakperp)
    deltayb = (kzpyb - kzmyb)/(2*deltakperp)   
    Dxa = (kzpxa + kzmxa - 2*kz0a)/(2*deltakperp**2)
    Dxb = (kzpxb + kzmxb - 2*kz0b)/(2*deltakperp**2)
    Dya = (kzpya + kzmya - 2*kz0a)/(2*deltakperp**2)
    Dyb = (kzpyb + kzmyb - 2*kz0b)/(2*deltakperp**2)
    kzpxpya, kzpxpyb = kz(deltakperp,deltakperp,epsilon,theta,phi)
    kzmxpya, kzmxpyb = kz(-deltakperp,deltakperp,epsilon,theta,phi)
    kzpxmya, kzpxmyb = kz(deltakperp,-deltakperp,epsilon,theta,phi)
    kzmxmya, kzmxmyb = kz(-deltakperp,-deltakperp,epsilon,theta,phi)
    Dxya = (kzpxpya - kzmxpya - kzpxmya + kzmxmya)/(4*deltakperp**2)
    Dxyb = (kzpxpyb - kzmxpyb - kzpxmyb + kzmxmyb)/(4*deltakperp**2)
    return np.array([kz0a,deltaxa,deltaya,Dxa,Dya,Dxya]) , np.array([kz0b,deltaxb,deltayb,Dxb,Dyb,Dxyb])

def n(uk,epsilon): # computes surface distance to origin of index ellipsoide (vectorial input possible)
    return np.sqrt(1/(uk[0]**2/epsilon[0]+uk[1]**2/epsilon[1]+uk[2]**2/epsilon[2]))

def ie(uk,epsilon): # produces coordinate vector to plot index ellipsoid (vectorial input possible)
    return n(uk,epsilon)*uk

def D(theta,phi,epsilon): # unit vectors of D field (scalar input only)
    # compute matrix M with MD=1/n^2D
    M11 = (uk(theta,phi)[1]**2+uk(theta,phi)[2]**2)/epsilon[0]
    M21 = - uk(theta,phi)[1]*uk(theta,phi)[0]/epsilon[0]
    M31 = - uk(theta,phi)[2]*uk(theta,phi)[0]/epsilon[0]
    M12 = - uk(theta,phi)[0]*uk(theta,phi)[1]/epsilon[1]
    M22 = (uk(theta,phi)[0]**2+uk(theta,phi)[2]**2)/epsilon[1]
    M32 = - uk(theta,phi)[2]*uk(theta,phi)[1]/epsilon[1]
    M13 = - uk(theta,phi)[0]*uk(theta,phi)[2]/epsilon[2]
    M23 = - uk(theta,phi)[1]*uk(theta,phi)[2]/epsilon[2]
    M33 = (uk(theta,phi)[0]**2+uk(theta,phi)[1]**2)/epsilon[2]
    M = np.array([[M11[0], M12[0], M13[0]], [M21[0], M22[0], M23[0]], [M31[0], M32[0], M33[0]]])
    #
    eigenvalues, eigenvectors = LA.eig(M)
    for index in range(3):
        if index == np.argmax(eigenvalues):
            Da = eigenvectors[:,index]
        elif index == np.argmin(eigenvalues):
            pass
        else:
            Db = eigenvectors[:,index]
    if LA.det(np.array([Da,Db,uk(theta,phi)[:,0]])) < 0:
        Da = -Da
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
    
def plot_ns(ax,theta0,phi0,epsilon,show_E,show_S,small):
    setattr(Axes3D, 'arrow3D', _arrow3D)
    
    if theta0 > np.pi/2:
        theta = np.linspace(np.pi/2, np.pi, 100)
    else:
        theta = np.linspace(0, np.pi/2, 100)
    phi = np.linspace((phi0//(np.pi/2))*np.pi/2, (phi0//(np.pi/2)+1)*np.pi/2, 100)

    nsxza,nsxzb = ns(uk(theta,(np.sign(np.cos(phi[len(phi)//2]))-1)*np.pi/2),epsilon)
    nsyza,nsyzb = ns(uk(theta,np.sign(np.sin(phi[len(phi)//2]))*np.pi/2),epsilon)
    nsxya,nsxyb = ns(uk(np.pi/2,phi),epsilon)

    ax.plot(nsyza[0], nsyza[1], nsyza[2], label=r'$n_a$', color='b')
    ax.plot(nsyzb[0], nsyzb[1], nsyzb[2], label=r'$n_b$', color='r')
    ax.plot(nsxza[0], nsxza[1], nsxza[2], color='b')
    ax.plot(nsxzb[0], nsxzb[1], nsxzb[2], color='r')
    ax.plot(nsxya[0], nsxya[1], nsxya[2], color='b')
    ax.plot(nsxyb[0], nsxyb[1], nsxyb[2], color='r')
    ax.legend()
    
    nsxyza,nsxyzb = ns(uk(theta,phi),epsilon)
    nsxyza = nsxyza.reshape((3, 100, 100))
    nsxyzb = nsxyzb.reshape((3, 100, 100))

    ax.plot_surface(nsxyza[0], nsxyza[1], nsxyza[2], color='b', alpha = .2)
    ax.plot_surface(nsxyzb[0], nsxyzb[1], nsxyzb[2], color='r', alpha = .2)
    
    if small:
        ax.get_xaxis().set_ticklabels([])
        ax.get_yaxis().set_ticklabels([])
        ax.get_zaxis().set_ticklabels([])
        ax.xaxis.labelpad=-10
        ax.yaxis.labelpad=-10
        ax.zaxis.labelpad=-10

    ax.set_xlabel(r'$k_1c/\omega$')
    ax.set_ylabel(r'$k_2c/\omega$')
    ax.set_zlabel(r'$k_3c/\omega$')
    
    ax.arrow3D(0,0,0,uk(theta0,phi0)[0][0],uk(theta0,phi0)[1][0],uk(theta0,phi0)[2][0], mutation_scale=10, arrowstyle="-|>", color = 'k', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0, clip_on = False)
    if small:
        ax.text(uk(theta0,phi0)[0][0],uk(theta0,phi0)[1][0],uk(theta0,phi0)[2][0],r'$z$',color='k', clip_on = False)
    else:
        ax.text(uk(theta0,phi0)[0][0],uk(theta0,phi0)[1][0],uk(theta0,phi0)[2][0],r'$\mathbf{u}^{\rm k}$',color='k', clip_on = False)
    ax.arrow3D(0,0,0,np.sqrt(np.amax(epsilon))*uk(theta0,phi0)[0][0],np.sqrt(np.amax(epsilon))*uk(theta0,phi0)[1][0],np.sqrt(np.amax(epsilon))*uk(theta0,phi0)[2][0], mutation_scale=10, arrowstyle="-", color = 'k', linestyle="dotted", lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0, clip_on = False)
    
    nsa,nsb = ns(uk(theta0,phi0),epsilon)
    if not small:
        ax.text(nsa[0][0],nsa[1][0],nsa[2][0],r'$n^a$',color='b', clip_on = False)
        ax.text(nsb[0][0],nsb[1][0],nsb[2][0],r'$n^b$',color='r', clip_on = False)
    ax.plot(nsa[0][0],nsa[1][0],nsa[2][0], color='b', marker = '.')
    ax.plot(nsb[0][0],nsb[1][0],nsb[2][0], color='r', marker = '.')

    if np.abs(np.sin(phi0))>1.e-3 and np.abs(np.cos(phi0))>1.e-3 and np.abs(np.sin(theta0))>1.e-3:
        nsthetaa,nsthetab = ns(uk(theta,phi0),epsilon)
        ax.plot(nsthetaa[0], nsthetaa[1], nsthetaa[2], color='b', linestyle = ':')
        ax.plot(nsthetab[0], nsthetab[1], nsthetab[2], color='r', linestyle = ':')

        nsphia,nsphib = ns(uk(theta0,phi),epsilon)
        ax.plot(nsphia[0], nsphia[1], nsphia[2], color='b', linestyle = ':')
        ax.plot(nsphib[0], nsphib[1], nsphib[2], color='r', linestyle = ':')

    ax.arrow3D(0,0,0,1.1*ax.get_xlim()[np.argmax(np.abs(ax.get_xlim()))],0,0, mutation_scale=10, arrowstyle="-|>", color = 'k',  shrinkA=0,  shrinkB=0, clip_on = False)
    ax.arrow3D(0,0,0,0,1.1*ax.get_ylim()[np.argmax(np.abs(ax.get_ylim()))],0, mutation_scale=10, arrowstyle="-|>", color = 'k',  shrinkA=0,  shrinkB=0, clip_on = False)
    ax.arrow3D(0,0,0,0,0,1.1*ax.get_zlim()[np.argmax(np.abs(ax.get_zlim()))], mutation_scale=10, arrowstyle="-|>", color = 'k',  shrinkA=0,  shrinkB=0, clip_on = False)
    
    Da,Db = D(theta0,phi0,epsilon)
    
    Da=Da*0.4*np.sqrt(np.amax(epsilon))
    Db=Db*0.4*np.sqrt(np.amax(epsilon))
    
    if small:
        ax.arrow3D(0,0,0,Da[0],Da[1],Da[2], mutation_scale=10, arrowstyle="-|>", color = 'k', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0, clip_on = False)
        ax.text(Da[0],Da[1],Da[2],r'$x$',color='k', clip_on = False)
        ax.arrow3D(0,0,0,Db[0],Db[1],Db[2], mutation_scale=10, arrowstyle="-|>", color = 'k', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0, clip_on = False)
        ax.text(Db[0],Db[1],Db[2],r'$y$',color='k', clip_on = False)        
    else:
        ax.arrow3D(0,0,0,Da[0],Da[1],Da[2], mutation_scale=10, arrowstyle="-|>", color = 'b', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0, clip_on = False)
        ax.text(Da[0],Da[1],Da[2],r'$\mathbf{D}^a$',color='b', clip_on = False)
        ax.arrow3D(0,0,0,Db[0],Db[1],Db[2], mutation_scale=10, arrowstyle="-|>", color = 'r', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0, clip_on = False)
        ax.text(Db[0],Db[1],Db[2],r'$\mathbf{D}^b$',color='r', clip_on = False)
    
    if show_E == 'show':
        Ea,Eb = E(theta0,phi0,epsilon)
    
        Ea=Ea*0.4*np.sqrt(np.amax(epsilon))
        Eb=Eb*0.4*np.sqrt(np.amax(epsilon))
        
        ax.arrow3D(nsa[0][0],nsa[1][0],nsa[2][0],Ea[0],Ea[1],Ea[2], mutation_scale=10, arrowstyle="-|>", color = 'b', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0, clip_on = False)
        ax.text(Ea[0]+nsa[0][0],Ea[1]+nsa[1][0],Ea[2]+nsa[2][0],r'$\mathbf{E}^a$',color='b', clip_on = False)
        ax.arrow3D(nsb[0][0],nsb[1][0],nsb[2][0],Eb[0],Eb[1],Eb[2], mutation_scale=10, arrowstyle="-|>", color = 'r', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0, clip_on = False)
        ax.text(Eb[0]+nsb[0][0],Eb[1]+nsb[1][0],Eb[2]+nsb[2][0],r'$\mathbf{E}^b$',color='r', clip_on = False)

    if show_S == 'show':
        Sa,Sb = S(theta0,phi0,epsilon)
    
        Sa=Sa*0.4*np.sqrt(np.amax(epsilon))
        Sb=Sb*0.4*np.sqrt(np.amax(epsilon))
    
        ax.arrow3D(nsa[0][0],nsa[1][0],nsa[2][0],Sa[0],Sa[1],Sa[2], mutation_scale=10, arrowstyle="-|>", color = 'b', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0, clip_on = False)
        ax.text(Sa[0]+nsa[0][0],Sa[1]+nsa[1][0],Sa[2]+nsa[2][0],r'$\mathbf{S}^a$',color='b', clip_on = False)
        ax.arrow3D(nsb[0][0],nsb[1][0],nsb[2][0],Sb[0],Sb[1],Sb[2], mutation_scale=10, arrowstyle="-|>", color = 'r', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0, clip_on = False)
        ax.text(Sb[0]+nsb[0][0],Sb[1]+nsb[1][0],Sb[2]+nsb[2][0],r'$\mathbf{S}^b$',color='r', clip_on = False)

def plot_kz(ax,theta0,phi0,epsilon,mode):
    setattr(Axes3D, 'arrow3D', _arrow3D)
    
    na,nb = dr(uk(theta0,phi0),epsilon)
    kperp = np.linspace(0, 0.8*np.minimum(na,nb), 10)
    phi = np.linspace(-np.pi, np.pi/2, 100)
    KPERP,PHI = np.meshgrid(kperp,phi)
    KX = KPERP*np.cos(PHI)
    KY = KPERP*np.sin(PHI)

    vkz = np.vectorize(kz, excluded=[2,3,4])
    kza,kzb = vkz(KX,KY,epsilon,theta0,phi0)
    kza = kza.reshape((100, 10))
    kzb = kzb.reshape((100, 10))

    taylora, taylorb = kz_taylor(epsilon,theta0,phi0)

    if (mode == 'a'):
        ax.plot_surface(KX, KY, kza, color='b', alpha = .2)
        ax.plot_wireframe(KX, KY, kza, color='k', linewidths = .5 ,rstride=5,cstride=1)
        ax.plot_surface(KX, KY, taylora[0] + taylora[1]*KX + taylora[2]*KY + taylora[3]*KX**2 + taylora[4]*KY**2 + taylora[5]*KX*KY, color='g', alpha = .2)
        ax.plot_wireframe(KX, KY, taylora[0] + taylora[1]*KX + taylora[2]*KY + taylora[3]*KX**2 + taylora[4]*KY**2 + taylora[5]*KX*KY, color='k', linewidths = .5 ,rstride=5,cstride=1)
        ax.set_title(r'$k^{\rm a}_zc/\omega$')
    else:
        ax.plot_surface(KX, KY, kzb, color='r', alpha = .2)
        ax.plot_wireframe(KX, KY, kzb, color='k', linewidths = .5 ,rstride=5,cstride=1)
        ax.plot_surface(KX, KY, taylorb[0] + taylorb[1]*KX + taylorb[2]*KY + taylorb[3]*KX**2 + taylorb[4]*KY**2 + taylorb[5]*KX*KY, color='g', alpha = .2)
        ax.plot_wireframe(KX, KY, taylorb[0] + taylorb[1]*KX + taylorb[2]*KY + taylorb[3]*KX**2 + taylorb[4]*KY**2 + taylorb[5]*KX*KY, color='k', linewidths = .5 ,rstride=5,cstride=1)
        ax.set_title(r'$k^{\rm b}_zc/\omega$')
    zlimits = ax.get_zlim()
    ax.arrow3D(0,0,zlimits[0],0,0,1.25*(zlimits[1]-zlimits[0]), mutation_scale=10, arrowstyle="-|>", color = 'k', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0, clip_on = False)
    ax.text(0,0,1.25*zlimits[1]-0.25*zlimits[0],r'$z$',color='k', clip_on = False)
    ax.set_xlabel(r'$k_xc/\omega$')
    ax.set_ylabel(r'$k_yc/\omega$')

def plot_ns_uniaxial(ax,theta0,epsilon,show_E,show_S):
    
    theta = np.linspace(-np.pi, np.pi, 250)

    if epsilon[2] > epsilon[0]: # check if positive uniaxial
        nsyzor,nsyze = ns(uk(theta,np.pi/2),epsilon)
        nsor,nse = ns(uk(theta0,np.pi/2),epsilon)
        Dor,De = D(theta0,np.pi/2,epsilon)
        Eor,Ee = E(theta0,np.pi/2,epsilon)
        Sor,Se = S(theta0,np.pi/2,epsilon)
    else:
        nsyze,nsyzor = ns(uk(theta,np.pi/2),epsilon)
        nse,nsor = ns(uk(theta0,np.pi/2),epsilon)
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

    ax.annotate("", xy=(0, 0), xytext=(uk(theta0,np.pi/2)[1][0],uk(theta0,np.pi/2)[2][0]), arrowprops=dict(arrowstyle="<-", color = 'k', lw = 2, alpha = 0.8, shrinkA=0, shrinkB=0))  
    ax.annotate(r'$\mathbf{u}^{\rm k}$', xy=(0, 0), xytext=(uk(theta0,np.pi/2)[1][0]+0.15*De[1],uk(theta0,np.pi/2)[2][0]+0.15*De[2]), color = 'k',horizontalalignment='center', verticalalignment='center')
    ax.annotate("", xy=(0, 0), xytext=(np.sqrt(np.amax(epsilon))*uk(theta0,np.pi/2)[1][0],np.sqrt(np.amax(epsilon))*uk(theta0,np.pi/2)[2][0]), arrowprops=dict(arrowstyle="-", color = 'k', linestyle="dotted", lw = 2, alpha = 0.8, shrinkA=0, shrinkB=0))  
    
    ax.plot(nsor[1][0],nsor[2][0], color='b', marker = '.')
    ax.plot(nse[1][0],nse[2][0], color='r', marker = '.')

    theta_circ = np.linspace(np.pi/2, theta0, 100)
    ax.plot(0.15*np.sqrt(epsilon[0])*np.sin(theta_circ), 0.15*np.sqrt(epsilon[0])*np.cos(theta_circ), color='k')
    ax.annotate(r"$\theta$", xy=(0, 0), xytext=(0.075*np.sqrt(epsilon[0])*np.sin(theta_circ[50]), 0.075*np.sqrt(epsilon[0])*np.cos(theta_circ[50])),horizontalalignment='center', verticalalignment='center')
    
    De=De*0.4*np.sqrt(np.amax(epsilon))
    ax.annotate("", xy=(0, 0), xytext=(De[1],De[2]), arrowprops=dict(arrowstyle="<-", color = 'r', lw = 2, alpha = 0.8, shrinkA=0, shrinkB=0))    
    ax.annotate(r'$\mathbf{D}^{\rm e}$', xy=(0, 0), xytext=(1.1*De[1],1.1*De[2]), color = 'r', horizontalalignment='center', verticalalignment='center')
        
    if show_E == 'show':
        Ee=Ee*0.4*np.sqrt(np.amax(epsilon))

        ax.annotate("", xy=(nse[1][0],nse[2][0]), xytext=(nse[1][0]+Ee[1],nse[2][0]+Ee[2]), arrowprops=dict(arrowstyle="<-", color = 'r', lw = 2, alpha = 0.8, shrinkA=0, shrinkB=0))    
        ax.annotate(r'$\mathbf{E}^{\rm e}$', xy=(nse[1][0],nse[2][0]), xytext=(nse[1][0]+1.1*Ee[1],nse[2][0]+1.1*Ee[2]), color = 'r', horizontalalignment='center', verticalalignment='center')

    if show_S == 'show':
        Se=Se*0.4*np.sqrt(np.amax(epsilon))

        ax.annotate("", xy=(nse[1][0],nse[2][0]), xytext=(nse[1][0]+Se[1],nse[2][0]+Se[2]), arrowprops=dict(arrowstyle="<-", color = 'r', lw = 2, alpha = 0.8, shrinkA=0, shrinkB=0))    
        ax.annotate(r'$\mathbf{S}^{\rm e}$', xy=(nse[1][0],nse[2][0]), xytext=(nse[1][0]+1.1*Se[1],nse[2][0]+1.1*Se[2]), color = 'r', horizontalalignment='center', verticalalignment='center') 

def plot_ie(ax,theta0,phi0,epsilon,show_E,show_H):
    setattr(Axes3D, 'arrow3D', _arrow3D)
    
    theta = np.linspace(0, np.pi, 200)
    phi = np.linspace(0, 2*np.pi, 200)

    iexy = ie(uk(np.pi/2,phi),epsilon)

    ax.plot(iexy[0], iexy[1], iexy[2], label=r'index ellipsoid', color='grey')
    
    iexyz = ie(uk(theta,phi),epsilon).reshape((3, 200, 200))
    
    ax.plot_surface(iexyz[0], iexyz[1], iexyz[2],alpha=0.1)

    ax.set_xlabel(r'$k_1c/\omega$')
    ax.set_ylabel(r'$k_2c/\omega$')
    ax.set_zlabel(r'$k_3c/\omega$')
    
    ax.arrow3D(0,0,0,uk(theta0,phi0)[0][0],uk(theta0,phi0)[1][0],uk(theta0,phi0)[2][0], mutation_scale=10, arrowstyle="-|>", color = 'k', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0, clip_on = False)
    ax.text(uk(theta0,phi0)[0][0],uk(theta0,phi0)[1][0],uk(theta0,phi0)[2][0],r'$\mathbf{u}^{\rm k}$',color='k', clip_on = False)

    ax.arrow3D(0,0,0,1.3*ax.get_xlim()[np.argmax(np.abs(ax.get_xlim()))],0,0, mutation_scale=10, arrowstyle="-|>", color = 'k',  shrinkA=0,  shrinkB=0, clip_on = False)
    ax.arrow3D(0,0,0,0,1.3*ax.get_ylim()[np.argmax(np.abs(ax.get_ylim()))],0, mutation_scale=10, arrowstyle="-|>", color = 'k',  shrinkA=0,  shrinkB=0, clip_on = False)
    ax.arrow3D(0,0,0,0,0,1.3*ax.get_zlim()[np.argmax(np.abs(ax.get_zlim()))], mutation_scale=10, arrowstyle="-|>", color = 'k',  shrinkA=0,  shrinkB=0, clip_on = False)
    ax.arrow3D(0,0,0,-1.3*ax.get_xlim()[np.argmax(np.abs(ax.get_xlim()))],0,0, mutation_scale=10, arrowstyle="-|>", color = 'k',  shrinkA=0,  shrinkB=0, clip_on = False)
    ax.arrow3D(0,0,0,0,-1.3*ax.get_ylim()[np.argmax(np.abs(ax.get_ylim()))],0, mutation_scale=10, arrowstyle="-|>", color = 'k',  shrinkA=0,  shrinkB=0, clip_on = False)
    ax.arrow3D(0,0,0,0,0,-1.3*ax.get_zlim()[np.argmax(np.abs(ax.get_zlim()))], mutation_scale=10, arrowstyle="-|>", color = 'k',  shrinkA=0,  shrinkB=0, clip_on = False)

    ax.plot(0,0,np.sqrt(epsilon[2]), color='k', marker = '.')
    ax.plot(0,0,-np.sqrt(epsilon[2]), color='k', marker = '.')
    ax.plot(0,np.sqrt(epsilon[1]),0, color='k', marker = '.')
    ax.plot(0,-np.sqrt(epsilon[1]),0, color='k', marker = '.')
    ax.plot(np.sqrt(epsilon[0]),0,0, color='k', marker = '.')
    ax.plot(-np.sqrt(epsilon[0]),0,0, color='k', marker = '.')
    
    Da,Db = D(theta0,phi0,epsilon)
    
    ax.arrow3D(0,0,0,Da[0],Da[1],Da[2], mutation_scale=10, arrowstyle="-|>", color = 'b', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0, clip_on = False)
    ax.text(Da[0],Da[1],Da[2],r'$\mathbf{D}^a$',color='b', clip_on = False)
    ax.arrow3D(0,0,0,Db[0],Db[1],Db[2], mutation_scale=10, arrowstyle="-|>", color = 'r', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0, clip_on = False)
    ax.text(Db[0],Db[1],Db[2],r'$\mathbf{D}^b$',color='r', clip_on = False)
    
    na=Da*n(Da,epsilon)
    nb=Db*n(Db,epsilon)
    
    iel = ie(np.outer(na,np.cos(phi)) + np.outer(nb,np.sin(phi)),epsilon)
    ax.plot(iel[0,:], iel[1,:], iel[2,:], label=r'index ellipse', color='k', linestyle = ':')
    ax.legend()
        
    ax.arrow3D(0,0,0,na[0],na[1],na[2], mutation_scale=10, arrowstyle="-", color = 'b', linestyle="dotted", lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0, clip_on = False)
    ax.text(na[0],na[1],na[2],r'$n^a$',color='b', clip_on = False)
    ax.plot(na[0],na[1],na[2], color='b', marker = '.')
    ax.arrow3D(0,0,0,nb[0],nb[1],nb[2], mutation_scale=10, arrowstyle="-", color = 'r', linestyle="dotted", lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0, clip_on = False)
    ax.text(nb[0],nb[1],nb[2],r'$n^b$',color='r', clip_on = False)
    ax.plot(nb[0],nb[1],nb[2], color='r', marker = '.')
    
    if show_E == 'show':
        Ea,Eb = E(theta0,phi0,epsilon)
        
        ax.arrow3D(na[0],na[1],na[2],Ea[0],Ea[1],Ea[2], mutation_scale=10, arrowstyle="-|>", color = 'b', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0, clip_on = False)
        ax.text(Ea[0]+na[0],Ea[1]+na[1],Ea[2]+na[2],r'$\mathbf{E}^a$',color='b', clip_on = False)
        ax.arrow3D(nb[0],nb[1],nb[2],Eb[0],Eb[1],Eb[2], mutation_scale=10, arrowstyle="-|>", color = 'r', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0, clip_on = False)
        ax.text(Eb[0]+nb[0],Eb[1]+nb[1],Eb[2]+nb[2],r'$\mathbf{E}^b$',color='r', clip_on = False)

    if show_H == 'show':
        Ha,Hb = H(theta0,phi0,epsilon)
    
        ax.arrow3D(na[0],na[1],na[2],Ha[0],Ha[1],Ha[2], mutation_scale=10, arrowstyle="-|>", color = 'b', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0, clip_on = False)
        ax.text(Ha[0]+na[0],Ha[1]+na[1],Ha[2]+na[2],r'$\mathbf{H}^a$',color='b', clip_on = False)
        ax.arrow3D(nb[0],nb[1],nb[2],Hb[0],Hb[1],Hb[2], mutation_scale=10, arrowstyle="-|>", color = 'r', lw = 2, alpha = 0.8,  shrinkA=0,  shrinkB=0, clip_on = False)
        ax.text(Hb[0]+nb[0],Hb[1]+nb[1],Hb[2]+nb[2],r'$\mathbf{H}^b$',color='r', clip_on = False)