import numpy as np

def mTE(kfz,z): # fundamental matrix TE polarization
    if np.abs(kfz) == 0:
        m = np.array([[1,z],[0,1]])
    else:
        m = np.array([[np.cos(kfz*z),np.sin(kfz*z)/kfz],[-np.sin(kfz*z)*kfz,np.cos(kfz*z)]])
    return m

def mTM(kfz,epsilon_f,z): # fundamental matrix TM polarization
    if np.abs(kfz) == 0:
        m = np.array([[1,z*epsilon_f],[0,1]])
    else:
        m = np.array([[np.cos(kfz*z),np.sin(kfz*z)*epsilon_f/kfz],[-np.sin(kfz*z)*kfz/epsilon_f,np.cos(kfz*z)]])
    return m  

def MP(kx,d1,epsilon_f1,d2,epsilon_f2,N): # normalized TE and TM matrices for N double layers with thicknesses d1,d2 and dielectric constants epsilon_f1,epsilon_f1
    kf1z = np.sqrt(epsilon_f1*(2*np.pi)**2+0j-kx**2)
    kf2z = np.sqrt(epsilon_f2*(2*np.pi)**2+0j-kx**2)
    m1TE = mTE(kf1z,d1)
    m2TE = mTE(kf2z,d2)
    m1TM = mTM(kf1z,epsilon_f1,d1)
    m2TM = mTM(kf2z,epsilon_f2,d2)
    MTE = np.matmul(m2TE,m1TE)
    MTM = np.matmul(m2TM,m1TM)
    MTEP = np.linalg.matrix_power(MTE,N)
    MTMP = np.linalg.matrix_power(MTM,N)
    return MTEP,MTMP

def MP_lambda(d1,epsilon_f1,d2,epsilon_f2,lambdav,N): # matrix for N double layers with thicknesses d1,d2 and dielectric constants epsilon_f1,epsilon_f1 at normal incidence
    kf1z = np.sqrt(epsilon_f1)*2*np.pi/lambdav
    kf2z = np.sqrt(epsilon_f2)*2*np.pi/lambdav
    m1 = mTE(kf1z,d1)
    m2 = mTE(kf2z,d2)
    M = np.linalg.matrix_power(np.matmul(m2,m1),N)
    return M

def MP_Bloch(kx,d1,epsilon_f1,d2,epsilon_f2,N): # same as MP above, but calculated with Bloch wave approach
    kf1z = np.sqrt(epsilon_f1*(2*np.pi)**2+0j-kx**2)
    kf2z = np.sqrt(epsilon_f2*(2*np.pi)**2+0j-kx**2)
    m1TE = mTE(kf1z,d1)
    m2TE = mTE(kf2z,d2)
    m1TM = mTM(kf1z,epsilon_f1,d1)
    m2TM = mTM(kf2z,epsilon_f2,d2)
    MTE = np.matmul(m2TE,m1TE)
    MTM = np.matmul(m2TM,m1TM)
    kzTE = np.arccos((MTE[1,1]+MTE[0,0])/2)/(d1+d2)
    kzTM = np.arccos((MTM[1,1]+MTM[0,0])/2)/(d1+d2)
    kappafTE = -1j*(np.exp(1j*kzTE*(d1+d2))-MTE[0,0])/MTE[0,1]
    if np.real(kappafTE) < 0:
        kzTE = -kzTE
        kappafTE = -1j*(np.exp(1j*kzTE*(d1+d2))-MTE[0,0])/MTE[0,1]
    kappafTM = -1j*(np.exp(1j*kzTM*(d1+d2))-MTM[0,0])/MTM[0,1]
    if np.real(kappafTM) < 0:
        kzTM = -kzTM
        kappafTM = -1j*(np.exp(1j*kzTM*(d1+d2))-MTM[0,0])/MTM[0,1]
    kappabTE = 1j*(np.exp(1j*kzTE*(d1+d2))-MTE[1,1])/MTE[0,1]
    kappabTM = 1j*(np.exp(1j*kzTM*(d1+d2))-MTM[1,1])/MTM[0,1]
    if kappabTE == kappafTE:
        MTEP = np.array([[0,1],[0,0]])
    else:
        MTEP = 1/(kappabTE-kappafTE)*np.array([[kappabTE*np.exp(1j*kzTE*N*(d1+d2))-kappafTE*np.exp(-1j*kzTE*N*(d1+d2)),-2*np.sin(kzTE*N*(d1+d2))],
                                              [-2*kappafTE*kappabTE*np.sin(kzTE*N*(d1+d2)),-kappafTE*np.exp(1j*kzTE*N*(d1+d2))+kappabTE*np.exp(-1j*kzTE*N*(d1+d2))]])
    if kappabTM == kappafTM:
        MTMP = np.array([[0,1],[0,0]])
    else:
        MTMP = 1/(kappabTM-kappafTM)*np.array([[kappabTM*np.exp(1j*kzTM*N*(d1+d2))-kappafTM*np.exp(-1j*kzTM*N*(d1+d2)),-2*np.sin(kzTM*N*(d1+d2))],
                                              [-2*kappafTM*kappabTM*np.sin(kzTM*N*(d1+d2)),-kappafTM*np.exp(1j*kzTM*N*(d1+d2))+kappabTM*np.exp(-1j*kzTM*N*(d1+d2))]])        
    return MTEP,MTMP

def KSC(epsilon_s,epsilon_c,phi): # normalized kx and kz in substrate and cladding
    kx = np.sqrt(epsilon_s)*np.sin(phi)*2*np.pi
    ksz = np.sqrt(epsilon_s*(2*np.pi)**2-kx**2)
    kcz = np.sqrt(epsilon_c*(2*np.pi)**2-kx**2)
    return kx,ksz,kcz

def KSC_lambda(epsilon_s,epsilon_c,lambdav): # kx and kz in substrate and cladding at normal incidence
    ksz = np.sqrt(epsilon_s)*2*np.pi/lambdav
    kcz = np.sqrt(epsilon_c)*2*np.pi/lambdav
    return ksz,kcz

def RTAU(ksz,kcz,epsilon_s,epsilon_c,MTE,MTM): # coefficients of reflection and transmission and transmissivity for system characterized by matrices MTE and MTM
    NTE = ksz*MTE[1,1]+kcz*MTE[0,0]+1j*MTE[1,0]-1j*ksz*kcz*MTE[0,1]
    RTE = (ksz*MTE[1,1]-kcz*MTE[0,0]-1j*MTE[1,0]-1j*ksz*kcz*MTE[0,1])/NTE
    NTM = ksz*MTM[1,1]/epsilon_s+kcz*MTM[0,0]/epsilon_c+1j*MTM[1,0]-1j*ksz*kcz*MTM[0,1]/(epsilon_s*epsilon_c)
    RTM = -(ksz*MTM[1,1]/epsilon_s-kcz*MTM[0,0]/epsilon_c-1j*MTM[1,0]-1j*ksz*kcz*MTM[0,1]/(epsilon_s*epsilon_c))/NTM # for x component fo electric field (negative of magnetic coeff.)
    TTE = 2*ksz/NTE
    tauTE = np.real(kcz)/ksz*np.abs(TTE)**2
    TTM = 2*ksz/epsilon_s/NTM
    tauTM = np.real(kcz/epsilon_c)*epsilon_s/ksz*np.abs(TTM)**2
    return RTE,RTM,TTE,TTM,tauTE,tauTM

def RTAU_lambda(ksz,kcz,epsilon_s,epsilon_c,M): # coefficients of reflection and transmission and transmissivity for system characterized by matrix M at normal incidence
    DENOM = ksz*M[1,1]+kcz*M[0,0]+1j*M[1,0]-1j*ksz*kcz*M[0,1]
    R = (ksz*M[1,1]-kcz*M[0,0]-1j*M[1,0]-1j*ksz*kcz*M[0,1])/DENOM
    T = 2*ksz/DENOM
    tau = np.real(kcz)/np.real(ksz)*np.abs(T)**2
    return R,T,tau

def ourangle(z): # angle of pi is replaced by -pi
    ourangle = np.angle(z)
    if ourangle == np.pi:
        ourangle = -np.pi
    return ourangle

def plot_curves_vs_angle(ax,phi,curves,labels,colors,phi_min, phi_max):
    if np.floor(8*phi_max/np.pi)-np.ceil(8*phi_min/np.pi) >= 1:
        for index in range(len(labels)):
            ax.plot(phi,curves[index],colors[index],label=labels[index])
        ax.set_xticks([0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2])
        ax.set_xticklabels([r'$0$', r'$\pi/8$', r'$\pi/4$', r'$3\pi/8$', r'$\pi/2$'])
        ax.set_xlabel(r'$\varphi_{\rm i}$')
        ax.set_xlim([phi_min, phi_max])
    else:
        for index in range(len(labels)):
            ax.plot(phi/np.pi,curves[index],colors[index],label=labels[index])
        ax.set_xlabel(r'$\varphi_{\rm i}/\pi$')
        ax.set_xlim([phi_min/np.pi, phi_max/np.pi])       
    ax.set_ylabel(','.join(labels))
    ax.legend()
    
def plot_amplitude(ax,M,F,G,z):
    Fplot = M[0,0]*F+M[0,1]*G
    Gplot = M[1,0]*F+M[1,1]*G
    FGabs = np.abs(Fplot)+np.abs(Gplot)
    if np.any(np.less(FGabs,1.e-6)):
        Fplot[np.where(FGabs<1.e-6)[0][0]:]=0
        Gplot[np.where(FGabs<1.e-6)[0][0]:]=0
    F=Fplot[-1]
    G=Gplot[-1]
    ax.plot(z,np.abs(Fplot),'b')
    return F,G