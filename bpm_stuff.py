import numpy as np

def propagation_init(Nx,u0,epsilon,delta_x,Nz,delta_z,Nabs,D,a):
    U0 = np.fft.fft(u0)
    u = np.zeros((Nz,Nx),dtype=np.complex128)
    u[0,:] = u0
    epsilon0 = np.average(epsilon)
    prop = np.exp(-1j*D*np.fft.fftfreq(Nx,delta_x)**2*np.pi/np.sqrt(epsilon0)*delta_z - 1j*a*np.fft.fftfreq(Nx,delta_x)**3*4*np.pi**2/np.sqrt(epsilon0)*delta_z)
    absorb = np.ones(Nx)
    if Nabs > 0:
        absorb[0:10*Nabs] = 1-1/np.cosh(5.3/Nabs*np.arange(10*Nabs)) 
        absorb[Nx-10*Nabs:Nx] = 1-1/np.cosh(5.3/Nabs*np.arange(10*Nabs-1,-1,-1))
    return U0,u,epsilon0,prop,absorb    

def propagation_parax(Nx,u0,delta_x,Nz,delta_z,Nabs): # solves du/dz = i/(4pi) d^2u/dx^2
    U0,u,epsilon0,prop,absorb = propagation_init(Nx,u0,1,delta_x,Nz,delta_z,Nabs,1,0)
    for index in range(1,Nz):
        u[index,:] = np.fft.ifft(U0*prop)*absorb
        U0 = np.fft.fft(u[index,:])
    return u

def propagation_wg(Nx,u0,epsilon,delta_x,Nz,delta_z,Nabs): # all lengths normalized to lambda
    U0,u,epsilon0,prop,absorb = propagation_init(Nx,u0,epsilon,delta_x,Nz,delta_z,Nabs,1,0)
    for index in range(1,Nz):
        u[index,:] = np.fft.ifft(U0*prop)*np.exp(1j*np.pi*(epsilon[index,:]-epsilon0)/np.sqrt(epsilon0)*delta_z)*absorb
        U0 = np.fft.fft(u[index,:])
    return u

def propagation_nls(Nx,u0,delta_x,Nz,delta_z,Nabs,N,D,Gamma,a,b): # solves du/dz = i D/(4pi) d^2u/dx^2 + i Gamma/(2pi) N^2  |u|^2 u + a/(2pi) d^3u/dx^3 - b/(4pi) u 
    U0,u,epsilon0,prop,absorb = propagation_init(Nx,u0,1,delta_x,Nz,delta_z,Nabs,D,a)
    for index in range(1,Nz):
        u0 = u[index-1,:]*np.exp((1j*Gamma*N**2*(u[index-1,:].real**2+u[index-1,:].imag**2)-b/2)*delta_z/(2*np.pi))*absorb
        u[index,:] = np.fft.ifft(np.fft.fft(u0)*prop)
    return u

def propagation_discrete(N,u0,Nz,delta_z):
    U0 = np.fft.fft(u0,2*N)
    u = np.zeros((Nz,N),dtype=np.complex128)
    u[0,:] = u0
    prop = np.fft.fftshift(np.exp(1j*np.cos(np.linspace(-np.pi,np.pi,2*N,endpoint=False))*delta_z))
    for index in range(1,Nz):
        U0 = U0*prop
        u[index,:] = np.fft.ifft(U0)[0:N]
    return u

def zeropadding_2D(u,x,y,Nx,Ny): # enlarge box by symmetric zero-padding (beam stays centered)
    u_new = np.zeros((Ny,Nx),dtype=np.complex128)
    Ny_old = np.shape(u)[0]
    Nx_old = np.shape(u)[1]
    u_new[0:Ny_old,0:Nx_old] = u
    u_new = np.roll(u_new, ((Nx-Nx_old)//2, (Ny-Ny_old)//2), axis=(1, 0))
    x_new = np.linspace(x[0]*Nx/Nx_old,-x[0]*Nx/Nx_old,Nx,endpoint=False)
    y_new = np.linspace(y[0]*Ny/Ny_old,-y[0]*Ny/Ny_old,Ny,endpoint=False)
    return u_new,x_new,y_new
    
def cropping_2D(u,x,y,Nx,Ny): # reduce box by symmetric cropping (beam stays centered)
    u_new = np.zeros((Ny,Nx),dtype=np.complex128)
    Ny_old = np.shape(u)[0]
    Nx_old = np.shape(u)[1]
    u_new = u[(Ny_old-Ny)//2:Ny_old-(Ny_old-Ny)//2,(Nx_old-Nx)//2:Nx_old-(Nx_old-Nx)//2]
    x_new = x[(Nx_old-Nx)//2:Nx_old-(Nx_old-Nx)//2]
    y_new = y[(Ny_old-Ny)//2:Ny_old-(Ny_old-Ny)//2]
    return u_new,x_new,y_new

def reshaping_2D(u,x,y,Nx,Ny): # change box by symmetric zero-padding or cropping (beam stays centered)
    Ny_old = np.shape(u)[0]
    Nx_old = np.shape(u)[1]
    if Nx >= Nx_old:
        if Ny >= Ny_old:
            u_new,x_new,y_new = zeropadding_2D(u,x,y,Nx,Ny)
        else:
            u_help,x_help,y_help = cropping_2D(u,x,y,Nx_old,Ny)
            u_new,x_new,y_new = zeropadding_2D(u_help,x_help,y_help,Nx,Ny)
    else:
        if Ny <= Ny_old:
            u_new,x_new,y_new = cropping_2D(u,x,y,Nx,Ny)
        else:
            u_help,x_help,y_help = cropping_2D(u,x,y,Nx,Ny_old)
            u_new,x_new,y_new = zeropadding_2D(u_help,x_help,y_help,Nx,Ny)
    return u_new,x_new,y_new

def init_2D_grid(k,z,profile='SG',alpha=1,beta=1):
    if profile ==  'SG':
        Lx = 8*(1+alpha*z/k.real*beta)
        Lkx = 4*(1+alpha)*beta
        Ly = Lx
        Lky = Lkx
        Nx0 = 128
        Ny0 = 128
        Nx = Nx0*(np.ceil(Lx*Lkx/512)).astype(int)
        Ny = Ny0*(np.ceil(Ly*Lky/512)).astype(int)
        Nkx = Nx
        Nky = Ny
    factor_x = (np.ceil(2*Nx/Lx/Lkx)).astype(int)
    factor_y = (np.ceil(2*Ny/Ly/Lky)).astype(int)
    x, delta_x = np.linspace(-factor_x*Lx/2,factor_x*Lx/2,factor_x*Nx,endpoint=False, retstep=True)
    y, delta_y = np.linspace(-factor_y*Ly/2,factor_y*Ly/2,factor_y*Ny,endpoint=False, retstep=True)
    kx = 2*np.pi*np.fft.fftshift(np.fft.fftfreq(factor_x*Nx,delta_x))
    ky = 2*np.pi*np.fft.fftshift(np.fft.fftfreq(factor_y*Ny,delta_y))
    return x,y,Nx,Ny,kx,ky,Nkx,Nky,Nx0,Ny0
    
def init_2D_beam(x,y,profile='SG',alpha=1,beta=1):
    X,Y = np.meshgrid(x,y, indexing='xy')
    if profile ==  'SG':
        u0 = np.exp(-(X**2+(beta*Y)**2)**alpha+0j)
    return u0

def init_prop_2D(kx,ky,k,delta_z):
    KX,KY = np.meshgrid(kx,ky, indexing='xy')
    KZ = np.sqrt(k**2-KX**2-KY**2+0j)
    prop = np.exp(1j*KZ*delta_z)
    return prop

def init_prop_2D_parax(kx,ky,k,delta_z):
    KX,KY = np.meshgrid(kx,ky, indexing='xy')
    KZ = k-(KX**2+KY**2)/(2*k)
    prop = np.exp(1j*KZ*delta_z)
    return prop

def propagation_2D(u0,prop):
    U0 = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(u0)))
    U = U0*prop
    u = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(U)))
    return u,U,U0

def init_prop_2D_aniso(KZa,KZb,delta_z):
    return np.exp(1j*KZa*delta_z), np.exp(1j*KZb*delta_z)

def propagation_2D_aniso(D0_x,D0_y,Dax,Day,Daz,Dbx,Dby,Dbz,prop_a,prop_b):
    FTD0_x = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(D0_x)))
    FTD0_y = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(D0_y)))
    U0_a = (Dby*FTD0_x-Dbx*FTD0_y)/(Dax*Dby-Day*Dbx)
    U0_b = (Dax*FTD0_y-Day*FTD0_x)/(Dax*Dby-Day*Dbx)
    U_a = U0_a*prop_a
    U_b = U0_b*prop_b 
    FTDa_x = U_a*Dax
    FTDa_y = U_a*Day
    FTDa_z = U_a*Daz
    FTDb_x = U_b*Dbx
    FTDb_y = U_b*Dby
    FTDb_z = U_b*Dbz
    FTD_x = FTDa_x + FTDb_x
    FTD_y = FTDa_y + FTDb_y
    FTD_z = FTDa_z + FTDb_z
    D_x = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(FTD_x)))
    D_y = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(FTD_y))) 
    D_z = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(FTD_z))) 
    return D_x,D_y,D_z,FTD_x,FTD_y,FTD_z,FTDa_x,FTDa_y,FTDa_z,FTDb_x,FTDb_y,FTDb_z