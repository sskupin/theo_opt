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