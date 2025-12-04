import numpy as np

def plot_spect(ax,xlim,ylabel,kx,k,V,foufac):
    ax.plot(np.fft.fftshift(kx)/k,(np.abs(np.fft.fftshift(V))*foufac)**2,'b')
    ax.set_xlim(xlim)
    ax.set_xlabel(r'$k_x/k$')
    ax.set_ylabel(ylabel)
    
def plot_conf(ax,xlim,ylabel,x,f,v):
    ax.plot(x/f,np.abs(v)**2,'b')
    ax.set_xlim(xlim)
    ax.set_ylim([-np.max((np.abs(v))**2)*0.05,np.max((np.abs(v))**2)*1.15])
    ax.set_xlabel(r'$x/f$')
    ax.set_ylabel(ylabel)
    
def lens(ax,xlim,ylabel,propagation,k,x,v,AL,f):
    if propagation == 'paraxial':
        v = np.where(np.abs(x)<AL,v,0)*np.exp(-1j*k*x**2/(2*f))
    else:
        v = np.where(np.abs(x)<AL,v,0)*np.exp(-1j*k*np.sqrt(f**2+x**2))
    ax.plot(x/f,np.abs(v)**2,'b')
    ax.set_xlim(xlim)
    ax.set_ylim([-np.max((np.abs(v))**2)*0.05,np.max((np.abs(v))**2)*1.15])
    if AL/(f)<xlim[1]:
        if AL/(f)>0.04*xlim[1]:
           ax.annotate(r'2A', xy=(0,np.max((np.abs(v))**2)*1.075),horizontalalignment='center', verticalalignment='center')
        if AL/(f)>0.1*xlim[1]:
           ax.annotate(r'', xy=(-AL/(f),np.max((np.abs(v))**2)*1.075), xytext=(-0.05*xlim[1],np.max((np.abs(v))**2)*1.075), arrowprops=dict(arrowstyle='->'))
           ax.annotate(r'', xy=(AL/(f),np.max((np.abs(v))**2)*1.075), xytext=(0.05*xlim[1],np.max((np.abs(v))**2)*1.075), arrowprops=dict(arrowstyle='->'))
    ax.plot([-AL/f,-AL/f],ax.get_ylim(),'k:',[AL/f,AL/f],ax.get_ylim(),'k:')
    ax.set_xlabel(r'$x/f$')
    ax.set_ylabel(ylabel)
    return v
    
def prop(propagation,kx,k,f,V,v):
    if propagation == 'paraxial':
        V = V*np.exp(-1j*kx**2/(2*k)*f)
    else:
        V = V*np.exp(1j*np.sqrt(k**2-kx**2+0j)*f) 
    v = np.fft.ifft(V) 
    return v,V