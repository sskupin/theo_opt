import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
root = Tk.Tk()
root.title("SPP at Vacuum Drude Interface")

def n_eff(epsilon_s,epsilon_c): # compute n_eff
    return np.sqrt(epsilon_s*epsilon_c/(epsilon_s+epsilon_c))

def mode_profile(epsilon_s,epsilon_c,n_eff): # computing mode profile
    v = np.sqrt(n_eff**2-epsilon_s)*2*np.pi
    w = np.sqrt(n_eff**2-epsilon_c)*2*np.pi
    xrange = 8/np.minimum(np.real(v),np.real(w))
    x = np.linspace(-xrange, xrange, num=4001, endpoint=True) # x in units of lambda
    Fs = np.exp(v*x[(x<0)])
    Fc = np.exp(-w*x[(x>=0)])
    F = np.concatenate((Fs,Fc))
    Gxs = n_eff*Fs/epsilon_s*2*np.pi
    Gxc = n_eff*Fc/epsilon_c*2*np.pi
    Gx = np.concatenate((Gxs,Gxc))    
    Gzs = 1j*v*Fs/epsilon_s
    Gzc = -1j*w*Fc/epsilon_c
    Gz = np.concatenate((Gzs,Gzc))      
        
    return F,Gx,Gz,x

def plot_mode(ax,x,F1,F2,label1,label2,unit):
    ax.plot(x,F1,'b',label=label1)
    ax.plot(x,F2,'r',label=label2)
    ax.set_ylim(ax.get_ylim())
    ax.plot([0,0],ax.get_ylim(),'k:')
    ax.set_xlabel(r'$x/\lambda$')
    ax.set_ylabel(label1+','+label2+' '+unit)
    ax.legend()
    ax.text(0.2, 0.55, r'substrate', verticalalignment='center', horizontalalignment='center', transform=ax.transAxes)
    ax.text(0.8, 0.55, r'cladding', verticalalignment='center', horizontalalignment='center', transform=ax.transAxes)

def initialize():
    log_omegaptau_double.set(np.log(90))
    omega0_double.set(.65)
    
    calculate()

def calculate():
    omegaptau = np.exp(float(log_omegaptau_double.get()))
    omega0 = float(omega0_double.get())
        
    f.clf()
    omega = np.linspace(1/400, 2, num=4000, endpoint=True) # omega in units of omega_p
    epsilon_s = 1-omegaptau/(omegaptau*omega**2+1j*omega)
    epsilon_s0 = 1-omegaptau/(omegaptau*omega0**2+1j*omega0)
    a1 = f.add_subplot(231)
    a1.plot(omega,np.real(n_eff(epsilon_s,1)),'b',omega0,np.real(n_eff(epsilon_s0,1)),'bo')
    a1.set_xlim([0,2])
    a1.set_xlabel(r'$\omega$')
    a1.set_xticks([0,1/np.sqrt(2),1])
    a1.set_xticklabels([r'$0$',r'$\omega_{\rm p}/\sqrt{2}$', r'$\omega_{\rm p}$'])
    a1.set_ylabel(r'$k^{\prime}c/\omega$')
    a1.plot([1/np.sqrt(2),1/np.sqrt(2)],a1.set_ylim(),'r:',[1,1],a1                                                                   .set_ylim(),'r:')
    a2 = f.add_subplot(232)
    a2.semilogy(omega,np.imag(n_eff(epsilon_s,1)),'b',omega0,np.imag(n_eff(epsilon_s0,1)),'bo')
    a2.set_xlim([0,2])
    a2.set_xlabel(r'$\omega$')
    a2.set_xticks([0,1/np.sqrt(2),1])
    a2.set_xticklabels([r'$0$',r'$\omega_{\rm p}/\sqrt{2}$', r'$\omega_{\rm p}$'])
    a2.set_ylabel(r'$k^{\prime\prime}c/\omega$')
    a2.plot([1/np.sqrt(2),1/np.sqrt(2)],a2.set_ylim(),'r:',[1,1],a2.set_ylim(),'r:')
          
    F,Gx,Gz,x = mode_profile(1,epsilon_s0,n_eff(1,epsilon_s0))
    a3 = f.add_subplot(233)
    Sx = np.real(-Gz*np.conj(F))
    Sz = np.real(Gx*np.conj(F))
    Smax = np.amax(np.sqrt(Sx**2+Sz**2))
    plot_mode(a3,x,Sz/Smax,Sx/Smax,r'$S_z$',r'$S_x$',r'$[|\mathbf{S}|_\mathrm{max}\equiv1]$')
            
    a4 = f.add_subplot(234)
    plot_mode(a4,x,np.real(F),np.imag(F),r'$\Re \mathcal{H}_y$',r'$\Im \mathcal{H}_y$',r'$[\mathcal{H}_y(x=0)\equiv1]$')
    a5 = f.add_subplot(235)
    plot_mode(a5,x,np.real(Gx),np.imag(Gx),r'$\Re \mathcal{E}_x$',r'$\Im \mathcal{E}_x$',r'$[\mathcal{H}_y(x=0)\equiv1]$')           
    a6 = f.add_subplot(236)
    plot_mode(a6,x,np.real(Gz),np.imag(Gz),r'$\Re \mathcal{E}_z$',r'$\Im \mathcal{E}_z$',r'$[\mathcal{H}_y(x=0)\equiv1]$')   
    plt.tight_layout()
            
#    plt.savefig('spp_drude.pdf',bbox_inches='tight',dpi=300, transparent=True)

    canvas.draw()

f = plt.figure(1,[10,5])
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

log_omegaptau_double = Tk.DoubleVar()
omega0_double = Tk.DoubleVar()

initialize()

row = 1
row = gui.create_title(mainframe,"Absorption parameter",row)
row = gui.create_logslider_with_latex(mainframe,r"cladding $\omega_{\rm p}\tau_{\rm f} =$",log_omegaptau_double,1,1000,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_title(mainframe,"Field parameter",row)
row = gui.create_slider_with_latex(mainframe,r"frequency $\omega/\omega_{\rm p} =$",omega0_double,0.01,1.99,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)