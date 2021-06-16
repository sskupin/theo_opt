import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.optimize as spo
import tkinter as Tk
import gui_stuff as gui

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rc('text', usetex=True)
mpl.rc('text.latex', preamble=r'\usepackage{cmbright}')
mpl.rcParams.update({'font.size': 10})

root = Tk.Tk()
root.title("TE Lossy Mode")

def mTE(kfx,x):
    if kfx == 0.:
        mTE01 = x
    else:
        mTE01 = np.sin(kfx*x)/kfx
    return np.array([[np.cos(kfx*x),mTE01],[-np.sin(kfx*x)*kfx,np.cos(kfx*x)]])

def RTE(epsilon_s,d1,epsilon_f1,epsilon_c,n_eff): 
    ksx = 1j*np.sqrt(n_eff**2-epsilon_s+0j)*2*np.pi
    kf1x = np.sqrt(epsilon_f1-n_eff**2+0j)*2*np.pi
    kcx = 1j*np.sqrt(n_eff**2-epsilon_c+0j)*2*np.pi
    MTE = mTE(kf1x,d1)
    
    return (ksx*MTE[1,1]-kcx*MTE[0,0]-1j*MTE[1,0]-1j*ksx*kcx*MTE[0,1])/(ksx*MTE[1,1]+kcx*MTE[0,0]+1j*MTE[1,0]-1j*ksx*kcx*MTE[0,1])

def mode_profile(epsilon_s,d1,epsilon_f1,epsilon_c,n_eff): # computing mode profile
    gs = np.sqrt(n_eff**2-epsilon_s+0j)*2*np.pi
    kf1x = np.sqrt(epsilon_f1-n_eff**2+0j)*2*np.pi
    gc = np.sqrt(n_eff**2-epsilon_c+0j)*2*np.pi
    
    xs = np.linspace(-2, 0, num=101, endpoint=True)
    xf1 = np.linspace(0, d1, num=101, endpoint=True)
    xc = np.linspace(d1, d1+2, num=101, endpoint=True)
    
    Fs = np.exp(gs*xs)
    Gxs = -n_eff*Fs*2*np.pi
    Gzs = gs*Fs
    
    MTE = mTE(kf1x,xf1)
    Ff1 = MTE[0,0]*Fs[-1]+MTE[0,1]*Gzs[-1]
    Gxf1 = -n_eff*Ff1*2*np.pi
    Gzf1 = MTE[1,0]*Fs[-1]+MTE[1,1]*Gzs[-1]
    
    Fc = Ff1[-1]*np.exp(-gc*(xc-d1))
    Gxc = -n_eff*Fc*2*np.pi
    Gzc = -gc*Fc
        
    return np.concatenate((Fs,Ff1,Fc)),np.concatenate((Gxs,Gxf1,Gxc)),np.concatenate((-1j*Gzs,-1j*Gzf1,-1j*Gzc)),np.concatenate((xs,xf1,xc))

def n_eff(epsilon_s,d1,epsilon_f1,epsilon_c,initial_guess):
    def func(params):
        n_eff_real, n_eff_imag = params
        return 1/np.abs(RTE(epsilon_s,d1,epsilon_f1,epsilon_c,n_eff_real+1j*n_eff_imag)) 
    result = spo.minimize(func, [np.real(initial_guess),np.imag(initial_guess)], bounds = ((np.sqrt(epsilon_s), None), (0, None)), tol = 1.e-8)    

    return result.x[0]+1j*result.x[1]  
      
def initialize():
    epsilon_f1_imag_double.set(.4)
    
    calculate()

def calculate():
    epsilon_f1_imag = epsilon_f1_imag_double.get()
    a1.cla()       
    a1.plot(epsilon_f_imag,np.real(n_eff_f),'b')
    a1.plot([epsilon_f_imag[0],epsilon_f_imag[-1]],[n_eff_mode,n_eff_mode],'k:')
    a1.set_xlim([epsilon_f_imag[0],epsilon_f_imag[-1]])
    a1.set_xlabel(r'$\varepsilon_{\rm f}^{\prime\prime}$')
    a1.set_ylabel(r'$n^{\prime}_{\rm eff}$')

    a2.cla()                   
    a2.semilogy(epsilon_f_imag,np.imag(n_eff_f),'b')
    a2.set_xlim([epsilon_f_imag[0],epsilon_f_imag[-1]])
    a2.set_xlabel(r'$\varepsilon_{\rm f}^{\prime\prime}$')
    a2.set_ylabel(r'$n^{\prime\prime}_{\rm eff}$')

    a3.cla()
    a3bis.cla()         
    lns1 = a3.plot(x_mode-d1/2,np.abs(F_mode),'k:')
    n_eff_lossy = np.interp(epsilon_f1_imag, epsilon_f_imag, n_eff_f)
    a1.plot(epsilon_f1_imag,np.real(n_eff_lossy),'bo')
    a2.plot(epsilon_f1_imag,np.imag(n_eff_lossy),'bo')
    F,Gx,Gz,x = mode_profile(epsilon_s,d1,epsilon_f1_real+1j*epsilon_f1_imag,epsilon_c,n_eff_lossy) # compute lossy mode profile
    lns2 = a3.plot(x-d1/2,np.abs(F),'b')
    a3.set_xlabel(r'$x/\lambda$')
    a3.set_ylabel(r'$|E_y|/|E_y(x=0)|$')    
    lns3 = a3bis.plot([x[0]-d1/2,-d1/2,-d1/2,d1/2,d1/2,x[-1]-d1/2],[epsilon_s,epsilon_s,epsilon_f1_real,epsilon_f1_real,epsilon_c,epsilon_c],'g')
    a3.axvspan(-d1/2, d1/2, color='0.875')
    a3bis.annotate(r'$\varepsilon_{\rm f}^{\prime\prime}=$ '+str(round(epsilon_f1_imag,4)), xy=(0,(epsilon_s+epsilon_c)/2),horizontalalignment='center', verticalalignment='center')
    a3bis.set_ylabel(r'$\varepsilon^{\prime}$')
    a3.set_xlim([x[0]-d1/2,x[-1]-d1/2])
    a3.set_ylim([0,2])
    a3.legend(lns1+lns2+lns3,[r'$|E_y|$ ideal mode',r'$|E_y|$ lossy mode',r'$\varepsilon^{\prime}$'])
     
    a4.cla()      
    Sx = np.real(F*np.conj(Gz))
    Sz = np.real(-F*np.conj(Gx))
    Smax = np.amax(np.sqrt(Sx**2+Sz**2))
    a4.plot(x-d1/2,Sx/Smax,'b')
    a4.set_xlabel(r'$x/\lambda$')
    a4.set_ylabel(r'$S_x/|\mathbf{S}|_\mathrm{max}$')
    a4.plot([x[0]-d1/2,x[-1]-d1/2],[0,0],'k:')
    a4.set_xlim([x[0]-d1/2,x[-1]-d1/2])
    a4.set_ylim([-.02,.08])
    
    plt.tight_layout()
      
#    plt.savefig('lossy_mode.pdf',bbox_inches='tight',dpi=300, transparent=True)

    canvas.draw()
            
f = plt.figure(1,[8,4.75])
gs = mpl.gridspec.GridSpec(2, 2, width_ratios=[1, 3], height_ratios=[1, 1])
a1 = f.add_subplot(gs[0])
a2 = f.add_subplot(gs[2])
a3 = f.add_subplot(gs[1]) 
a3bis = a3.twinx()
a4 = f.add_subplot(gs[3])  
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

epsilon_s = 2.
d1 = 1.
epsilon_f1_real = 2.25
epsilon_c = 1.

epsilon_f_imag = np.linspace(0,.5, num=100)
n_eff_mode = np.real((n_eff(epsilon_s,d1,epsilon_f1_real,epsilon_c,1.47))) # compute film waveguide mode
F_mode,Gx,Gz,x_mode = mode_profile(epsilon_s,d1,epsilon_f1_real,epsilon_c,n_eff_mode) # compute film waveguide mode profile
vn_eff = np.vectorize(n_eff)
n_eff_f = vn_eff(epsilon_s,d1,epsilon_f1_real+1j*epsilon_f_imag,epsilon_c,1.47+1j*0.001)

epsilon_f1_imag_double = Tk.DoubleVar()

initialize()

row = 1
row = gui.create_slider_with_latex(mainframe,r'absorption in film $\varepsilon_{\rm f}'' =$',epsilon_f1_imag_double,0,.5,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)