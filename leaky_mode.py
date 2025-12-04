import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.optimize as spo
import tkinter as Tk
import sys
sys.path.append('./aux')
import gui_stuff as gui
import strat_stuff as strat

gui.set_rcParams()
title = "Leaky Dielectric Slab Waveguides - TE"
root = Tk.Tk()
root.title(title)

def RTE(epsilon_s,d1,epsilon_f1,d2,epsilon_f2,epsilon_c,n_eff,sign): 
    ksx = 1j*np.sqrt(n_eff**2-epsilon_s+0j)*2*np.pi
    kf1x = np.sqrt(epsilon_f1-n_eff**2+0j)*2*np.pi
    kf2x = np.sqrt(epsilon_f2-n_eff**2+0j)*2*np.pi
    kcx = sign*1j*np.sqrt(n_eff**2-epsilon_c+0j)*2*np.pi
    m1TE = strat.mTE(kf1x,d1)
    m2TE = strat.mTE(kf2x,d2)
    MTE = np.matmul(m2TE,m1TE)
    return (ksx*MTE[1,1]-kcx*MTE[0,0]-1j*MTE[1,0]-1j*ksx*kcx*MTE[0,1])/(ksx*MTE[1,1]+kcx*MTE[0,0]+1j*MTE[1,0]-1j*ksx*kcx*MTE[0,1])

def mode_profile(epsilon_s,d1,epsilon_f1,d2,epsilon_f2,epsilon_c,n_eff,sign): # computing mode profile
    gs = np.sqrt(n_eff**2-epsilon_s+0j)*2*np.pi
    kf1x = np.sqrt(epsilon_f1-n_eff**2+0j)*2*np.pi
    kf2x = np.sqrt(epsilon_f2-n_eff**2+0j)*2*np.pi
    gc = sign*np.sqrt(n_eff**2-epsilon_c+0j)*2*np.pi
    
    xs = np.linspace(-2, 0, num=101, endpoint=True)
    xf1 = np.linspace(0, d1, num=101, endpoint=True)
    xf2 = np.linspace(d1, d1+d2, num=101, endpoint=True) 
    xc = np.linspace(d1+d2, 3, num=501, endpoint=True)
    
    Fs = np.exp(gs*xs)
    Gxs = -n_eff*Fs*2*np.pi
    Gzs = gs*Fs
    
    MTE = strat.mTE(kf1x,xf1)
    Ff1 = MTE[0,0]*Fs[-1]+MTE[0,1]*Gzs[-1]
    Gxf1 = -n_eff*Ff1*2*np.pi
    Gzf1 = MTE[1,0]*Fs[-1]+MTE[1,1]*Gzs[-1]
    
    MTE = strat.mTE(kf2x,xf2-d1)
    Ff2 = MTE[0,0]*Ff1[-1]+MTE[0,1]*Gzf1[-1]
    Gxf2= -n_eff*Ff2*2*np.pi
    Gzf2 = MTE[1,0]*Ff1[-1]+MTE[1,1]*Gzf1[-1]    
    
    Fc = Ff2[-1]*np.exp(-gc*(xc-d1-d2))
    Gxc = -n_eff*Fc*2*np.pi
    Gzc = -gc*Fc
        
    return np.concatenate((Fs,Ff1,Ff2,Fc)),np.concatenate((Gxs,Gxf1,Gxf2,Gxc)),np.concatenate((-1j*Gzs,-1j*Gzf1,-1j*Gzf2,-1j*Gzc)),np.concatenate((xs,xf1,xf2,xc))

def n_eff(epsilon_s,d1,epsilon_f1,d2,epsilon_f2,epsilon_c,sign,initial_guess):
    def func(params):
        n_eff_real, n_eff_imag = params
        return 1/np.abs(RTE(epsilon_s,d1,epsilon_f1,d2,epsilon_f2,epsilon_c,n_eff_real+1j*n_eff_imag,sign)) 
    result = spo.minimize(func, [np.real(initial_guess),np.imag(initial_guess)], bounds = ((np.sqrt(epsilon_s), None), (0, None)), tol = 1.e-8)    

    return result.x[0]+1j*result.x[1]  
      
def initialize():
    d2_double.set(.2)
    show_mode.set("noshow")
    calculate()
    
def show_manual():
    gui.show_manual("man/leaky_mode.png",title)

def calculate():
    gui.change_cursor(root,"trek")
    d2 = d2_double.get()
    mode = show_mode.get()
    f.clf()  
    a1 = f.add_subplot(gs[0])    
    a1.plot(d,np.real(n_eff_d),'b')
    if mode == 'show':  
        a1.plot([d[0],d[-1]],[n_eff_mode,n_eff_mode],'k:')
    a1.set_xlim([d[0],d[-1]])
    a1.set_xlabel(r'$d_{\rm b}/\lambda$')
    a1.set_ylabel(r'$k^{\prime}c/\omega$')

    a2 = f.add_subplot(gs[2])                  
    a2.semilogy(d,np.imag(n_eff_d),'b')
    a2.set_xlim([d[0],d[-1]])
    a2.set_xlabel(r'$d_{\rm b}/\lambda$')
    a2.set_ylabel(r'$k^{\prime\prime}c/\omega$')
    
    a3 = f.add_subplot(gs[1]) 
    a3bis = a3.twinx()  
    if mode == 'show':
        lns1 = a3.plot(x_mode-d1/2,np.abs(F_mode),'k:')
    n_eff_leaky = np.interp(d2, d, n_eff_d)
    a1.plot(d2,np.real(n_eff_leaky),'bo')
    a2.plot(d2,np.imag(n_eff_leaky),'bo')
    F,Gx,Gz,x = mode_profile(epsilon_s,d1,epsilon_f1,d2,epsilon_f2,epsilon_c,n_eff_leaky,-1) # compute leaky mode profile
    lns2 = a3.plot(x-d1/2,np.abs(F),'b')
    a3.set_xlabel(r'$x/\lambda$')
    a3.set_ylabel(r'$|E_y|$ $[|E_y(x=0)|\equiv1]$')    
    lns3 = a3bis.plot([x[0]-d1/2,-d1/2,-d1/2,d1-d1/2,d1-d1/2,d1/2+d2,d1/2+d2,x[-1]-d1/2],[epsilon_s,epsilon_s,epsilon_f1,epsilon_f1,epsilon_f2,epsilon_f2,epsilon_c,epsilon_c],'g')
    if mode == 'show':
        lns4 = a3bis.plot([x[0]-d1/2,-d1/2,-d1/2,d1-d1/2,d1-d1/2,x[-1]-d1/2],[epsilon_s,epsilon_s,epsilon_f1,epsilon_f1,epsilon_f2,epsilon_f2],'g:')
    a3bis.annotate(r'', xy=(d1/2,1.0), xytext=(d1/2-.5,1.0), arrowprops=dict(arrowstyle='->'))
    a3bis.annotate(r'', xy=(d1/2+d2,1.0), xytext=(d1/2+d2+.5,1.0), arrowprops=dict(arrowstyle='->'))
    a3bis.annotate(r'$d_{\rm b}$', xy=(d1/2+d2/2,1.0),horizontalalignment='center', verticalalignment='bottom')
    a3bis.annotate(r'', xy=(-d1/2,epsilon_f1), xytext=(-d1/2-.5,epsilon_f1), arrowprops=dict(arrowstyle='->'))
    a3bis.annotate(r'', xy=(d1/2,epsilon_f1), xytext=(d1/2+.5,epsilon_f1), arrowprops=dict(arrowstyle='->'))
    a3bis.annotate(r'$d_{\rm f}$', xy=(0,epsilon_f1),horizontalalignment='center', verticalalignment='bottom')
    a3bis.set_ylabel(r'$\varepsilon$')
    a3bis.set_ylim([a3bis.get_ylim()[0],a3bis.get_ylim()[1]*1.035])
    a3.set_xlim([x[0]-d1/2,x[-1]-d1/2])
    a3.set_ylim([0,2])
    if mode == 'show':
        a3.legend(lns2+lns3+lns1+lns4,[r'$|E_y|$ leaky mode',r'$\varepsilon$ leaky structure',r'$|E_y|$ guided mode',r'$\varepsilon$ guiding structure'],loc = 6)
    else:
        a3.legend(lns2+lns3,[r'$|E_y|$ leaky mode',r'$\varepsilon$ leaky structure'],loc = 6)
         
    a4 = f.add_subplot(gs[3])  
    Sx = np.real(F*np.conj(Gz))
    Sz = np.real(-F*np.conj(Gx))
    Smax = np.amax(np.sqrt(Sx**2+Sz**2))
    a4.plot(x-d1/2,Sx/Smax,'b')
    a4.set_xlabel(r'$x/\lambda$')
    a4.set_ylabel(r'$S_x/|\mathbf{S}|_\mathrm{max}$')
    a4.plot([x[0]-d1/2,x[-1]-d1/2],[0,0],'k:')
    a4.set_ylim([-.01,.12])
    a4.set_xlim([x[0]-d1/2,x[-1]-d1/2])
    
    plt.tight_layout() 
    
#    plt.savefig('leaky_mode.pdf',bbox_inches='tight',dpi=300, transparent=True)

    canvas.draw()
    gui.change_cursor(root,"arrow")
            
f = plt.figure(1,[8,4.75])
gs = mpl.gridspec.GridSpec(2, 2, width_ratios=[1, 3], height_ratios=[1, 1])

canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

epsilon_s = 2.
d1 = 1.
epsilon_f1 = 2.25
epsilon_f2 = 1.
epsilon_c = 2.2

d = np.linspace(0.05,.8, num=100)
n_eff_mode = np.real((n_eff(epsilon_s,d1,epsilon_f1,d1,epsilon_f2,epsilon_f2,1,1.47))) # compute film waveguide mode
F_mode,Gx,Gz,x_mode = mode_profile(epsilon_s,d1,epsilon_f1,d1,epsilon_f2,epsilon_f2,n_eff_mode,1) # compute film waveguide mode profile
vn_eff = np.vectorize(n_eff)
n_eff_d = vn_eff(epsilon_s,d1,epsilon_f1,d,epsilon_f2,epsilon_c,-1,1.465+1j*0.001)

d2_double = Tk.DoubleVar()
show_mode = Tk.StringVar()

initialize()

row = 1
row = gui.create_formula_with_latex(mainframe,r'substrate $\varepsilon_{\rm s} =$',r'$2.0$',row)
row = gui.create_formula_with_latex(mainframe,r'film $\varepsilon_{\rm f} =$',r'$2.25$',row)
row = gui.create_formula_with_latex(mainframe,r'film $\varepsilon_{\rm c_1} =$',r'$1.0$',row)
row = gui.create_formula_with_latex(mainframe,r'cladding $\varepsilon_{\rm c_2} =$',r'$2.2$',row)
row = gui.create_spacer(mainframe,row)
row = gui.create_formula_with_latex(mainframe,r'film thickness $d_{\rm f}/\lambda=$',r'$1.0$',row)
row = gui.create_slider_with_latex(mainframe,r'distance $d_{\rm b}/\lambda =$',d2_double,.05,.8,row,increment=.01)
row = gui.create_checkbutton(mainframe,"show guided mode",'noshow','show',show_mode,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)