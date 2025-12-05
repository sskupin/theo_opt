import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import sys
sys.path.append('./stuff')
import gui_stuff as gui
import bpm_stuff as bpm

gui.set_rcParams()
title = "Bright Solitons"
root = Tk.Tk()
root.title(title)

def initialize():
    LZ_double.set(2)
    N_double.set(1)
    kappa_double.set(0)
    a_double.set(0)
    b_double.set(0)
    var_string[0].set("+1")     # D
    var_string[1].set("+1")     # Gamma   
    var_string[2].set("sech")   # profile  
    calculate()
    
def show_manual():
    gui.show_manual("man/soliton.png",title) 
    
def calculate():
    gui.change_cursor(root,"trek")
    Neta = 4096
    Leta = 40 
    NZ = 10001
    LZ = LZ_double.get()
    Nabs = 128
    N = N_double.get()
    kappa = kappa_double.get()
    D = float(var_string[0].get())
    a = a_double.get()
    Gamma = float(var_string[1].get())
    profile = var_string[2].get()
    b = b_double.get()
        
    eta, delta_eta = np.linspace(-Leta/2,Leta/2,Neta,endpoint=False, retstep=True)
    if profile == 'sech':
        u0 =  1/np.cosh(eta)*np.exp(1j*kappa*eta)
    elif profile == 'Gaussian':
        u0 = np.exp(-eta**2)*np.exp(1j*kappa*eta)
    else:
        u0 = np.exp(-eta**8)*np.exp(1j*kappa*eta)
    Z, delta_Z = np.linspace(0,LZ,NZ,endpoint=True, retstep=True)
    
    u = bpm.propagation_nls(Neta,u0,delta_eta,NZ,delta_Z*2*np.pi**2,Nabs,N,D,Gamma,a,b) # multiplication of delta_z by 2pi^2 to take into account scaling of LZ
    
    f.clf()
    
    a1 = plt.subplot2grid((5, 4), (1, 1), rowspan=4, colspan=2)
    im = a1.imshow(np.abs(np.transpose(u))**2 ,extent=[Z[0], Z[-1], eta[0], eta[-1]] , aspect='auto', origin='lower', vmin=0, cmap='jet')
    a1.set_xlabel(r'$Z/\pi$')
    a1.set_ylabel(r'$\eta$')
    ac = plt.subplot2grid((5, 4), (0, 1), colspan=2)
    f.colorbar(im, cax=ac, orientation='horizontal')
            
    a2 = plt.subplot2grid((5, 4), (1, 0), rowspan=4)       
    keta = 2*np.pi*np.fft.fftshift(np.fft.fftfreq(Neta,delta_eta))
    U0 = np.fft.fftshift(np.fft.ifft(np.fft.fftshift(u0))) # ifft because we assume time dependent problem
    lns1 = a2.plot(np.abs(U0)/np.max(np.abs(U0)),keta,'r', label=r'$|U_0(\mu)|$')
    a2.set_xlabel(r'$|u_0|/N$, $|U_0|/\mathrm{max}|U_0|$')
    a2.set_ylabel(r'$\mu$')
    a2.set_ylim([-20,20])
    a2.invert_yaxis()
    a2bis = a2.twinx()
    lns2 = a2bis.plot(np.abs(u0)/np.max(np.abs(u0)),eta, label=r'$|u_0(\eta)|$',color='b')
    a2bis.set_ylim([-10,10])
    a2bis.invert_yaxis()
    a2bis.set(yticklabels=[]) 
    lns = lns2+lns1
    labs = [l.get_label() for l in lns]
    a2.legend(lns, labs, bbox_to_anchor=(0, 1.05, 1, 0), loc="lower right")
    
    a3 = plt.subplot2grid((5, 4), (1, 3), rowspan=4, sharey=a1)
    lns1 = a3.plot(np.abs(u[-1,:])/np.max(np.abs(u0)),eta,'b', label=r'$|u(\eta,Z_L)|$')
    a3.set_xlabel(r'$|u|/N, |U|/\mathrm{max}|U_0|$')
    a3.set_ylabel(r'$\eta$')
    a3.set_ylim([-10,10])
    a3.invert_yaxis()   
    a3bis = a3.twinx()
    U = np.fft.fftshift(np.fft.ifft(np.fft.fftshift(u[-1,:]))) # ifft because we assume time dependent problem
    lns2 = a3bis.plot(np.abs(U)/np.max(np.abs(U0)),keta,'r', label=r'$|U(\mu,Z_L)|$')
    a3bis.set_ylim([-20,20])
    a3bis.set_ylabel(r'$\mu$')
    a3bis.invert_yaxis()
    lns = lns1+lns2
    labs = [l.get_label() for l in lns]
    a3.legend(lns, labs, bbox_to_anchor=(0, 1.05, 1, 0), loc="lower right")
        
    plt.tight_layout()

    ac.set_position([0.35,0.825,0.3,0.025])  
    ac.xaxis.tick_top()
    ac.set_xlabel(r'$|u|^2$')
    ac.xaxis.set_label_position('top') 
     
#    plt.savefig('soliton.pdf',bbox_inches='tight',dpi=300, transparent=True)
            
    canvas.draw()
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[10,5])

canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

LZ_double = Tk.DoubleVar()
N_double = Tk.DoubleVar()
kappa_double = Tk.DoubleVar()
a_double = Tk.DoubleVar()
b_double = Tk.DoubleVar()
var_string = gui.create_stringvar_vector(3)

initialize()

row = 1
row = gui.create_slider_with_latex(mainframe,r'propagation length $Z_L/\pi=$',LZ_double,0.5,2,row,increment=.025)
row = gui.create_formula_with_latex(mainframe,r'$\mathrm{i}\partial_Z u +\frac{1}{2}\partial^2_\eta u + |u|^2 u = $',r'$\mathrm{i}a \partial^3_\eta u - \mathrm{i}\frac{b}{2}u$',row)
row = gui.create_radiobutton_single_column(mainframe,[u'input pulse profile:','sech','Gaussian','super-Gaussian'],var_string[2],3,row)
row = gui.create_slider_with_latex(mainframe,r'amplitude (Soliton order) $N=$',N_double,0.5,3.5,row,increment=.05)
row = gui.create_slider_with_latex(mainframe,r'frequency shift $\kappa=$',kappa_double,-0.5,0.5,row,increment=.025)
row = gui.create_slider_with_latex(mainframe,r'third-order dispersion $a=$',a_double,-0.2,0.2,row,increment=.01)
row = gui.create_slider_with_latex(mainframe,r'linear losses $b=$',b_double,0,.1,row,increment=.01)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)