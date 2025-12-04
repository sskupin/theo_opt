import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import sys
sys.path.append('./aux')
import gui_stuff as gui
import bpm_stuff as bpm

gui.set_rcParams()
title = "1D Waveguide Array"
root = Tk.Tk()
root.title(title)

def initialize():
    w0_double.set(1)
    tilt_double.set(0.25)
    calculate()
    
def show_manual():
    gui.show_manual("man/wg_array.png",title) 
    
def calculate():
    gui.change_cursor(root,"trek")
    w0 = w0_double.get()
    tilt = tilt_double.get()

    N = 32
    n = np.linspace(0,N,N,endpoint=False)
    a0 = np.exp(-(n-N/2)**2/w0**2+1j*n*tilt*np.pi)
    Nz = 201
    L = 10 # 2 \kappa z_max, normalized diffraction length for tilt = 0 is .5
    delta_z = L/(Nz-1)
        
    a = bpm.propagation_discrete(N,a0,Nz,delta_z)
    
    f.clf()
        
    a1 = f.add_subplot(221)
    a1.plot(n,np.abs(a[0,:]), 'b-o', label=r'input beam')
    a1.plot(n,np.abs(a[Nz-1,:]), 'r-o', label=r'output beam')
    a1.set_xlim([-0.5,N-.5])
    a1.set_xlabel(r'\# of waveguide $n$')
    a1.set_ylabel(r'$|\bar{a}_n|$ [arb.u.]') 
    a1.legend()

    a2 = f.add_subplot(223)
    U0 = np.fft.fftshift(np.fft.fft(a0,8*N))
    kxd = np.linspace(-np.pi,np.pi,8*N,endpoint=False)
    lns1 = a2.plot(kxd,np.abs(U0)/np.amax(np.abs(U0)),'b', label=r'input spect. $\hat{a}$')
    a2.set_xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi])
    a2.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$', r'$\pi$'])
    a2.set_xlim([-np.pi,np.pi])
    a2.set_xlabel(r'$k_x d$')
    a2.set_ylabel(r'$|\hat{a}|$ [arb.u.]') 
    a2bis = a2.twinx()
    lns2 = a2bis.plot(kxd,2*np.cos(kxd),'r', label=r'disp. rel. $\mathcal{K}_z$')
    a2bis.set_ylabel(r'$\mathcal{K}_z/\kappa$')
    lns = lns1+lns2
    labs = [l.get_label() for l in lns]
    a2.legend(lns, labs, loc=4)

    a3 = f.add_subplot(122)
    a3.imshow(np.abs(a),extent=[-0.5, N-.5, L, 0], aspect=6, origin='upper', cmap='jet')
    for index in range(0,N):
        a3.plot([.5+index,.5+index],[0,L],'k',linewidth=2)
    a3.set_xlim([-0.5,N-.5])
    a3.set_xlabel(r'\# of waveguide $n$')
    a3.set_ylim([L,0])
    a3.set_ylabel(r'$2 \kappa z$') 
            
    plt.tight_layout()
                
#    plt.savefig('waveguide_array.pdf',bbox_inches='tight',dpi=300, transparent=True)

    canvas.draw()   
    gui.change_cursor(root,"arrow")    

f = plt.figure(1,[7,5])

canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

w0_double = Tk.DoubleVar()
tilt_double = Tk.DoubleVar()

initialize()

row = 1
row = gui.create_title(mainframe,'input beam profile',row)
row = gui.create_formula_with_latex(mainframe,r'$\bar{a}_n(z=0)\propto$',r'$\exp\!\left[ -\frac{(n-16)^2}{W_0^2} +\mathrm{i}\, \pi n\, \alpha_0 \right]$',row)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r"beam width $W_0 =$",w0_double,0.1,5,row,increment=.1)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r"beam tilt $\alpha_0 =$",tilt_double,-1,1,row,increment=.05)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)