import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import sys
sys.path.append('./stuff')
import gui_stuff as gui
import bpm_stuff as bpm

gui.set_rcParams()
title = "Scalar Beam Propagation in Homogeneous and Isotropic Media"
root = Tk.Tk()
root.title(title)

def initialize():
    var_string[0].set("1")   # \Re \varepsilon
    var_string[1].set("0")   # \Im \varepsilon
    var_string[2].set("exact")       
    var_double[0].set(10)    # w_0/\lambda
    var_double[1].set(1)     # \alpha
    var_double[2].set(1)     # \beta
    var_double[3].set(2)     # z/L_F  
    calculate()
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    calculate() # because sliders may have changed
    
def show_manual():
    gui.show_manual("man/scalar_beam_prop.png",title)
    
def calculate():
    gui.change_cursor(root,"trek")
    try:  
        epsilon = float(var_string[0].get()) + 1j*float(var_string[1].get())    
        k = 2*np.pi*var_double[0].get()*np.sqrt(epsilon) # in units of 1/w_0
        alpha = var_double[1].get()
        beta = var_double[2].get()
        z = var_double[3].get()*2*np.pi*var_double[0].get() # in units of w_0
        propagation = var_string[2].get()
        
        if epsilon.real < 1: 
            gui.input_error("Real part of relative permittivity has to be larger than one. Re-initializing ...",reinitialize)  
        elif epsilon.imag < 0: 
            gui.input_error("Imaginary part of relative permittivity must not be negative. Re-initializing ...",reinitialize)  
        else:

            f.clf()
            
            x_calc,y_calc,Nx,Ny,kx_calc,ky_calc,Nkx,Nky,Nx0,Ny0 = bpm.init_2D_grid(k,z,alpha=alpha,beta=beta)
            
            u0_calc = bpm.init_2D_beam(x_calc,y_calc,alpha=alpha,beta=beta)
            
            if propagation == 'paraxial':
                prop = bpm.init_prop_2D_parax(kx_calc,ky_calc,k,z)
            else:
                prop = bpm.init_prop_2D(kx_calc,ky_calc,k,z)
            
            u_calc,U_calc,U0_calc = bpm.propagation_2D(u0_calc,prop)
            
            u0,x0,y0 = bpm.reshaping_2D(u0_calc,x_calc,y_calc,Nx0,Ny0)
            U0,kx,ky = bpm.reshaping_2D(U0_calc,kx_calc,ky_calc,Nkx,Nky)
            u,x,y = bpm.reshaping_2D(u_calc,x_calc,y_calc,Nx,Ny)
            U,kx,ky = bpm.reshaping_2D(U_calc,kx_calc,ky_calc,Nkx,Nky)
    
            a1 = f.add_subplot(221)
            im1 = a1.imshow(np.abs(u0) ,extent=[x0[0], x0[-1], y0[0], y0[-1]] , aspect='equal', origin='lower', vmin=0, cmap='jet')
            a1.annotate(r'$|u_0|/|u_0|_{\rm max}$', xy=(0.9*x0[0],0.75*y0[-1]),horizontalalignment='left', verticalalignment='bottom', color='w', size=14)
            a1.set_xlabel(r'$x/w_0$')
            a1.set_ylabel(r'$y/w_0$')
            plt.colorbar(im1,location='top',shrink=0.75)
            
            a2 = f.add_subplot(222)
            im2 = a2.imshow(np.abs(U0)/np.amax(np.abs(U0)) ,extent=[kx[0], kx[-1], ky[0], ky[-1]] , aspect='equal', origin='lower', vmin=0, cmap='jet')
            a2.annotate(r'$|U_0|/|U_0|_{\rm max}$', xy=(0.9*kx[0],0.75*ky[-1]),horizontalalignment='left', verticalalignment='bottom', color='w', size=14)
            a2.set_xlabel(r'$k_x w_0$')
            a2.set_ylabel(r'$k_y w_0$')
            plt.colorbar(im2,location='top',shrink=0.75)
            
            a3 = f.add_subplot(223)
            im3 = a3.imshow(np.abs(u) ,extent=[x[0], x[-1], y[0], y[-1]] , aspect='equal', origin='lower', vmin=0, cmap='jet')
            a3.annotate(r'$|u(z)|/|u_0|_{\rm max}$', xy=(0.9*x[0],0.75*y[-1]),horizontalalignment='left', verticalalignment='bottom', color='w', size=14)
            a3.set_xlabel(r'$x/w_0$')
            a3.set_ylabel(r'$y/w_0$')
            plt.colorbar(im3,location='top',shrink=0.75)
            
            a4 = f.add_subplot(224)
            im4 = a4.imshow(np.abs(U)/np.amax(np.abs(U0)) ,extent=[kx[0], kx[-1], ky[0], ky[-1]] , aspect='equal', origin='lower', vmin=0, cmap='jet')
            a4.annotate(r'$|U(z)|/|U_0|_{\rm max}$', xy=(0.9*kx[0],0.75*ky[-1]),horizontalalignment='left', verticalalignment='bottom', color='w', size=14)
            a4.set_xlabel(r'$k_x w_0$')
            a4.set_ylabel(r'$k_y w_0$')
            plt.colorbar(im4,location='top',shrink=0.75)
            
            plt.tight_layout()
                
#            plt.savefig('scalar_beam_prop.pdf',bbox_inches='tight',dpi=300, transparent=True)
            
            gui.copy_stringvar_vector(var_string,var_save)

            canvas.draw()       
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[6,6])

canvas = gui.create_canvas(root,f)
canvas.draw() # for faster feedback to user on startup
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(3)
var_save = gui.create_stringvar_vector(3)
var_double = gui.create_doublevar_vector(4)

initialize()

row = 1
row = gui.create_formula_with_latex(mainframe,r'$u_0 \propto $',r'$\exp\!\left\{-\left[\left(\frac{x}{w_0}\right)^{2} + \left(\beta\frac{y}{w_0}\right)^{2} \right]^{\alpha}\right\}$',row)
row = gui.create_slider_with_latex(mainframe,r'beam width $w_0/\lambda=$',var_double[0],1,20,row,increment=.5)
row = gui.create_slider_with_latex(mainframe,r'degree of super-Gaussian $\alpha=$',var_double[1],1,4,row,increment=.1)
row = gui.create_slider_with_latex(mainframe,r'beam ellipticity parameter $\beta=$',var_double[2],1,2,row,increment=.1)
row = gui.create_slider_with_latex(mainframe,r'propagation distance $z/L_{\rm F}=z\lambda/(2 \pi w_0^2)=$',var_double[3],0,5,row,increment=.25)
row = gui.create_double_entry_with_latex(mainframe,r"relative permittivity: $\varepsilon'$ =",var_string[0],r"$\varepsilon''$ =",var_string[1],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_radiobutton(mainframe,['propagation:','paraxial','exact'],var_string[2],2,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)