import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import bpm_stuff as bpm

gui.set_rcParams()
root = Tk.Tk()
root.title("Scalar Beam Propagation in Homogeneous and Isotropic Media")

def initialize():
    var_string[0].set("128")   # Nx
    var_string[1].set("25")    # Lx in w0
    var_string[2].set("100")   # Nz
    var_string[3].set("2")     # Lz in LF 
    var_string[4].set("0")     # Nabs
    var_string[5].set("1")     # alpha
    gui.copy_stringvar_vector(var_string,var_save)
    calculate()
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    
def calculate():
    gui.change_cursor(root,"trek")
    try:
        Nz = int(float(var_string[2].get())) + 1 # Number of points = Number of steps + 1
        var_string[2].set(Nz-1)
        z = float(var_string[3].get())
        Nabs = int(float(var_string[4].get()))
        var_string[4].set(Nabs)
        alpha = float(var_string[5].get())    
        
        if alpha < 1: 
            gui.input_error("Order of Super-Gaussian should be larger or equal one. Re-initializing ...",reinitialize)  
        else:

            f.clf()
            
            k = 10
            
            x,y,kx,ky = bpm.init_2D_grid(k,z,alpha=alpha)
            
            u0 = bpm.init_2D_beam(x,y,alpha=alpha)
            
            prop = bpm.init_prop_2D(kx,ky,k,z)
            
            u,U,U0 = bpm.propagation_2D(u0,prop)
    
            a1 = f.add_subplot(221)
            a1.imshow(np.abs(u0) ,extent=[x[0], x[-1], y[0], y[-1]] , aspect='equal', origin='lower', cmap='jet')
            a1.set_xlabel(r'$x/w_0$')
            a1.set_ylabel(r'$y/w_0$')
            
            a2 = f.add_subplot(222)
            a2.imshow(np.abs(U0) ,extent=[kx[0], kx[-1], ky[0], ky[-1]] , aspect='equal', origin='lower', cmap='jet')
            a2.set_xlabel(r'$k_x w_0$')
            a2.set_ylabel(r'$k_y w_0$')
            
            a3 = f.add_subplot(223)
            a3.imshow(np.abs(u) ,extent=[x[0], x[-1], y[0], y[-1]] , aspect='equal', origin='lower', cmap='jet')
            a3.set_xlabel(r'$x/w_0$')
            a3.set_ylabel(r'$y/w_0$')
            
            a4 = f.add_subplot(224)
            a4.imshow(np.abs(U) ,extent=[kx[0], kx[-1], ky[0], ky[-1]] , aspect='equal', origin='lower', cmap='jet')
            a4.set_xlabel(r'$k_x w_0$')
            a4.set_ylabel(r'$k_y w_0$')
            
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

var_string = gui.create_stringvar_vector(6)
var_save = gui.create_stringvar_vector(6)

initialize()

row = 1
row = gui.create_entry_with_latex(mainframe,r"Transverse number of points $N_x=$",var_string[0],row)
row = gui.create_entry_with_latex(mainframe,r"Transverse box length $L_x/w_0=$",var_string[1],row)
row = gui.create_entry_with_latex(mainframe,r"Longitudinal number of steps $N_z=$",var_string[2],row)
row = gui.create_entry_with_latex(mainframe,r"Longitudinal box length $L_z/L_{\rm F}=$",var_string[3],row)
row = gui.create_entry_with_latex(mainframe,r"Number of absorber points $N_{\rm abs}=$",var_string[4],row)
row = gui.create_formula_with_latex(mainframe,r'$\partial_z v = $',r'$\mathrm{i}\frac{1}{2k}\partial^2_x v $',row)
row = gui.create_formula_with_latex(mainframe,r'$v_0=$',r'$\exp\!\left[ -\left(\frac{x}{w_0}\right)^{2\alpha} \right]$',row)
row = gui.create_entry_with_latex(mainframe,r"Degree of super-Gaussian $\alpha=$",var_string[5],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)