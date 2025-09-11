import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import aniso_stuff as ani

gui.set_rcParams()
title = "Normal Surfaces Uniaxial"
root = Tk.Tk()
root.title(title)

def initialize():
    var_string[0].set("2")   # epsilon_or
    var_string[1].set("3")   # epsilon_e
    var_string[2].set("no_show")   # show E
    var_string[3].set("no_show")   # show S
    var_string[6].set("0.33")  # theta/pi, defined by the normal to the optical axis (k3) and uk (Sect. 2.2.2)
    gui.copy_stringvar_vector(var_string,var_stringsave)
    calculate()
    
def reinitialize():
    gui.copy_stringvar_vector(var_stringsave,var_string)

def show_manual():
    gui.show_manual("taylor_series.png",title)
    
def calculate():
    gui.change_cursor(root,"trek")
    try:
        epsilon = np.array([float(var_string[0].get()),float(var_string[0].get()),float(var_string[1].get())])
        theta = np.mod(float(var_string[6].get())+1,2)-1 # shift everything in the interval [-1,1]
        var_string[6].set(str(round(theta,len(var_string[6].get()))))
        theta = theta*np.pi
        theta0 = np.pi/2 - theta # get proper polar angle defined by k3 and uk 
 
        if (epsilon < 1).any() or (epsilon > 12).any(): 
            gui.input_error("Tensor elements must be between 1 and 12. Re-initializing ...",reinitialize)
        else:
            
            f.clf()
        
            ax = f.add_subplot(1,1,1)
            
            ani.plot_ns_uniaxial(ax,theta0,epsilon,var_string[2].get(),var_string[3].get())
            
            na,nb = ani.dr(ani.uk(theta0,0),epsilon)
            if epsilon[2] > epsilon[0]: # check if positive uniaxial
                var_string[4].set(round(na[0],4))
                var_string[5].set(round(nb[0],4))
            else:
                var_string[4].set(round(nb[0],4))
                var_string[5].set(round(na[0],4))  
                
#            plt.savefig('normal_surface_uniaxial.pdf',bbox_inches='tight',dpi=300, transparent=True)
            
            gui.copy_stringvar_vector(var_string,var_stringsave)

            canvas.draw()       
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")
        
f = plt.figure(1,[5,4])

canvas = gui.create_canvas(root,f)
canvas.draw() # for faster feedback to user on startup
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(7)
var_stringsave = gui.create_stringvar_vector(7)

initialize()

row = 1
row = gui.create_entry_with_latex(mainframe,r"Dielectric tensor element $\varepsilon_{\rm or}=$",var_string[0],row)
row = gui.create_entry_with_latex(mainframe,r"Dielectric tensor element $\varepsilon_{\rm e}=$",var_string[1],row)
row = gui.create_entry_with_latex(mainframe,r'Angle of \textbf{u}$^{\rm k}$, $\theta/\pi=$',var_string[6],row)
row = gui.create_label_with_latex(mainframe,r'index $n_{\rm or}=$',var_string[4],row)
row = gui.create_label_with_latex(mainframe,r'index $n_{\rm e}=$',var_string[5],row)
row = gui.create_double_checkbutton_with_latex(mainframe,r'show $\mathbf{E}^{\rm e}$','no_show','show',var_string[2],r'show $\mathbf{S}^{\rm e}$','no_show','show',var_string[3],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)