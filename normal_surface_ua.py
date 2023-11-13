import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import aniso_stuff as ani

gui.set_rcParams()
root = Tk.Tk()
root.title("Normal Surfaces Uniaxial")

def initialize():
    var_string[0].set("1")   # epsilon_or
    var_string[1].set("2")   # epsilon_e
    var_string[2].set("no_show")   # show E
    var_string[3].set("no_show")   # show S
    var_double[0].set(1/3)  # theta/pi, defined by the normal to the optical axis (k3) and uk (Sect. 2.2.2)
    gui.copy_stringvar_vector(var_string,var_stringsave)
    calculate()
    
def reinitialize():
    gui.copy_stringvar_vector(var_stringsave,var_string)
    calculate() # because slider may have changed
    
def calculate():
    gui.change_cursor(root,"trek")
    try:
        epsilon = np.array([float(var_string[0].get()),float(var_string[0].get()),float(var_string[1].get())])
        theta = var_double[0].get()*np.pi
        theta0 = np.pi/2 - theta # get proper elevation angle defined by k3 and uk 
 
        if epsilon[0] <= 0 or epsilon[1] <= 0  or epsilon[2] <= 0: 
            gui.input_error("Tensor elements have to be positive. Re-initializing ...",reinitialize)
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
        
f = plt.figure(1,[5,5])

canvas = gui.create_canvas(root,f)
canvas.draw() # for faster feedback to user on startup
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(6)
var_stringsave = gui.create_stringvar_vector(6)
var_double = gui.create_doublevar_vector(1)

initialize()

row = 1
row = gui.create_entry_with_latex(mainframe,r"Dielectric tensor element $\varepsilon_{\rm or}=$",var_string[0],row)
row = gui.create_entry_with_latex(mainframe,r"Dielectric tensor element $\varepsilon_{\rm e}=$",var_string[1],row)
row = gui.create_slider_with_latex(mainframe,r'Angle of propagation direction $\theta/\pi=$',var_double[0],-1,1,row)
row = gui.create_label_with_latex(mainframe,r'index $n_{\rm or}=$',var_string[4],row)
row = gui.create_label_with_latex(mainframe,r'index $n_{\rm e}=$',var_string[5],row)
row = gui.create_double_checkbutton_with_latex(mainframe,r'show $\mathbf{E}^{\rm e}$','no_show','show',var_string[2],r'show $\mathbf{S}^{\rm e}$','no_show','show',var_string[3],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)