import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import aniso_stuff as ani

gui.set_rcParams()
title = "Index Ellipsoid in Anisotropic Crystals"
root = Tk.Tk()
root.title(title)

def initialize():
    var_string[0].set("2")   # epsilon1
    var_string[1].set("3")   # epsilon2
    var_string[2].set("4")   # epsilon3
    var_string[3].set("no_show")   # show E
    var_string[4].set("no_show")   # show H
    var_string[7].set("0.25")  # theta0/pi
    var_string[8].set("1.85")  # phi0/pi
    var_string[9].set("0.4")  # theta_view/180
    var_string[10].set("1.4") # phi_view/180
    calculate()
    
def show_manual():
    gui.show_manual("taylor_series.png",title)
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    calculate() # because sliders may have changed
    
def calculate():
    gui.change_cursor(root,"trek")
    try:
        epsilon = np.array([float(var_string[0].get()),float(var_string[1].get()),float(var_string[2].get())])
        theta0 = float(var_string[7].get())*np.pi
        phi0 = float(var_string[8].get())*np.pi
        theta_view = float(var_string[9].get())*180
        phi_view = float(var_string[10].get())*180
 
        if (epsilon < 1).any() or (epsilon > 12).any(): 
            gui.input_error("Tensor elements must be between 1 and 12. Re-initializing ...",reinitialize)
        else:
            
            f.clf()
        
            ax = f.add_subplot(1,1,1, projection='3d')
            ax.view_init(azim=phi_view, elev=90-theta_view)
            
            ani.plot_ie(ax,theta0,phi0,epsilon,var_string[3].get(),var_string[4].get())
 
            limit = np.sqrt(np.amax(epsilon))
            ax.set_xlim([-limit,limit])
            ax.set_ylim([-limit,limit])
            ax.set_zlim([-limit,limit])
            limits = np.array([getattr(ax, f'get_{axis}lim')() for axis in 'xyz'])
            ax.set_box_aspect(np.ptp(limits, axis = 1))
            ax.locator_params(nbins=4)
            
            na,nb = ani.dr(ani.uk(theta0,phi0),epsilon)
            var_string[5].set(round(na[0],4))
            var_string[6].set(round(nb[0],4))
                
#            plt.savefig('index_ellipsoid.pdf',bbox_inches='tight',dpi=300, transparent=True)
            
            gui.copy_stringvar_vector(var_string,var_save)

            canvas.draw()       
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")
       
f = plt.figure(1,[8,8])

canvas = gui.create_canvas(root,f)
canvas.draw() # for faster feedback to user on startup
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(11)
var_save = gui.create_stringvar_vector(11)

initialize()

row = 1
row = gui.create_entry_with_latex(mainframe,r"Dielectric tensor element $\varepsilon_1=$",var_string[0],row)
row = gui.create_entry_with_latex(mainframe,r"Dielectric tensor element $\varepsilon_2=$",var_string[1],row)
row = gui.create_entry_with_latex(mainframe,r"Dielectric tensor element $\varepsilon_3=$",var_string[2],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry_with_latex(mainframe,r'Azimuthal angle of \textbf{u}$^{\rm k}$, $\varphi/\pi=$',var_string[8],row)
row = gui.create_entry_with_latex(mainframe,r'Polar angle of \textbf{u}$^{\rm k}$, $\theta/\pi=$',var_string[7],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_label_with_latex(mainframe,r'index $n_a=$',var_string[5],row)
row = gui.create_label_with_latex(mainframe,r'index $n_b=$',var_string[6],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry_with_latex(mainframe,r'Azimuthal angle of view, $\varphi_{\rm view}/\pi=$',var_string[10],row)
row = gui.create_entry_with_latex(mainframe,r'Polar angle of view, $\theta_{\rm view}/\pi=$',var_string[9],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_checkbutton_with_latex(mainframe,r'show $\mathbf{E}^{a,b}$','no_show','show',var_string[3],r'show $\mathbf{H}^{a,b}$','no_show','show',var_string[4],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)