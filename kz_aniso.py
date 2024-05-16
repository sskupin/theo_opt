import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import aniso_stuff as ani

gui.set_rcParams()
root = Tk.Tk()
root.title("Longitdinal wavevector components in anisotropic media")

def initialize():
    var_string[0].set("2")   # epsilon1
    var_string[1].set("3")   # epsilon2
    var_string[2].set("4")   # epsilon3
    var_double[0].set(0.3)  # theta0/pi
    var_double[1].set(0.4)  # phi0/pi
    calculate()
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    calculate() # because sliders may have changed
    
def calculate():
    gui.change_cursor(root,"trek")
    try:
        epsilon = np.array([float(var_string[0].get()),float(var_string[1].get()),float(var_string[2].get())])
        theta0 = var_double[0].get()*np.pi
        phi0 = var_double[1].get()*np.pi
        theta_view = 0.3*180
        phi_view = 0.8*180
 
        if epsilon[0] <= 0 or epsilon[1] <= 0  or epsilon[2] <= 0: 
            gui.input_error("Tensor elements have to be positive. Re-initializing ...",reinitialize)
        else:
            
            f.clf()
                    
            ax1 = f.add_subplot(1,2,1, projection='3d')
            ax1.view_init(azim=phi_view, elev=90-theta_view)
            ax2 = f.add_subplot(1,2,2, projection='3d')
            ax2.view_init(azim=phi_view, elev=90-theta_view)
            ax3 = f.add_subplot(1,2,2, projection='3d')
            ax3.view_init(azim=phi_view, elev=90-theta_view)
            ax3.set_position([0.35,0.7,0.3,0.3])
            
            ani.plot_kz(ax1,theta0,phi0,epsilon,'a')
            ani.plot_kz(ax2,theta0,phi0,epsilon,'b') 
            ani.plot_ns(ax3,theta0,phi0,epsilon,'no_show','no_show',True)
                
            plt.savefig('kz_aniso.pdf',bbox_inches='tight',dpi=300, transparent=True)
            
            gui.copy_stringvar_vector(var_string,var_save)

            canvas.draw()       
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")
       
f = plt.figure(1,[12,6])

canvas = gui.create_canvas(root,f)
canvas.draw() # for faster feedback to user on startup
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(3)
var_save = gui.create_stringvar_vector(7)
var_double = gui.create_doublevar_vector(2)

initialize()

row = 1
row = gui.create_entry_with_latex(mainframe,r"Dielectric tensor element $\varepsilon_1=$",var_string[0],row)
row = gui.create_entry_with_latex(mainframe,r"Dielectric tensor element $\varepsilon_2=$",var_string[1],row)
row = gui.create_entry_with_latex(mainframe,r"Dielectric tensor element $\varepsilon_3=$",var_string[2],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r'Azimuth of propagation direction $\varphi/\pi=$',var_double[1],0,2,row)
row = gui.create_slider_with_latex(mainframe,r'Elevation of propagation direction $\theta/\pi=$',var_double[0],0,1,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)