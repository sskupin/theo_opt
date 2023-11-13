import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import aniso_stuff as ani

gui.set_rcParams()
root = Tk.Tk()
root.title("Normal Surfaces")

def initialize():
    var_string[0].set("1")   # epsilon1
    var_string[1].set("2")   # epsilon2
    var_string[2].set("3")   # epsilon3
    var_string[3].set("no_show")   # show E
    var_string[4].set("no_show")   # show S
    var_double[0].set(1/4)  # theta0/pi
    var_double[1].set(1/4)  # phi0/pi
    var_double[2].set(1/3)  # theta_view/180
    var_double[3].set(0.41) # phi_view/180
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
        theta_view = var_double[2].get()*180
        phi_view = var_double[3].get()*180
 
        if epsilon[0] <= 0 or epsilon[1] <= 0  or epsilon[2] <= 0: 
            gui.input_error("Tensor elements have to be positive. Re-initializing ...",reinitialize)
        else:
            
            f.clf()
        
            ax = f.add_subplot(1,1,1, projection='3d')
            ax.set_box_aspect(aspect = (1,1,1))
            ax.view_init(azim=phi_view, elev=90-theta_view)
            
            ani.plot_ns(ax,theta0,phi0,epsilon,var_string[3].get(),var_string[4].get())
            
            na,nb = ani.dr(ani.uk(theta0,phi0),epsilon)
            var_string[5].set(round(na[0],4))
            var_string[6].set(round(nb[0],4))
                
#            plt.savefig('normal_surface.pdf',bbox_inches='tight',dpi=300, transparent=True)
            
            gui.copy_stringvar_vector(var_string,var_save)

            canvas.draw()       
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")
       
f = plt.figure(1,[6,6])

canvas = gui.create_canvas(root,f)
canvas.draw() # for faster feedback to user on startup
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(7)
var_save = gui.create_stringvar_vector(7)
var_double = gui.create_doublevar_vector(4)

initialize()

row = 1
row = gui.create_entry_with_latex(mainframe,r"Dielectric tensor element $\varepsilon_1=$",var_string[0],row)
row = gui.create_entry_with_latex(mainframe,r"Dielectric tensor element $\varepsilon_2=$",var_string[1],row)
row = gui.create_entry_with_latex(mainframe,r"Dielectric tensor element $\varepsilon_3=$",var_string[2],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r'Azimuth of propagation direction $\varphi/\pi=$',var_double[1],0,2,row)
row = gui.create_slider_with_latex(mainframe,r'Elevation of propagation direction $\theta/\pi=$',var_double[0],0,1,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_label_with_latex(mainframe,r'index $n_a=$',var_string[5],row)
row = gui.create_label_with_latex(mainframe,r'index $n_b=$',var_string[6],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r'Azimuth of view $\varphi_{\rm view}/\pi=$',var_double[3],0,2,row)
row = gui.create_slider_with_latex(mainframe,r'Elevation of view $\theta_{\rm view}/\pi=$',var_double[2],0,1,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_checkbutton_with_latex(mainframe,r'show $\mathbf{E}^{a,b}$','no_show','show',var_string[3],r'show $\mathbf{S}^{a,b}$','no_show','show',var_string[4],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)