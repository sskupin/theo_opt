import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import aniso_stuff as ani

gui.set_rcParams()
root = Tk.Tk()
root.title("Normal Surfaces")

def initialize():
    global var_save
    var_string[0].set("1")   # epsilon1
    var_string[1].set("2")   # epsilon2
    var_string[2].set("3")   # epsilon3
    var_string[3].set("no_show")   # show E
    var_string[4].set("no_show")   # show S
    theta0_double.set(1/4)
    phi0_double.set(1/4)
    theta_view_double.set(1/3)
    phi_view_double.set(0.41)
    calculate()
    
def reinitialize():
    global var_string
    gui.copy_stringvar_vector(var_save,var_string)
    calculate()
    
def calculate():
    global var_save
    try:
        epsilon = np.array([float(var_string[0].get()),float(var_string[1].get()),float(var_string[2].get())])
        theta0 = theta0_double.get()*np.pi
        phi0 = phi0_double.get()*np.pi
        theta_view = theta_view_double.get()*180
        phi_view = phi_view_double.get()*180
 
        if epsilon[0] <= 0 or epsilon[1] <= 0  or epsilon[2] <= 0: 
            gui.input_error("Tensor elements have to be positive. Re-initializing ...",reinitialize)
        else:
            
            f.clf()
        
            ax = f.add_subplot(1,1,1, projection='3d')
            ax.set_box_aspect(aspect = (1,1,1))
            ax.view_init(azim=phi_view, elev=90-theta_view)
            
            ani.plot_ns(ax,theta0,phi0,epsilon,var_string[3].get(),var_string[4].get())
            
            na,nb = ani.dr(ani.uk(theta0,phi0),epsilon)
            na_string.set(round(na[0],4))
            nb_string.set(round(nb[0],4))
                
            plt.savefig('normal_surface.pdf',bbox_inches='tight',dpi=300, transparent=True)
            
            gui.copy_stringvar_vector(var_string,var_save)

            canvas.draw()       
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", initialize)
        
f = plt.figure(1,[8,8])

canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(5)
var_save = gui.create_stringvar_vector(5)
theta0_double = Tk.DoubleVar()
phi0_double = Tk.DoubleVar()
theta_view_double = Tk.DoubleVar()
phi_view_double = Tk.DoubleVar()
na_string = Tk.StringVar()
nb_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_entry_with_latex(mainframe,r"Dielectric tensor element $\varepsilon_1=$",var_string[0],row)
row = gui.create_entry_with_latex(mainframe,r"Dielectric tensor element $\varepsilon_2=$",var_string[1],row)
row = gui.create_entry_with_latex(mainframe,r"Dielectric tensor element $\varepsilon_3=$",var_string[2],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r'Azimuth of propagation direction $\varphi_0/\pi=$',phi0_double,0,2,row)
row = gui.create_slider_with_latex(mainframe,r'Elevation of propagation direction $\theta_0/\pi=$',theta0_double,0,1,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_label_with_latex(mainframe,r'index $n_a=$',na_string,row)
row = gui.create_label_with_latex(mainframe,r'index $n_b=$',nb_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_slider_with_latex(mainframe,r'Azimuth of view $\varphi_{\rm view}/\pi=$',phi_view_double,0,2,row)
row = gui.create_slider_with_latex(mainframe,r'Elevation of view $\theta_{\rm view}/\pi=$',theta_view_double,0,1,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_checkbutton_with_latex(mainframe,r'show $\mathbf{E}^{a,b}$','no_show','show',var_string[3],r'show $\mathbf{S}^{a,b}$','no_show','show',var_string[4],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)