import numpy as np
import EMpy
import matplotlib.pyplot as plt
import tkinter as Tk
import sys
sys.path.append('./stuff')
import gui_stuff as gui

gui.set_rcParams()
title = "Dielectric Strip Waveguide - Finite-Difference Method"
root = Tk.Tk()
root.title(title)

all_indices = np.zeros((24,2)) # mutable 

def stretch(z,zl,zr,factor):
    zstretch = np.copy(z)
    i, = np.where(z < zl)
    a = factor*z[0]/(z[i[-1]]-z[i[0]])**2
    for index in np.nditer(i):
        zstretch[index] = z[index] + a*(z[i[-1]]-z[index])**2
    j, = np.where(z > zr)
    a = factor*z[-1]/(z[j[-1]]-z[j[0]])**2
    for index in np.nditer(j):
        zstretch[index] = z[index] + a*(z[index]-z[j[0]])**2
    return zstretch

def initialize():
    var_string[0].set("11.825")
    var_string[1].set("11.56")
    var_string[2].set("1.")
    var_string[3].set("1")
    var_string[4].set("4")
    var_string[5].set(".5")
    var_string[6].set("2")
    gui.copy_stringvar_vector(var_string,var_save)
    calculate()
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    calculate() 
    
def show_manual():
    gui.show_manual("man/strip_wg_fd.png",title) 

def draw_interfaces(ax,epsilon_s,epsilon_c,h,d,w,y,color='w'):    
    if epsilon_s == epsilon_c and d == 0:
        ax.plot([-w/2,w/2],[0,0],color,[-w/2,w/2],[h-d,h-d],color,[-w/2,-w/2],[0,h-d],color,[w/2,w/2],[0,h-d],color)
    else:
        ax.plot([y[1],y[-2]],[-d,-d],color,[y[1],-w/2],[0,0],color,[w/2,y[-2]],[0,0],color,[-w/2,w/2],[h-d,h-d],color)
        ax.plot([-w/2,-w/2],[0,h-d],color,[w/2,w/2],[0,h-d],color)

def show_indices():
    top = Tk.Toplevel()
    top.title("eff. indices of all computed modes")
    scrollbar = Tk.Scrollbar(top, orient="vertical")
    listbox = Tk.Listbox(top, width=60, height=20, yscrollcommand=scrollbar.set)
    scrollbar.config(command=listbox.yview)
    scrollbar.pack(side="right", fill="y")
    listbox.pack(side="left",fill="both", expand=True)
    for index in range(all_indices.shape[0]):
        if all_indices[index,0] != 0:
            if all_indices[index,1] == 0:
                listbox.insert(Tk.END, "mode number "+str(index)+": eff. index = "+str(np.around(all_indices[index,0],decimals=5)))
            else:
                listbox.insert(Tk.END, "mode number "+str(index)+": eff. index = "+str(np.around(all_indices[index,0],decimals=5))+" (spurious mode)")
    top.mainloop() 

def calculate():
    gui.change_cursor(root,"trek")
    try:        
        epsilon_f = float(var_string[0].get())
        epsilon_s = float(var_string[1].get())
        epsilon_c = float(var_string[2].get())
        h = float(var_string[3].get())
        w = float(var_string[4].get())
        d = float(var_string[5].get())
        mu = int(float(var_string[6].get()))
        var_string[6].set(mu)
        
        def epsfunc(x, y):
            xx, yy = np.meshgrid(x, y)
            eps = np.where((xx.T <= h-d) * (np.abs(yy.T) <= w/2), epsilon_f, epsilon_c)
            eps = np.where((xx.T <= 0), epsilon_f, eps)
            eps = np.where((xx.T <= -d), epsilon_s, eps)
            return eps

        if mu < 0 or mu > 23:
            gui.input_error("Mode index must lie between 0 and 10. Re-initializing with previous parameters...",reinitialize)
        elif w <= 0 or h <= 0 or d < 0:
            gui.input_error("Some waveguide dimensions are invalid. Re-initializing with previous parameters...",reinitialize)
        elif epsilon_c <= 0 or epsilon_s <= 0 or epsilon_f <= epsilon_s or epsilon_f <= epsilon_c: 
            gui.input_error("All dielectric constants have to be positive, with the film one being the largest. Re-initializing with previous parameters...",reinitialize)
        else:  
            if d >= 0.9*h:
                gui.input_error("Film thickness d must be smaller than 0.9h. Adjusting to d=0.8h ...")
                d = 0.8*h
                var_string[5].set(d)
            xlim = 80*(h-d)/(4*np.floor(10*(h-d)/h))
            x = np.linspace(-xlim, xlim, 193) 
            xc = EMpy.utils.centered1d(x)
            y = np.linspace(-w, w, 193)
            yc = EMpy.utils.centered1d(y)
            xstretch = stretch(x,-3*h/2,3*h/2,2)
            ystretch = stretch(y,-3*w/4,3*w/4,2)
            ystretchc = EMpy.utils.centered1d(ystretch)
            neigs = mu+1
            tol = 1e-6
            boundary = '0000'
            
            solver = EMpy.modesolvers.FD.VFDModeSolver(1, xstretch, ystretch, epsfunc, boundary).solve(neigs, tol)
            mask_x = (abs(x[1:-1])>1/100) * (abs(x[1:-1]-h+d)>1/100)
            mask_xc = (abs(xc)>1/100) * (abs(xc-h+d)>1/100)
            mask_y = (abs(y[1:-1]-w/2)>1/100) * (abs(y[1:-1]+w/2)>1/100)
            mask_yc = (abs(yc-w/2)>1/100) * (abs(yc+w/2)>1/100)
            mask_ystretch = (ystretch[1:-1]+w<0) | (ystretch[1:-1]-w>0)
            mask_ystretchc = (ystretchc+w<0) | (ystretchc-w>0)
            all_indices[:,0] = 0
            for index in range(neigs):
                Ex,Ey,Ez,Hx,Hy,Hz = solver.modes[index].get_fields_for_FDTD(x,y)
                Exstretch,Eystretch,Ezstretch,Hxstretch,Hystretch,Hzstretch = solver.modes[index].get_fields_for_FDTD(xstretch,ystretch)
                max_E = np.amax([np.amax(abs(Ex[mask_xc,:][:,mask_y])),np.amax(abs(Ey[mask_x,:][:,mask_yc])),np.amax(abs(Ez))])
                max_E_boundary = np.amax([np.amax(abs(Exstretch[:,mask_ystretch])),np.amax(abs(Eystretch[:,mask_ystretchc])),np.amax(abs(Ezstretch[:,mask_ystretch]))])
                if (solver.modes[index].neff > np.maximum(np.sqrt(epsilon_s),np.sqrt(epsilon_c))):
                    all_indices[index,0] = np.real(solver.modes[index].neff)
                    if max_E_boundary < 0.9*max_E: # real mode                        
                        all_indices[index,1] = 0
                    else: # spurious mode
                        all_indices[index,1] = 1
                else:
                    gui.input_error("Mode index too high. Reducing index...")
                    mu = index-1
                    var_string[6].set(mu)
                    break
            if mu<0:
                gui.input_error("No mode found. Re-initializing with previous parameters...",reinitialize)
            else:
                var_string[7].set(np.around(np.real(solver.modes[mu].neff),decimals=5))
                if all_indices[mu,1] == 0:
                    var_string[8].set("")
                else:
                    var_string[8].set("(spurious mode)")
                   
                f.clf()  

                Ex,Ey,Ez,Hx,Hy,Hz = solver.modes[mu].get_fields_for_FDTD(x,y)
                max_E = np.amax([np.amax(abs(Ex[mask_xc,:][:,mask_y])),np.amax(abs(Ey[mask_x,:][:,mask_yc])),np.amax(abs(Ez))])
                max_H = np.amax([np.amax(abs(Hx)),np.amax(abs(Hy)),np.amax(abs(Hz))])
                a1 = f.add_subplot(3, 2, 1)
                a1.imshow(abs(Ex), origin='lower', extent=[y[1], y[-2], xc[0], xc[-1]], vmin=0, vmax=max_E, aspect=1, cmap='jet')
                draw_interfaces(a1,epsilon_s,epsilon_c,h,d,w,y)
                a1.text(.95,.95,r'$|E_x|$', color='w', verticalalignment='top', horizontalalignment='right', transform=a1.transAxes)
                a1.set_ylabel(r'$x/\lambda$')
                a3 = f.add_subplot(3, 2, 3)
                a3.imshow(abs(Ey), origin='lower', extent=[yc[0], yc[-1], x[1], x[-2]], vmin=0, vmax=max_E, aspect=1, cmap='jet')
                draw_interfaces(a3,epsilon_s,epsilon_c,h,d,w,y)
                a3.text(.95,.95,r'$|E_y|$', color='w', verticalalignment='top', horizontalalignment='right', transform=a3.transAxes)
                a3.set_ylabel(r'$x/\lambda$')
                a5 = f.add_subplot(3, 2, 5)
                a5.imshow(abs(Ez), origin='lower', extent=[y[1], y[-2], x[1], x[-2]], vmin=0, vmax=max_E, aspect=1, cmap='jet')
                draw_interfaces(a5,epsilon_s,epsilon_c,h,d,w,y)
                a5.text(.95,.95,r'$|E_z|$', color='w', verticalalignment='top', horizontalalignment='right', transform=a5.transAxes)
                a5.set_xlabel(r'$y/\lambda$')
                a5.set_ylabel(r'$x/\lambda$')
                a2 = f.add_subplot(3, 2, 2)
                a2.imshow(abs(Hx), origin='lower', extent=[yc[0], yc[-1], x[1], x[-2]], vmin=0, vmax=max_H, aspect=1, cmap='jet')
                draw_interfaces(a2,epsilon_s,epsilon_c,h,d,w,y)
                a2.text(.95,.95,r'$|H_x|$', color='w', verticalalignment='top', horizontalalignment='right', transform=a2.transAxes)
                a4 = f.add_subplot(3, 2, 4)
                a4.imshow(abs(Hy), origin='lower', extent=[y[1], y[-2], xc[0], xc[-1]], vmin=0, vmax=max_H, aspect=1, cmap='jet')
                draw_interfaces(a4,epsilon_s,epsilon_c,h,d,w,y)
                a4.text(.95,.95,r'$|H_y|$', color='w', verticalalignment='top', horizontalalignment='right', transform=a4.transAxes)
                a6 = f.add_subplot(3, 2, 6)
                a6.imshow(abs(Hz), origin='lower', extent=[yc[0], yc[-1], xc[0], xc[-1]], vmin=0, vmax=max_H, aspect=1, cmap='jet')
                draw_interfaces(a6,epsilon_s,epsilon_c,h,d,w,y)
                a6.text(.95,.95,r'$|H_z|$', color='w', verticalalignment='top', horizontalalignment='right', transform=a6.transAxes)
                a6.set_xlabel(r'$y/\lambda$')
                
                plt.tight_layout()
            
#                plt.savefig('strip_waveguide_finite_difference.pdf',bbox_inches='tight',dpi=300, transparent=True)

                gui.copy_stringvar_vector(var_string,var_save)

                canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing with previous parameters...", reinitialize)
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[6.5,6.5])
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(9)
var_save = gui.create_stringvar_vector(9)

initialize()

row = 1
row = gui.create_entry_with_latex(mainframe,r'film $\varepsilon_{\rm f} =$',var_string[0],row)
row = gui.create_entry_with_latex(mainframe,r'substrate $\varepsilon_{\rm s} =$',var_string[1],row)
row = gui.create_entry_with_latex(mainframe,r'cladding $\varepsilon_{\rm c} =$',var_string[2],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry_with_latex(mainframe,r'structure height $h/\lambda =$',var_string[3],row)
row = gui.create_entry_with_latex(mainframe,r'structure width $w/\lambda =$',var_string[4],row)
row = gui.create_entry_with_latex(mainframe,r'film thickness $d/\lambda =$',var_string[5],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry_with_latex(mainframe,r'mode index $\mu=$',var_string[6],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_label_with_latex(mainframe,r'eff.\ index $n_{\rm eff}=$',var_string[7],row)
row = gui.create_label(mainframe,"",var_string[8],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"other eff. indices",show_indices,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)