import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import sys
sys.path.append('./aux')
import gui_stuff as gui
import film_stuff as film

gui.set_rcParams()
title = "Dielectric Strip Waveguide - Effective-Index Method"
root = Tk.Tk()
root.title(title)

def n_eff_d_cladding(epsilon_f,epsilon_med,epsilon_min,d,polarization_1,n_eff_1,p):
    n_eff_1c_d = np.sqrt(epsilon_med)
    if film.number_of_modes(epsilon_f,epsilon_med,epsilon_min,d,polarization_1) > p:
        n_eff_1c_d = film.n_eff_d(epsilon_f,epsilon_med,epsilon_min,n_eff_1,d,p,polarization_1)[0]
    return n_eff_1c_d

def initialize():
    var_string[0].set("TE")
    var_string[1].set("11.825")
    var_string[2].set("11.56")
    var_string[3].set("1.")
    var_string[4].set("1")
    var_string[5].set("4")
    var_string[6].set(".5")
    gui.copy_stringvar_vector(var_string,var_save)   
    calculate()
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    calculate() 

def show_manual():
    gui.show_manual("man/strip_wg_ei.png",title) 

def calculate():
    gui.change_cursor(root,"trek")
    try:
        polarization_1 = var_string[0].get()
        epsilon_f = float(var_string[1].get())
        epsilon_s = float(var_string[2].get())
        epsilon_c = float(var_string[3].get())
        h = float(var_string[4].get())
        w = float(var_string[5].get())
        d = float(var_string[6].get())
        if polarization_1 == 'TE':
            polarization_2 = 'TM'
            title = 'quasi-TE'
        else:
            polarization_2 = 'TE'
            title = 'quasi-TM'
            
        if w <= 0 or h <= 0 or d < 0 or d >= 0.9*h:
            gui.input_error("Some waveguide dimensions are invalid. Re-initializing with previous parameters...",reinitialize)
        elif epsilon_c <= 0 or epsilon_s <= 0 or epsilon_f <= epsilon_s or epsilon_f <= epsilon_c: 
            gui.input_error("All dielectric constants have to be positive, with the film one being the largest. Re-initializing with previous parameters...",reinitialize)
        else:  
            epsilon_min = np.minimum(epsilon_s,epsilon_c)
            epsilon_med = np.maximum(epsilon_s,epsilon_c)
            n_eff_1 = np.linspace(np.sqrt(epsilon_med), np.sqrt(epsilon_f), num=200, endpoint=False)
            vinvdisp = np.vectorize(film.invdisp)
            pmax = film.number_of_modes(epsilon_f,epsilon_med,epsilon_min,h,polarization_1)
            number_of_2D_modes = 0
            if pmax <= 12:
                for p in range(pmax):
                    n_eff_1f_h = film.n_eff_d(epsilon_f,epsilon_med,epsilon_min,n_eff_1,h,p,polarization_1)[0]
                    n_eff_1c_d = n_eff_d_cladding(epsilon_f,epsilon_med,epsilon_min,d,polarization_1,n_eff_1,p)
                    qmax = film.number_of_modes(n_eff_1f_h**2,n_eff_1c_d**2,n_eff_1c_d**2,2*w,polarization_2)
                    number_of_2D_modes = number_of_2D_modes + qmax
            if number_of_2D_modes <= 0:
                gui.input_error("No mode found in the plotting window. Re-initializing with previous parameters...",reinitialize)    
            elif number_of_2D_modes > 12:
                gui.input_error("Too many modes found for readable plotting. Re-initializing with previous parameters...",reinitialize)
            else:
                f.clf()
                a1 = f.add_subplot(111)
                a1bis = a1.twinx()
                sqrtepsilon = [np.sqrt(epsilon_med),np.sqrt(epsilon_f)]
                if epsilon_med > epsilon_s:
                    labels = [r'$\sqrt{\varepsilon_{\rm c}}$', r'$\sqrt{\varepsilon_{\rm f}}$']  
                else:
                    labels = [r'$\sqrt{\varepsilon_{\rm s}}$', r'$\sqrt{\varepsilon_{\rm f}}$']                    
                a1bis.set_yticks(sqrtepsilon) 
                a1bis.set_yticklabels(labels)
                n_eff_max = np.sqrt(epsilon_f)
                n_eff_2 = np.zeros(len(n_eff_1))
                for p in range(pmax):
                    n_eff_1f_h = film.n_eff_d(epsilon_f,epsilon_med,epsilon_min,n_eff_1,h,p,polarization_1)[0]
                    n_eff_1c_d = n_eff_d_cladding(epsilon_f,epsilon_med,epsilon_min,d,polarization_1,n_eff_1,p)
                    n_eff_2 = np.linspace(n_eff_1c_d, n_eff_1f_h, num=200, endpoint=False)
                    qmax = film.number_of_modes(n_eff_1f_h**2,n_eff_1c_d**2,n_eff_1c_d**2,2*w,polarization_2) 
                    a1.plot([0, 2*w], [n_eff_1c_d, n_eff_1c_d], 'k--', label=r'cut-off $p=%2d$'%p)
                    for q in range(qmax):
                        a1.plot(np.real(q/(2*np.sqrt(n_eff_1f_h**2-n_eff_2**2+0j))+vinvdisp(n_eff_1f_h**2,n_eff_1c_d**2,n_eff_1c_d**2,n_eff_2,w,polarization_2)), n_eff_2, label=r'$p=%2d$'%p+r', $q=%2d$'%q)
                        n_eff_2_d = film.n_eff_d(n_eff_1f_h**2,n_eff_1c_d**2,n_eff_1c_d**2,n_eff_2,w,q,polarization_2)
                        if len(n_eff_2_d)>0:
                                a1.plot(w, n_eff_2_d[0], 'o',[w,2*w],[n_eff_2_d[0],n_eff_2_d[0]],':',color='gray')
                                labels = list(a1bis.get_yticklabels()) + [str(round(n_eff_2_d[0],4))]
                                a1bis.set_yticks(list(a1bis.get_yticks()) + [n_eff_2_d[0]])
                                a1bis.set_yticklabels(labels)
                a1.axis('tight')
                a1.set_xlim([0,2*w])
                a1.set_ylim([n_eff_1[0],n_eff_max])
                a1.set_xlabel(r'$w/\lambda$')
                a1.set_ylabel(r'$n_{\rm eff}$')
                a1.text(.5,.925,title,horizontalalignment='center',transform=a1.transAxes)
                a1.legend()
                a1bis.set_ylim([n_eff_1[0],n_eff_max])
                
#                plt.savefig('strip_waveguide_effective_index.pdf',bbox_inches='tight',dpi=300, transparent=True)
                
                gui.copy_stringvar_vector(var_string,var_save)

                canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing with previous parameters...", reinitialize)
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[6,5])
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(7)
var_save = gui.create_stringvar_vector(7)

initialize()

row = 1
row = gui.create_radiobutton(mainframe,['quasi-polarization:','TE','TM'],var_string[0],2,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry_with_latex(mainframe,r'film $\varepsilon_{\rm f} =$',var_string[1],row)
row = gui.create_entry_with_latex(mainframe,r'substrate $\varepsilon_{\rm s} =$',var_string[2],row)
row = gui.create_entry_with_latex(mainframe,r'cladding $\varepsilon_{\rm c} =$',var_string[3],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry_with_latex(mainframe,r'structure height $h/\lambda =$',var_string[4],row)
row = gui.create_entry_with_latex(mainframe,r'structure width $w/\lambda =$',var_string[5],row)
row = gui.create_entry_with_latex(mainframe,r'film thickness $d/\lambda =$',var_string[6],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)