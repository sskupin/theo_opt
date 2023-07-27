import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import film_stuff as film

gui.set_rcParams()
root = Tk.Tk()
root.title("Dielectric Strip Waveguide by Effective Index Method")

def n_eff_d_cladding(epsilon_f,epsilon_med,epsilon_min,d,polarization_1,n_eff_1,p):
    n_eff_1c_d = np.sqrt(epsilon_med)
    if film.number_of_modes(epsilon_f,epsilon_med,epsilon_min,d,polarization_1) > p:
        n_eff_1c_d = film.n_eff_d(epsilon_f,epsilon_med,epsilon_min,n_eff_1,d,p,polarization_1)[0]
    return n_eff_1c_d

def initialize():
    global polarization_save,epsilon_f_save,epsilon_s_save,epsilon_c_save,h_string_save,w_string_save,d_string_save
    polarization_string.set("TE")
    epsilon_f_string.set("11.825")
    epsilon_s_string.set("11.56")
    epsilon_c_string.set("1.")
    h_string.set("1")
    w_string.set("4")
    d_string.set(".5")
    
    polarization_save = polarization_string.get()
    epsilon_f_save = epsilon_f_string.get()
    epsilon_s_save = epsilon_s_string.get()
    epsilon_c_save = epsilon_c_string.get()
    h_string_save = h_string.get()
    w_string_save = w_string.get()
    d_string_save = d_string.get()
    
    calculate()
    
def reinitialize():
    global polarization_save,epsilon_f_save,epsilon_s_save,epsilon_c_save,h_string_save,w_string_save,d_string_save
    polarization_string.set(polarization_save)
    epsilon_f_string.set(epsilon_f_save)
    epsilon_s_string.set(epsilon_s_save)
    epsilon_c_string.set(epsilon_c_save)    
    h_string.set(h_string_save)
    w_string.set(w_string_save)
    d_string.set(d_string_save)    

def calculate():
    global polarization_save,epsilon_f_save,epsilon_s_save,epsilon_c_save,h_string_save,w_string_save,d_string_save
    try:
        polarization_1 = polarization_string.get()
        epsilon_f = float(epsilon_f_string.get())
        epsilon_s = float(epsilon_s_string.get())
        epsilon_c = float(epsilon_c_string.get())
        h = float(h_string.get())
        w = float(w_string.get())
        d = float(d_string.get())
        if polarization_1 == 'TE':
            polarization_2 = 'TM'
            title = 'quasi-TE'
        else:
            polarization_2 = 'TE'
            title = 'quasi-TM'
            
        if w <= 0 or h <= 0 or d < 0 or d >= h:
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
                
                polarization_save = polarization_string.get()
                epsilon_f_save = epsilon_f_string.get()
                epsilon_s_save = epsilon_s_string.get()
                epsilon_c_save = epsilon_c_string.get()
                h_string_save = h_string.get()
                w_string_save = w_string.get()
                d_string_save = d_string.get()

                canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing with previous parameters...", reinitialize)

f = plt.figure(1,[6,5])
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

polarization_string = Tk.StringVar()
epsilon_f_string = Tk.StringVar()
epsilon_s_string = Tk.StringVar()
epsilon_c_string = Tk.StringVar()
h_string = Tk.StringVar()
w_string = Tk.StringVar()
d_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_radiobutton(mainframe,['quasi-polarization:','TE','TM'],polarization_string,2,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry_with_latex(mainframe,r'film $\varepsilon_{\rm f} =$',epsilon_f_string,row)
row = gui.create_entry_with_latex(mainframe,r'substrate $\varepsilon_{\rm s} =$',epsilon_s_string,row)
row = gui.create_entry_with_latex(mainframe,r'cladding $\varepsilon_{\rm c} =$',epsilon_c_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry_with_latex(mainframe,r'structure height $h/\lambda =$',h_string,row)
row = gui.create_entry_with_latex(mainframe,r'structure width $w/\lambda =$',w_string,row)
row = gui.create_entry_with_latex(mainframe,r'film thickness $d/\lambda =$',d_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)