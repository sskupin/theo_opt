import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import tkinter as Tk
import gui_stuff as gui
import film_stuff as film

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rc('text', usetex=True)
mpl.rc('text.latex', preamble=r'\usepackage{cmbright}')
mpl.rcParams.update({'font.size': 10})

root = Tk.Tk()
root.title("Slab or Film Waveguide")

def plot_mode(ax,x,d,F):
    ax.plot(x*d,F,'b')
    ax.set_xlim([x[0]*d,x[-1]*d])
    ylimits = ax.get_ylim()
    ax.set_ylim([ylimits[0],ylimits[1]+0.05*(ylimits[1]-ylimits[0])])
    ax.plot([-0.5*d,-0.5*d],ax.get_ylim(),'k:',[0.5*d,0.5*d],ax.get_ylim(),'k:')
    ax.set_xlabel(r'$x/\lambda$')
    ax.text(0.5, 0.95, r'film', verticalalignment='center', horizontalalignment='center', transform=ax.transAxes)
    ax.text(0.2, 0.95, r'substrate', verticalalignment='center', horizontalalignment='center', transform=ax.transAxes)
    ax.text(0.8, 0.95, r'cladding', verticalalignment='center', horizontalalignment='center', transform=ax.transAxes)
    
def initialize():
    global polarization_save,epsilon_f_save,epsilon_s_save,epsilon_c_save,d_string_save,mu_string_save
    polarization_string.set("TE")
    epsilon_f_string.set("2.25")
    epsilon_s_string.set("2.")
    epsilon_c_string.set("1.")
    d_string.set("3.")
    mu_string.set("0")
    
    polarization_save = polarization_string.get()
    epsilon_f_save = epsilon_f_string.get()
    epsilon_s_save = epsilon_s_string.get()
    epsilon_c_save = epsilon_c_string.get()
    d_string_save = d_string.get()
    mu_string_save = mu_string.get()
    
    calculate()
    
def reinitialize():
    global polarization_save,epsilon_f_save,epsilon_s_save,epsilon_c_save,d_string_save,mu_string_save
    polarization_string.set(polarization_save)
    epsilon_f_string.set(epsilon_f_save)
    epsilon_s_string.set(epsilon_s_save)
    epsilon_c_string.set(epsilon_c_save)
    d_string.set(d_string_save)
    mu_string.set(mu_string_save)

#    calculate()

def calculate():
    global polarization_save,epsilon_f_save,epsilon_s_save,epsilon_c_save,d_string_save,mu_string_save
    try:
        polarization = polarization_string.get()
        epsilon_f = float(epsilon_f_string.get())
        epsilon_s = float(epsilon_s_string.get())
        epsilon_c = float(epsilon_c_string.get())
        d = float(d_string.get())
        mu = int(float(mu_string.get()))
        mu_string.set(mu)
        
        if epsilon_s < epsilon_c:
            gui.input_error("Substrate epsilon smaller than cladding epsilon. Switching values...")
            epsilon_s = float(epsilon_c_string.get())
            epsilon_c = float(epsilon_s_string.get())
            epsilon_s_string.set(epsilon_s)
            epsilon_c_string.set(epsilon_c)
        if mu < 0:
            gui.input_error("Mode index negative. Re-initializing with previous parameters...",reinitialize)
        elif d < 0:
            gui.input_error("Film thickness negative. Re-initializing with previous parameters...",reinitialize)
        elif epsilon_f == 0 or epsilon_s == 0 or epsilon_c == 0:
            gui.input_error("Some dielectric constants are zero. Re-initializing with previous parameters...",reinitialize)
        elif (epsilon_f > 0 and epsilon_f < epsilon_s) or film.number_of_modes(epsilon_f,epsilon_s,epsilon_c,2*d,polarization) <= 0:
            gui.input_error("No mode found in the plotting window. Re-initializing with previous parameters...",reinitialize)
        elif film.number_of_modes(epsilon_f,epsilon_s,epsilon_c,2*d,polarization) > 10:
            gui.input_error("More than 10 modes found, that is too many for readable plotting of DR. Re-initializing with previous parameters...",reinitialize)
        else:
            if epsilon_f > 0:
                n_eff,step = np.linspace(np.sqrt(np.maximum(0,epsilon_s)), np.sqrt(epsilon_f), num=500, endpoint=False, retstep=True)
            elif -epsilon_f > epsilon_c: # connect to film-cladding SPP at d -> infinity
                n_eff = np.linspace(np.sqrt(epsilon_s), np.sqrt(np.maximum(epsilon_s,epsilon_f*epsilon_c/(epsilon_f+epsilon_c))), num=500, endpoint=True)
            else:
                n_eff = np.linspace(np.sqrt(epsilon_s),5*np.sqrt(epsilon_s), num=500, endpoint=True)
            vinvdisp = np.vectorize(film.invdisp)
            f.clf()
            a1 = f.add_subplot(221)
            n_eff_max = n_eff[-1]
            for index in range(film.number_of_modes(epsilon_f,epsilon_s,epsilon_c,2*d,polarization)):
                if ((epsilon_f < 0 and -epsilon_f > epsilon_s) or (polarization == 'TM' and -epsilon_s > epsilon_f)) and index == 0: # connect to film-substrate SPP at d -> infinity
                    n_eff_extra = np.linspace(np.sqrt(epsilon_f*epsilon_s/(epsilon_f+epsilon_s)), 1.5*np.sqrt(epsilon_f*epsilon_s/(epsilon_f+epsilon_s)), num=500, endpoint=True)
                    n_eff_max = np.maximum(n_eff_max, n_eff_extra[-1])
                    a1.plot(np.real(vinvdisp(epsilon_f,epsilon_s,epsilon_c,n_eff_extra,d,polarization)), n_eff_extra, label=r'$\mu=%2d$'%index)
                elif epsilon_f < 0: # connect to film-cladding SPP at d -> infinity
                    n_eff_max = np.maximum(n_eff_max, n_eff[-1] + 0.01*(n_eff[-1]-n_eff[0]))
                    a1.plot(np.real(vinvdisp(epsilon_f,epsilon_s,epsilon_c,n_eff,d,polarization)), n_eff, label=r'$\mu=%2d$'%index)
                elif polarization == 'TM' and ((-epsilon_s > epsilon_f and index == 1) or (epsilon_s < 0 and -epsilon_c > epsilon_f  and index == 0)): # connect to film-cladding SPP at d -> infinity
                    n_eff_joint = np.concatenate((n_eff, np.linspace(n_eff[-1]+step, np.sqrt(epsilon_f*epsilon_c/(epsilon_f+epsilon_c)), num=500, endpoint=True)))
                    n_eff_max = np.maximum(n_eff_max, n_eff_joint[-1] + 0.01*(n_eff_joint[-1]-n_eff_joint[0]))
                    a1.plot(np.real(vinvdisp(epsilon_f,epsilon_s,epsilon_c,n_eff_joint,d,polarization)), n_eff_joint, label=r'$\mu=%2d$'%index)
                elif polarization == 'TM' and epsilon_s < 0 and index == 0: # account for groundstate range above sqrt(epsilon_f)
                    n_eff_joint = np.concatenate((n_eff, np.linspace(n_eff[-1]+step, 9*n_eff[-1]-8*n_eff[0], num=500, endpoint=True)))
                    n_eff_max = np.maximum(n_eff_max, n_eff_joint[-1] + 0.01*(n_eff_joint[-1]-n_eff_joint[0]))
                    a1.plot(np.real(vinvdisp(epsilon_f,epsilon_s,epsilon_c,n_eff_joint,d,polarization)), n_eff_joint, label=r'$\mu=%2d$'%index)
                elif polarization == 'TM' and -epsilon_c > epsilon_s and index == 0: # connect to SPP for d -> 0
                    if -epsilon_c > epsilon_f: # connect to film-cladding SPP at d -> infinity
                        n_eff_joint = np.concatenate((n_eff, np.linspace(np.maximum(np.sqrt(epsilon_s*epsilon_c/(epsilon_s+epsilon_c)),n_eff[-1]+step), 
                                                                         np.maximum(np.sqrt(epsilon_f*epsilon_c/(epsilon_f+epsilon_c)),n_eff[-1]), num=500, endpoint=True))) 
                    else:
                        n_eff_joint = np.concatenate((n_eff, np.linspace(np.maximum(np.sqrt(epsilon_s*epsilon_c/(epsilon_s+epsilon_c)),n_eff[-1]+step), 
                                                                         np.maximum(1.5*np.sqrt(epsilon_s*epsilon_c/(epsilon_s+epsilon_c)),n_eff[-1]), num=500, endpoint=True)))
                    n_eff_max = np.maximum(n_eff_max, n_eff_joint[-1] + 0.01*(n_eff_joint[-1]-n_eff_joint[0]))
                    a1.plot(np.real(vinvdisp(epsilon_f,epsilon_s,epsilon_c,n_eff_joint,d,polarization)), n_eff_joint, label=r'$\mu=%2d$'%index)
                elif polarization == 'TM' and -epsilon_s > epsilon_f: # usual modes when film-substrate SPP exists
                    a1.plot(np.real((index-1)/(2*np.sqrt(epsilon_f-n_eff**2+0j))+vinvdisp(epsilon_f,epsilon_s,epsilon_c,n_eff,d,polarization)), n_eff, label=r'$\mu=%2d$'%index)
                else:
                    a1.plot(np.real(index/(2*np.sqrt(epsilon_f-n_eff**2+0j))+vinvdisp(epsilon_f,epsilon_s,epsilon_c,n_eff,d,polarization)), n_eff, label=r'$\mu=%2d$'%index)
            plt.axis('tight')
            xmin, xmax = a1.get_xlim()
            if 2*d <= xmin or 0 >= xmax:
                gui.input_error("No mode found in the plotting window. Re-initializing with previous parameters...",reinitialize)            
            else:
                plt.xlim([0,2*d])
                plt.ylim([n_eff[0],np.max([n_eff_max,np.real(np.sqrt(epsilon_f+0j))])])
                plt.xlabel(r'$d/\lambda$')
                plt.ylabel(r'$n_{\rm eff}$')
                plt.legend()
                if epsilon_f > 0 and n_eff_max > 10*n_eff[-1]-9*n_eff[0]:
                    gui.input_warning("DR plot may not be readable for chosen parameters. Good luck...")
                else:
                    a1.twinx()
                    sqrtepsilon = [np.sqrt(np.maximum(0,epsilon_s)),np.real(np.sqrt(epsilon_f+0j))]
                    if epsilon_s<0:
                        labels = [r'$0$', r'$\sqrt{\varepsilon_{\rm f}}$']
                    else:
                        labels = [r'$\sqrt{\varepsilon_{\rm s}}$', r'$\sqrt{\varepsilon_{\rm f}}$']                    
                    plt.yticks(sqrtepsilon, labels)
                    plt.ylim(a1.get_ylim())                
                nm = film.number_of_modes(epsilon_f,epsilon_s,epsilon_c,d,polarization)
                if nm > 0 and nm < mu+1:
                    gui.input_error("Mode index too high. Reducing index...")
                    mu = nm-1
                    mu_string.set(mu)
                if ((epsilon_f < 0 and -epsilon_f > epsilon_s) or (polarization == 'TM' and -epsilon_s > epsilon_f)) and mu == 0: 
                    n_eff_selected = film.n_eff_d(epsilon_f,epsilon_s,epsilon_c,n_eff_extra,d,0,polarization)
                elif epsilon_f < 0:
                    n_eff_selected = film.n_eff_d(epsilon_f,epsilon_s,epsilon_c,n_eff,d,0,polarization)
                elif polarization == 'TM' and ((-epsilon_s > epsilon_f and mu == 1) or ((epsilon_s < 0 or -epsilon_c > epsilon_s)  and mu == 0)): 
                    n_eff_selected = film.n_eff_d(epsilon_f,epsilon_s,epsilon_c,n_eff_joint,d,0,polarization)
                elif polarization == 'TM' and -epsilon_s > epsilon_f:
                    n_eff_selected = film.n_eff_d(epsilon_f,epsilon_s,epsilon_c,n_eff,d,mu-1,polarization)
                else:
                    n_eff_selected = film.n_eff_d(epsilon_f,epsilon_s,epsilon_c,n_eff,d,mu,polarization)
                cs = ['b','r','g']
                if n_eff_selected.size > 0:
                    x = np.linspace(-2., 2, num=401, endpoint=True) # x in units of d
                    for index in range(n_eff_selected.size):
                        if index == 0:
                            neff_string.set(np.real(n_eff_selected[0]))
                            a1.plot(d,n_eff_selected[0],cs[0]+'o')
                            F,Gx,Gz = film.mode_profile(epsilon_f,epsilon_s,epsilon_c,n_eff_selected[0],d,x,polarization)
                            a2 = f.add_subplot(222)
                            a3 = f.add_subplot(223)
                            a4 = f.add_subplot(224)               
                        if index > 0:
                            neff_string.set(np.real(n_eff_selected[index]))
                            a1.plot(d,n_eff_selected[index],cs[index]+'o')
                            Fbis,Gxbis,Gzbis = film.mode_profile(epsilon_f,epsilon_s,epsilon_c,n_eff_selected[index],d,x,polarization)
                            if np.isnan(Fbis[0]) == False:
                                a2.plot(x*d,Fbis,cs[index])
                                a3.plot(x*d,Gxbis,cs[index])
                                a4.plot(x*d,Gzbis,cs[index])
                    if np.isnan(F[0]) == False:
                        plot_mode(a2,x,d,F)
                        plot_mode(a3,x,d,Gx)
                        plot_mode(a4,x,d,Gz)
                        if polarization == 'TE':
                            a2.set_ylabel(r'$E_y/E_y(x=-d/2)$')
                            a3.set_ylabel(r'$H_x/H_x(x=-d/2)$')
                            a4.set_ylabel(r'$H_z/H_z(x=-d/2)$')
                        elif polarization == 'TM':
                            a2.set_ylabel(r'$H_y/H_y(x=-d/2)$')
                            a3.set_ylabel(r'$E_x/E_x(x=-d/2)$')
                            a4.set_ylabel(r'$E_z/E_z(x=-d/2)$')      
                else:
                    gui.input_warning("Mode does not exist or is out of plotting range for given parameters. Plotting DR only...")
                plt.tight_layout()
                
#                plt.savefig('film_waveguide.pdf',bbox_inches='tight',dpi=300, transparent=True)

                polarization_save = polarization_string.get()
                epsilon_f_save = epsilon_f_string.get()
                epsilon_s_save = epsilon_s_string.get()
                epsilon_c_save = epsilon_c_string.get()
                d_string_save = d_string.get()
                mu_string_save = mu_string.get()

                canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing with previous parameters...", reinitialize)

f = plt.figure(1,[8,6])
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

polarization_string = Tk.StringVar()
epsilon_f_string = Tk.StringVar()
epsilon_s_string = Tk.StringVar()
epsilon_c_string = Tk.StringVar()
d_string = Tk.StringVar()
mu_string = Tk.StringVar()
neff_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_radiobutton(mainframe,['polarization:','TE','TM'],polarization_string,2,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry_with_latex(mainframe,r'film $\varepsilon_{\rm f} =$',epsilon_f_string,row)
row = gui.create_entry_with_latex(mainframe,r'substrate $\varepsilon_{\rm s} =$',epsilon_s_string,row)
row = gui.create_entry_with_latex(mainframe,r'cladding $\varepsilon_{\rm c} =$',epsilon_c_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry_with_latex(mainframe,r'film thickness $d/\lambda=$',d_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry_with_latex(mainframe,r'mode index $\mu=$',mu_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_label_with_latex(mainframe,r'eff.\ index $n_{\rm eff}=$',neff_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)