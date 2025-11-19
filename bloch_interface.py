import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import strat_stuff as strat
import gui_stuff as gui

gui.set_rcParams()
title = "Reflection and Transmission at Infinite Bloch Media"
root = Tk.Tk()
root.title(title)

def reflection_transmission(epsilon_s,epsilon_f1,epsilon_f2,d1,d2,phi): # computing coefficients of reflection and transmission
    kx,ksz,kczTE,kappaTE,dummy = strat.KSC_Bloch(epsilon_s,d1,epsilon_f1,d2,epsilon_f2,'TE',phi)
    kx,ksz,kczTM,kappaTM,dummy = strat.KSC_Bloch(epsilon_s,d1,epsilon_f1,d2,epsilon_f2,'TM',phi)       
    RTE = (ksz-kappaTE)/(ksz+kappaTE)
    RTM = (kappaTM-ksz/epsilon_s)/(ksz/epsilon_s+kappaTM) # for electric field (negative of magnetic coeff.)
    return RTE,RTM,1-np.abs(RTE)**2,1-np.abs(RTM)**2

def initialize():
    var_string[0].set("1") # epsilon_s
    var_string[1].set("0.16") # d1 in units of lambda
    var_string[2].set("2.122") # Re epsilon_f1
    var_string[3].set("0") # Im epsilon_f1
    var_string[4].set("0.08") # d2 in units of lambda
    var_string[5].set("6.675") # Re epsilon_f2
    var_string[6].set("0") # Im epsilon_f2
    gui.copy_stringvar_vector(var_string,var_save)
    calculate()
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    calculate()  

def show_manual():
    gui.show_manual("taylor_series.png",title)

def calculate():
    gui.change_cursor(root,"trek")
    try:
        epsilon_s = float(var_string[0].get())
        d1 = float(var_string[1].get())
        epsilon_f1_real = float(var_string[2].get())
        epsilon_f1_imag = float(var_string[3].get())
        d2 = float(var_string[4].get())
        epsilon_f2_real = float(var_string[5].get())
        epsilon_f2_imag = float(var_string[6].get())  
        
        if epsilon_s <= 0 or d1 <= 0 or d2 <= 0 or (epsilon_f1_real == 0 and epsilon_f1_imag == 0) or (epsilon_f2_real == 0 and epsilon_f2_imag == 0): 
            gui.input_error("Values out of range. Re-initializing ...", reinitialize)
        else:
            f.clf()
            phi = np.linspace(0, np.pi/2, num=401, endpoint=False) # angle of incidence
            epsilon_f1 = epsilon_f1_real + 1j*epsilon_f1_imag
            epsilon_f2 = epsilon_f2_real + 1j*epsilon_f2_imag
            vreflection_transmission = np.vectorize(reflection_transmission)
            RTE,RTM,tauTE,tauTM = vreflection_transmission(epsilon_s,epsilon_f1,epsilon_f2,d1,d2,phi)
            a1 = f.add_subplot(131)
            strat.plot_curves_vs_angle(a1,phi,[np.abs(RTE)**2,tauTE],[r'$\rho_{\rm TE}$',r'$\tau_{\rm TE}$'],['b','r'],phi[0],phi[-1])
            left, right = a1.get_ylim()
            a1.set_ylim([min(left,-0.025),max(right,1.025)])
            a2 = f.add_subplot(132)
            strat.plot_curves_vs_angle(a2,phi,[np.abs(RTM)**2,tauTM],[r'$\rho_{\rm TM}$',r'$\tau_{\rm TM}$'],['b','r'],phi[0],phi[-1])
            left, right = a2.get_ylim()
            a2.set_ylim([min(left,-0.025),max(right,1.025)])
            a3 = f.add_subplot(133)
            strat.plot_curves_vs_angle(a3,phi,[np.angle(RTE),np.angle(RTM)],[r'$\theta_{\rm TE}$',r'$\theta_{\rm TM}$'],['b','r'],phi[0],phi[-1])  
            a3.set_yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi])
            a3.set_yticklabels([r'$-\pi$', r'$-\frac{\pi}{2}$', r'$0$', r'$\frac{\pi}{2}$', r'$\pi$'])
            a3.set_ylim([-1.05*np.pi,1.05*np.pi])
            
            plt.tight_layout()
            
#            plt.savefig('Bloch_interface.pdf',bbox_inches='tight',dpi=300, transparent=True)

            gui.copy_stringvar_vector(var_string,var_save)

            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[7,2])
canvas = gui.create_canvas(root,f)
canvas.draw()
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(7)
var_save = gui.create_stringvar_vector(7)

initialize()

row = 1
row = gui.create_entry(mainframe,u"substrate: \u03B5 =",var_string[0],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry(mainframe,u"film 1 thickness: d/\u03BB =",var_string[1],row)
row = gui.create_double_entry(mainframe,u"film 1: \u03B5' =",var_string[2],u"\u03B5'' =",var_string[3],row)
row = gui.create_entry(mainframe,u"film 2 thickness: d/\u03BB =",var_string[4],row)
row = gui.create_double_entry(mainframe,u"film 2: \u03B5' =",var_string[5],u"\u03B5'' =",var_string[6],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)