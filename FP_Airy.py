import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
title = "Transmission of a Fabry-Perot Resonator (Airy Formula)"
root = Tk.Tk()
root.title(title)

def transmission(rho0,tau0,rhoD,tauD,delta): # computing transmission
    taumax = tau0*tauD/(1-np.sqrt(rho0*rhoD))**2
    F = 4*np.sqrt(rho0*rhoD)/(1-np.sqrt(rho0*rhoD))**2
    tau = taumax/(1+F*np.sin(delta/2)**2)
    
    return taumax,F,tau

def initialize():
    var_string[0].set("0.99")
    var_string[1].set("0.99")
    var_string[2].set("lossless")
    var_string[5].set(-1)
    var_string[6].set(3)
    gui.copy_stringvar_vector(var_string,var_save) 
    calculate()

def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    calculate()  


def calculate():
    gui.change_cursor(root,"trek")
    try:
        rho0 = float(var_string[0].get())
        rhoD = float(var_string[1].get())
        absorption = var_string[2].get()
        delta_min = float(var_string[5].get())*np.pi
        delta_max = float(var_string[6].get())*np.pi
        if absorption == 'lossless':
            tau0 = 1-rho0
            tauD = 1-rhoD
            var_string[3].set(str(round(tau0,len(var_string[0].get()))))
            var_string[4].set(str(round(tauD,len(var_string[1].get()))))
        else:
            tau0 = float(var_string[3].get())
            tauD = float(var_string[4].get())    
            
        if rho0 <= 0 or rho0 >= 1 or rhoD <= 0 or rhoD >= 1 or tau0 <= 0 or tauD <= 0 or delta_min >= delta_max:
            gui.input_error("Values out of range. Re-initializing ...", reinitialize)
        elif rho0+tau0 > 1 or rhoD+tauD > 1:
            gui.input_error("Sum of mirror reflectance and transmittance must not exceed one. Re-initializing ...", reinitialize)
        else:
            f.clf()
            delta,delta_delta = np.linspace(delta_min, delta_max , num=20001, endpoint=False, retstep=True) 
            delta_FWHM = 2*(1-np.sqrt(rho0*rhoD))/(rho0*rhoD)**0.25
            if delta_FWHM < 10*delta_delta:
                gui.input_error("Resonances too sharp for being displayed. Reduce range. Re-initializing ...", reinitialize)
            vtransmission = np.vectorize(transmission)
            taumax,F,tau = vtransmission(rho0,tau0,rhoD,tauD,delta)
            plt.plot(delta/np.pi,tau,'r')
            ax = f.gca()
            ax.set_xlabel(r'$\delta/\pi$')
            ax.autoscale(enable=True, axis='x', tight=True)
            ax.set_ylabel(r'$\tau$')
            ax.set_ylim([-0.025,1.025])
            ax.set_title(r'Finesse $\Phi=$ '+str(round(np.pi*np.sqrt(F[0])/2 ,5)))
            plt.tight_layout()  
            
#            plt.savefig('FP_airy.pdf',bbox_inches='tight',dpi=300, transparent=True)

            gui.copy_stringvar_vector(var_string,var_save)

            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[4,2])
canvas = gui.create_canvas(root,f)
canvas.draw()
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(7)
var_save = gui.create_stringvar_vector(7)

initialize()

row = 1
row = gui.create_double_entry_with_latex(mainframe,r'$\rho_0=$',var_string[0],r'$\rho_\mathrm{D}=$',var_string[1],row)
row = gui.create_checkbutton(mainframe,"lossless mirrors",'lossy','lossless',var_string[2],row)
row = gui.create_double_entry_with_latex(mainframe,r'$\tau_0=$',var_string[3],r'$\tau_\mathrm{D}=$',var_string[4],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_entry_with_latex(mainframe,r'$\delta/\pi>$',var_string[5],r'$\delta/\pi<$',var_string[6],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)