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
    rho0_string.set("0.99")
    rhoD_string.set("0.99")
    absorption_string.set("lossless")
    delta_min_string.set(-1)
    delta_max_string.set(3)
        
    calculate()
    
def calculate():
    try:
        rho0 = float(rho0_string.get())
        rhoD = float(rhoD_string.get())
        absorption = absorption_string.get()
        delta_min = float(delta_min_string.get())*np.pi
        delta_max = float(delta_max_string.get())*np.pi
        if absorption == 'lossless':
            tau0 = 1-rho0
            tauD = 1-rhoD
            tau0_string.set('')
            tauD_string.set('')
        else:
            tau0 = float(tau0_string.get())
            tauD = float(tauD_string.get())    
            
        if rho0 <= 0 or rho0 >= 1 or rhoD <= 0 or rhoD >= 1 or tau0 <= 0 or rho0+tau0 > 1 or tauD <= 0 or rhoD+tauD > 1 or delta_min >= delta_max:
            gui.input_error(initialize)
        else:
            f.clf()
            delta = np.linspace(delta_min, delta_max , num=10001, endpoint=False) 
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

            canvas.draw()
    except ValueError: gui.input_error(initialize)

f = plt.figure(1,[4,2])
canvas = gui.create_canvas(root,f)
canvas.draw()
mainframe = gui.create_mainframe(root)

rho0_string = Tk.StringVar()
rhoD_string = Tk.StringVar()
tau0_string = Tk.StringVar()
tauD_string = Tk.StringVar()
absorption_string = Tk.StringVar()
delta_min_string = Tk.StringVar()
delta_max_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_double_entry_with_latex(mainframe,r'$\rho_0=$',rho0_string,r'$\rho_\mathrm{D}=$',rhoD_string,row)
row = gui.create_checkbutton(mainframe,"lossy mirrors",'lossless','lossy',absorption_string,row)
row = gui.create_double_entry_with_latex(mainframe,r'$\tau_0=$',tau0_string,r'$\tau_\mathrm{D}=$',tauD_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_entry_with_latex(mainframe,r'$\delta/\pi>$',delta_min_string,r'$\delta/\pi<$',delta_max_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)