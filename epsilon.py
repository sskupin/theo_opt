import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as spc
import tkinter as Tk
import gui_stuff as gui
import media as media

gui.set_rcParams()
title = "Dielectric Function for Selected Media"
root = Tk.Tk()
root.title(title)

def initialize():
    var_string[0].set("400") # lambda_min
    var_string[1].set("2400") # lambda_max
    var_string[2].set("GaAs") # medium
    gui.copy_stringvar_vector(var_string,var_save)
    calculate()
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    
def show_manual():
    gui.show_manual("taylor_series.png",title)
    
def calculate():
    gui.change_cursor(root,"trek")
    try:
        lambda_min = float(var_string[0].get())
        lambda_max = float(var_string[1].get())
        medium = var_string[2].get()
        
        if lambda_min < 400 or lambda_max > 2400 or lambda_min >= lambda_max:
            gui.input_error("Wavelength range between 400 and 2400 nm. Re-initializing ...", reinitialize)
        else:
            f.clf()
            a = f.add_subplot(111)
            omega_max = media.omega2lambda(lambda_min) # in 10^{15} s^{-1}
            omega_min = media.omega2lambda(lambda_max) # in 10^{15} s^{-1}
            omega = np.linspace(omega_min, omega_max, num=10001, endpoint=False) # angular frequency in 10^{15} s^{-1}
            lambdav = media.omega2lambda(omega) # vacuum wavelength in nm
            lns1 = a.plot(omega,np.real(media.epsilon(medium,lambdav)),'b',label=r'$\varepsilon^{\prime}$')
            a.set_xlim([omega_min, omega_max])
            a.set_xlabel(r'$\omega$ [$10^{15}$ s$^{-1}$]')
            a.set_ylabel(r'$\varepsilon^{\prime}$')
            if np.amax(np.imag(media.epsilon(medium,lambdav)))>0:
                abis = a.twinx()
                abis.yaxis.get_major_formatter().set_powerlimits((0, 3))
                lns2 = abis.semilogy(omega,np.imag(media.epsilon(medium,lambdav)),'r',label=r'$\varepsilon^{\prime\prime}$')
                ylim = abis.get_ylim()
                ylim = [np.amax([ylim[0],1e-6]),ylim[1]]
                abis.set_ylim(ylim)
                abis.set_ylabel(r'$\varepsilon^{\prime\prime}$')
                a.legend(lns1+lns2,[r'$\varepsilon^{\prime}$',r'$\varepsilon^{\prime\prime}$'])
            abis2 = a.secondary_xaxis("top",functions=(media.omega2lambda, media.omega2lambda))
            abis2.set_xticks(2.e-6*np.pi*spc.c / a.get_xticks())
            abis2.set_xlabel(r'$\lambda$ [nm]')
            plt.tight_layout()
            
#            plt.savefig('epsilon.pdf',bbox_inches='tight',dpi=300, transparent=True)
            
            gui.copy_stringvar_vector(var_string,var_save)

            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[7,3.5])
canvas = gui.create_canvas(root,f)
canvas.draw() # for faster feedback to user on startup
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(3)
var_save = gui.create_stringvar_vector(3)

initialize()

row = 1
row = gui.create_radiobutton_with_latex(mainframe,[r'medium:',r'Al',r'Si',r'Al$_{0.7}$Ga$_{0.3}$GAs',r'AlAs',r'Al$_{0.315}$Ga$_{0.685}$As',r'GaAs',r'TiO$_2$ ($\varepsilon_{\rm or}$)',r'Ag',r'fused silica',r'BaSF'],['medium:','Al','Si','AlGaAs (70% Al)','AlAs','AlGaAs (31.5% Al)','GaAs','TiO2','Ag','fused silica','BaSF'],var_string[2],10,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_entry_with_latex(mainframe,r"$\lambda$ [nm] $>$",var_string[0],u"$\lambda$ [nm] $<$",var_string[1],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)