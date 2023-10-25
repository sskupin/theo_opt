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

def epsilon(medium,lambdav):
    if medium == "Al":
        epsilon_medium = media.Al(lambdav)
    elif medium == "Si":
        epsilon_medium = media.Si(lambdav)
    elif medium == "AlAs":
        epsilon_medium = media.AlAs(lambdav)
    elif medium == "GaAs":
        epsilon_medium = media.GaAs(lambdav)
    elif medium == "AlGaAs (70% Al)":
        epsilon_medium = media.AlGaAs70(lambdav)
    elif medium == "AlGaAs (31.5% Al)":
        epsilon_medium = media.AlGaAs31(lambdav)
    elif medium == "TiO2":
        epsilon_medium = media.TiO2(lambdav)
    elif medium == "Ag":
        epsilon_medium = media.Ag(lambdav)
    elif medium == "fused silica":
        epsilon_medium = media.silica(lambdav)
    elif medium == "BaSF":
        epsilon_medium = media.BaSF(lambdav)
    else:
        print("Oops! Medium not known")
        
    return epsilon_medium
    
def omega2lambda(x):
    return 2.e-6*np.pi*spc.c*np.reciprocal(x , out=np.zeros_like(x), where=x!=0)

def initialize():
    global lambda_min_save,lambda_max_save
    lambda_min_string.set("400")
    lambda_max_string.set("2400")
    medium_string.set("GaAs")
    
    lambda_min_save = lambda_min_string.get()
    lambda_max_save = lambda_max_string.get()
    
    calculate()
    
def reinitialize():
    global lambda_min_save,lambda_max_save
    lambda_min_string.set(lambda_min_save)
    lambda_max_string.set(lambda_max_save)
    
def show_manual():
    top = Tk.Toplevel()
    img = gui.read_image("taylor_series.png")
    gui.show_image(top,title,img)
    gui.mainloop_safe_for_mac(top)
    
def calculate():
    global lambda_min_save,lambda_max_save
    try:
        lambda_min = float(lambda_min_string.get())
        lambda_max = float(lambda_max_string.get())
        medium = medium_string.get()
        
        if lambda_min < 400 or lambda_max > 2400 or lambda_min >= lambda_max:
            gui.input_error("Wavelength range between 400 and 2400 nm. Re-initializing ...", reinitialize)
        else:
            f.clf()
            a = f.add_subplot(111)
            omega_max = omega2lambda(lambda_min) # in 10^{15} s^{-1}
            omega_min = omega2lambda(lambda_max) # in 10^{15} s^{-1}
            omega = np.linspace(omega_min, omega_max, num=10001, endpoint=False) # angular frequency in 10^{15} s^{-1}
            lambdav = omega2lambda(omega) # vacuum wavelength in nm
            lns1 = a.plot(omega,np.real(epsilon(medium,lambdav)),'b',label=r'$\varepsilon^{\prime}$')
            a.set_xlim([omega_min, omega_max])
            a.set_xlabel(r'$\omega$ [$10^{15}$ s$^{-1}$]')
            a.set_ylabel(r'$\varepsilon^{\prime}$')
            if np.amax(np.imag(epsilon(medium,lambdav)))>0:
                abis = a.twinx()
                abis.yaxis.get_major_formatter().set_powerlimits((0, 3))
                lns2 = abis.semilogy(omega,np.imag(epsilon(medium,lambdav)),'r',label=r'$\varepsilon^{\prime\prime}$')
                ylim = abis.get_ylim()
                ylim = [np.amax([ylim[0],1e-6]),ylim[1]]
                abis.set_ylim(ylim)
                abis.set_ylabel(r'$\varepsilon^{\prime\prime}$')
                a.legend(lns1+lns2,[r'$\varepsilon^{\prime}$',r'$\varepsilon^{\prime\prime}$'])
            abis2 = a.secondary_xaxis("top",functions=(omega2lambda, omega2lambda))
            abis2.set_xticks(2.e-6*np.pi*spc.c / a.get_xticks())
            abis2.set_xlabel(r'$\lambda$ [nm]')
            plt.tight_layout()
            
#            plt.savefig('epsilon.pdf',bbox_inches='tight',dpi=300, transparent=True)
            
            lambda_min_save = lambda_min_string.get()
            lambda_max_save = lambda_max_string.get()

            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)

f = plt.figure(1,[7,3.5])
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

lambda_min_string = Tk.StringVar()
lambda_max_string = Tk.StringVar()
medium_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_radiobutton(mainframe,['medium:','Al','Si','AlGaAs (70% Al)','AlAs','AlGaAs (31.5% Al)','GaAs','TiO2','Ag','fused silica','BaSF'],medium_string,10,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_entry(mainframe,u"\u03bb [nm] >",lambda_min_string,u"\u03bb [nm] <",lambda_max_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)