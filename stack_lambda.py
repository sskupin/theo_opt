import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import strat_stuff as strat
import media as media

gui.set_rcParams()
title = "Reflection and Transmission at Stacks -- Spectral Dependence"
root = Tk.Tk()
root.title(title)

def reflection_transmission(epsilon_s,d1,epsilon_f1,d2,epsilon_f2,N,epsilon_c,lambdav): # computing coefficients of reflection and transmission
    ksz,kcz = strat.KSC_lambda(epsilon_s,epsilon_c,lambdav)    
    M = strat.MP_lambda(d1,epsilon_f1,d2,epsilon_f2,lambdav,N)
    R,T,tau = strat.RTAU_lambda(ksz,kcz,epsilon_s,epsilon_c,M)
    return R,tau

def plot_subplot_twinx(ax,lambdav,curves,labels,colors):
    axbis = ax.twinx()
    axbis.yaxis.get_major_formatter().set_powerlimits((0, 3))
    lns1 = ax.plot(lambdav,curves[0],colors[0],label=labels[0])
    ax.set_xlabel(r'$\lambda$ [nm]')
    ax.set_ylabel(labels[0])
    lns2 = axbis.plot(lambdav,curves[1],colors[1],label=labels[1])
    axbis.set_ylabel(labels[1])
    ax.legend(lns1+lns2,labels)
    return axbis
    
def plot_subplot(ax,lambdav,curves,labels,colors):
    for index in range(len(labels)):
        ax.plot(lambdav,curves[index],colors[index],label=labels[index])
    ax.set_xlabel(r'$\lambda$ [nm]')
    ax.set_ylabel(','.join(labels))
    ax.legend()
    
def initialize():
    var_string[0].set("Vacuum") # substrate medium
    var_string[1].set("543") # thickness layer 1 in nm
    var_string[2].set("fused silica") # layer 1 medium
    var_string[3].set("0") # thickness layer 2 in nm
    var_string[4].set("fused silica") # layer 2 medium
    var_string[5].set("1") # Number of periods
    var_string[6].set("Si") # cladding medium
    var_string[7].set("400") # mininmum vacuum wavelenght in nm
    var_string[8].set("800") # maximum vacuum wavelenght in nm
    var_string[9].set("false") # use smae axis scale for dielectric funcions
    gui.copy_stringvar_vector(var_string,var_save)
    calculate()
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    calculate()  
    
def calculate():
    gui.change_cursor(root,"trek")
    try:
        substrate = var_string[0].get()
        d1 = float(var_string[1].get())
        film1 = var_string[2].get()
        d2 = float(var_string[3].get())
        film2 = var_string[4].get()
        N = int(var_string[5].get())
        cladding = var_string[6].get()
        lambda_min = float(var_string[7].get())
        lambda_max = float(var_string[8].get())
        epsilonscale = var_string[9].get()
        
        if d1 < 0 or d2 < 0 or N < 0:
            gui.input_error("Values out of range. Re-initializing ...", reinitialize)
        elif lambda_min < 400 or lambda_max > 2400 or lambda_min >= lambda_max:
            gui.input_error("Wavelength range between 400 and 2400 nm. Re-initializing ...", reinitialize)
        else:
            f.clf()
            lambdav = np.linspace(lambda_min, lambda_max, num=10001, endpoint=False) # vacuum wavelength in nm
            epsilon_s = media.epsilon(substrate,lambdav)
            epsilon_c = media.epsilon(cladding,lambdav)
            epsilon_f1 = media.epsilon(film1,lambdav)
            epsilon_f2 = media.epsilon(film2,lambdav)
            epsilonp_min = min(np.amin(np.real(epsilon_f1)),np.amin(np.real(epsilon_f2)),np.amin(np.real(epsilon_c)))
            epsilonp_min = epsilonp_min - .1*np.abs(epsilonp_min)
            epsilonp_max = max(np.amax(np.real(epsilon_f1)),np.amax(np.real(epsilon_f2)),np.amax(np.real(epsilon_c)))
            epsilonp_max = epsilonp_max + .1*np.abs(epsilonp_max)
            epsilonpp_max = 1.1*max(np.amax(np.imag(epsilon_f1)),np.amax(np.imag(epsilon_f2)),np.amax(np.imag(epsilon_c)))
            vreflection_transmission = np.vectorize(reflection_transmission)
            R,tau = vreflection_transmission(epsilon_s,d1,epsilon_f1,d2,epsilon_f2,N,epsilon_c,lambdav)
            a1 = f.add_subplot(221)
            plot_subplot(a1,lambdav,[np.abs(R)**2,tau,1-np.abs(R)**2-tau],[r'$\rho$',r'$\tau$',r'$a$'],['b','r','g'])
            a1.set_xlim([lambda_min, lambda_max])
            a1.set_ylim([-0.025,1.025])
            a2 = f.add_subplot(222)
            a2bis = plot_subplot_twinx(a2,lambdav,[np.real(epsilon_c),np.imag(epsilon_c)],[r'$\varepsilon_{\rm c}^{\prime}$',r'$\varepsilon_{\rm c}^{\prime\prime}$'],['b','r'])
            a2.set_xlim([lambda_min, lambda_max])
            if epsilonscale == 'true':
                a2.set_ylim([epsilonp_min,epsilonp_max])
                a2bis.set_ylim([0,max(1.e-3,epsilonpp_max)])
            else:
                a2.set_ylim([np.amin(np.real(epsilon_c)) - .1*np.abs(np.amin(np.real(epsilon_c))),np.amax(np.real(epsilon_c)) + .1*np.abs(np.amax(np.real(epsilon_c)))])
                a2bis.set_ylim([0,max(1.e-3,1.1*np.amax(np.imag(epsilon_c)))])
            a2.text(0.5, 0.925, r'cladding', verticalalignment='center', horizontalalignment='center', transform=a2.transAxes)
            a3 = f.add_subplot(223)
            a3bis = plot_subplot_twinx(a3,lambdav,[np.real(epsilon_f1),np.imag(epsilon_f1)],[r'$\varepsilon_{\rm f1}^{\prime}$',r'$\varepsilon_{\rm f1}^{\prime\prime}$'],['b','r'])
            a3.set_xlim([lambda_min, lambda_max])
            if epsilonscale == 'true':
                a3.set_ylim([epsilonp_min,epsilonp_max])
                a3bis.set_ylim([0,max(1.e-3,epsilonpp_max)])
            else:
                a3.set_ylim([np.amin(np.real(epsilon_f1)) - .1*np.abs(np.amin(np.real(epsilon_f1))),np.amax(np.real(epsilon_f1)) + .1*np.abs(np.amax(np.real(epsilon_f1)))])
                a3bis.set_ylim([0,max(1.e-3,1.1*np.amax(np.imag(epsilon_f1)))])
            a3.text(0.5, 0.925, r'film 1', verticalalignment='center', horizontalalignment='center', transform=a3.transAxes)
            a4 = f.add_subplot(224)
            a4bis = plot_subplot_twinx(a4,lambdav,[np.real(epsilon_f2),np.imag(epsilon_f2)],[r'$\varepsilon_{\rm f2}^{\prime}$',r'$\varepsilon_{\rm f2}^{\prime\prime}$'],['b','r'])
            a4.set_xlim([lambda_min, lambda_max])
            if epsilonscale == 'true':
                a4.set_ylim([epsilonp_min,epsilonp_max])
                a4bis.set_ylim([0,max(1.e-3,epsilonpp_max)])
            else:
                a4.set_ylim([np.amin(np.real(epsilon_f2)) - .1*np.abs(np.amin(np.real(epsilon_f2))),np.amax(np.real(epsilon_f2)) + .1*np.abs(np.amax(np.real(epsilon_f2)))])            
                a4bis.set_ylim([0,max(1.e-3,1.1*np.amax(np.imag(epsilon_f2)))])
            a4.text(0.5, 0.925, r'film 2', verticalalignment='center', horizontalalignment='center', transform=a4.transAxes)

            plt.tight_layout()  
            
 #           plt.savefig('stack_lambda.pdf',bbox_inches='tight',dpi=300, transparent=True)

            gui.copy_stringvar_vector(var_string,var_save)

            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[8,6])
canvas = gui.create_canvas(root,f)
canvas.draw()
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(10)
var_save = gui.create_stringvar_vector(10)

initialize()

row = 1
row = gui.create_radiobutton(mainframe,['substrate medium:','Vacuum','fused silica'],var_string[0],2,row)
row = gui.create_entry(mainframe,u"film 1 thickness: d [nm] =",var_string[1],row)
row = gui.create_radiobutton(mainframe,['film 1 medium:','Vacuum','Si','AlGaAs (70% Al)','AlAs','AlGaAs (31.5% Al)','GaAs','TiO2','Ag','fused silica','BaSF'],var_string[2],10,row)
row = gui.create_entry(mainframe,u"film 2 thickness: d [nm] =",var_string[3],row)
row = gui.create_radiobutton(mainframe,['film 2 medium:','Vacuum','Si','AlGaAs (70% Al)','AlAs','AlGaAs (31.5% Al)','GaAs','TiO2','Ag','fused silica','BaSF'],var_string[4],10,row)
row = gui.create_entry(mainframe,u"number of periods =",var_string[5],row)
row = gui.create_radiobutton(mainframe,['cladding medium:','Vacuum','Si','AlGaAs (70% Al)','AlAs','AlGaAs (31.5% Al)','GaAs','TiO2','Ag','fused silica','BaSF'],var_string[6],10,row)
row = gui.create_double_entry(mainframe,u"\u03bb [nm] >",var_string[7],u"\u03bb [nm] <",var_string[8],row)
row = gui.create_checkbutton(mainframe,'show dielectric functions on same scale','false','true',var_string[9],row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)