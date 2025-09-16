import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as spc
import tkinter as Tk
import strat_stuff as strat
import gui_stuff as gui
import media as media


gui.set_rcParams()
title = "Reflection and Transmission at Stacks of Bloch Media - Normal Incidence"
root = Tk.Tk()
root.title(title)

def reflection_transmission(epsilon_s,da,epsilon_fa,db,epsilon_fb,N,epsilon_c,lambdav): # computing coefficients of reflection and transmission
    ksz,kcz = strat.KSC_lambda(epsilon_s,epsilon_c,lambdav)
    M = np.array([[1,0],[0,1]])
    for index in range(len(N)):
        MP = strat.MP_Bloch_lambda(da[index],epsilon_fa,db[index],epsilon_fb,lambdav,N[index])
        M = np.matmul(MP,M)    
    R,T,tau = strat.RTAU_lambda(ksz,kcz,epsilon_s,epsilon_c,M)    
    return R,tau
        
def initialize():
    var_string[0].set("Vacuum") # substrate medium
    var_string[1].set("fused silica") # layer 1 medium
    var_string[2].set("TiO2") # layer 2 medium
    var_string[3].set("20") # Number of periods PC 1
    var_string[4].set("140") # d1 of medium 1 in PC 1 in nm
    var_string[5].set("130") # d2 of medium 2 in PC 1 in nm
    var_string[6].set("20") # Number of periods PC 2
    var_string[7].set("140") # d1 of medium 1 in PC 2 in nm
    var_string[8].set("180") # d2 of medium 2 in PC 2 in nm    
    var_string[9].set("20") # Number of periods PC 3
    var_string[10].set("140") # d1 of medium 1 in PC 3 in nm
    var_string[11].set("80") # d2 of medium 2 in PC 3 in nm
    var_string[12].set("fused silica") # cladding medium
    var_string[13].set("600")
    var_string[14].set("1600")
    
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
        lambda_min = float(var_string[13].get())
        lambda_max = float(var_string[14].get())
        substrate = var_string[0].get()
        da = np.zeros(3)
        db = np.zeros(3)
        filma = var_string[1].get()
        filmb = var_string[2].get()
        cladding = var_string[12].get()
        N = np.zeros(3,dtype=int)
        for index in range(3):
            da[index] = float(var_string[4+index*3].get())
            db[index] = float(var_string[5+index*3].get())
            N[index] = int(var_string[3+index*3].get())
        
        if (da <= 0).any() or (db < 0).any() or lambda_min < 400 or lambda_max > 2400 or lambda_min >= lambda_max:
            gui.input_error("Values out of range. Re-initializing ...", reinitialize)
        elif (N < 0).any() or (N > 50).any():
            gui.input_error("Number of periods must be between 0 and 50. Re-initializing ...", reinitialize)
        else:
            f.clf()
            lambdav = np.linspace(lambda_min, lambda_max, num=10001, endpoint=False) # vacuum wavelength in nm
            epsilon_s = media.epsilon(substrate,lambdav)
            epsilon_c = media.epsilon(cladding,lambdav)
            epsilon_fa = media.epsilon(filma,lambdav)
            epsilon_fb = media.epsilon(filmb,lambdav)

            vreflection_transmission = np.vectorize(reflection_transmission, excluded=[1,3,5])
            R,tau = vreflection_transmission(epsilon_s,da, epsilon_fa,db,epsilon_fb,N,epsilon_c,lambdav)
            a1 = f.add_subplot(211)
            if any(1-np.abs(R)**2-tau>10.e-8):
                strat.plot_subplot_lambda(a1,lambdav,[np.abs(R)**2,tau,1-np.abs(R)**2-tau],[r'$\rho$',r'$\tau$',r'$a$'],['b','r','g'])
            else:
                strat.plot_subplot_lambda(a1,lambdav,[np.abs(R)**2,tau,],[r'$\rho$',r'$\tau$'],['b','r'])
            a1.set_xlim([lambda_min, lambda_max])
            a1.set_ylim([-0.025,1.025])
            a2 = f.add_subplot(212)
            theta = np.angle(R)
            a2.plot(lambdav,theta,'b')
            a2.set_xlabel(r'$\lambda$ [nm]')
            a2.set_xlim([lambda_min, lambda_max])
            a2.set_ylabel(r'$\theta$')
                     
            plt.tight_layout()  
            
#            plt.savefig('Bloch_stack_lambda.pdf',bbox_inches='tight',dpi=300, transparent=True)

            gui.copy_stringvar_vector(var_string,var_save)

            canvas.draw()
    except ValueError: gui.input_error()
    gui.change_cursor(root,"arrow")
    
f = plt.figure(1,[4,4])
canvas = gui.create_canvas(root,f)
canvas.draw()
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(15)
var_save = gui.create_stringvar_vector(15)

initialize()

row = 1
row = gui.create_radiobutton(mainframe,['substrate medium:','Vacuum','fused silica'],var_string[0],2,row)
row = gui.create_radiobutton_with_latex(mainframe,[r'film 1 medium:',r'Vacuum',r'Al$_{0.7}$Ga$_{0.3}$GAs',r'AlAs',r'Al$_{0.315}$Ga$_{0.685}$As',r'TiO$_2$ ($\varepsilon_{\rm or}$)',r'fused silica'],['film 1 medium:','Vacuum','AlGaAs (70% Al)','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica'],var_string[1],6,row)
row = gui.create_radiobutton_with_latex(mainframe,[r'film 2 medium:',r'Vacuum',r'Al$_{0.7}$Ga$_{0.3}$GAs',r'AlAs',r'Al$_{0.315}$Ga$_{0.685}$As',r'TiO$_2$ ($\varepsilon_{\rm or}$)',r'fused silica'],['film 2 medium:','Vacuum','AlGaAs (70% Al)','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica'],var_string[2],6,row)
row = gui.create_triple_entry(mainframe,u"PC 1: N =",var_string[3],u"film 1: d [nm] =",var_string[4],u"film 2: d [nm] =",var_string[5],row)
row = gui.create_triple_entry(mainframe,u"PC 2: N =",var_string[6],u"film 1: d [nm] =",var_string[7],u"film 2: d [nm] =",var_string[8],row)
row = gui.create_triple_entry(mainframe,u"PC 3: N =",var_string[9],u"film 1: d [nm] =",var_string[10],u"film 2: d [nm] =",var_string[11],row)
row = gui.create_radiobutton_with_latex(mainframe,[r'cladding medium:',r'Vacuum',r'Al$_{0.7}$Ga$_{0.3}$GAs',r'AlAs',r'Al$_{0.315}$Ga$_{0.685}$As',r'TiO$_2$ ($\varepsilon_{\rm or}$)',r'fused silica'],['film 1 medium:','Vacuum','AlGaAs (70% Al)','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica'],var_string[12],6,row)
row = gui.create_double_entry(mainframe,u"\u03bb [nm] >",var_string[13],u"\u03bb [nm] <",var_string[14],row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)