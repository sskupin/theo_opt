import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import sys
sys.path.append('./stuff')
import gui_stuff as gui
import strat_stuff as strat
import media as media

gui.set_rcParams()
title = "Transmission of Fabry-Perot Resonators - Normal Incidence"
root = Tk.Tk()
root.title(title)

def reflection_transmission(epsilon_s,da,epsilon_fa0,epsilon_fa1,db,epsilon_fb0,epsilon_fb1,N,D,epsilon_f,epsilon_c,lambdav): # computing coefficients of reflection and transmission
    ksz,kcz = strat.KSC_lambda(epsilon_s,epsilon_c,lambdav)    
    M0 = strat.MP_lambda(da[0],epsilon_fa0,db[0],epsilon_fb0,lambdav,N[0])
    Mf = strat.MP_lambda(D,epsilon_f,0,epsilon_f,lambdav,1)
    M = np.matmul(Mf,M0)
    M1 = strat.MP_lambda(da[1],epsilon_fa1,db[1],epsilon_fb1,lambdav,N[1])
    M = np.matmul(M1,M)
    R,T,tau = strat.RTAU_lambda(ksz,kcz,epsilon_s,epsilon_c,M)   
    return R,tau

def initialize():
    var_string[0].set("Vacuum") # substrate medium
    var_string[1].set("61") # mirror 1: thickness layer 1 in nm
    var_string[2].set("TiO2") # mirror 1: layer 1 medium
    var_string[3].set("109") # mirror 1: thickness layer 2 in nm
    var_string[4].set("fused silica") # mirror 1: layer 2 medium
    var_string[5].set("5") # mirror 1: Number of periods
    var_string[6].set("10000") # cavity thickness D in nm
    var_string[7].set("Vacuum") # cavity medium
    var_string[8].set("109") # mirror 2: thickness layer 1 in nm
    var_string[9].set("fused silica") # mirror 2: layer 1 medium
    var_string[10].set("61") # mirror 2: thickness layer 2 in nm
    var_string[11].set("TiO2") # mirror 2: layer 2 medium
    var_string[12].set("5") # mirror 2: Number of periods
    var_string[13].set("Vacuum") # cladding medium
    var_string[14].set("600") # mininmum vacuum wavelenght in nm
    var_string[15].set("700") # maximum vacuum wavelenght in nm
    var_string[16].set("no_show_all")
    gui.copy_stringvar_vector(var_string,var_save) 
    calculate()
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    calculate()  

def show_manual():
    gui.show_manual("man/fp_lambda.png",title)
    
def calculate():
    gui.change_cursor(root,"trek")
    try:
        lambda_min = float(var_string[14].get())
        lambda_max = float(var_string[15].get())
        cladding = var_string[13].get()
        substrate = var_string[0].get()
        da = np.zeros(2)
        db = np.zeros(2)
        filma = np.empty(2, dtype=object)
        filmb = np.empty(2, dtype=object)
        N = np.zeros(2,dtype=int)
        da[0] = float(var_string[1].get())
        da[1] = float(var_string[8].get())
        filma[0] = var_string[2].get()
        filma[1] = var_string[9].get()
        db[0] = float(var_string[3].get())
        db[1] = float(var_string[10].get())
        filmb[0] = var_string[4].get()
        filmb[1] = var_string[11].get()
        N[0] = int(var_string[5].get())
        N[1] = int(var_string[12].get())
        D = float(var_string[6].get())
        film = var_string[7].get()

        if (da <= 0).any() or (db < 0).any() or (D<0):
            gui.input_error("Values out of range. Re-initializing ...", reinitialize)
        elif lambda_min < 400 or lambda_max > 2400 or lambda_min >= lambda_max:
            gui.input_error("Wavelength range between 400 and 2400 nm. Re-initializing ...", reinitialize)
        elif (N < 0).any() or (N > 50).any():
            gui.input_error("Number of periods must be between 0 and 50. Re-initializing ...", reinitialize)
        else:
            f.clf()
            lambdav = np.linspace(lambda_min, lambda_max, num=10001, endpoint=False) # vacuum wavelength in nm
            epsilon_s = media.epsilon(substrate,lambdav)
            epsilon_c = media.epsilon(cladding,lambdav)
            epsilon_fa = np.zeros([len(N),len(lambdav)])+0j
            epsilon_fb = np.zeros([len(N),len(lambdav)])+0j
            for index in range(len(N)):
                epsilon_fa[index,:] = media.epsilon(filma[index],lambdav)
                epsilon_fb[index,:] = media.epsilon(filmb[index],lambdav)
            epsilon_f = media.epsilon(film,lambdav)
            vreflection_transmission = np.vectorize(reflection_transmission, excluded=[1,4,7])
            R,tau = vreflection_transmission(epsilon_s,da,epsilon_fa[0,:],epsilon_fa[1,:],db,epsilon_fb[0,:],epsilon_fb[1,:],N,D,epsilon_f,epsilon_c,lambdav)
            a1 = f.add_subplot(111)
            a1.plot(lambdav,tau,'r',label=r'$\tau$')
            a1.set_xlabel(r'$\lambda$ [nm]')
            a1.set_ylabel(r'$\tau$')
            a1.set_xlim([lambda_min, lambda_max])
            a1.set_ylim([-0.025,1.025])
            if var_string[16].get()=="show_all":
                a1.plot(lambdav,np.abs(R)**2,'b',label=r'$\rho$')
                a1.plot(lambdav,1-np.abs(R)**2-tau,'g',label=r'$a$')
                a1.set_ylabel(r'$\rho$, $\tau$, $a$')
                a1.legend()                
            plt.tight_layout()  
            
#            plt.savefig('FP_lambda.pdf',bbox_inches='tight',dpi=300, transparent=True)

            gui.copy_stringvar_vector(var_string,var_save)

            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")
    
f = plt.figure(1,[5,2.5])
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(17)
var_save = gui.create_stringvar_vector(17)

initialize()

row = 1
row = gui.create_radiobutton(mainframe,['substrate medium:','Vacuum','fused silica'],var_string[0],2,row)
row = gui.create_triple_entry(mainframe,u"mirror 1: N =",var_string[5],u"d\u2081 [nm] =",var_string[1],u"d\u2082 [nm] =",var_string[3],row)
row = gui.create_radiobutton(mainframe,['film 1 medium:','Ag','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica','BaSF'],var_string[2],6,row)
row = gui.create_radiobutton(mainframe,['film 2 medium:','Ag','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica','BaSF'],var_string[4],6,row)
row = gui.create_radiobutton_with_entry(mainframe,u"cavity: D [nm] =",var_string[6],['Vacuum','fused silica'],var_string[7],2,row)
row = gui.create_triple_entry(mainframe,u"mirror 2: N =",var_string[12],u"d\u2081 [nm] =",var_string[8],u"d\u2082 [nm] =",var_string[10],row)
row = gui.create_radiobutton(mainframe,['film 1 medium:','Ag','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica','BaSF'],var_string[9],6,row)
row = gui.create_radiobutton(mainframe,['film 2 medium:','Ag','AlAs','AlGaAs (31.5% Al)','TiO2','fused silica','BaSF'],var_string[11],6,row)
row = gui.create_radiobutton(mainframe,['cladding medium:','Vacuum','fused silica'],var_string[13],2,row)
row = gui.create_double_entry(mainframe,u"\u03bb [nm] >",var_string[14],u"\u03bb [nm] <",var_string[15],row)
row = gui.create_checkbutton_with_latex(mainframe,r'show transmittance only','show_all','no_show_all',var_string[16],row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)