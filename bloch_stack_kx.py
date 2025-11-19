import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import strat_stuff as strat
import gui_stuff as gui

gui.set_rcParams()
title = "Reflection and Transmission at Stacks of Bloch Media - Fixed Frequency"
root = Tk.Tk()
root.title(title)

def reflection_transmission(epsilon_s,da,epsilon_fa,db,epsilon_fb,N,epsilon_c,phi): # computing coefficients of reflection and transmission
    kx,ksz,kcz = strat.KSC(epsilon_s,epsilon_c,phi)
    MTE = np.array([[1,0],[0,1]])
    MTM = np.array([[1,0],[0,1]])
    for index in range(len(N)):
        MTEP,MTMP = strat.MP_Bloch(kx,da[index],epsilon_fa[index],db[index],epsilon_fb[index],N[index])
        MTE = np.matmul(MTEP,MTE)
        MTM = np.matmul(MTMP,MTM)
    RTE,RTM,TTE,TTM,tauTE,tauTM = strat.RTAU(ksz,kcz,epsilon_s,epsilon_c,MTE,MTM) 
    return RTE,RTM,tauTE,tauTM
       
def initialize():
    var_string[0].set("1") # epsilon_s
    var_string[1].set("20") # Number of periods medium 1
    var_string[2].set("0") # Number of periods medium 2
    var_string[3].set("0") # Number of periods medium 3
    var_string[4].set("0.16") # d1 of medium 1 in units of lambda
    var_string[5].set("1") # d1 of medium 2 in units of lambda
    var_string[6].set("1") # d1 of medium 3 in units of lambda
    var_string[7].set("2.122") # Re epsilon_f1 of medium 1
    var_string[8].set("1") # Re epsilon_f1 of medium 2
    var_string[9].set("1") # Re epsilon_f1 of medium 3
    var_string[10].set("0") # Im epsilon_f1 of medium 1
    var_string[11].set("0") # Im epsilon_f1 of medium 2
    var_string[12].set("0") # Im epsilon_f1 of medium 3
    var_string[13].set("0.08") # d2 of medium 1 in units of lambda
    var_string[14].set("1") # d2 of medium 2 in units of lambda
    var_string[15].set("1") # d2 of medium 3 in units of lambda
    var_string[16].set("6.675") # Re epsilon_f2 of medium 1
    var_string[17].set("1") # Re epsilon_f2 of medium 2
    var_string[18].set("1") # Re epsilon_f2 of medium 3
    var_string[19].set("0") # Im epsilon_f2 of medium 1
    var_string[20].set("0") # Im epsilon_f2 of medium 2
    var_string[21].set("0") # Im epsilon_f2 of medium 3    
    var_string[22].set("2.122")
    var_string[23].set("0")
    var_string[24].set("0")
    var_string[25].set(".5")
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
        da = np.zeros(3)
        db = np.zeros(3)
        epsilon_fa_real = np.zeros(3)
        epsilon_fa_imag = np.zeros(3)
        epsilon_fb_real = np.zeros(3)
        epsilon_fb_imag = np.zeros(3)
        N = np.zeros(3,dtype=int)
        for index in range(3):
            N[index] = int(var_string[1+index].get())
            da[index] = float(var_string[4+index].get())
            epsilon_fa_real[index] = float(var_string[7+index].get())
            epsilon_fa_imag[index] = float(var_string[10+index].get())
            db[index] = float(var_string[13+index].get())
            epsilon_fb_real[index] = float(var_string[16+index].get())       
            epsilon_fb_imag[index] = float(var_string[19+index].get())                    
        epsilon_c_real = float(var_string[22].get())
        epsilon_c_imag = float(var_string[23].get())
        phi_min = float(var_string[24].get())*np.pi
        phi_max = float(var_string[25].get())*np.pi
        
        if (epsilon_c_real == 0 and epsilon_c_imag == 0) or (epsilon_fa_real +1j* epsilon_fa_imag == 0+0j).any() or (epsilon_fb_real +1j* epsilon_fb_imag == 0+0j).any() or epsilon_s <= 0\
            or (da <= 0).any() or (db < 0).any() or phi_max > np.pi/2 or phi_min < 0 or phi_min >= phi_max:
            gui.input_error("Values out of range. Re-initializing ...", reinitialize)
        elif (N < 0).any() or (N > 50).any():
            gui.input_error("Number of periods must be between 0 and 50. Re-initializing ...", reinitialize)
        else:
            f.clf()
            phi = np.linspace(phi_min, phi_max, num=10001, endpoint=False) # angle of incidence
            epsilon_c = epsilon_c_real+1j*epsilon_c_imag
            epsilon_fa = epsilon_fa_real+1j*epsilon_fa_imag
            epsilon_fb = epsilon_fb_real+1j*epsilon_fb_imag
            vreflection_transmission = np.vectorize(reflection_transmission, excluded=[1,2,3,4,5])
            RTE,RTM,tauTE,tauTM = vreflection_transmission(epsilon_s,da,epsilon_fa,db,epsilon_fb,N,epsilon_c+0j,phi)
            a1 = f.add_subplot(131)
            a2 = f.add_subplot(132)
            if (epsilon_fa_imag != 0).any() or (epsilon_fb_imag != 0).any():
                strat.plot_curves_vs_angle(a1,phi,[np.abs(RTE)**2,tauTE,1-np.abs(RTE)**2-tauTE],[r'$\rho_{\rm TE}$',r'$\tau_{\rm TE}$',r'$a_{\rm TE}$'],['b','r','g'],phi_min, phi_max)
                strat.plot_curves_vs_angle(a2,phi,[np.abs(RTM)**2,tauTM,1-np.abs(RTM)**2-tauTM],[r'$\rho_{\rm TM}$',r'$\tau_{\rm TM}$',r'$a_{\rm TM}$'],['b','r','g'],phi_min, phi_max)                
            else:
                strat.plot_curves_vs_angle(a1,phi,[np.abs(RTE)**2,tauTE],[r'$\rho_{\rm TE}$',r'$\tau_{\rm TE}$'],['b','r'],phi_min, phi_max)
                strat.plot_curves_vs_angle(a2,phi,[np.abs(RTM)**2,tauTM],[r'$\rho_{\rm TM}$',r'$\tau_{\rm TM}$'],['b','r'],phi_min, phi_max)
            a1.set_ylim([-0.025,1.025])
            a2.set_ylim([-0.025,1.025])
            a3 = f.add_subplot(133)
            strat.plot_curves_vs_angle(a3,phi,[np.angle(RTE),np.angle(RTM)],[r'$\theta_{\rm TE}$',r'$\theta_{\rm TM}$'],['b','r'],phi_min, phi_max)  
            a3.set_yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi])
            a3.set_yticklabels([r'$-\pi$', r'$-\frac{\pi}{2}$', r'$0$', r'$\frac{\pi}{2}$', r'$\pi$'])
            a3.set_ylim([-1.05*np.pi,1.05*np.pi])
                     
            plt.tight_layout()  
            
#            plt.savefig('Bloch_stack.pdf',bbox_inches='tight',dpi=300, transparent=True)
            
            gui.copy_stringvar_vector(var_string,var_save)

            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[7,2])
canvas = gui.create_canvas(root,f)
canvas.draw()
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(26)
var_save = gui.create_stringvar_vector(26)

initialize()

row = 1
row = gui.create_entry(mainframe,u"substrate: \u03B5 =",var_string[0],row)
row = gui.create_entry(mainframe,u"PC 1: N =",var_string[1],row)
row = gui.create_triple_entry(mainframe,u"d\u2081/\u03BB =",var_string[4],u"\u03B5\u2081' =",var_string[7],u"\u03B5\u2081'' =",var_string[10],row)
row = gui.create_triple_entry(mainframe,u"d\u2082/\u03BB =",var_string[13],u"\u03B5\u2082' =",var_string[16],u"\u03B5\u2082'' =",var_string[19],row)
row = gui.create_entry(mainframe,u"PC 2: N =",var_string[2],row)
row = gui.create_triple_entry(mainframe,u"d\u2081/\u03BB =",var_string[5],u"\u03B5\u2081' =",var_string[8],u"\u03B5\u2081'' =",var_string[11],row)
row = gui.create_triple_entry(mainframe,u"d\u2082/\u03BB =",var_string[14],u"\u03B5\u2082' =",var_string[17],u"\u03B5\u2082'' =",var_string[20],row)
row = gui.create_entry(mainframe,u"PC 3: N =",var_string[3],row)
row = gui.create_triple_entry(mainframe,u"d\u2081/\u03BB =",var_string[6],u"\u03B5\u2081' =",var_string[9],u"\u03B5\u2081'' =",var_string[12],row)
row = gui.create_triple_entry(mainframe,u"d\u2082/\u03BB =",var_string[15],u"\u03B5\u2082' =",var_string[18],u"\u03B5\u2082'' =",var_string[21],row)
row = gui.create_entry(mainframe,u"cladding: \u03B5' =",var_string[22],row)
row = gui.create_entry(mainframe,u"cladding: \u03B5'' =",var_string[23],row)
row = gui.create_double_entry_with_latex(mainframe,r'$\varphi_\mathrm{i}/\pi>$',var_string[24],r'$\varphi_\mathrm{i}/\pi<$',var_string[25],row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)