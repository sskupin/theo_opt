import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import sys
sys.path.append('./aux')
import gui_stuff as gui
import strat_stuff as strat

gui.set_rcParams()
title = "Transmission of Fabry-Perot Resonators - Fixed Frequency"
root = Tk.Tk()
root.title(title)

def reflection_transmission(epsilon_s,da1,epsilon_fa1,db1,epsilon_fb1,N1,D,epsilon_f,da2,epsilon_fa2,db2,epsilon_fb2,N2,epsilon_c,phi): # computing coefficients of reflection and transmission
    kx,ksz,kcz = strat.KSC(epsilon_s,epsilon_c,phi)
    MTE1,MTM1 = strat.MP(kx,da1,epsilon_fa1,db1,epsilon_fb1,N1)
    MTEf,MTMf = strat.MP(kx,D,epsilon_f,0,epsilon_f,1)
    MTE = np.matmul(MTEf,MTE1)
    MTM = np.matmul(MTMf,MTM1)
    MTE2,MTM2 = strat.MP(kx,da2,epsilon_fa2,db2,epsilon_fb2,N2)
    MTE = np.matmul(MTE2,MTE)
    MTM = np.matmul(MTM2,MTM)
    RTE,RTM,TTE,TTM,tauTE,tauTM = strat.RTAU(ksz,kcz,epsilon_s,epsilon_c,MTE,MTM)     
    return RTE,RTM,tauTE,tauTM

def initialize():
    var_string[0].set("1") # epsilon_s
    var_string[1].set("0.147") # da1 in units of lambda
    var_string[2].set("2.13") # Re epsilon_fa1
    var_string[3].set("0") # Im epsilon_fa1
    var_string[4].set("0.13") # db1 in units of lambda
    var_string[5].set("6.68") # Re epsilon_fb1
    var_string[6].set("0") # Im epsilon_fb1 
    var_string[7].set("7") # N1
    var_string[8].set("5") # D in units of lambda
    var_string[9].set("2.5") # Re epsilon_f
    var_string[10].set("0") # Im epsilon_f
    var_string[11].set("0.13") # da2 in units of lambda
    var_string[12].set("6.68") # Re epsilon_fa2
    var_string[13].set("0") # Im epsilon_fa2
    var_string[14].set("0.147") # db2 in units of lambda
    var_string[15].set("2.13") # Re epsilon_fb2
    var_string[16].set("0") # Im epsilon_fb2 
    var_string[17].set("7") # N2
    var_string[18].set("1") # epsilon_c
    var_string[19].set("0")
    var_string[20].set("0.7")
    var_string[21].set("no_show_all")
    gui.copy_stringvar_vector(var_string,var_save)
    calculate()
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    calculate()  
    
def show_manual():
    gui.show_manual("man/fp_kx.png",title)

def calculate():
    gui.change_cursor(root,"trek")
    try:
        epsilon_s = float(var_string[0].get())
        da1 = float(var_string[1].get())
        epsilon_fa1_real = float(var_string[2].get())
        epsilon_fa1_imag = float(var_string[3].get())
        db1 = float(var_string[4].get())
        epsilon_fb1_real = float(var_string[5].get())
        epsilon_fb1_imag = float(var_string[6].get())
        N1 = int(var_string[7].get())
        D = float(var_string[8].get())
        epsilon_f_real = float(var_string[9].get())
        epsilon_f_imag = float(var_string[10].get())        
        da2 = float(var_string[11].get())
        epsilon_fa2_real = float(var_string[12].get())
        epsilon_fa2_imag = float(var_string[13].get())
        db2 = float(var_string[14].get())
        epsilon_fb2_real = float(var_string[15].get())
        epsilon_fb2_imag = float(var_string[16].get())        
        N2 = int(var_string[17].get())
        epsilon_c = float(var_string[18].get())   
        sinphi_min = float(var_string[19].get())
        sinphi_max = float(var_string[20].get())

        if epsilon_c <= 0 or (epsilon_fa1_real == 0 and epsilon_fa1_imag == 0) or (epsilon_fa2_real == 0 and epsilon_fa2_imag == 0)\
           or (epsilon_fb1_real == 0 and epsilon_fb1_imag == 0) or (epsilon_fb2_real == 0 and epsilon_fb2_imag == 0) or epsilon_s <= 0 or sinphi_min<0 or sinphi_max>1\
           or sinphi_min>=sinphi_max or da1 < 0 or da2 < 0 or db1 < 0 or db2 < 0 or D < 0 or (epsilon_f_real == 0 and epsilon_f_imag == 0):
           gui.input_error("Values out of range. Re-initializing ...", reinitialize)
        elif N1 < 0 or N1 > 50 or N2 < 0 or N2 > 50:
           gui.input_error("Number of periods must be between 0 and 50. Re-initializing ...", reinitialize)
        else:
            f.clf()
            phi = np.linspace(np.arcsin(sinphi_min), np.arcsin(sinphi_max), num=10001, endpoint=False) # angle of incidence
            epsilon_fa1 = epsilon_fa1_real + 1j*epsilon_fa1_imag
            epsilon_fb1 = epsilon_fb1_real + 1j*epsilon_fb1_imag
            epsilon_f = epsilon_f_real + 1j*epsilon_f_imag
            epsilon_fa2 = epsilon_fa2_real + 1j*epsilon_fa2_imag
            epsilon_fb2 = epsilon_fb2_real + 1j*epsilon_fb2_imag
            vreflection_transmission = np.vectorize(reflection_transmission)
            RTE,RTM,tauTE,tauTM = vreflection_transmission(epsilon_s,da1,epsilon_fa1,db1,epsilon_fb1,N1,D,epsilon_f,da2,epsilon_fa2,db2,epsilon_fb2,N2,epsilon_c,phi)
            a1 = f.add_subplot(211)
            a1.plot(np.sqrt(epsilon_s)*np.sin(phi),tauTE,'r',label=r'$\tau_{\rm TE}$')
            a1.set_xlabel(r'$\sqrt{\varepsilon_{\rm s}}\sin\varphi_{\rm i}=k_x \lambda/(2\pi)$')
            a1.set_xlim([np.sqrt(epsilon_s)*sinphi_min,np.sqrt(epsilon_s)*sinphi_max])
            a1.set_ylabel(r'$\tau_{\rm TE}$')
            a1.set_ylim([-0.025,0.025 + np.maximum(1,np.amax(tauTE))])
            if var_string[21].get()=="show_all":
                a1.plot(np.sqrt(epsilon_s)*np.sin(phi),np.abs(RTE)**2,'b',label=r'$\rho_{\rm TE}$')
                a1.set_ylabel(r'$\rho_{\rm TE}$, $\tau_{\rm TE}$')
                if epsilon_fa1_imag != 0 or epsilon_fa2_imag != 0 or epsilon_f_imag != 0 or epsilon_fb1_imag != 0 or epsilon_fb2_imag != 0:
                    a1.plot(np.sqrt(epsilon_s)*np.sin(phi),1-np.abs(RTE)**2-tauTE,'g',label=r'$a_{\rm TE}$')
                    a1.set_ylabel(r'$\rho_{\rm TE}$, $\tau_{\rm TE}$, $a_{\rm TE}$')
                    a1.set_ylim([np.minimum(-0.025,-0.025 + np.amin(1-np.abs(RTE)**2-tauTE)),0.025 + np.amax([1,np.amax(tauTE),np.amax(np.abs(RTE)**2)])])
                a1.legend()
            
            a2 = f.add_subplot(212)
            a2.plot(np.sqrt(epsilon_s)*np.sin(phi),tauTM,'r',label=r'$\tau_{\rm TM}$')
            a2.set_xlabel(r'$\sqrt{\varepsilon_{\rm s}}\sin\varphi_{\rm i}=k_x \lambda/(2\pi)$')
            a2.set_xlim([np.sqrt(epsilon_s)*sinphi_min,np.sqrt(epsilon_s)*sinphi_max])
            a2.set_ylabel(r'$\tau_{\rm TM}$') 
            a2.set_ylim([-0.025,1.025*np.maximum(1,np.amax(tauTM))])
            if var_string[21].get()=="show_all":
                a2.plot(np.sqrt(epsilon_s)*np.sin(phi),np.abs(RTM)**2,'b',label=r'$\rho_{\rm TM}$')
                a2.set_ylabel(r'$\rho_{\rm TM}$, $\tau_{\rm TM}$')
                if epsilon_fa1_imag != 0 or epsilon_fa2_imag != 0 or epsilon_f_imag != 0 or epsilon_fb1_imag != 0 or epsilon_fb2_imag != 0:
                    a2.plot(np.sqrt(epsilon_s)*np.sin(phi),1-np.abs(RTM)**2-tauTM,'g',label=r'$a_{\rm TM}$')
                    a2.set_ylabel(r'$\rho_{\rm TM}$, $\tau_{\rm TM}$, $a_{\rm TM}$')
                    a2.set_ylim([np.minimum(-0.025,-0.025 + np.amin(1-np.abs(RTM)**2-tauTM)),0.025 + np.amax([1,np.amax(tauTM),np.amax(np.abs(RTM)**2)])])
                a2.legend()
                     
            plt.tight_layout()  
            
#            plt.savefig('FP_kx.pdf',bbox_inches='tight',dpi=300, transparent=True)

            gui.copy_stringvar_vector(var_string,var_save)

            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[5,5])
canvas = gui.create_canvas(root,f)
canvas.draw()
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(22)
var_save = gui.create_stringvar_vector(22)

initialize()

row = 1
row = gui.create_entry(mainframe,u"substrate: \u03B5 =",var_string[0],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry(mainframe,u"mirror 1: N =",var_string[7],row)
row = gui.create_triple_entry(mainframe,u"d\u2081/\u03BB =",var_string[1],u"\u03B5\u2081' =",var_string[2],u"\u03B5\u2081'' =",var_string[3],row)
row = gui.create_triple_entry(mainframe,u"d\u2082/\u03BB =",var_string[4],u"\u03B5\u2082' =",var_string[5],u"\u03B5\u2082'' =",var_string[6],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_triple_entry(mainframe,u"cavity: D/\u03BB =",var_string[8],u"\u03B5' =",var_string[9],u"\u03B5'' =",var_string[10],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry(mainframe,u"mirror 2: N =",var_string[17],row)
row = gui.create_triple_entry(mainframe,u"d\u2081/\u03BB =",var_string[11],u"\u03B5\u2081' =",var_string[12],u"\u03B5\u2081'' =",var_string[13],row)
row = gui.create_triple_entry(mainframe,u"d\u2082/\u03BB =",var_string[14],u"\u03B5\u2082' =",var_string[15],u"\u03B5\u2082'' =",var_string[16],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry(mainframe,u"cladding: \u03B5 =",var_string[18],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_entry_with_latex(mainframe,r'$\sin\varphi_\mathrm{i}>$',var_string[19],r'$\sin\varphi_\mathrm{i}<$',var_string[20],row)
row = gui.create_checkbutton_with_latex(mainframe,r'show transmittance only','show_all','no_show_all',var_string[21],row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)