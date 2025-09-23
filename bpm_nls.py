import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import bpm_stuff as bpm

gui.set_rcParams()
title = "Beam Propagation Method - Nonlinear Schroedinger Equation "
root = Tk.Tk()
root.title(title)

def initialize():
    var_string[0].set("128")   # Neta
    var_string[1].set("25")    # Leta
    var_string[2].set("200")   # NZ
    var_string[3].set("1")     # LZ in units of pi
    var_string[4].set("0")     # Nabs
    var_string[5].set("2")     # soliton order N
    var_string[6].set("+1")     # D
    var_string[7].set("+1")     # Gamma    
    var_string[8].set("sech")   # profile    
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
        Neta = int(float(var_string[0].get()))
        var_string[0].set(Neta)
        Leta = float(var_string[1].get())
        NZ = int(float(var_string[2].get())) + 1 # Number of points = Number of steps + 1
        var_string[2].set(NZ-1)
        LZ = float(var_string[3].get()) 
        Nabs = int(float(var_string[4].get()))
        var_string[4].set(Nabs)    
        N = float(var_string[5].get())
        D = float(var_string[6].get())    
        Gamma = float(var_string[7].get())
        profile = var_string[8].get()
        
        if Neta <= 1 or Leta <= 0 or NZ <= 1 or LZ <= 0: 
            gui.input_error("Box lengths or numbers of points invalid. Re-initializing ...",reinitialize)
        elif Nabs < 0 or 20*Nabs > Neta: 
            gui.input_error("Number of absober points must be nonzero and must not exceed Neta/20. Re-initializing ...",reinitialize)
        elif N <= 0: 
            gui.input_error("Soliton order must be positive. Re-initializing ...",reinitialize)    
        else:

            eta, delta_eta = np.linspace(-Leta/2,Leta/2,Neta,endpoint=False, retstep=True)
            if profile == 'sech':
                u0 =  1/np.cosh(eta)
            elif profile == 'Gaussian':
                u0 = np.exp(-eta**2)
            elif profile == 'super-Gaussian':
                u0 = np.exp(-eta**8)
            else:
                rng = np.random.default_rng()
                u0 = np.ones(Neta) + 0.001 * rng.standard_normal(Neta)
            Z, delta_Z = np.linspace(0,LZ,NZ,endpoint=True, retstep=True)
    
            u = bpm.propagation_nls(Neta,u0,delta_eta,NZ,delta_Z*2*np.pi**2,Nabs,N,D,Gamma,0,0) # multiplication of delta_z by 2pi^2 to take into account scaling of LZ
    
            f.clf()
    
            a1 = plt.subplot2grid((1, 4), (0, 1), colspan=2)
            a1.imshow(np.abs(np.transpose(u))**2 ,extent=[Z[0], Z[-1], eta[0], eta[-1]] , aspect='auto', origin='lower', vmin=0, cmap='jet')
            a1.set_xlabel(r'$Z/\pi$')
            a1.set_ylabel(r'$\eta$')
            
            a2 = plt.subplot2grid((1, 4), (0, 0))       
            keta = 2*np.pi*np.fft.fftshift(np.fft.fftfreq(Neta,delta_eta))
            U0 = np.fft.fftshift(np.fft.fft(np.fft.fftshift(u0)))
            lns1 = a2.plot(np.abs(U0)/np.max(np.abs(U0)),keta,'r--', label=r'$|U_0(K_\eta)|$')
            a2.set_xlabel(r'$|u_0|/N$, $|U_0|/\mathrm{max}|U_0|$')
            a2.set_ylabel(r'$K_\eta$')
            a2.set_ylim([keta[0],keta[-1]])
            a2.invert_yaxis()
            a2bis = a2.twinx()
            lns2 = a2bis.plot(np.abs(u0)/np.max(np.abs(u0)),eta, label=r'$|u_0(\eta)|$',color='b')
            a2bis.set_ylim([eta[0],eta[-1]])
            a2bis.invert_yaxis()
            a2bis.set(yticklabels=[]) 
            lns = lns2+lns1
            labs = [l.get_label() for l in lns]
            a2.legend(lns, labs, loc=1)
    
            a3 = plt.subplot2grid((1, 4), (0, 3), sharey=a1)
            lns1 = a3.plot(np.abs(u[-1,:])/np.max(np.abs(u0)),eta,'b', label=r'$|u(\eta,Z=L_Z)|$')
            a3.set_xlabel(r'$|u|/N, |U|/\mathrm{max}|U_0|$')
            a3.set_ylabel(r'$\eta$')
            a3.invert_yaxis()   
            a3bis = a3.twinx()
            U = np.fft.fftshift(np.fft.fft(np.fft.fftshift(u[-1,:])))
            lns2 = a3bis.plot(np.abs(U)/np.max(np.abs(U0)),keta,'r--', label=r'$|U(K_\eta,Z=L_Z)|$')
            a3bis.set_ylim([keta[0],keta[-1]])
            a3bis.set_ylabel(r'$K_\eta$')
            a3bis.invert_yaxis()
            lns = lns1+lns2
            labs = [l.get_label() for l in lns]
            a3.legend(lns, labs, loc=1)
            
            plt.tight_layout()
                
#            plt.savefig('bpm_nls.pdf',bbox_inches='tight',dpi=300, transparent=True)
            
            gui.copy_stringvar_vector(var_string,var_save)

            canvas.draw()       
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", initialize)
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[10,4])

canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(9)
var_save = gui.create_stringvar_vector(9)

initialize()

row = 1
row = gui.create_entry_with_latex(mainframe,r"transverse number of points $N_\eta=$",var_string[0],row)
row = gui.create_entry_with_latex(mainframe,r"transverse box width $L_\eta=$",var_string[1],row)
row = gui.create_entry_with_latex(mainframe,r"longitudinal number of steps $N_Z=$",var_string[2],row)
row = gui.create_entry_with_latex(mainframe,r"longitudinal box length $L_Z/\pi=$",var_string[3],row)
row = gui.create_entry_with_latex(mainframe,r"number of absorber points $N_{\rm abs}=$",var_string[4],row)
row = gui.create_formula_with_latex(mainframe,r'$\partial_Z u - \mathrm{i}\frac{D}{2}\partial^2_\eta u = $',r'$\mathrm{i} \Gamma |u|^2 u$',row)
row = gui.create_entry_with_latex(mainframe,r"amplitude (Soliton order) $N=$",var_string[5],row)
row = gui.create_radiobutton(mainframe,['sign of D:','+1','-1'],var_string[6],2,row)
row = gui.create_radiobutton(mainframe,[u'sign of \u0393:','+1','-1'],var_string[7],2,row)
row = gui.create_radiobutton_single_column(mainframe,[u'input beam profile:','sech','Gaussian','super-Gaussian','noisy plane wave'],var_string[8],4,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)