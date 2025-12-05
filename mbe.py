import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
import tkinter as Tk
import sys
sys.path.append('./stuff')
import gui_stuff as gui

gui.set_rcParams()
title = "Maxwell-Bloch Equations"
root = Tk.Tk()
root.title(title)
  
def initialize():
    var_string[0].set("2")
    var_string[1].set("0")
    var_string[2].set("0")
    var_string[3].set("0")
    var_string[4].set("0")
    var_string[5].set("-1")
    var_string[6].set("4")
    var_string[7].set("sech")
    
    gui.copy_stringvar_vector(var_string,var_save) 
    calculate()
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    calculate()
    
def show_manual():
    gui.show_manual("man/mbe.png",title)

def calculate():
    gui.change_cursor(root,"trek")
    try:
        E0 = float(var_string[0].get())*np.sqrt(np.pi)
        C0 = float(var_string[1].get())
        Delta = float(var_string[2].get())
        T1inv = float(var_string[3].get())
        T2inv = float(var_string[4].get())
        GE = float(var_string[5].get())
        Lz = float(var_string[6].get())
        profile = var_string[7].get()
    
        if T1inv < 0: 
            gui.input_error("Energy relaxation time have to be positive. Re-initializing with previous parameters ...",reinitialize)
        elif T2inv < 0:
            gui.input_error("Dephasing time has to be positive. Re-initializing with previous parameters ...",reinitialize)
        elif np.abs(GE) > 1:
            gui.input_error("Equilibrium inversion has to be in the interval [-N,N]. Re-initializing with previous parameters ...",reinitialize)
        elif np.abs(E0) > 10*np.sqrt(np.pi):
            gui.input_error("Normalized pulse area too large. Re-initializing with previous parameters ...",reinitialize)
        elif np.abs(C0) > 10:
            gui.input_error("Chirp parameter too large. Re-initializing with previous parameters ...",reinitialize)
        elif np.abs(Delta) > 100:
            gui.input_error("Detuning parameter too large. Re-initializing with previous parameters ...",reinitialize)
        elif Lz <= 0:
            gui.input_error("Propagation length myst be positive. Re-initializing with previous parameters ...",reinitialize)
        elif Lz > 30:
            gui.input_error("Propagation length too large. Re-initializing with previous parameters ...",reinitialize)
        else:                 
            Nt = 512
            Nz = 1024
            Lt = 6
            T = np.linspace(-Lt/2,Lt/2,Nt,endpoint=True)
            Z = np.linspace(0,Lz,Nz,endpoint=True)
            E = np.zeros((Nt,Nz),dtype=np.complex128)
            if profile == 'Gaussian':
                E[:,0] = E0*np.exp(-(1+C0*1j)*T**2)
            else:
                E[:,0] = E0/np.cosh(np.sqrt(np.pi)*T)*np.exp(-C0*1j*T**2) # soliton with same amplitude and pulse area as unchirped Gaussian

            def compute_gamma(E_local):
                def compute_E(T_local):
                    return np.interp(T_local, T, E_local)
    
                def compute_rhs_obe(T_local,A): 
                    rhs1 = -T2inv*A[0]-Delta*A[1]+A[2]*np.imag(compute_E(T_local))
                    rhs2 = -T2inv*A[1]+Delta*A[0]-A[2]*np.real(compute_E(T_local))
                    rhs3 = -T1inv*(A[2]+1)+A[1]*np.real(compute_E(T_local))-A[0]*np.imag(compute_E(T_local))
                    return np.array([rhs1, rhs2, rhs3])
                
                A = np.zeros(3)
                A[0] = 0 # P1
                A[1] = 0 # P2
                A[2] = GE # gamma      
             
                solobe = spi.solve_ivp(compute_rhs_obe, [T[0], T[-1]], A, t_eval=T)
    
                return solobe.y[2,:]       
            
            def compute_P(B):
                def compute_E(T_local):
                    return np.interp(T_local, T, B[0:Nt]+1j*B[Nt:2*Nt])
    
                def compute_rhs_obe(T_local,A): 
                    rhs1 = -T2inv*A[0]-Delta*A[1]+A[2]*np.imag(compute_E(T_local))
                    rhs2 = -T2inv*A[1]+Delta*A[0]-A[2]*np.real(compute_E(T_local))
                    rhs3 = -T1inv*(A[2]+1)+A[1]*np.real(compute_E(T_local))-A[0]*np.imag(compute_E(T_local))
                    return np.array([rhs1, rhs2, rhs3])
                
                A = np.zeros(3)
                A[0] = 0 # P1
                A[1] = 0 # P2
                A[2] = GE # gamma      
             
                solobe = spi.solve_ivp(compute_rhs_obe, [T[0], T[-1]], A, t_eval=T)
    
                return solobe.y[0,:]+1j*solobe.y[1,:]             
            
            def compute_rhs_mbe(Z,B): 
                rhs1 = -np.imag(compute_P(B))
                rhs2 = np.real(compute_P(B))
                return np.concatenate((rhs1, rhs2))
        
            f.clf()
            B = np.zeros(2*Nt)
            B[0:Nt] = np.real(E[:,0]) # E1
            B[Nt:2*Nt] = np.imag(E[:,0]) # E2
            
            solmbe = spi.solve_ivp(compute_rhs_mbe, [Z[0], Z[-1]], B, t_eval=Z)
    
            E = solmbe.y[0:Nt,:]+1j*solmbe.y[Nt:2*Nt,:]
            
            a1 = plt.subplot2grid((5, 4), (1, 1), rowspan=4, colspan=2)
            im = a1.imshow(np.abs(E/E0)**2 ,extent=[Z[0], Z[-1], T[0], T[-1]] , aspect='auto', origin='lower', vmin=0, cmap='jet')
            a1.set_xlabel(r'$T_{\rm p} \kappa z$')
            a1.set_ylabel(r'$\tau/T_{\rm p}$')
            ac = plt.subplot2grid((5, 4), (0, 1), colspan=2)
            f.colorbar(im, cax=ac, orientation='horizontal')

            a2 = plt.subplot2grid((5, 4), (1, 0), rowspan=4)       
            a2.plot(np.abs(E[:,0])/E0,T,'r', label=r'$|\tilde E(z=0)|/E_0$')
            a2.plot(compute_gamma(E[:,0]),T,'g', label=r'$\gamma^{\rm I}(z=0)/N$')
            a2.set_xlabel(r'$|\tilde E|/E_0$, $\gamma^{\rm I}/N$')
            a2.set_ylabel(r'$\tau/T_{\rm p}$')
            a2.set_ylim([-Lt/2,Lt/2])
            a2.legend(bbox_to_anchor=(0, 1.05, 1, 0), loc="lower right")
            
            a3 = plt.subplot2grid((5, 4), (1, 3), rowspan=4, sharey=a1)
            a3.plot(np.abs(E[:,-1])/E0,T,'r', label=r'$|\tilde E(z=z_{\rm L})|/E_0$')
            a3.plot(compute_gamma(E[:,-1]),T,'g', label=r'$\gamma^{\rm I}(z=z_{\rm L})/N$')
            a3.set_xlabel(r'$|\tilde E|/E_0$, $\gamma^{\rm I}/N$')
            a3.set_ylabel(r'$\tau/T_{\rm p}$')
            a3.set_ylim([-Lt/2,Lt/2])
            a3.legend(bbox_to_anchor=(0, 1.05, 1, 0), loc="lower right")
              
            plt.tight_layout()
            
            ac.set_position([0.35,0.825,0.3,0.025])  
            ac.xaxis.tick_top()
            ac.set_xlabel(r'$|\tilde E/E_0|^2$')
            ac.xaxis.set_label_position('top') 
           
#            plt.savefig('mbe.pdf',bbox_inches='tight',dpi=300, transparent=True)
            
            gui.copy_stringvar_vector(var_string,var_save)
    
            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing with previous parameters ...", reinitialize)
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[10,5])
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(8)
var_save = gui.create_stringvar_vector(8)

initialize()

row = 1
row = gui.create_entry_with_latex(mainframe,r"normalized pulse area $\alpha=\frac{d}{\sqrt{\pi}\hbar}T_{\rm p}E_0$",var_string[0],row)
row = gui.create_entry_with_latex(mainframe,r"chirp parameter $C_0=$",var_string[1],row)
row = gui.create_entry_with_latex(mainframe,r"detuning $T_{\rm p}\Delta=$",var_string[2],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry_with_latex(mainframe,r"dephasing $T_{\rm p}/T_2=$",var_string[4],row)
row = gui.create_entry_with_latex(mainframe,r"inversion relaxation $T_{\rm p}/T_1=$",var_string[3],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry_with_latex(mainframe,r"initial inversion $\gamma^{\rm I}(\tau=-3T_{\rm p})/N=$",var_string[5],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry_with_latex(mainframe,r"propagation length $Z_L=T_{\rm p} \kappa z_{\rm L}=$",var_string[6],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_radiobutton_single_column_with_latex(mainframe,[r'input pulse:',r'$E_0\,\mathrm{sech}\!\left(\sqrt{\pi}\tau/T_\mathrm{p}\right)\exp\!\left(-iC_0\tau^2/T^2_\mathrm{p}\right)$',r'$E_0\exp\!\left[-(1+iC_0)\tau^2/T^2_\mathrm{p}\right]$'],[u'sech','Gaussian'],var_string[7],2,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)

