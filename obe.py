# TO DO: use gui.create_stringvar_vector etc. and remove global variables

import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
title = "Optical Bloch Equations"
root = Tk.Tk()
root.title(title)
  
def initialize():
    var_string[0].set("3")
    var_string[1].set("0")
    var_string[2].set("0")
    var_string[3].set("0")
    var_string[4].set("0")
    var_string[5].set("-1")
    var_string[6].set("PHASE")
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
        E0 = float(var_string[0].get())*np.sqrt(np.pi)
        C0 = float(var_string[1].get())
        Delta = float(var_string[2].get())
        T1inv = float(var_string[3].get())
        T2inv = float(var_string[4].get())
        GE = float(var_string[5].get())
        RI = var_string[6].get()

        if T1inv < 0: 
            gui.input_error("Energy relaxation time must be positive. Re-initializing with previous parameters ...",reinitialize)
        elif T2inv < 0:
            gui.input_error("Dephasing time must be positive. Re-initializing with previous parameters ...",reinitialize)
        elif np.abs(GE) > 1:
            gui.input_error("Equilibrium inversion has to be in the interval [-N,N]. Re-initializing with previous parameters ...",reinitialize)
        elif E0 <= 0:
            gui.input_error("Normalized pulse area must be positive. Re-initializing with previous parameters ...",reinitialize)
        elif E0 > 10*np.sqrt(np.pi):
            gui.input_error("Normalized pulse area too large. Re-initializing with previous parameters ...",reinitialize)
        elif np.abs(C0) > 10:
            gui.input_error("Chirp parameter too large. Re-initializing with previous parameters ...",reinitialize)
        elif np.abs(Delta) > 100:
            gui.input_error("Detuning parameter too large. Re-initializing with previous parameters ...",reinitialize)
        else:         
            def compute_E(T):
                return E0*np.exp(-(1+C0*1j)*T**2)

            def compute_rhs(T,A): 
                rhs1 = -T2inv*A[0]-Delta*A[1]+A[2]*np.imag(compute_E(T))
                rhs2 = -T2inv*A[1]+Delta*A[0]-A[2]*np.real(compute_E(T))
                rhs3 = -T1inv*(A[2]+1)+A[1]*np.real(compute_E(T))-A[0]*np.imag(compute_E(T))
                return np.array([rhs1, rhs2, rhs3])
        
            f.clf()
            A = np.zeros(3)
            A[0] = 0 # P1
            A[1] = 0 # P2
            A[2] = GE # gamma        

            sol = spi.solve_ivp(compute_rhs, [-3, 3], A, max_step = 1.e-2)

            P = sol.y[0,:]+1j*sol.y[1,:]
            E = compute_E(sol.t)  
            
            a1 = f.add_subplot(121)
            if RI == 'RI':
                a1.plot(sol.t, np.real(E)/E0, 'b', label=r'$\Re\tilde{E}/E_0$')
                a1.plot(sol.t, np.imag(E)/E0, 'r', label=r'$\Im\tilde{E}/E_0$')
                a1.set_ylabel(r'$\Re\tilde{E}/E_0$,$\Im\tilde{E}/E_0$, $\gamma^{\rm I}/N$') 
            else:
                a1.plot(sol.t, np.abs(E)/E0, 'b', label=r'$|\tilde{E}|/E_0$')
                a1.set_ylabel(r'$|\tilde{E}|/E_0$, $\gamma^{\rm I}/N$')
            a1.plot(sol.t, sol.y[2,:], 'g', label=r'$\gamma^{\rm I}/N$')
            a1.legend()
            a1.set_xlabel(r'$t/T_{\rm p}$')
            
        
            a2 = f.add_subplot(122)
            if RI == 'RI':
                lns1 = a2.plot(sol.t, np.real(P), 'b', label=r'$\Re\tilde{P}/(dN)$')
                lns2 = a2.plot(sol.t, np.imag(P), 'r', label=r'$\Im\tilde{P}/(dN)$')
                a2.set_ylabel(r'$\Re\tilde{P}/(dN)$, $\Im\tilde{P}/(dN)$')
            else:
                lns1 = a2.plot(sol.t, np.abs(P), 'b', label=r'$|\tilde{P}|/(dN)$')
                a2.set_ylabel(r'$|\tilde{P}|/(dN)$')
                a2bis = a2.twinx()
                lns2 = a2bis.plot(sol.t[1:], np.angle(P[1:]/E[1:]), 'r', label=r'$\mathrm{arg}\tilde{P}-\mathrm{arg}\tilde{E}$')
                a2bis.set_ylabel(r'arg$\tilde{P}$ - arg$\tilde{E}$')
                a2bis.set_yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi])
                a2bis.set_yticklabels([r'$-\pi$', r'$-\frac{\pi}{2}$', r'$0$', r'$\frac{\pi}{2}$', r'$\pi$'])
                a2bis.set_ylim([-1.05*np.pi,1.05*np.pi])
            lns = lns1+lns2
            labs = [l.get_label() for l in lns]
            a2.legend(lns, labs)
            a2.set_xlabel(r'$t/T_{\rm p}$')
              
            plt.tight_layout()
           
#            plt.savefig('obe.pdf',bbox_inches='tight',dpi=300, transparent=True)
            
            gui.copy_stringvar_vector(var_string,var_save)

            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing with previous parameters ...", reinitialize)
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[9,3])
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(7)
var_save = gui.create_stringvar_vector(7)

initialize()

row = 1
row = gui.create_entry_with_latex(mainframe,r"normalized pulse area $\alpha=\frac{d}{\sqrt{\pi}\hbar}T_{\rm p}E_0$",var_string[0],row)
row = gui.create_entry_with_latex(mainframe,r"chirp parameter $C_0=$",var_string[1],row)
row = gui.create_entry_with_latex(mainframe,r"detuning $T_{\rm p}\Delta=$",var_string[2],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry_with_latex(mainframe,r"dephasing $T_{\rm p}/T_2=$",var_string[4],row)
row = gui.create_entry_with_latex(mainframe,r"inversion relaxation $T_{\rm p}/T_1=$",var_string[3],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry_with_latex(mainframe,r"initial inversion $\gamma^{\rm I}(t=-3T_{\rm p})/N=$",var_string[5],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_checkbutton(mainframe,"show Re and Im",'PHASE','RI',var_string[6],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)

