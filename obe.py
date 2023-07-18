import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
root = Tk.Tk()
root.title("Optical Bloch Equations")
  
def initialize():
    global E0_save,C0_save,Delta_save,T1inv_save,T2inv_save,GE_save,RI_save
    E0_string.set("3")
    C0_string.set("0")
    Delta_string.set("0")
    T1inv_string.set("0")
    T2inv_string.set("0")
    GE_string.set("-1")
    RI_string.set("PHASE")
    
    E0_save = float(E0_string.get())
    C0_save = float(C0_string.get())
    Delta_save = float(Delta_string.get())
    T1inv_save = float(T1inv_string.get())
    T2inv_save = float(T2inv_string.get())
    GE_save = float(GE_string.get())
    RI_save = RI_string.get()
    
    calculate()
    
def reinitialize():
    global E0_save,C0_save,Delta_save,T1inv_save,T2inv_save,GE_save,RI_save
    E0_string.set(E0_save)
    C0_string.set(C0_save)
    Delta_string.set(Delta_save)
    T1inv_string.set(T1inv_save)
    T2inv_string.set(T2inv_save)
    GE_string.set(GE_save)
    RI_string.set(RI_save)
    
    calculate()

def calculate():
    global E0_save,C0_save,Delta_save,T1inv_save,T2inv_save,GE_save,RI_save
    try:
        E0 = float(E0_string.get())*np.sqrt(np.pi)
        C0 = float(C0_string.get())
        Delta = float(Delta_string.get())
        T1inv = float(T1inv_string.get())
        T2inv = float(T2inv_string.get())
        GE = float(GE_string.get())
        RI = RI_string.get()

        if T1inv < 0: 
            gui.input_error("Energy relaxation time have to be positive. Re-initializing with previous parameters ...",reinitialize)
        elif T2inv < 0:
            gui.input_error("Dephasing time has to be positive. Re-initializing with previous parameters ...",reinitialize)
        elif np.abs(GE) > 1:
            gui.input_error("Equilibrium inversion has to be in the interval [-N,N]. Re-initializing with previous parameters ...",reinitialize)
        elif np.abs(E0) > 10*np.sqrt(np.pi):
            gui.input_error("Normalized amplitude too large. Re-initializing with previous parameters ...",reinitialize)
        elif np.abs(C0) > 10:
            gui.input_error("Chirp parameter too large. Re-initializing with previous parameters ...",reinitialize)
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
            
            E0_save = E0_string.get()
            C0_save = C0_string.get()
            Delta_save = Delta_string.get()
            T1inv_save = T1inv_string.get()
            T2inv_save = T2inv_string.get()
            GE_save = GE_string.get() 
            RI_save = RI_string.get()

            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing with previous parameters ...", reinitialize)

f = plt.figure(1,[9,3])
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

E0_string = Tk.StringVar()
C0_string = Tk.StringVar()
Delta_string = Tk.StringVar()
T1inv_string = Tk.StringVar()
T2inv_string = Tk.StringVar()
GE_string = Tk.StringVar()
RI_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_entry_with_latex(mainframe,r"Normalized amplitude $\alpha=\frac{T_{\rm p}d}{\sqrt{\pi}\hbar}E_0$",E0_string,row)
row = gui.create_entry_with_latex(mainframe,r"Chirp parameter $C_0=$",C0_string,row)
row = gui.create_entry_with_latex(mainframe,r"Detuning $T_{\rm p}\Delta=$",Delta_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry_with_latex(mainframe,r"Dephasing $T_{\rm p}/T_2=$",T2inv_string,row)
row = gui.create_entry_with_latex(mainframe,r"Inversion relaxation $T_{\rm p}/T_1=$",T1inv_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry_with_latex(mainframe,r"Initial inversion $\gamma^{\rm I}(t=-3T_{\rm p})/N=$",GE_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_checkbutton(mainframe,"show Re and Im",'PHASE','RI',RI_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)

