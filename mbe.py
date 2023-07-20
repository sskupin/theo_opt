import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
root = Tk.Tk()
root.title("Maxwell-Bloch Equations")
  
def initialize():
    global E0_save,C0_save,Delta_save,T1inv_save,T2inv_save,GE_save,ZL_save,Pulse_save
    E0_string.set("2")
    C0_string.set("0")
    Delta_string.set("0")
    T1inv_string.set("0")
    T2inv_string.set("0")
    GE_string.set("-1")
    ZL_string.set("6")
    Pulse_string.set("sech")
    
    E0_save = float(E0_string.get())
    C0_save = float(C0_string.get())
    Delta_save = float(Delta_string.get())
    T1inv_save = float(T1inv_string.get())
    T2inv_save = float(T2inv_string.get())
    GE_save = float(GE_string.get())
    ZL_save = float(ZL_string.get())
    Pulse_save = Pulse_string.get()
    
    calculate()
    
def reinitialize():
    global E0_save,C0_save,Delta_save,T1inv_save,T2inv_save,GE_save,ZL_save,Pulse_save
    E0_string.set(E0_save)
    C0_string.set(C0_save)
    Delta_string.set(Delta_save)
    T1inv_string.set(T1inv_save)
    T2inv_string.set(T2inv_save)
    GE_string.set(GE_save)
    ZL_string.set(ZL_save)
    Pulse_string.set(Pulse_save)
    
    calculate()

def calculate():
    global E0_save,C0_save,Delta_save,T1inv_save,T2inv_save,GE_save,ZL_save,Pulse_save
    try:
        E0 = float(E0_string.get())*np.sqrt(np.pi)
        C0 = float(C0_string.get())
        Delta = float(Delta_string.get())
        T1inv = float(T1inv_string.get())
        T2inv = float(T2inv_string.get())
        GE = float(GE_string.get())
        Lz = float(ZL_string.get())
        profile = Pulse_string.get()
    
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
        elif Lz <= 0:
            gui.input_error("Propagation length myst be positive. Re-initializing with previous parameters ...",reinitialize)
        elif Lz > 10:
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
            
            E0_save = E0_string.get()
            C0_save = C0_string.get()
            Delta_save = Delta_string.get()
            T1inv_save = T1inv_string.get()
            T2inv_save = T2inv_string.get()
            GE_save = GE_string.get() 
            ZL_save = ZL_string.get()
            Pulse_save = Pulse_string.get()
    
            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing with previous parameters ...", reinitialize)

f = plt.figure(1,[10,5])
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

E0_string = Tk.StringVar()
C0_string = Tk.StringVar()
Delta_string = Tk.StringVar()
T1inv_string = Tk.StringVar()
T2inv_string = Tk.StringVar()
GE_string = Tk.StringVar()
ZL_string = Tk.StringVar()
Pulse_string = Tk.StringVar()

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
row = gui.create_entry_with_latex(mainframe,r"Propagation length $Z_L=T_{\rm p} \kappa z_{\rm L}=$",ZL_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_radiobutton_single_column_with_latex(mainframe,[r'Input pulse:',r'$E_0\,\mathrm{sech}\!\left(\sqrt{\pi}\tau/T_\mathrm{p}\right)\exp\!\left(-iC_0\tau^2/T^2_\mathrm{p}\right)$',r'$E_0\exp\!\left[-(1+iC_0)\tau^2/T^2_\mathrm{p}\right]$'],[u'sech','Gaussian'],Pulse_string,2,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)

