import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import scipy.integrate as spi
import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
root = Tk.Tk()
root.title("3-wave mixing -- SFG in UDPA")

def initialize():
    global var_save
    var_string[0].set("1") # L/L_c
    var_string[1].set("1") # L/L_nl
    var_string[2].set("noshow") # show exact solution
    var_string[3].set("noshow") # show I_P
    var_string[4].set("0.01") # I_1/I_P
    var_string[5].set("1.1") # omega_{\rm 2}/omega_1
    gui.copy_stringvar_vector(var_string,var_save)    
    calculate()
    
def reinitialize():
    global var_string
    gui.copy_stringvar_vector(var_save,var_string)
    calculate()

def calculate():
    global var_save
    try:
        LLc = float(var_string[0].get())
        LLnl = float(var_string[1].get())
        I1IP = float(var_string[4].get())
        omega2omega1 = float(var_string[5].get())
        
        if LLnl <= 0:
            gui.input_error("Pump amplitude must be positive. Re-initializing with previous parameters...",reinitialize)
        elif I1IP <= 0: 
            gui.input_error("Intensity ratio must be positive. Re-initializing with previous parameters...",reinitialize)
        elif omega2omega1 <= 0:
            gui.input_error("Frequency ratio must be positive. Re-initializing with previous parameters...",reinitialize)
        else: 
            if LLnl > 100:
                gui.input_error("Pump amplitude too large. Reducing...")
                LLnl = 100
                var_string[1].set(LLnl)
                
            Z = np.linspace(0, LLnl, num=1000)
            Ldeltakh = LLc*np.pi/2
            LGamma = np.sqrt(LLnl**2 + Ldeltakh**2)
            I1N = np.abs(np.cos(LGamma * Z / LLnl) - 1j * Ldeltakh / LGamma * np.sin(LGamma * Z / LLnl))**2
            I3N = np.abs(LLnl / LGamma * np.sin(LGamma * Z / LLnl))**2

            f.clf() 
        
            a1 = f.add_subplot(gs[1:, 0])
            lns = a1.plot(Z, I1N, 'r', label=r'$I^{\rm UDPA}_1/I_{10}$')
            lns1 = a1.plot(Z, I3N, 'b', label=r'$\omega_1 I^{\rm UDPA}_3/(\omega_3I_{10})$')
            lns = lns + lns1
            if var_string[3].get() == 'showIP':
                lns1 = a1.plot([Z[0],Z[-1]], [1/I1IP,1/I1IP], 'g', label=r'$I^{\rm UDPA}_{\rm P}/I_{10}$')
                lns = lns + lns1     

            if var_string[2].get() == 'showexact':
                omega3omega2 = 1 + 1 / omega2omega1
                omega3omega1 = 1 + omega2omega1
                LLnl_exact = LLnl * np.sqrt((1 + I1IP) / omega3omega2)
                s = LLc/LLnl_exact * np.pi/2
                A = np.zeros(3) + 0j 
                I1I = I1IP/ (1 + I1IP)
                IPI = 1 / (1 + I1IP)
                A[0] = np.sqrt(I1I * omega3omega1)
                A[1] = np.sqrt(IPI * omega3omega2)
        
                def compute_rhs(Z,A): # computes rhs of ode system, A[j-1] = A_j
                    rhs1 = 1j * A[2] * np.conj(A[1])
                    rhs2 = 1j * A[2] * np.conj(A[0])
                    rhs3 = 1j * (A[0] * A[1] - 2 * s * A[2])
                    return np.array([rhs1, rhs2, rhs3])
        
                sol = spi.solve_ivp(compute_rhs, [0, LLnl_exact], A, max_step = 1.e-3*LLnl_exact)

                lns1 = a1.plot(sol.t * LLnl / LLnl_exact, np.abs(sol.y[0,:])**2 / omega3omega1 / I1I, 'r:', label=r'$I_1/I_{10}$')
                lns = lns + lns1
                lns1 = a1.plot(sol.t * LLnl / LLnl_exact, np.abs(sol.y[2,:])**2 / omega3omega1 / I1I, 'b:', label=r'$\omega_1I_3/(\omega_3I_{10})$')
                lns = lns + lns1      
                if var_string[3].get() == 'showIP':
                    lns1 = a1.plot(sol.t * LLnl / LLnl_exact, np.abs(sol.y[1,:])**2 / omega3omega2 / I1I, 'g:', label=r'$I_{\rm P}/I_{10}$')
                    lns = lns + lns1                      
            
            a1.set_xlim([0,LLnl])
            a1.set_xlabel(r'$Z = z/L_{\rm nl}$')
            a1.set_ylabel(r'Normalized Intensities')
            labs = [l.get_label() for l in lns]
            a1.legend(lns, labs, bbox_to_anchor=(0, 1.1, 1, 0), loc="lower left", mode="expand", ncol=2)
        
#        plt.savefig('3wm_sfg.pdf',bbox_inches='tight',dpi=300, transparent=True)
        
        gui.copy_stringvar_vector(var_string,var_save)

        canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)

f = plt.figure(1,[7,5])
gs = GridSpec(4, 1, figure=f)
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(6)
var_save = gui.create_stringvar_vector(6)

initialize()

row = 1
row = gui.create_description(mainframe,'phase mismatch:',row)
row = gui.create_entry_with_latex(mainframe,r'$\textrm{sgn}(\Delta k) L / L_{\rm c} = L \Delta k / \pi =$',var_string[0],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_description(mainframe,'nonlinear interaction strength:',row)
row = gui.create_entry_with_latex(mainframe,r'$L / L_{\rm nl} = L \chi \sqrt{\omega_1\omega_3I_{\rm P0}} = $',var_string[1],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_checkbutton_with_latex(mainframe,r'show exact solution','noshow','showexact',var_string[2],r'show $I_{\rm P}$','noshow','showIP',var_string[3],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_description(mainframe,'initial condition:',row)
row = gui.create_entry_with_latex(mainframe,r'$I_{10} / I_{\rm P0} = $',var_string[4],row)
row = gui.create_description(mainframe,'frequency ratio:',row)
row = gui.create_entry_with_latex(mainframe,r'$\omega_2 / \omega_1 = $',var_string[5],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)