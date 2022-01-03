import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import scipy.integrate as spi
import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
root = Tk.Tk()
root.title("3-wave mixing -- DFG in UDPA")

def initialize():
    global var_save
    var_string[0].set("1") # L/L_c
    var_string[1].set("1") # L/L_nl
    var_string[2].set("noshow") # show exact solution
    var_string[3].set("noshow") # show I_P
    var_string[4].set("0.01") # I_S/I_P
    var_string[5].set("1.1") # omega_{\rm I}/omega_{\rm S}
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
        ISIP = float(var_string[4].get())
        omegaIomegaS = float(var_string[5].get())
        
        if LLnl <= 0:
            gui.input_error("Pump amplitude must be positive. Re-initializing with previous parameters...",reinitialize)
        elif ISIP <= 0: 
            gui.input_error("Intensity ratio must be positive. Re-initializing with previous parameters...",reinitialize)
        elif omegaIomegaS <= 0:
            gui.input_error("Frequency ratio must be positive. Re-initializing with previous parameters...",reinitialize)
        else: 
            if LLnl > 100:
                gui.input_error("Pump amplitude too large. Reducing...")
                LLnl = 100
                var_string[1].set(LLnl)
                
            Z = np.linspace(0, LLnl, num=1000)
            Ldeltakh = LLc*np.pi/2
            Lg = np.sqrt(LLnl**2 - Ldeltakh**2 + 0j)
            ISN = np.abs(np.cosh(Lg * Z / LLnl) + 1j * Ldeltakh / Lg * np.sinh(Lg * Z / LLnl))**2
            IIN = np.abs(LLnl / Lg * np.sinh(Lg * Z / LLnl))**2

            f.clf() 
        
            a1 = f.add_subplot(gs[1:, 0])
            lns = a1.plot(Z, ISN, 'r', label=r'$I^{\rm UDPA}_{\rm S}/I_{\rm S0}$')
            lns1 = a1.plot(Z, IIN, 'g', label=r'$\omega_{\rm S}I^{\rm UDPA}_{\rm I}/(\omega_{\rm I}I_{\rm S0})$')
            lns = lns + lns1
            if var_string[3].get() == 'showIP':
                lns1 = a1.plot([Z[0],Z[-1]], [1/ISIP,1/ISIP], 'b', label=r'$I^{\rm UDPA}_{\rm P}/I_{\rm S0}$')
                lns = lns + lns1     

            if var_string[2].get() == 'showexact':
                LLnl_exact = LLnl * np.sqrt(1 + ISIP)
                s = LLc/LLnl_exact * np.pi/2
                A = np.zeros(3) + 0j 
                omegaSomegaP = 1 / (1 + omegaIomegaS)
                ISI = ISIP/ (1 + ISIP)
                IPI = 1 / (1 + ISIP)
                A[0] = np.sqrt(ISI / omegaSomegaP)
                A[2] = np.sqrt(IPI)
        
                def compute_rhs(Z,A): # computes rhs of ode system, A[j-1] = A_j
                    rhs1 = 1j * A[2] * np.conj(A[1])
                    rhs2 = 1j * A[2] * np.conj(A[0])
                    rhs3 = 1j * (A[0] * A[1] - 2 * s * A[2])
                    return np.array([rhs1, rhs2, rhs3])
        
                sol = spi.solve_ivp(compute_rhs, [0, LLnl_exact], A, max_step = 1.e-3*LLnl_exact)

                lns1 = a1.plot(sol.t * LLnl / LLnl_exact, np.abs(sol.y[0,:])**2 * omegaSomegaP / ISI, 'r:', label=r'$I_{\rm S}/I_{\rm S0}$')
                lns = lns + lns1
                lns1 = a1.plot(sol.t * LLnl / LLnl_exact, np.abs(sol.y[1,:])**2 * omegaSomegaP / ISI, 'g:', label=r'$\omega_{\rm S}I_{\rm I}/(\omega_{\rm I}I_{\rm S0})$')
                lns = lns + lns1      
                if var_string[3].get() == 'showIP':
                    lns1 = a1.plot(sol.t * LLnl / LLnl_exact, np.abs(sol.y[2,:])**2 / ISI, 'b:', label=r'$I_{\rm P}/I_{\rm S0}$')
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
row = gui.create_entry_with_latex(mainframe,r'$L / L_{\rm nl} = L \chi \sqrt{\omega_{\rm S}\omega_{\rm I}I_{\rm P0}} = $',var_string[1],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_checkbutton_with_latex(mainframe,r'show exact solution','noshow','showexact',var_string[2],r'show $I_{\rm P}$','noshow','showIP',var_string[3],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_description(mainframe,'initial condition:',row)
row = gui.create_entry_with_latex(mainframe,r'$I_{\rm S0} / I_{\rm P0} = $',var_string[4],row)
row = gui.create_description(mainframe,'frequency ratio:',row)
row = gui.create_entry_with_latex(mainframe,r'$\omega_{\rm I} / \omega_{\rm S} = $',var_string[5],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)