import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import scipy.integrate as spi
import tkinter as Tk
import sys
sys.path.append('./aux')
import gui_stuff as gui

gui.set_rcParams()
title = "Three-Wave Mixing in UDPA -- Type I SHG"
root = Tk.Tk()
root.title(title)

def initialize():
    var_string[0].set("5") # L/L_c
    var_string[1].set("1") # L/L_nl
    var_string[2].set("noshow") # show exact solution
    var_string[3].set("noshow") # show I_P
    gui.copy_stringvar_vector(var_string,var_save)    
    calculate()
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    calculate()
    
def show_manual():
    gui.show_manual("man/3wm_shg.png",title) 

def calculate():
    gui.change_cursor(root,"trek")
    try:
        LLc = float(var_string[0].get())
        LLnl = float(var_string[1].get())
        
        if LLnl <= 0:
            gui.input_error("Nonlinear interaction strength must be positive. Re-initializing with previous parameters...",reinitialize)
        else: 
            if LLnl > 100:
                gui.input_error("Pump amplitude too large. Reducing...")
                LLnl = 100
                var_string[1].set(LLnl)
                
            Z = np.linspace(0, LLnl, num=1000)
            Ldeltakh = LLc*np.pi/2
            I2N = (Z * np.sinc(Ldeltakh * Z / LLnl / np.pi))**2

            f.clf() 
        
            a1 = f.add_subplot(gs[1:, 0])
            lns = a1.plot(Z, I2N, 'b', label=r'$I^{\rm UDPA}_{2 \omega}/I_{\rm P0}$')
            if var_string[3].get() == 'showIP':
                lns1 = a1.plot([Z[0],Z[-1]], [1,1], 'r', label=r'$I^{\rm UDPA}_{\rm P}/I_{\rm P0}$')
                lns = lns + lns1     

            if var_string[2].get() == 'showexact':
                s = LLc/LLnl * np.pi/2
                A = np.zeros(2) + 0j 
                A[0] = 1
        
                def compute_rhs(Z,A): # computes rhs of ode system, A[j-1] = A_j
                    rhs1 = 1j * A[1] * np.conj(A[0])
                    rhs2 = 1j * (A[0]**2 - 2 * s * A[1])
                    return np.array([rhs1, rhs2])
        
                sol = spi.solve_ivp(compute_rhs, [0, LLnl], A, max_step = 1.e-3*LLnl)

                lns1 = a1.plot(sol.t, np.abs(sol.y[1,:])**2, 'b:', label=r'$I_{2\omega}/I_{\rm P0}$')
                lns = lns + lns1      
                if var_string[3].get() == 'showIP':
                    lns1 = a1.plot(sol.t, np.abs(sol.y[0,:])**2, 'r:', label=r'$I_{\rm P}/I_{\rm P0}$')
                    lns = lns + lns1                      
                    
            a1.set_xlim([0,LLnl])
            ylim = a1.get_ylim()
            ymax = np.amin([2,ylim[1]])
            a1.set_ylim([-ymax/20,ymax])
            a1.set_xlabel(r'$Z = z/L_{\rm nl}$')
            a1.set_ylabel(r'Normalized Intensities')
            labs = [l.get_label() for l in lns]
            a1.legend(lns, labs, bbox_to_anchor=(0, 1.1, 1, 0), loc="lower left", mode="expand", ncol=2)
        
#        plt.savefig('3wm_shg.pdf',bbox_inches='tight',dpi=300, transparent=True)
        
        gui.copy_stringvar_vector(var_string,var_save)

        canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[7,5])
gs = GridSpec(4, 1, figure=f)
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(4)
var_save = gui.create_stringvar_vector(4)

initialize()

row = 1
row = gui.create_description(mainframe,'phase mismatch:',row)
row = gui.create_entry_with_latex(mainframe,r'$\textrm{sgn}(\Delta k)L / L_{\rm c} = L \Delta k / \pi =$',var_string[0],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_description(mainframe,'nonlinear interaction strength:',row)
row = gui.create_entry_with_latex(mainframe,r'$L / L_{\rm nl} = L |\chi| \omega \sqrt{I_{\rm P0}} = $',var_string[1],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_checkbutton_with_latex(mainframe,r'show exact solution','noshow','showexact',var_string[2],r'show $I_{\rm P}$','noshow','showIP',var_string[3],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)