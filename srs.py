import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import scipy.integrate as spi
import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
root = Tk.Tk()
root.title("SRS -- Stokes Anti-Stokes coupling in UDPA")

def initialize():
    global var_save
    var_string[0].set(".1") # \Delta k_{nl} / g_R
    var_string[1].set("10") # L*g_R
    var_string[2].set("1.2") # \delta
    var_string[3].set("0.01") # \kappa_L/g_R
    var_string[4].set("0") # I_{AS0}/I_{S0}
    var_string[5].set("0") # \varphi_{S0} + \varphi_{AS0}
    gui.copy_stringvar_vector(var_string,var_save)    
    calculate()
    
def reinitialize():
    global var_string
    gui.copy_stringvar_vector(var_save,var_string)
    calculate()

def calculate():
    global var_save
    gui.change_cursor(root,"trek")
    try:
        s = float(var_string[0].get())
        LLnl = float(var_string[1].get())
        delta = float(var_string[2].get())
        kappaL = float(var_string[3].get())
        IASIS = float(var_string[4].get())
        phiSAS = float(var_string[5].get())
        
        if LLnl <= 0:
            gui.input_error("Propagation range must be positive. Re-initializing with previous parameters...",reinitialize)
        elif delta <= 1: 
            gui.input_error("Frequency ratio must be larger than one. Re-initializing with previous parameters...",reinitialize)
        elif IASIS < 0: 
            gui.input_error("Intensity ratio must be non-negative. Re-initializing with previous parameters...",reinitialize)
        else: 
            if LLnl > 1000:
                gui.input_error("Propagation range too large. Reducing...")
                LLnl = 1000
                var_string[1].set(LLnl)
                
            kappa = np.sqrt(delta)/2 + 1j*kappaL

            f.clf() 
        
            a1 = f.add_subplot(gs[1:, 0])
            
            ### analytical solution for Stokes input only
            # Z = np.linspace(0, LLnl, num=1000)
            # sq = np.sqrt(((1-delta)/4)**2-s**2/4+kappaL**2-1j*((1+delta)*s/4+np.sqrt(delta)*kappaL))
            # lambda1 = ((1-delta)/4 + sq)/1j
            # lambda2 = ((1-delta)/4 - sq)/1j
            # print(lambda1,lambda2)
            # bs = 1/(lambda1-lambda2)*((lambda1-s/2-1j*delta/2)*np.exp(1j*lambda1*Z) - (lambda2-s/2-1j*delta/2)*np.exp(1j*lambda2*Z))*np.exp(1j*s*Z/2)
            # bas = 1j*kappa/(lambda1-lambda2)*(np.exp(1j*lambda1*Z)-np.exp(1j*lambda2*Z))*np.exp(-1j*s*Z/2)
            # a1.plot(Z, np.abs(bs)**2)
            # a1.plot(Z, np.abs(bas)**2)

            A = np.zeros(2) + 0j 
            A[0] = 1
            A[1] = np.sqrt(IASIS)*np.exp(1j*phiSAS)
        
            def compute_rhs(Z,A): # computes rhs of ode system, A[j-1] = A_j
                rhs1 = kappa * A[1] + (1/2 - 1j*s/2) * A[0]
                rhs2 = -kappa * A[0] - (delta/2 - 1j*s/2) * A[1]
                return np.array([rhs1, rhs2])
        
            sol = spi.solve_ivp(compute_rhs, [0, LLnl], A, max_step = 1.e-3*LLnl)

            lns = a1.plot(sol.t, np.abs(sol.y[0,:])**2, 'r', label=r'$I_{\rm S}/I_{\rm S0}$')
            lns1 = a1.plot(sol.t, np.abs(sol.y[1,:])**2, 'b', label=r'$\omega_{\rm S}I_{\rm AS}/(\omega_{\rm AS}I_{\rm S0})$')
            lns = lns + lns1                      
            
            a1.set_xlim([0,LLnl])
            a1.set_xlabel(r'$Z = g_{\rm R}z$')
            a1.set_ylabel(r'Normalized Intensities')
            labs = [l.get_label() for l in lns]
            a1.legend(lns, labs, bbox_to_anchor=(0, 1.1, 1, 0), loc="lower left", mode="expand", ncol=2)
        
        # plt.savefig('srs.pdf',bbox_inches='tight',dpi=300, transparent=True)
        
        gui.copy_stringvar_vector(var_string,var_save)

        canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[7,5])
gs = GridSpec(4, 1, figure=f)
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(6)
var_save = gui.create_stringvar_vector(6)

initialize()

row = 1
row = gui.create_description(mainframe,'nonlinear phase mismatch:',row)
row = gui.create_entry_with_latex(mainframe,r'$\Delta k_{\rm nl} / g_{\rm R} =$',var_string[0],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_description(mainframe,'interaction length:',row)
row = gui.create_entry_with_latex(mainframe,r'$g_{\rm R}L = $',var_string[1],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_description(mainframe,'frequency ratio:',row)
row = gui.create_entry_with_latex(mainframe,r'$\delta = $',var_string[2],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_description(mainframe,'laser SM:',row)
row = gui.create_entry_with_latex(mainframe,r'$\kappa_{\rm L}/ g_{\rm R}= $',var_string[3],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_description(mainframe,'initial condition:',row)
row = gui.create_entry_with_latex(mainframe,r'$\omega_{\rm S}I_{\rm AS0} / (\omega_{\rm AS}I_{\rm S0} )= $',var_string[4],row)
row = gui.create_entry_with_latex(mainframe,r'$\varphi_{\rm S0} +  \varphi_{\rm AS0} - 2 \varphi_{\rm L0} = $',var_string[5],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)