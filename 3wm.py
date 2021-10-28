import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
import scipy.integrate as spi
import tkinter as Tk
import gui_stuff as gui

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rc('text', usetex=True)
mpl.rc('text.latex', preamble=r'\usepackage{cmbright}')
mpl.rcParams.update({'font.size': 10})

root = Tk.Tk()
root.title("3-wave mixing -- rigorous solution")

def initialize():
    global var_save
    var_string[0].set("1") # L/L_c
    var_string[1].set("1") # L/L_nl
    var_string[2].set("2") # rho_2/rho_1
    var_string[3].set("0") # rho_3/rho_1
    var_string[4].set("0") # theta
    var_string[5].set("2") # omega_2/omega_1
    var_string[6].set("showrho2") # show rho_2
    var_string[7].set("showrho3") # show rho_3
    var_string[8].set("noshow") # show theta
    var_string[9].set("noshow") # show Manley Rowe
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
        rho2rho1 = float(var_string[2].get())
        rho3rho1 = float(var_string[3].get())
        theta = float(var_string[4].get())
        omega2omega1 = float(var_string[5].get())
        
        if LLnl <= 0:
            gui.input_error("Nonlinear length must be positive. Re-initializing with previous parameters...",reinitialize)
        elif rho2rho1 < 0 or rho3rho1 < 0: 
            gui.input_error("Amplitude ratios must not be negative. Re-initializing with previous parameters...",reinitialize)
        elif omega2omega1 < 0:
            gui.input_error("Frequency ratio must not be negative. Re-initializing with previous parameters...",reinitialize)
        else: 
            if LLnl > 100:
                gui.input_error("Total amplitude too large. Reducing...")
                LLnl = 100
                var_string[1].set(LLnl)
                
            omega1omega3 = 1/(1+omega2omega1)
            omega2omega3 = 1/(1+1/omega2omega1)
            rho1 = np.sqrt(1/(omega1omega3 + omega2omega3 * rho2rho1**2 + rho3rho1**2))
            s = LLc/LLnl * np.pi/2
            A = np.zeros(3) + 0j 
            A[0] = rho1 * np.exp(1j * theta)
            A[1] = rho2rho1 * rho1
            A[2] = rho3rho1 * rho1
        
            def compute_rhs(Z,A): # computes rhs of ode system, A[j-1] = A_j
                rhs1 = 1j * A[2] * np.conj(A[1])
                rhs2 = 1j * A[2] * np.conj(A[0])
                rhs3 = 1j * (A[0] * A[1] - 2 * s * A[2])
                return np.array([rhs1, rhs2, rhs3])
        
            sol = spi.solve_ivp(compute_rhs, [0, LLnl], A, max_step = 1.e-3*LLnl)

            f.clf() 
        
            a1 = f.add_subplot(gs[1:, 0])
            lns = a1.plot(sol.t, np.abs(sol.y[0,:])**2, 'r', label=r'$\rho_1^2=\omega_3 I_1/(\omega_1 I)$')
            if var_string[6].get() == 'showrho2':
                lns1 = a1.plot(sol.t, np.abs(sol.y[1,:])**2, 'g', label=r'$\rho_2^2=\omega_3 I_2/(\omega_2 I)$')
                lns = lns + lns1
            if var_string[7].get() == 'showrho3':            
                lns1 = a1.plot(sol.t, np.abs(sol.y[2,:])**2, 'b', label=r'$\rho_3^2=I_3/I$')
                lns = lns + lns1
            if var_string[8].get() == 'showtheta':     
                a1bis = a1.twinx()
                lns1 = a1bis.plot(sol.t, np.angle(sol.y[0,:]) + np.angle(sol.y[1,:]) - np.angle(sol.y[2,:]), 'm', label=r'$\theta = \phi_1 + \phi_2 - \phi_3$')
                lns = lns + lns1
                a1bis.set_ylabel(r'$\theta = \phi_1 + \phi_2 - \phi_3$', color='m')
                a1bis.tick_params(axis='y', labelcolor='m') 
                a1bis.set_yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi])
                a1bis.set_yticklabels([r'$-\pi$', r'$-\frac{\pi}{2}$', r'$0$', r'$\frac{\pi}{2}$', r'$\pi$'])
                a1bis.set_ylim([-1.1*np.pi,1.1*np.pi])
            if var_string[9].get() == 'showMR': 
                if var_string[6].get() == 'showrho2': 
                    if np.abs(A[0]) >= np.abs(A[1]):
                        lns1 = a1.plot(sol.t, np.abs(sol.y[0,:])**2-np.abs(sol.y[1,:])**2, 'k:', label=r'$\rho_1^2-\rho_2^2$')
                        lns = lns + lns1
                    else:
                        lns1 = a1.plot(sol.t, np.abs(sol.y[1,:])**2-np.abs(sol.y[0,:])**2, 'k:', label=r'$\rho_2^2-\rho_1^2$')
                        lns = lns + lns1
                    if var_string[7].get() == 'showrho3':
                        lns1 = a1.plot(sol.t, np.abs(sol.y[1,:])**2+np.abs(sol.y[2,:])**2, 'k:', label=r'$\rho_2^2+\rho_3^2$')
                        lns = lns + lns1
                if var_string[7].get() == 'showrho3': 
                    lns1 = a1.plot(sol.t, np.abs(sol.y[0,:])**2+np.abs(sol.y[2,:])**2, 'k:', label=r'$\rho_1^2+\rho_3^2$')
                    lns = lns + lns1
                    
            a1.set_xlabel(r'$Z = z/L_{\rm nl}$')
            a1.set_ylabel(r'$\rho_1^2,~\rho_2^2,~\rho_3^2$')
            a1.set_title(r'$s = \pi L_{\rm nl} / (2 L_{\rm c}) =$ '+str(round(s,4)))
            labs = [l.get_label() for l in lns]
            a1.legend(lns, labs, bbox_to_anchor=(0, 1.1, 1, 0), loc="lower left", mode="expand", ncol=2)
        
#        plt.savefig('3wm.pdf',bbox_inches='tight',dpi=300, transparent=True)
        
        gui.copy_stringvar_vector(var_string,var_save)

        canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)

f = plt.figure(1,[7,5])
gs = GridSpec(4, 1, figure=f)
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(10)
var_save = gui.create_stringvar_vector(10)

initialize()

row = 1
row = gui.create_description(mainframe,'phase mismatch:',row)
row = gui.create_entry_with_latex(mainframe,r'$L / L_{\rm c} = L \Delta k / \pi =$',var_string[0],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_description(mainframe,'total amplitude:',row)
row = gui.create_entry_with_latex(mainframe,r'$L / L_{\rm nl} = L \chi \sqrt{\omega_1\omega_2I} = $',var_string[1],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_description(mainframe,'initial conditions:',row)
row = gui.create_entry_with_latex(mainframe,r'$\rho_2(Z=0) / \rho_1(Z=0) = $',var_string[2],row)
row = gui.create_entry_with_latex(mainframe,r'$\rho_3(Z=0) / \rho_1(Z=0) = $',var_string[3],row)
row = gui.create_entry_with_latex(mainframe,r'$\theta(Z=0) = $',var_string[4],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_description(mainframe,'frequency ratio:',row)
row = gui.create_entry_with_latex(mainframe,r'$\omega_2 / \omega_1 = $',var_string[5],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_checkbutton_with_latex(mainframe,r'show $\rho_2^2$','noshow','showrho2',var_string[6],r'show $\theta$','noshow','showtheta',var_string[8],row)
row = gui.create_double_checkbutton_with_latex(mainframe,r'show $\rho_3^2$','noshow','showrho3',var_string[7],r'Manley-Rowe','noshow','showMR',var_string[9],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)