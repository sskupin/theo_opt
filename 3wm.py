import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
import scipy.integrate as spi
import scipy.special as sps
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
    var_string[0].set(".1") # s
    var_string[1].set("5") # L/L_nl
    var_string[2].set("2") # rho_1/rho_2
    var_string[3].set("1") # rho_3/rho_2
    var_string[4].set("1") # theta
    var_string[5].set("2") # omega_2/omega_1
    var_string[6].set("showI1") # show I_1
    var_string[7].set("showI2") # show I_2
    var_string[8].set("showI3") # show I_3
    var_string[9].set("noshow") # show theta
    gui.copy_stringvar_vector(var_string,var_save)    
    calculate()
    
def reinitialize():
    global var_string
    gui.copy_stringvar_vector(var_save,var_string)
    calculate()

def calculate():
    global var_save
    try:
        s = float(var_string[0].get())
        LLnl = float(var_string[1].get())
        I1I2 = float(var_string[2].get())
        I3I2 = float(var_string[3].get())
        theta = float(var_string[4].get())
        omega2omega1 = float(var_string[5].get())
        
        if LLnl <= 0:
            gui.input_error("Propagation range must be positive. Re-initializing with previous parameters...",reinitialize)
        elif I1I2 < 0 or I3I2 < 0: 
            gui.input_error("Intensity ratios must not be negative. Re-initializing with previous parameters...",reinitialize)
        elif omega2omega1 <= 0:
            gui.input_error("Frequency ratio must be positive. Re-initializing with previous parameters...",reinitialize)
        else: 
                
            omega3omega1 = 1 + omega2omega1
            omega3omega2 = 1 + 1 / omega2omega1
            rho2 = np.sqrt(omega3omega2/(1 + I1I2 + I3I2))
            rho1 = np.sqrt(omega2omega1 * I1I2) * rho2
            rho3 = np.sqrt(I3I2 / omega3omega2) * rho2
            
            A = np.zeros(3) + 0j 
            A[0] = rho1
            A[1] = rho2 * np.exp(1j * theta)
            A[2] = rho3
            
            H = rho1 * rho2 * rho3 * np.cos(theta) - s * rho3**2
            C2 = rho1**2 + rho3**2
            C3 = rho2**2 + rho3**2
            U = np.sort(np.roots([1, - C2 - C3 - s**2, C2 * C3 - 2 * s * H, - H**2]))
            m = (U[1] - U[0]) / (U[2] - U[0])
            Km = sps.ellipk(m)
            ZP = 4 * Km / np.sqrt(U[2] - U[0])
            
            if LLnl > 2*ZP:
                gui.input_error("Propagation range too large. Reducing ...")
                LLnl = 2*ZP
                var_string[1].set(LLnl)
        
            def compute_rhs(Z,A): # computes rhs of ode system, A[j-1] = A_j
                rhs1 = 1j * A[2] * np.conj(A[1])
                rhs2 = 1j * A[2] * np.conj(A[0])
                rhs3 = 1j * (A[0] * A[1] - 2 * s * A[2])
                return np.array([rhs1, rhs2, rhs3])
        
            sol = spi.solve_ivp(compute_rhs, [0, LLnl], A, max_step = 1.e-3*LLnl)

            f.clf() 
        
            a1 = f.add_subplot(gs[1:, 0])
# =============================================================================\
#             Plot analytical solution
#             Z = np.linspace(0, LLnl, num=1000)
#             alpha = -np.arcsin(np.sqrt((rho3**2-U[0])/(U[1]-U[0])))*np.sign(np.sin(theta))
#             F = sps.ellipkinc(alpha, m)
#             solution = U[0] + (U[1]-U[0])*sps.ellipj(np.sqrt(U[2]-U[0])*Z+F,m)[0]**2
#             a1.plot(Z, solution)
#             a1.plot(Z, (C3 - solution)/ omega3omega2)
#             a1.plot(Z, (C2 - solution)/ omega3omega1)
# =============================================================================
            lns = []
            if var_string[6].get() == 'showI1':
                lns1 = a1.plot(sol.t, np.abs(sol.y[0,:])**2 / omega3omega1, 'r', label=r'$I_1/I$')
                lns = lns + lns1
            if var_string[7].get() == 'showI2':
                lns1 = a1.plot(sol.t, np.abs(sol.y[1,:])**2 / omega3omega2, 'g', label=r'$I_2/I$')
                lns = lns + lns1                
            if var_string[8].get() == 'showI3':            
                lns1 = a1.plot(sol.t, np.abs(sol.y[2,:])**2, 'b', label=r'$I_3/I$')
                lns = lns + lns1
            ylim = a1.get_ylim()
            a1.plot([ZP/2,ZP/2],ylim,'k:')
            a1.plot([ZP,ZP],ylim,'k:')
            a1.plot([3*ZP/2,3*ZP/2],ylim,'k:')
            a1.set_ylim(ylim)
            a1.set_xlim([0,LLnl])
            if var_string[9].get() == 'showtheta':     
                a1bis = a1.twinx()
                theta = np.angle(sol.y[0,:]) + np.angle(sol.y[1,:]) - np.angle(sol.y[2,:])
                theta = (theta + np.pi) % (2 * np.pi) - np.pi
                lns1 = a1bis.plot(sol.t[1:], theta[1:], 'm', label=r'$\theta = \phi_1 + \phi_2 - \phi_3$')
                lns = lns + lns1
                a1bis.set_ylabel(r'$\theta = \phi_1 + \phi_2 - \phi_3$', color='m')
                a1bis.tick_params(axis='y', labelcolor='m') 
                a1bis.set_yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi])
                a1bis.set_yticklabels([r'$-\pi$', r'$-\frac{\pi}{2}$', r'$0$', r'$\frac{\pi}{2}$', r'$\pi$'])
                a1bis.set_ylim([-1.1*np.pi,1.1*np.pi])
            
            a1.set_xlabel(r'$Z = z/L_{\rm nl}$')
            a1.set_ylabel(r'$I_1/I,~I_2/I,~I_3/I$')
            a1.set_title(r'$Z_{\rm P}/2 =$ '+str(round(ZP/2,4)))
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
row = gui.create_entry_with_latex(mainframe,r'$s = \frac{\pi}{2}\textrm{sgn}(\Delta k)L_{\rm nl}/L_{\rm c}  = \frac{\Delta k}{2 \chi \sqrt{\omega_1\omega_2I}}   =$',var_string[0],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_description(mainframe,'propagation range:',row)
row = gui.create_entry_with_latex(mainframe,r'$L / L_{\rm nl} = L \chi \sqrt{\omega_1\omega_2I} = $',var_string[1],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_description(mainframe,'initial conditions:',row)
row = gui.create_entry_with_latex(mainframe,r'$I_{10} / I_{20} = $',var_string[2],row)
row = gui.create_entry_with_latex(mainframe,r'$I_{30} / I_{20} = $',var_string[3],row)
row = gui.create_entry_with_latex(mainframe,r'$\theta_0 = $',var_string[4],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_description(mainframe,'frequency ratio:',row)
row = gui.create_entry_with_latex(mainframe,r'$\omega_2 / \omega_1 = $',var_string[5],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_checkbutton_with_latex(mainframe,r'show $I_1$','noshow','showI1',var_string[6],r'show $I_2$','noshow','showI2',var_string[7],row)
row = gui.create_double_checkbutton_with_latex(mainframe,r'show $I_3$','noshow','showI3',var_string[8],r'show $\theta$','noshow','showtheta',var_string[9],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)