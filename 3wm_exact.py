import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import scipy.integrate as spi
import scipy.special as sps
import tkinter as Tk
import sys
sys.path.append('./stuff')
import gui_stuff as gui

gui.set_rcParams()
title = "Three-Wave Mixing - Exact Numerical Solution"
root = Tk.Tk()
root.title(title)

def initialize():
    var_string[0].set(".36") # s
    var_string[1].set("10") # \kappa L
    var_string[2].set("0") # rho_1/rho_2
    var_string[3].set("9") # rho_3/rho_2
    var_string[4].set("1") # theta
    var_string[5].set("0.66") # omega_2/omega_1
    var_string[6].set("showI1") # show I_1
    var_string[7].set("showI2") # show I_2
    var_string[8].set("showI3") # show I_3
    var_string[9].set("noshow") # show theta
    var_string[10].set("noTypeI") # switch to Type I SHG
    gui.copy_stringvar_vector(var_string,var_save)    
    calculate()
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    calculate()
    
def show_manual():
    gui.show_manual("man/3wm_exact.png",title) 

def calculate():
    gui.change_cursor(root,"trek")
    try:
        TypeI = False
        if var_string[10].get() == 'TypeI':
            TypeI = True
            var_string[2].set("1")
            var_string[5].set("1")
            var_string[6].set("noshow")
        s = float(var_string[0].get())
        LLnl = float(var_string[1].get())
        I1I2 = float(var_string[2].get())
        I3I2 = float(var_string[3].get())
        theta = float(var_string[4].get())
        omega2omega1 = float(var_string[5].get())
        
        if LLnl == 0:
            gui.input_error("Nonlinear interaction strength must not be zero. Re-initializing with previous parameters...",reinitialize)
        elif np.abs(LLnl) > 20:
            gui.input_error("Nonlinear interaction strength too large. Re-initializing with previous parameters...",reinitialize)
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
            U = np.real(np.sort(np.roots([1, - C2 - C3 - s**2, C2 * C3 - 2 * s * H, - H**2])))
            m = (U[1] - U[0]) / (U[2] - U[0])
            Km = sps.ellipk(m)
            ZP = 4 * Km / np.sqrt(U[2] - U[0])
            
            if np.abs(LLnl) > 2*np.abs(ZP):
                gui.input_error("Propagation range too large. Reducing ...")
                LLnl = 2*ZP
                var_string[1].set(LLnl)
        
            def compute_rhs(Z,A): # computes rhs of ode system, A[j-1] = A_j
                rhs1 = 1j * A[2] * np.conj(A[1])
                rhs2 = 1j * A[2] * np.conj(A[0])
                rhs3 = 1j * (A[0] * A[1] - 2 * s * A[2])
                return np.array([rhs1, rhs2, rhs3])
        
            sol = spi.solve_ivp(compute_rhs, [0, LLnl], A, max_step = 1.e-3*np.abs(LLnl))

            f.clf() 
        
            a1 = f.add_subplot(gs[5:, 0:])
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
                lns1 = a1.plot(sol.t, np.abs(sol.y[0,:])**2 / omega3omega1, 'r', label=r'$I_1/I_0$')
                lns = lns + lns1
            if var_string[7].get() == 'showI2':
                if TypeI:
                    lns1 = a1.plot(sol.t, 2*np.abs(sol.y[1,:])**2 / omega3omega2, 'g', label=r'$I_{\omega}/I_0$')
                else:
                    lns1 = a1.plot(sol.t, np.abs(sol.y[1,:])**2 / omega3omega2, 'g', label=r'$I_2/I_0$')
                lns = lns + lns1                
            if var_string[8].get() == 'showI3':  
                if TypeI:
                    lns1 = a1.plot(sol.t, np.abs(sol.y[2,:])**2, 'b', label=r'$I_{2\omega}/I_0$')
                else:
                    lns1 = a1.plot(sol.t, np.abs(sol.y[2,:])**2, 'b', label=r'$I_3/I_0$')
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
                if TypeI:
                    lns1 = a1bis.plot(sol.t[1:], theta[1:], 'm', label=r'$\theta = 2\varphi_\omega - \varphi_{2\omega}$')
                    a1bis.set_ylabel(r'$\theta = 2\varphi_\omega - \varphi_{2\omega}$', color='m')
                else:
                    lns1 = a1bis.plot(sol.t[1:], theta[1:], 'm', label=r'$\theta = \varphi_1 + \varphi_2 - \varphi_3$')
                    a1bis.set_ylabel(r'$\theta = \varphi_1 + \varphi_2 - \varphi_3$', color='m')
                lns = lns + lns1
                a1bis.tick_params(axis='y', labelcolor='m') 
                a1bis.set_yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi])
                a1bis.set_yticklabels([r'$-\pi$', r'$-\frac{\pi}{2}$', r'$0$', r'$\frac{\pi}{2}$', r'$\pi$'])
                a1bis.set_ylim([-1.1*np.pi,1.1*np.pi])
            
            a1.set_xlabel(r'$Z = \kappa z$')
            if TypeI:
                a1.set_ylabel(r'$I_\omega/I_0,~I_{2\omega}/I_0$')
            else:
                a1.set_ylabel(r'$I_1/I_0,~I_2/I_0,~I_3/I_0$')
            a1.set_title(r'$Z_{\rm P}/2 =$ '+str(round(ZP/2,4)))
            labs = [l.get_label() for l in lns]
            a1.legend(lns, labs, bbox_to_anchor=(0., 1.12, .5, 0), loc="lower left", mode="expand", ncol=1)
            
            def V(X):
                return -2 * (X*(C2-X)*(C3-X) - (H+s*X)**2)
            
            a2 = f.add_subplot(gs[0:3, 2]) # plot potential
            if U[1] > U[0]:
                X = np.linspace(U[0]-0.1*(U[1]-U[0]),U[1]+0.1*(U[1]-U[0]), num=100, endpoint=True)
            elif U[2] > U[0] and V((U[2]-U[0])/2) < 0:
                X = np.linspace(U[0]-0.1*(U[2]-U[0]),U[2]+0.1*(U[2]-U[0]), num=100, endpoint=True)
            else:
                X = np.linspace(U[0]-0.1*(U[2]-U[0]),U[0]+0.1*(U[2]-U[0]), num=100, endpoint=True)
            a2.plot(X,V(X),'k')
            a2.plot(X,0*X,'k:')
            a2.autoscale(enable=True, axis='x', tight=True)
            if TypeI:
                a2.set_xlabel(r'$U=I_{2\omega}/I_0$')
            else:
                a2.set_xlabel(r'$U=I_3/I_0$')
            a2.set_ylabel(r'$V(U)$')
        
#        plt.savefig('3wm.pdf',bbox_inches='tight',dpi=300, transparent=True)
        
        gui.copy_stringvar_vector(var_string,var_save)

        canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[7,5])
gs = GridSpec(16, 3, figure=f)
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(11)
var_save = gui.create_stringvar_vector(11)

initialize()

row = 1
row = gui.create_description(mainframe,'phase mismatch:',row)
row = gui.create_entry_with_latex(mainframe,r'$s = \frac{\pi}{2}\textrm{sgn}(\Delta k)/(\kappa L_{\rm c})  = \frac{\Delta k}{2 \chi \sqrt{\omega_1\omega_2I_0}}   =$',var_string[0],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_description(mainframe,'nonlinear interaction strength:',row)
row = gui.create_entry_with_latex(mainframe,r'$\kappa L = L \chi \sqrt{\omega_1\omega_2I_0} = $',var_string[1],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_description(mainframe,'initial conditions:',row)
row = gui.create_entry_with_latex(mainframe,r'$I_{10} / I_{20} = $',var_string[2],row)
row = gui.create_entry_with_latex(mainframe,r'$I_{30} / I_{20} = $',var_string[3],row)
row = gui.create_entry_with_latex(mainframe,r'(meaningful only if $I_{10} \ne 0$ and $I_{30} \ne 0$) $\theta_0 = $',var_string[4],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_description(mainframe,'frequency ratio:',row)
row = gui.create_entry_with_latex(mainframe,r'$\omega_2 / \omega_1 = $',var_string[5],row)
row = gui.create_checkbutton_with_latex(mainframe,r'Type I SHG with $I_{\omega} = 2I_2$ and $I_{2\omega} = I_3$','noTypeI','TypeI',var_string[10],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_checkbutton_with_latex(mainframe,r'show $I_1$','noshow','showI1',var_string[6],r'show $I_2$','noshow','showI2',var_string[7],row)
row = gui.create_double_checkbutton_with_latex(mainframe,r'show $I_3$','noshow','showI3',var_string[8],r'show $\theta$','noshow','showtheta',var_string[9],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)