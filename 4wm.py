import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import scipy.integrate as spi
import scipy.special as sps
import tkinter as Tk
import sys
sys.path.append('./aux')
import gui_stuff as gui

gui.set_rcParams()
title = "Four-Wave Mixing - Exact Numerical Solution"
root = Tk.Tk()
root.title(title)

def initialize():
    var_string[0].set(".2") # s
    var_string[1].set("15") # L/L_nl
    var_string[2].set("1.4") # rho_1/rho_2
    var_string[3].set("0.01") # rho_3/rho_2
    var_string[4].set("0") # rho_4/rho_2    
    var_string[5].set("0") # theta
    var_string[6].set("0.67") # omega_2/omega_1
    var_string[7].set("0.65") # omega_3/omega_1    
    var_string[8].set("showP1") # show P_1
    var_string[9].set("showP2") # show P_2
    var_string[10].set("showP3") # show P_3
    var_string[11].set("showP4") # show P_4
    var_string[12].set("noshow") # show theta
    var_string[13].set("-0.8") # Q_1
    var_string[14].set("-0.4") # Q_2
    var_string[15].set("0.4") # Q_3   
    var_string[16].set("0.8") # Q_4
    var_string[17].set("nodeg") # 
    gui.copy_stringvar_vector(var_string,var_save)    
    calculate()
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
    calculate()
    
def show_manual():
    gui.show_manual("man/4wm.png",title) 

def calculate():
    gui.change_cursor(root,"trek")
    try:
        deg = False
        if var_string[17].get() == 'deg':
            deg = True
            var_string[2].set("1")
            var_string[6].set("1")
            var_string[8].set("noshow")
        s = float(var_string[0].get())
        LLnl = float(var_string[1].get())
        P1P2 = float(var_string[2].get())
        P3P2 = float(var_string[3].get())
        P4P2 = float(var_string[4].get())
        theta = float(var_string[5].get())
        omega2omega1 = float(var_string[6].get())
        omega3omega1 = float(var_string[7].get())
        Q1 = float(var_string[13].get())
        Q2 = float(var_string[14].get())
        Q3 = float(var_string[15].get())
        Q4 = float(var_string[16].get())
        
        if LLnl == 0:
            gui.input_error("Nonlinear interaction strength must not be zero. Re-initializing with previous parameters...",reinitialize)
        elif np.abs(LLnl) > 20:
            gui.input_error("Nonlinear interaction strength too large. Re-initializing with previous parameters...",reinitialize)
        elif P1P2 < 0 or P3P2 < 0 or P4P2 < 0: 
            gui.input_error("Power ratios must not be negative. Re-initializing with previous parameters...",reinitialize)
        elif omega2omega1 <= 0 or omega3omega1 <= 0:
            gui.input_error("Frequency ratio must be positive. Re-initializing with previous parameters...",reinitialize)
        else: 
                
            omega4omega1 = 1 + omega2omega1 - omega3omega1
            omega4omega2 = omega4omega1/omega2omega1
            omega4omega3 = omega4omega1/omega3omega1
        if omega4omega3 <= 0:
            gui.input_error("4th frequency must be positive. Re-initializing with previous parameters...",reinitialize)
        else:
            omega3omega2 = omega3omega1/omega2omega1
            rho2 = np.sqrt(omega4omega2/(1 + P1P2 + P3P2 + P4P2))
            rho1 = np.sqrt(omega2omega1 * P1P2) * rho2
            rho3 = np.sqrt(P3P2 / omega3omega2) * rho2
            rho4 = np.sqrt(P4P2 / omega4omega2) * rho2
            
            A = np.zeros(4) + 0j 
            A[0] = rho1
            A[1] = rho2 * np.exp(1j * theta)
            A[2] = rho3
            A[3] = rho4
            
            H = rho1 * rho2 * rho3 * rho4 * np.cos(theta) - s * rho4**2 + (Q1*rho1**4+Q2*rho2**4-Q3*rho3**4-Q4*rho4**4)/4
            C1 = rho1**2 + rho4**2
            C2 = rho2**2 + rho4**2
            C3 = rho3**2 - rho4**2
            
            C02 = 1 - (Q1+Q2-Q3-Q4)**2/16
            if C02 == 0: # here we cheat a bit, to avoid analytical expressions for ZP and C02 = 0
                Q4 = Q4 + 1.e-3
                C02 = 1 - (Q1+Q2-Q3-Q4)**2/16
                H = rho1 * rho2 * rho3 * rho4 * np.cos(theta) - s * rho4**2 + (Q1*rho1**4+Q2*rho2**4-Q3*rho3**4-Q4*rho4**4)/4
            PF3 = -C1-C2+C3-(2*(s+(1/2)*Q1*C1+(1/2)*Q2*C2+(1/2)*Q3*C3))*(-(1/4)*Q1-(1/4)*Q2+(1/4)*Q3+(1/4)*Q4)
            PF2 = C1*C2+(-C1-C2)*C3-(2*(H-(1/4)*Q1*C1**2-(1/4)*Q2*C2**2+(1/4)*Q3*C3**2))*(-(1/4)*Q1-(1/4)*Q2+(1/4)*Q3+(1/4)*Q4)-(s+(1/2)*Q1*C1+(1/2)*Q2*C2+(1/2)*Q3*C3)**2
            PF1 = C1*C2*C3-(2*(H-(1/4)*Q1*C1**2-(1/4)*Q2*C2**2+(1/4)*Q3*C3**2))*(s+(1/2)*Q1*C1+(1/2)*Q2*C2+(1/2)*Q3*C3)
            PF0 = -(H-(1/4)*Q1*C1**2-(1/4)*Q2*C2**2+(1/4)*Q3*C3**2)**2
            U = np.real(np.sort(np.roots([C02, PF3, PF2, PF1, PF0])))
            if C02 > 0:
                m = (U[2] - U[1]) * (U[3] - U[0]) / ( (U[2] - U[0]) * (U[3] - U[1]) )
                Km = sps.ellipk(m)
                ZP = 4 * Km / np.sqrt(C02 * (U[2] - U[0]) * (U[3] - U[1]))
            elif C02 < 0:
                m = (U[3] - U[2]) * (U[1] - U[0]) / ( (U[2] - U[0]) * (U[3] - U[1]) )
                Km = sps.ellipk(m)
                ZP = 4 * Km / np.sqrt(- C02 * (U[2] - U[0]) * (U[3] - U[1]))
            
            if LLnl > 2*ZP:
                gui.input_error("Propagation range too large. Reducing ...")
                LLnl = 2*ZP
                var_string[1].set(LLnl)
        
            def compute_rhs(Z,A): # computes rhs of ode system, A[j-1] = A_j
                rhs1 = 1j * A[3] * A[2] * np.conj(A[1])
                rhs2 = 1j * A[3] * A[2] * np.conj(A[0]) + 1j*(Q1*np.abs(A[0])**2+Q2*np.abs(A[1])**2+Q3*np.abs(A[2])**2+Q4*np.abs(A[3])**2)*A[1] # put all Q_sigma here, only theta counts
                rhs3 = 1j * A[1] * A[0] * np.conj(A[3])
                rhs4 = 1j * (A[0] * A[1] * np.conj(A[2]) - 2 * s * A[3])
                return np.array([rhs1, rhs2, rhs3, rhs4])
        
            sol = spi.solve_ivp(compute_rhs, [0, LLnl], A, max_step = 1.e-3*np.abs(LLnl))

            f.clf() 
        
            a1 = f.add_subplot(gs[5:, 0:])
            lns = []
            if var_string[8].get() == 'showP1':
                lns1 = a1.plot(sol.t, np.abs(sol.y[0,:])**2 / omega4omega1, 'r', label=r'$P_1/P_0$')
                lns = lns + lns1
            if var_string[9].get() == 'showP2':
                if deg:
                    lns1 = a1.plot(sol.t, 2*np.abs(sol.y[1,:])**2 / omega4omega2, 'r', label=r'$P_{\rm P}/P_0$')
                else:
                    lns1 = a1.plot(sol.t, np.abs(sol.y[1,:])**2 / omega4omega2, 'g', label=r'$P_2/P_0$')
                lns = lns + lns1                
            if var_string[10].get() == 'showP3':  
                lns1 = a1.plot(sol.t, np.abs(sol.y[2,:])**2 / omega4omega3, 'y', label=r'$P_3/P_0$')
                lns = lns + lns1
            if var_string[11].get() == 'showP4':  
                lns1 = a1.plot(sol.t, np.abs(sol.y[3,:])**2, 'b', label=r'$P_4/P_0$')
                lns = lns + lns1
            ylim = a1.get_ylim()
            a1.plot([ZP/2,ZP/2],ylim,'k:')
            a1.plot([ZP,ZP],ylim,'k:')
            a1.plot([3*ZP/2,3*ZP/2],ylim,'k:')
            a1.set_ylim(ylim)
            a1.set_xlim([0,LLnl])
            if var_string[12].get() == 'showtheta':     
                a1bis = a1.twinx()
                theta = np.angle(sol.y[0,:]) + np.angle(sol.y[1,:]) - np.angle(sol.y[2,:]) - np.angle(sol.y[3,:])
                theta = (theta + np.pi) % (2 * np.pi) - np.pi
                if deg:
                    lns1 = a1bis.plot(sol.t[1:], theta[1:], 'm', label=r'$\theta = 2\varphi_{\rm P} - \varphi_3 - \varphi_4$')
                    a1bis.set_ylabel(r'$\theta = 2\varphi_{\rm P} - \varphi_3 - \varphi_4$', color='m')                   
                else:
                    lns1 = a1bis.plot(sol.t[1:], theta[1:], 'm', label=r'$\theta = \varphi_1 + \varphi_2 - \varphi_3 - \varphi_4$')
                    a1bis.set_ylabel(r'$\theta = \varphi_1 + \varphi_2 - \varphi_3 - \varphi_4$', color='m')
                lns = lns + lns1
                a1bis.tick_params(axis='y', labelcolor='m') 
                a1bis.set_yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi])
                a1bis.set_yticklabels([r'$-\pi$', r'$-\frac{\pi}{2}$', r'$0$', r'$\frac{\pi}{2}$', r'$\pi$'])
                a1bis.set_ylim([-1.1*np.pi,1.1*np.pi])
            
            a1.set_xlabel(r'$Z = \kappa z$')
            if deg:
                a1.set_ylabel(r'$P_{\rm P}/P_0,~P_3/P_0,~P_4/P_0$')
            else:
                a1.set_ylabel(r'$P_1/P_0,~P_2/P_0,~P_3/P_0,~P_4/P_0$')
            a1.set_title(r'$Z_{\rm P}/2 =$ '+str(round(ZP/2,4)))
            labs = [l.get_label() for l in lns]
            a1.legend(lns, labs, bbox_to_anchor=(0., 1.1, .5, 0), loc="lower left", mode="expand", ncol=1)
            
            def V(X):
                return -2 * (X*(C1-X)*(C2-X)*(C3+X) - (H+s*X-Q1*(C1-X)**2/4-Q2*(C2-X)**2/4+Q3*(C3+X)**2/4+Q4*X**2/4)**2)
            
            a2 = f.add_subplot(gs[0:3, 2]) # plot potential
            if U[1] > U[0] and V((U[1]+U[0])/2) < 0 and (U[1]+U[0])/2 > 0:
                minindex = 0
            elif U[2] > U[1] and V((U[2]+U[1])/2) < 0 and (U[2]+U[1])/2 > 0:
                minindex = 1
            elif U[3] > U[2] and V((U[3]+U[2])/2) < 0 and (U[3]+U[2])/2 > 0:
                minindex = 2
            else:
                minindex = 3
            if U[3] > U[2] and V((U[3]+U[2])/2) < 0 and (U[3]+U[2])/2 < 1:
                maxindex = 3
            elif U[2] > U[1] and V((U[2]+U[1])/2) < 0 and (U[2]+U[1])/2 < 1:
                maxindex=2
            elif U[1] > U[0] and V((U[1]+U[0])/2) < 0 and (U[1]+U[0])/2 < 1:
                maxindex = 1
            else:
                maxindex = 0
            if  maxindex > minindex:
                X = np.linspace(U[minindex]-0.1*(U[maxindex]-U[minindex]),U[maxindex]+0.1*(U[maxindex]-U[minindex]), num=100, endpoint=True)
            else:
                X = np.linspace(U[minindex]-0.1,U[minindex]+0.1, num=100, endpoint=True)
            a2.plot(X,V(X),'k')
            a2.plot(X,0*X,'k:')
            a2.autoscale(enable=True, axis='x', tight=True)
            a2.set_xlabel(r'$U=P_4/P_0$')
            a2.set_ylabel(r'$V(U)$')
        
#        plt.savefig('4wm.pdf',bbox_inches='tight',dpi=300, transparent=True)
        
        gui.copy_stringvar_vector(var_string,var_save)
        gui.change_cursor(root,"arrow")

        canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)

f = plt.figure(1,[7,5])
gs = GridSpec(16, 3, figure=f)
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(18)
var_save = gui.create_stringvar_vector(18)

initialize()

row = 1
row = gui.create_description(mainframe,'phase mismatch:',row)
row = gui.create_entry_with_latex(mainframe,r'$s = \frac{\Delta k}{2 \kappa} = \frac{\Delta k\sqrt{\omega_4}}{2\chi^{\rm F}P\sqrt{\omega_1\omega_2\omega_3}}= $',var_string[0],row)
row = gui.create_description(mainframe,'nonlinear interaction strength:',row)
row = gui.create_entry_with_latex(mainframe,r'$\kappa L = L \chi^{\rm F}P_0 \sqrt{\frac{\omega_1\omega_2\omega_3}{\omega_4}}= $',var_string[1],row)
row = gui.create_description(mainframe,'initial conditions:',row)
row = gui.create_double_entry_with_latex(mainframe,r'$P_{10} / P_{20} = $',var_string[2],r'$P_{30} / P_{20} = $',var_string[3],row)
row = gui.create_double_entry_with_latex(mainframe,r'$P_{40} / P_{20} = $',var_string[4],r'$\theta_0 = $',var_string[5],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_description(mainframe,'frequency ratios:',row)
row = gui.create_double_entry_with_latex(mainframe,r'$\omega_2 / \omega_1 = $',var_string[6],r'$\omega_3 / \omega_1 = $',var_string[7],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_description(mainframe,'phase modulation coefficients:',row)
row = gui.create_double_entry_with_latex(mainframe,r'$Q_1 = $',var_string[13],r'$Q_2 = $',var_string[14],row)
row = gui.create_double_entry_with_latex(mainframe,r'$Q_3 = $',var_string[15],r'$Q_4 = $',var_string[16],row)
row = gui.create_checkbutton_with_latex(mainframe,r'Degenerate FWM with $P_{\rm P} = 2P_2$','nodeg','deg',var_string[17],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_checkbutton_with_latex(mainframe,r'show $P_1$','noshow','showP1',var_string[8],r'show $P_2$','noshow','showP2',var_string[9],row)
row = gui.create_double_checkbutton_with_latex(mainframe,r'show $P_3$','noshow','showP3',var_string[10],r'show $P_4$','noshow','showP4',var_string[11],row)
row = gui.create_checkbutton_with_latex(mainframe,r'show $\theta$','noshow','showtheta',var_string[12],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_double_button(mainframe,"Manual",show_manual,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)