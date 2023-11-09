import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import strat_stuff as strat
                   
gui.set_rcParams()
root = Tk.Tk()
root.title("Reflection and Transmission at Interface")

def reflection_transmission(epsilon_s,epsilon_c,phi): # computing coefficients of reflection and transmission
    kx,ksz,kcz = strat.KSC(epsilon_s,epsilon_c,phi)
    RTE,RTM,TTE,TTM,tauTE,tauTM = strat.RTAU(ksz,kcz,epsilon_s,epsilon_c,np.identity(2),np.identity(2))
        
    return RTE,RTM,tauTE,tauTM

def plot_subplot(ax,phi,curve_1,curve_2,label_1,label_2):
    ax.plot(phi,curve_1,'b',label=label_1)
    ax.plot(phi,curve_2,'r',label=label_2)
    ax.set_xticks([0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2])
    ax.set_xticklabels([r'$0$', r'$\pi/8$', r'$\pi/4$', r'$3\pi/8$', r'$\pi/2$'])
    ax.set_xlabel(r'$\varphi_{\rm i}$')
    ax.set_xlim([0,np.pi/2])
    ax.set_ylabel(label_1+','+label_2)
    ax.legend()
    
def ourangle(z): # angle of pi is replaced by -pi
    ourangle = np.angle(z)
    if ourangle == np.pi:
        ourangle = -np.pi
    return ourangle
        
def initialize():
    var_string[0].set("1") # epsilon_s
    var_string[1].set("15.1") # epsilon_c_real
    var_string[2].set("0") # epsilon_c_imag
    gui.copy_stringvar_vector(var_string,var_save)
    calculate()    
    
def reinitialize():
    gui.copy_stringvar_vector(var_save,var_string)
        
def calculate():
    gui.change_cursor(root,"trek")
    try:
        epsilon_s = float(var_string[0].get())
        epsilon_c_real = float(var_string[1].get())
        epsilon_c_imag = float(var_string[2].get())
        
        if epsilon_s <= 0:
            gui.input_error("Substrate epsilon must be positive. Re-initializing with previous parameters...",reinitialize)
        elif epsilon_c_real == 0 and epsilon_c_imag == 0: 
            gui.input_error("Cladding epsilon must not be negative. Re-initializing with previous parameters...",reinitialize)
        else:
            f.clf()
            phi = np.linspace(0, np.pi/2, num=401, endpoint=False) # angle of incidence
            epsilon_c = epsilon_c_real + 1j*epsilon_c_imag
            RTE,RTM,tauTE,tauTM = reflection_transmission(epsilon_s,epsilon_c,phi)
            a1 = f.add_subplot(131)
            plot_subplot(a1,phi,np.abs(RTE)**2,tauTE,r'$\rho_{\rm TE}$',r'$\tau_{\rm TE}$')
            a1.set_ylim([-0.025,1.025])
            a2 = f.add_subplot(132)
            plot_subplot(a2,phi,np.abs(RTM)**2,tauTM,r'$\rho_{\rm TM}$',r'$\tau_{\rm TM}$')
            a2.set_ylim([-0.025,1.025])
            a3 = f.add_subplot(133)
            vourangle = np.vectorize(ourangle)
            if epsilon_c_real > epsilon_s:
                plot_subplot(a3,phi,vourangle(RTE),vourangle(RTM),r'$\theta_{\rm TE}$',r'$\theta_{\rm TM}$')  
            else:
                plot_subplot(a3,phi,np.angle(RTE),np.angle(RTM),r'$\theta_{\rm TE}$',r'$\theta_{\rm TM}$')  
            a3.set_yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi])
            a3.set_yticklabels([r'$-\pi$', r'$-\frac{\pi}{2}$', r'$0$', r'$\frac{\pi}{2}$', r'$\pi$'])
            a3.set_ylim([-1.05*np.pi,1.05*np.pi])
            
            if epsilon_c_imag <= 0.01*np.abs(epsilon_c_real) and epsilon_c_real > 0:
                phiB = np.arctan(np.sqrt(epsilon_c_real/epsilon_s))
                a2.plot([phiB,phiB],[-1.05,1.05],'k:')
                a2.text(2*phiB/np.pi, 1.03, r'$\varphi_{\rm iB}$', verticalalignment='bottom', horizontalalignment='right', transform=a2.transAxes)                          
                a3.plot([phiB,phiB],[-1.05*np.pi,1.05*np.pi],'k:')
                a3.text(2*phiB/np.pi, 1.03, r'$\varphi_{\rm iB}$', verticalalignment='bottom', horizontalalignment='right', transform=a3.transAxes)           
                if epsilon_c_real < epsilon_s:
                    phiC = np.arcsin(np.sqrt(epsilon_c_real/epsilon_s))
                    a1.plot([phiC,phiC],[-1.05,1.05],'k:')
                    a1.text(2*phiC/np.pi, 1.03, r'$\varphi_{\rm iC}$', verticalalignment='bottom', horizontalalignment='left', transform=a1.transAxes)                          
                    a2.plot([phiC,phiC],[-1.05,1.05],'k:')
                    a2.text(2*phiC/np.pi, 1.03, r'$\varphi_{\rm iC}$', verticalalignment='bottom', horizontalalignment='left', transform=a2.transAxes)                          
                    a3.plot([phiC,phiC],[-1.05*np.pi,1.05*np.pi],'k:')
                    a3.text(2*phiC/np.pi, 1.03, r'$\varphi_{\rm iC}$', verticalalignment='bottom', horizontalalignment='left', transform=a3.transAxes)           
            
            plt.tight_layout()
            
#            plt.savefig('interface.pdf',bbox_inches='tight',dpi=300, transparent=True)

            gui.copy_stringvar_vector(var_string,var_save)

            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", reinitialize)
    gui.change_cursor(root,"arrow")

f = plt.figure(1,[7,2])
canvas = gui.create_canvas(root,f)
canvas.draw() # for faster feedback to user on startup
mainframe = gui.create_mainframe(root)

var_string = gui.create_stringvar_vector(3)
var_save = gui.create_stringvar_vector(3)

initialize()

row = 1
row = gui.create_entry_with_latex(mainframe,r"substrate $\varepsilon_{\rm s} =$",var_string[0],row)
row = gui.create_entry_with_latex(mainframe,r"cladding $\varepsilon_{\rm c}' =$",var_string[1],row)
row = gui.create_entry_with_latex(mainframe,r"cladding $\varepsilon_{\rm c}'' =$",var_string[2],row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)