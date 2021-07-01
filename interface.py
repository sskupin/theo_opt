import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import tkinter as Tk
import gui_stuff as gui

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rc('text', usetex=True)
mpl.rc('text.latex', preamble=r'\usepackage{cmbright}')
mpl.rcParams.update({'font.size': 10})

root = Tk.Tk()
root.title("Reflection and Transmission at Interface")

def reflection_transmission(epsilon_s,epsilon_c,phi): # computing coefficients of reflection and transmission
    ksx = np.sqrt(epsilon_s)*np.sin(phi)
    ksz = np.sqrt(epsilon_s-ksx**2)
    kcz = np.sqrt(epsilon_c-ksx**2)
    RTE = (ksz-kcz)/(ksz+kcz)
    RTM = (epsilon_s*kcz-epsilon_c*ksz)/(epsilon_c*ksz+epsilon_s*kcz) # for electric field (negative of magnetic coeff.)
    tauTE = np.real(kcz)/ksz*np.abs(2*ksz/(ksz+kcz))**2
    tauTM = np.real(kcz/epsilon_c)*epsilon_s/ksz*(np.abs(2*epsilon_c*ksz/(epsilon_c*ksz+epsilon_s*kcz)))**2
        
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
    global epsilon_s_save,epsilon_c_real_save,epsilon_c_imag_save
    epsilon_s_string.set("1")
    epsilon_c_real_string.set("15.1")
    epsilon_c_imag_string.set("0")
    
    epsilon_s_save = epsilon_s_string.get()
    epsilon_c_real_save = epsilon_c_real_string.get()
    epsilon_c_imag_save = epsilon_c_imag_string.get()
        
    calculate()
    
def reinitialize():
    global epsilon_s_save,epsilon_c_real_save,epsilon_c_imag_save
    epsilon_s_string.set(epsilon_s_save)
    epsilon_c_real_string.set(epsilon_c_real_save)
    epsilon_c_imag_string.set(epsilon_c_imag_save)
        
def calculate():
    global epsilon_s_save,epsilon_c_real_save,epsilon_c_imag_save
    try:
        epsilon_s = float(epsilon_s_string.get())
        epsilon_c_real = float(epsilon_c_real_string.get())
        epsilon_c_imag = float(epsilon_c_imag_string.get())
        
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

            epsilon_s_save = epsilon_s_string.get()
            epsilon_c_real_save = epsilon_c_real_string.get()
            epsilon_c_imag_save = epsilon_c_imag_string.get()

            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing ...", initialize)

f = plt.figure(1,[7,2])
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

epsilon_s_string = Tk.StringVar()
epsilon_c_real_string = Tk.StringVar()
epsilon_c_imag_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_entry_with_latex(mainframe,r"substrate $\varepsilon_{\rm s} =$",epsilon_s_string,row)
row = gui.create_entry_with_latex(mainframe,r"cladding $\varepsilon_{\rm c}' =$",epsilon_c_real_string,row)
row = gui.create_entry_with_latex(mainframe,r"cladding $\varepsilon_{\rm c}'' =$",epsilon_c_imag_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)