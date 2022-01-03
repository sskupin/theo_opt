import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui

gui.set_rcParams()
root = Tk.Tk()
root.title("SPP")

def n_eff(epsilon_s,epsilon_c): # compute n_eff
    return np.sqrt(epsilon_s*epsilon_c/(epsilon_s+epsilon_c))

def mode_profile(epsilon_s,epsilon_c,n_eff,x): # computing mode profile
    v = np.sqrt(n_eff**2-epsilon_s)*2*np.pi
    w = np.sqrt(n_eff**2-epsilon_c)*2*np.pi
    Fs = np.exp(v*x[(x<0)])
    Fc = np.exp(-w*x[(x>=0)])
    F = np.concatenate((Fs,Fc))
    Gxs = Fs
    Gxc = Fc*epsilon_s/epsilon_c
    Gx = np.concatenate((Gxs,Gxc))    
    Gzs = Fs
    Gzc = -epsilon_s/epsilon_c*w/v*Fc
    Gz = np.concatenate((Gzs,Gzc))      
        
    return F,Gx,Gz

def plot_mode(ax,x,F):
    ax.plot(x,F,'b')
    ax.set_xlim([x[0],x[-1]])
    ylimits = ax.get_ylim()
    ax.set_ylim([ylimits[0],ylimits[1]+0.05*(ylimits[1]-ylimits[0])])
    ax.set_ylim(ax.get_ylim())
    ax.plot([0,0],ax.get_ylim(),'k:')
    ax.set_xlabel(r'$x/\lambda$')
    ax.text(0.2, 0.95, r'substrate', verticalalignment='center', horizontalalignment='center', transform=ax.transAxes)
    ax.text(0.8, 0.95, r'cladding', verticalalignment='center', horizontalalignment='center', transform=ax.transAxes)
    
def initialize():
    global epsilon_s_save,epsilon_c_save
    epsilon_s_string.set("1")
    epsilon_c_string.set("-5")
    
    epsilon_s_save = epsilon_s_string.get()
    epsilon_c_save = epsilon_c_string.get()
    
    calculate()
    
def reinitialize():
    global epsilon_s_save,epsilon_c_save
    epsilon_s_string.set(epsilon_s_save)
    epsilon_c_string.set(epsilon_c_save)

#    calculate()

def calculate():
    global epsilon_s_save,epsilon_c_save
    try:
        epsilon_s = float(epsilon_s_string.get())
        epsilon_c = float(epsilon_c_string.get())

        if epsilon_s < epsilon_c:
            gui.input_error("Substrate epsilon smaller than cladding epsilon. Switching values...")
            epsilon_s = float(epsilon_c_string.get())
            epsilon_c = float(epsilon_s_string.get())
            epsilon_s_string.set(epsilon_s)
            epsilon_c_string.set(epsilon_c)
        
        if epsilon_c > 0 or epsilon_s < 0 or -epsilon_s <= epsilon_c: 
            gui.input_error("No SPP exists for these parameters. Re-initializing with previous parameters...",reinitialize)
        else:
            n_eff_string.set(np.around(n_eff(epsilon_s,epsilon_c),decimals=4))
            f.clf()
            x = np.linspace(-2., 2, num=401, endpoint=True) # x in units of lambda
            F,Gx,Gz = mode_profile(epsilon_s,epsilon_c,n_eff(epsilon_s,epsilon_c),x)
            a1 = f.add_subplot(131)
            plot_mode(a1,x,F)
            a2 = f.add_subplot(132)
            plot_mode(a2,x,Gx)
            a3 = f.add_subplot(133)
            plot_mode(a3,x,Gz)
            a1.set_ylabel(r'$H_y/H_y(x=0)$')
            a2.set_ylabel(r'$E_x/E_x(x=0)$')
            a3.set_ylabel(r'$E_z/E_z(x=0)$')                
            plt.tight_layout()
            
#            plt.savefig('spp.pdf',bbox_inches='tight',dpi=300, transparent=True)

            epsilon_s_save = epsilon_s_string.get()
            epsilon_c_save = epsilon_c_string.get()

            canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing with previous parameters...",reinitialize)

f = plt.figure(1,[10,2.5])
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

epsilon_s_string = Tk.StringVar()
epsilon_c_string = Tk.StringVar()
n_eff_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_entry_with_latex(mainframe,r'substrate $\varepsilon_{\rm s} =$',epsilon_s_string,row)
row = gui.create_entry_with_latex(mainframe,r'cladding $\varepsilon_{\rm c} =$',epsilon_c_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_label_with_latex(mainframe,r'eff.\ index $n_{\rm eff}=$',n_eff_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)