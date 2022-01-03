import numpy as np
import matplotlib.pyplot as plt
import tkinter as Tk
import gui_stuff as gui
import scipy.special as sps
import scipy.optimize as spo

gui.set_rcParams()
root = Tk.Tk()
root.title("Fiber in weak guiding approximation")

def number_of_modes(V,m): # works for up to 10 modes
    return np.sum(sps.jn_zeros(m-1, 11) < V) + 1 - np.sign(m)
    
def initialize():
    global epsilon_co_save,epsilon_c_save,a_string_save,m_string_save,mu_string_save
    epsilon_co_string.set("2.19")
    epsilon_c_string.set("2.085")
    a_string.set("5.")
    m_string.set("1")
    mu_string.set("1")
    
    epsilon_co_save = epsilon_co_string.get()
    epsilon_c_save = epsilon_c_string.get()
    a_string_save = a_string.get()
    m_string_save = m_string.get()
    mu_string_save = mu_string.get()
    
    calculate()
    
def reinitialize():
    global epsilon_co_save,epsilon_c_save,a_string_save,m_string_save,mu_string_save
    epsilon_co_string.set(epsilon_co_save)
    epsilon_c_string.set(epsilon_c_save)
    a_string.set(a_string_save)
    m_string.set(m_string_save)
    mu_string.set(mu_string_save)
    
def calculate():
    global epsilon_co_save,epsilon_c_save,a_string_save,m_string_save,mu_string_save
    try:
        epsilon_co = float(epsilon_co_string.get())
        epsilon_c = float(epsilon_c_string.get())
        a = float(a_string.get())
        m = int(np.abs(float(m_string.get())))
        m_string.set(m)
        mu = int(float(mu_string.get()))
        mu_string.set(mu)
        
        if epsilon_co < epsilon_c:
            gui.input_error("Substrate epsilon smaller than cladding epsilon. Switching values...")
            epsilon_co = float(epsilon_c_string.get())
            epsilon_c = float(epsilon_co_string.get())
            epsilon_co_string.set(epsilon_co)
            epsilon_c_string.set(epsilon_c)
        if mu < 0:
            gui.input_error("Mode index negative. Re-initializing with previous parameters...",reinitialize)
        elif a <= 0:
            gui.input_error("Core radius negative. Re-initializing with previous parameters...",reinitialize)
        elif epsilon_co <= 0 or epsilon_c <= 0: 
            gui.input_error("All dielectric constants have to be positive. Re-initializing with previous parameters...",reinitialize)
        else:
            V0 = 2*np.pi*a*np.sqrt(epsilon_co-epsilon_c)
            Vmax = 2*V0
            if number_of_modes(Vmax,m) <= 0:
                gui.input_error("No mode found in the plotting window. Re-initializing with previous parameters...",reinitialize)
            elif number_of_modes(Vmax,m) > 10:
                gui.input_error("Too many modes found for readable plotting. Re-initializing with previous parameters...",reinitialize)
            else:
                besselzero = sps.jn_zeros(m-1, 10)
                Vdiscr = np.sort(np.concatenate((np.linspace(0, Vmax, num=500, endpoint=True),besselzero[besselzero < Vmax] + 1.e-12)));
                Btable = np.zeros([number_of_modes(Vmax,m),Vdiscr.size]);
                Btable[:,0] = np.nan;
                if m==0: # no cut-off for fundamental mode
                    Btable[0,0] = 0
                    
                for indexV in range(1,Vdiscr.size):
                    def disp(B): # dispersion relation in the form disp=0
                        return np.sqrt(1-B)*sps.jv(m-1,Vdiscr[indexV]*np.sqrt(1-B))*sps.kv(m,Vdiscr[indexV]*np.sqrt(B))+np.sqrt(B)*sps.jv(m,Vdiscr[indexV]*np.sqrt(1-B))*sps.kv(m-1,Vdiscr[indexV]*np.sqrt(B))
                    for indexM in range(number_of_modes(Vdiscr[indexV],m)-1,-1,-1):
                        if indexM == number_of_modes(Vdiscr[indexV],m)-1: # check if we are dealing with highest mode (lowest B)
                            lb=1.e-12
                        else:
                            lb=Btable[indexM+1,indexV]+1.e-12
                        if indexM == 0: # check if we are dealing with groundstate
                            ub = 1-1.e-12
                        else:
                            ub = Btable[indexM-1,indexV-1]
                        if disp(lb)*disp(ub) > 0:
                            Btable[indexM,indexV] = 0
                        else:
                            Btable[indexM,indexV] = spo.root_scalar(disp, bracket=[lb, ub], method='brentq').root
                    Btable[number_of_modes(Vdiscr[indexV],m):number_of_modes(Vmax,m),indexV] = np.nan
     
                f.clf()
                a1 = f.add_subplot(221)
                for indexM in range(number_of_modes(Vmax,m)):
                    a1.plot(Vdiscr, Btable[indexM,:], label=r'$\mu=%2d$'%indexM)
                plt.xlim([0,Vmax])
                plt.ylim([0,1])
                plt.xlabel(r'$V$')
                plt.ylabel(r'$B$')
                plt.legend()
                a1bisy = a1.twinx()
                a1bisy.set_ylim(a1.get_ylim())
                a1bisy.set_ylabel(r'$n_{\rm eff}$')
                a1bisy.set_yticks([0, 0.25, 0.5, 0.75, 1])
                a1bisy.set_yticklabels([str(round(np.sqrt(epsilon_c),4)),str(round(np.sqrt(.25*(epsilon_co-epsilon_c)+epsilon_c),4)),str(round(np.sqrt(.5*(epsilon_co-epsilon_c)+epsilon_c),4)),str(round(np.sqrt(.75*(epsilon_co-epsilon_c)+epsilon_c),4)),str(round(np.sqrt(epsilon_co),4))])
                a1bis = a1.twiny()
                a1bis.set_xlim(a1.get_xlim())
                a1bis.set_xlabel(r'$a/\lambda$', labelpad=10)
                a1bis.set_xticks([0,.5*V0,V0,1.5*V0,Vmax])
                a1bis.set_xticklabels([str(round(0,3)),str(round(.5*V0/(2*np.pi*np.sqrt(epsilon_co-epsilon_c)),3)),str(round(V0/(2*np.pi*np.sqrt(epsilon_co-epsilon_c)),3)),str(round(1.5*V0/(2*np.pi*np.sqrt(epsilon_co-epsilon_c)),3)),str(round(Vmax/(2*np.pi*np.sqrt(epsilon_co-epsilon_c)),3))])
                a1bis.tick_params(direction='out', pad=0)
    
                nm = number_of_modes(V0,m)
                if nm > 0 and nm < mu+1:
                    gui.input_error("Mode index too high. Reducing index...")
                    mu = nm-1
                    mu_string.set(mu)
                Bselect = Btable[mu,(np.abs(Vdiscr - V0)).argmin()]       
                def disp(B): # dispersion relation in the form disp=0
                    return np.sqrt(1-B)*sps.jv(m-1,V0*np.sqrt(1-B))*sps.kv(m,V0*np.sqrt(B))+np.sqrt(B)*sps.jv(m,V0*np.sqrt(1-B))*sps.kv(m-1,V0*np.sqrt(B))           
                Bselect = spo.root_scalar(disp, bracket=[Btable[mu,(np.abs(Vdiscr - V0)).argmin()-1], Btable[mu,(np.abs(Vdiscr - V0)).argmin()+1]], method='brentq').root
                a1.plot(V0,Bselect,'bo')
                neff_select = np.sqrt(Bselect*(epsilon_co-epsilon_c)+epsilon_c)
                neff_string.set(neff_select)
               
                rcore = np.linspace(0, 1, num=101, endpoint=True) # r in core in units of a
                rcladding = np.linspace(1, 3, num=201, endpoint=True) # r in cladding in units of a
                U = 2*np.pi*a*np.sqrt(epsilon_co-neff_select**2)
                W = 2*np.pi*a*np.sqrt(neff_select**2-epsilon_c)
                Rcore = sps.jv(m,U*rcore)
                if number_of_modes(V0,m) <= 0:
                    gui.input_warning("Mode does not exist or is out of plotting range for given parameters. Plotting DR only...")
                elif np.isinf(sps.kv(m,W)):
                    gui.input_warning("Mode too close to cut-off for plotting. Plotting DR only...")
                else:
                    a2 = f.add_subplot(222)
                    Rcladding = sps.kv(m,W*rcladding)*sps.jv(m,U)/sps.kv(m,W)
                    r = np.concatenate((rcore,rcladding))
                    R = np.concatenate((Rcore,Rcladding))
                    a2.plot(r, R)
                    plt.xlabel(r'$R$')
                    plt.ylabel(r'$\psi_R$')
                    plt.xlim([r[0],r[-1]])
                    a2bis = a2.twiny()
                    a2bis.set_xlim(a2.get_xlim())
                    a2bis.set_xticks([1])
                    a2bis.set_xticklabels([r'$a$'])
                    a2bis.tick_params(direction='out', pad=0)
                    a2.set_title(r'$m=$ '+str(m)+r', $\mu=$ '+str(mu)+r', $B=$ '+str(round(Bselect,4)))
                    a2.text(0.9, 0.9, r'LP$_{}$'.format(m)+r'$_{}$'.format(mu+1), verticalalignment='center', horizontalalignment='center', transform=a2.transAxes)
                
                    a3 = f.add_subplot(223, projection='polar')
                    azimuths = np.radians(np.linspace(0, 360, 200))
                    zeniths = r
                    rho, theta = np.meshgrid(zeniths, azimuths)
                    values = np.zeros((azimuths.size, zeniths.size))
                    for index in range(azimuths.size):
                        values[index,:] = R*np.cos(m*azimuths[index])
                    a3.contourf(theta, rho, values, levels=100, cmap="jet")
                    a3.set_rlabel_position(0)
                    a3.set_thetagrids([])
                    a3.set_rgrids([1],[r'$a$'],color='grey')
                    plt.yticks(fontsize=20)
                
                    if m != 0:
                        a4 = f.add_subplot(224, projection='polar')           
                        for index in range(azimuths.size):
                            values[index,:] = R*np.sin(m*azimuths[index])
                        a4.contourf(theta, rho, values, levels=100, cmap="jet")
                        a4.set_rlabel_position(0)
                        a4.set_thetagrids([])
                        a4.set_rgrids([1],[r'$a$'],color='grey')
                        plt.yticks(fontsize=20)           
    
                plt.tight_layout()
                    
    #            plt.savefig('fiber_wga.pdf',bbox_inches='tight',dpi=300, transparent=True)
    
                epsilon_co_save = epsilon_co_string.get()
                epsilon_c_save = epsilon_c_string.get()
                a_string_save = a_string.get()
                m_string_save = m_string.get()
                mu_string_save = mu_string.get()
    
                canvas.draw()
    except ValueError: gui.input_error("Unknown error. Re-initializing with previous parameters...", reinitialize)

f = plt.figure(1,[8,6])
canvas = gui.create_canvas(root,f)
mainframe = gui.create_mainframe(root)

epsilon_co_string = Tk.StringVar()
epsilon_c_string = Tk.StringVar()
a_string = Tk.StringVar()
m_string = Tk.StringVar()
mu_string = Tk.StringVar()
neff_string = Tk.StringVar()

initialize()

row = 1
row = gui.create_entry_with_latex(mainframe,r'core $\varepsilon_{\rm co} =$',epsilon_co_string,row)
row = gui.create_entry_with_latex(mainframe,r'cladding $\varepsilon_{\rm c} =$',epsilon_c_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry_with_latex(mainframe,r'core radius $a/\lambda=$',a_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry_with_latex(mainframe,r'azimuthal index $m=$',m_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_entry_with_latex(mainframe,r'mode index $\mu=$',mu_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_label_with_latex(mainframe,r'eff.\ index $n_{\rm eff}=$',neff_string,row)
row = gui.create_spacer(mainframe,row)
row = gui.create_button(mainframe,"Calculate",calculate,row)

gui.mainloop_safe_for_mac(root)